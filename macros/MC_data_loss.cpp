#include "../src/InjectionStudy_utils.h"
#include "../src/InjectionStudy_utils.C"
#include "../src/MC_data_loss_utils.h"
#include "../src/MC_data_loss_utils.cpp"

#include <TEnv.h>
#include <TFile.h>
#include <TSystem.h>
#include <TRandom3.h>

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

//This macro takes an MC and data file and compares them. It injects an asymmetry and computes a loss function based on how well the injected MC and Data 
// match in the given phi bins. The goal of this macro is to compute a loss function for a given Asymemtry and phi binning, to be used in the optimization of the asymmetry.

static bool is_RhoSig(const inj_event& event) {
	bool pipParentIsRho = event.truepipparent_pid == 213;
	bool phoParentsParentIsRho = event.truepho1parentparent_pid == 213 && event.truepho2parentparent_pid == 213;
	bool phoParentsSame = event.truepho1parent_id == event.truepho2parent_id;
	bool phoParentsPi0 = event.truepho1parent_pid == 111 && event.truepho2parent_pid == 111;
	return pipParentIsRho && phoParentsParentIsRho && phoParentsSame && phoParentsPi0;
}

static bool has_branch(TTree* tree, const char* name) {
	return tree && tree->GetBranch(name);
}

static std::string make_bin_tag(double low, double high) {
	std::ostringstream oss;
	oss << "phi_" << std::fixed << std::setprecision(2) << low << "_" << high;
	std::string tag = oss.str();
	for (char& c : tag) {
		if (c == '.') {
			c = 'p';
		}
		if (c == '-') {
			c = 'm';
		}
	}
	return tag;
}

void MC_data_loss(const char* mc_file,
				  const char* data_file,
				  const char* signal_asym = "A10",
				  const char* config_file = "config/pippi0_RGAinbending_tbinning.txt",
				  bool rewrite_injected = true) {
	const char* tree_name = "pippi0";
	const char* phi_var = "phi";
	const char* mh_var = "Mh";
	const char* bkg_asym = "A01";

	std::string combo_name = BuildComboName(mc_file, data_file);
	std::string out_base = "out/" + combo_name + "/";
	std::string injected_dir = "/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/" + combo_name + "/injected/";
	std::string plots_dir = out_base + "plots/";
	EnsureDir(injected_dir);
	EnsureDir(plots_dir);

	std::string mc_base = StripExtension(GetBaseName(mc_file));
	std::string injected_file = injected_dir + mc_base + "_injAsym_" + signal_asym + ".root";

	if (rewrite_injected || !std::filesystem::exists(injected_file)) {
		TFile* inFile = TFile::Open(mc_file, "READ");
		if (!inFile || inFile->IsZombie()) {
			std::cerr << "Failed to open MC file: " << mc_file << std::endl;
			return;
		}

		TTree* inTree = dynamic_cast<TTree*>(inFile->Get(tree_name));
		if (!inTree) {
			std::cerr << "Could not find tree '" << tree_name << "' in MC file." << std::endl;
			inFile->Close();
			return;
		}

		const char* phi_branch = has_branch(inTree, "phi_true") ? "phi_true" : phi_var;
		const char* t_branch = has_branch(inTree, "t_elec_true") ? "t_elec_true" : "t_elec";

		if (!has_branch(inTree, phi_branch) || !has_branch(inTree, t_branch)) {
			std::cerr << "Required MC branches are missing: " << phi_branch << " or " << t_branch << std::endl;
			inFile->Close();
			return;
		}

		double phi = 0.0;
		double Pol = 0.0;
		double eps = 0.0;
		double t = 0.0;
		double truepipparent_pid = 0.0;
		double truepho1parentparent_pid = 0.0;
		double truepho2parentparent_pid = 0.0;
		double truepho1parent_id = 0.0;
		double truepho2parent_id = 0.0;
		double truepho1parent_pid = 0.0;
		double truepho2parent_pid = 0.0;

		inTree->SetBranchAddress(phi_branch, &phi);
		inTree->SetBranchAddress("Pol", &Pol);
		inTree->SetBranchAddress("eps", &eps);
		inTree->SetBranchAddress(t_branch, &t);
		inTree->SetBranchAddress("truepipparent_pid", &truepipparent_pid);
		inTree->SetBranchAddress("truepho1parentparent_pid", &truepho1parentparent_pid);
		inTree->SetBranchAddress("truepho2parentparent_pid", &truepho2parentparent_pid);
		inTree->SetBranchAddress("truepho1parent_id", &truepho1parent_id);
		inTree->SetBranchAddress("truepho2parent_id", &truepho2parent_id);
		inTree->SetBranchAddress("truepho1parent_pid", &truepho1parent_pid);
		inTree->SetBranchAddress("truepho2parent_pid", &truepho2parent_pid);

		TFile* outFile = TFile::Open(injected_file.c_str(), "RECREATE");
		if (!outFile || outFile->IsZombie()) {
			std::cerr << "Failed to create injected MC file: " << injected_file << std::endl;
			inFile->Close();
			return;
		}

		outFile->cd();
		TTree* outTree = inTree->CloneTree(0);
		if (!outTree) {
			std::cerr << "Failed to clone MC tree structure." << std::endl;
			outFile->Close();
			inFile->Close();
			return;
		}

		int hel = 0;
		TBranch* helBranch = outTree->GetBranch("hel");
		if (!helBranch) {
			std::cerr << "Output tree does not contain branch 'hel'." << std::endl;
			outFile->Close();
			inFile->Close();
			return;
		}
		helBranch->SetAddress(&hel);

		TRandom3 seededRng(12345);
		gRandom = &seededRng;

		cross_section sig_sigma(signal_asym);
		cross_section bkg_sigma(bkg_asym);

		const Long64_t nEntries = inTree->GetEntries();
		std::cout << "Injecting "<<signal_asym << " asymmetry into MC" << std::endl;
		for (Long64_t ev = 0; ev < nEntries; ++ev) {
			if (nEntries > 0 && ev % (nEntries >= 100 ? nEntries / 100 : 1) == 0) {
				int progress = static_cast<int>((100 * ev) / nEntries);
				std::cout << progress << "% complete.\r" << std::flush;
        	}

			inTree->GetEntry(ev);
			inj_event event(phi, Pol, eps, t,
							truepipparent_pid,
							truepho1parentparent_pid,
							truepho2parentparent_pid,
							truepho1parent_id,
							truepho2parent_id,
							truepho1parent_pid,
							truepho2parent_pid);

			if (is_RhoSig(event)) {
				hel = static_cast<int>(std::round(sig_sigma.AssignmentRule(event)));
			} else {
				hel = static_cast<int>(std::round(bkg_sigma.AssignmentRule(event)));
			}

			outTree->Fill();
		}

		outTree->Write();
		outFile->Close();
		inFile->Close();
		std::cout << "Injected MC file saved: " << injected_file << std::endl;
	}

	TEnv env;
	env.ReadFile(config_file, kEnvLocal);
	// Swap this config file to change the default phi binning.
	std::vector<double> phi_edges = parse_csv_to_doubles(env.GetValue("phibn_edges", ""));
	if (phi_edges.size() < 2) {
		std::cerr << "Invalid phi binning in config: " << config_file << std::endl;
		return;
	}

	TFile* dataFile = TFile::Open(data_file, "READ");
	if (!dataFile || dataFile->IsZombie()) {
		std::cerr << "Failed to open data file: " << data_file << std::endl;
		return;
	}
	TTree* dataTree = dynamic_cast<TTree*>(dataFile->Get(tree_name));
	if (!dataTree) {
		std::cerr << "Could not find tree '" << tree_name << "' in data file." << std::endl;
		dataFile->Close();
		return;
	}

	TFile* mcInjFile = TFile::Open(injected_file.c_str(), "READ");
	if (!mcInjFile || mcInjFile->IsZombie()) {
		std::cerr << "Failed to open injected MC file: " << injected_file << std::endl;
		dataFile->Close();
		return;
	}
	TTree* mcTree = dynamic_cast<TTree*>(mcInjFile->Get(tree_name));
	if (!mcTree) {
		std::cerr << "Could not find tree '" << tree_name << "' in injected MC file." << std::endl;
		mcInjFile->Close();
		dataFile->Close();
		return;
	}

	std::string out_root = out_base + "mh_compare_phi.root";
	TFile* outFile = TFile::Open(out_root.c_str(), "RECREATE");
	if (!outFile || outFile->IsZombie()) {
		std::cerr << "Failed to create output file: " << out_root << std::endl;
		mcInjFile->Close();
		dataFile->Close();
		return;
	}

	filter data_filter(phi_var, mh_var);
	data_filter.SetMhRange(0.5, 1.3);
	data_filter.SetExtraCut("0.104 < Mdiphoton&&Mdiphoton < 0.164 && 0.85 < Mx&&Mx < 1.05 &&pho1_E>0.2&&pho2_E>0.2");  // Pi0 cut, Mx cut, Min E pho cut

	plotter mh_plotter(mh_var, 100, 0.5, 1.3);

	std::vector<TH1F*> data_hists;
	std::vector<TH1F*> mc_hists;
	std::vector<Chi2Result> chi2_vals;

	TTree* chi2_tree = new TTree("chi2_summary", "chi2_summary");
	double phi_min = 0.0;
	double phi_max = 0.0;
	double chi2_val = 0.0;
	double chi2_ndf = 0.0;
	double mc_scale = 1.0;
	int ndf = 0;
	chi2_tree->Branch("phi_min", &phi_min, "phi_min/D");
	chi2_tree->Branch("phi_max", &phi_max, "phi_max/D");
	chi2_tree->Branch("chi2", &chi2_val, "chi2/D");
	chi2_tree->Branch("ndf", &ndf, "ndf/I");
	chi2_tree->Branch("chi2_ndf", &chi2_ndf, "chi2_ndf/D");
	chi2_tree->Branch("mc_scale", &mc_scale, "mc_scale/D");

	for (size_t i = 0; i + 1 < phi_edges.size(); ++i) {
		std::cout << "Calculating loss for phi bin " << i << std::endl;
		phi_min = phi_edges[i];
		phi_max = phi_edges[i + 1];
		std::string bin_tag = make_bin_tag(phi_min, phi_max);
		std::string data_name = "data_Mh_" + bin_tag;
		std::string mc_name = "mc_Mh_" + bin_tag;

		outFile->cd();
		TH1F* h_data = mh_plotter.MakeMhHist(data_name, "RGA Rho+ Invariant Mass");
		TH1F* h_mc = mh_plotter.MakeMhHist(mc_name, "MC Mh");
		h_data->SetDirectory(outFile);
		h_mc->SetDirectory(outFile);

		TCut cut = data_filter.BuildBinCut(phi_min, phi_max);
		mh_plotter.FillHist(dataTree, h_data, cut);
		mh_plotter.FillHist(mcTree, h_mc, cut);

		mc_scale = mh_plotter.NormalizeToData(h_mc, h_data);
		Chi2Result chi2 = mh_plotter.ComputeChi2(h_data, h_mc);
		chi2.mc_scale = mc_scale;

		chi2_val = chi2.chi2;
		ndf = chi2.ndf;
		chi2_ndf = chi2.chi2_ndf;
		chi2_tree->Fill();

		h_data->Write();
		h_mc->Write();

		data_hists.push_back(h_data);
		mc_hists.push_back(h_mc);
		chi2_vals.push_back(chi2);
	}

	chi2_tree->Write();

	TCanvas* canvas = mh_plotter.MakePostageStamp(data_hists, mc_hists, phi_edges, chi2_vals,
												  "MC vs Data Mh per phi bin", "Mh_postage_stamp");

	
	TH1F* chi2_hist = new TH1F("chi2_ndf_phi", "#chi^{2}/ndf vs #phi", static_cast<int>(phi_edges.size() - 1), phi_edges.data());
	chi2_hist->SetDirectory(outFile);
	chi2_hist->SetMarkerStyle(20);
	chi2_hist->SetMarkerSize(0.9);
	chi2_hist->SetLineColor(kBlack);
	chi2_hist->SetMarkerColor(kBlack);
	for (size_t i = 0; i < chi2_vals.size(); ++i) {
		chi2_hist->SetBinContent(static_cast<int>(i + 1), chi2_vals[i].chi2_ndf);
	}
	chi2_hist->Write();

	TCanvas* chi2_canvas = new TCanvas("chi2_phi_canvas", "Chi2 per phi bin", 800, 600);
	chi2_hist->GetXaxis()->SetTitle("#phi");
	chi2_hist->GetYaxis()->SetTitle("#chi^{2}/ndf");
	chi2_hist->SetStats(0);
	chi2_hist->Draw("E1");
	chi2_canvas->Write();
	chi2_canvas->SaveAs((plots_dir + "Chi2_vs_phi.png").c_str());


	if (canvas) {
		outFile->cd();
		canvas->Write();
		canvas->SaveAs((plots_dir + "Mh_compare_phi.png").c_str());
	}

	outFile->Close();
	mcInjFile->Close();
	dataFile->Close();

	std::cout << "Outputs written to: " << out_base << std::endl;
}
