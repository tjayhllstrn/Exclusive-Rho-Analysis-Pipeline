// Implementation of MC/data loss study utilities.

#include "MC_data_loss_utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <sstream>

#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TString.h>

filter::filter(const char* phi_var, const char* mh_var)
	: phi_var_(phi_var),
	  mh_var_(mh_var),
	  extra_cut_(""),
	  use_mh_range_(false),
	  mh_min_(0.0),
	  mh_max_(0.0) {}

void filter::SetPhiVar(const char* name) {
	phi_var_ = name ? name : "phi";
}

void filter::SetMhVar(const char* name) {
	mh_var_ = name ? name : "Mh";
}

void filter::SetMhRange(double min_val, double max_val) {
	use_mh_range_ = true;
	mh_min_ = min_val;
	mh_max_ = max_val;
}

void filter::ClearMhRange() {
	use_mh_range_ = false;
}

void filter::SetExtraCut(const std::string& cut) {
	extra_cut_ = cut;
}

TCut filter::BuildBaseCut() const {
	std::ostringstream oss;
	bool has_clause = false;

	if (use_mh_range_) {
		oss << mh_var_ << ">" << mh_min_ << " && " << mh_var_ << "<" << mh_max_;
		has_clause = true;
	}

	if (!extra_cut_.empty()) {
		if (has_clause) {
			oss << " && ";
		}
		oss << "(" << extra_cut_ << ")";
	}

	return TCut(oss.str().c_str());
}

TCut filter::BuildPhiCut(double phi_min, double phi_max) const {
	std::ostringstream oss;
	oss << phi_var_ << ">" << phi_min << " && " << phi_var_ << "<" << phi_max;
	return TCut(oss.str().c_str());
}

TCut filter::BuildBinCut(double phi_min, double phi_max) const {
	TCut base_cut = BuildBaseCut();
	TCut phi_cut = BuildPhiCut(phi_min, phi_max);
	return base_cut && phi_cut;
}

const std::string& filter::PhiVar() const { return phi_var_; }
const std::string& filter::MhVar() const { return mh_var_; }
bool filter::HasMhRange() const { return use_mh_range_; }
double filter::MhMin() const { return mh_min_; }
double filter::MhMax() const { return mh_max_; }

plotter::plotter(const std::string& mh_var, int mh_bins, double mh_min, double mh_max)
	: mh_var_(mh_var),
	  mh_bins_(mh_bins),
	  mh_min_(mh_min),
	  mh_max_(mh_max) {}

TH1F* plotter::MakeMhHist(const std::string& name, const std::string& title) const {
	TH1F* hist = new TH1F(name.c_str(), title.c_str(), mh_bins_, mh_min_, mh_max_);
	hist->Sumw2();
	return hist;
}

void plotter::FillHist(TTree* tree, TH1F* hist, const TCut& cut) const {
	if (!tree || !hist) {
		return;
	}
	hist->Reset();
	std::string draw_expr = mh_var_ + ">>" + std::string(hist->GetName());
	tree->Draw(draw_expr.c_str(), cut, "goff");
}

double plotter::NormalizeToData(TH1F* mc_hist, const TH1F* data_hist) const {
	if (!mc_hist || !data_hist) {
		return 1.0;
	}
	double mc_int = mc_hist->Integral();
	double data_int = data_hist->Integral();
	if (mc_int <= 0.0 || data_int <= 0.0) {
		return 1.0;
	}
	double scale = data_int / mc_int;
	mc_hist->Scale(scale);
	return scale;
}

Chi2Result plotter::ComputeChi2(const TH1F* data_hist, const TH1F* mc_hist) const {
	Chi2Result result{0.0, 0, 0.0, 1.0};
	if (!data_hist || !mc_hist) {
		return result;
	}

	int nbins = data_hist->GetNbinsX();
	int used_bins = 0;
	double chi2 = 0.0;

	for (int i = 1; i <= nbins; ++i) {
		double data_val = data_hist->GetBinContent(i);
		double mc_val = mc_hist->GetBinContent(i);
		double data_err = data_hist->GetBinError(i);
		double mc_err = mc_hist->GetBinError(i);

		if (data_err <= 0.0) {
			data_err = std::sqrt(std::max(0.0, data_val));
		}
		if (mc_err <= 0.0) {
			mc_err = std::sqrt(std::max(0.0, mc_val));
		}

		double var = data_err * data_err + mc_err * mc_err;
		if (var <= 0.0) {
			continue;
		}
		double diff = data_val - mc_val;
		chi2 += diff * diff / var;
		++used_bins;
	}

	int ndf = (used_bins > 0) ? (used_bins - 1) : 0;
	result.chi2 = chi2;
	result.ndf = ndf;
	result.chi2_ndf = (ndf > 0) ? (chi2 / ndf) : 0.0;
	return result;
}

TCanvas* plotter::MakePostageStamp(const std::vector<TH1F*>& data_hists,
								   const std::vector<TH1F*>& mc_hists,
								   const std::vector<double>& phi_edges,
								   const std::vector<Chi2Result>& chi2_vals,
								   const std::string& title,
								   const std::string& canvas_name) const {
	size_t n_plots = data_hists.size();
	if (n_plots == 0) {
		return nullptr;
	}

	int cols = static_cast<int>(TMath::Ceil(TMath::Sqrt(static_cast<double>(n_plots))));
	int rows = static_cast<int>(TMath::Ceil(static_cast<double>(n_plots) / cols));

	TCanvas* canvas = new TCanvas(canvas_name.c_str(), title.c_str(), 1200, 800);
	canvas->Divide(cols, rows);

	for (size_t i = 0; i < n_plots; ++i) {
		canvas->cd(static_cast<int>(i + 1));
		if (!data_hists[i] || !mc_hists[i]) {
			continue;
		}

		StyleDataHist(data_hists[i]);
		StyleMcHist(mc_hists[i]);

		double max_val = std::max(data_hists[i]->GetMaximum(), mc_hists[i]->GetMaximum());
		data_hists[i]->SetMaximum(max_val * 1.2);
		data_hists[i]->Draw("E1");
		mc_hists[i]->Draw("HIST SAME");
		data_hists[i]->Draw("E1 SAME");

		TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
		leg->AddEntry(data_hists[i], "Data", "lep");
		leg->AddEntry(mc_hists[i], "MC (scaled)", "l");
		leg->SetBorderSize(0);
		leg->Draw();

		TLatex latex;
		latex.SetNDC();
		latex.SetTextSize(0.045);
		if (phi_edges.size() >= i + 2) {
			latex.DrawLatex(0.12, 0.85,
							Form("phi in [%.2f, %.2f]", phi_edges[i], phi_edges[i + 1]));
		}
		if (i < chi2_vals.size()) {
			latex.DrawLatex(0.12, 0.78,
							Form("#chi^{2}/ndf = %.2f", chi2_vals[i].chi2_ndf));
		}
	}

	canvas->Update();
	return canvas;
}

const std::string& plotter::MhVar() const { return mh_var_; }
int plotter::MhBins() const { return mh_bins_; }
double plotter::MhMin() const { return mh_min_; }
double plotter::MhMax() const { return mh_max_; }

void plotter::StyleDataHist(TH1F* hist) const {
	hist->SetStats(0);
	hist->GetXaxis()->SetTitle(mh_var_.c_str());
	hist->SetMarkerStyle(20);
	hist->SetMarkerSize(0.8);
	hist->SetLineColor(kBlack);
	hist->SetMarkerColor(kBlack);
}

void plotter::StyleMcHist(TH1F* hist) const {
	hist->SetLineColor(kRed + 1);
	hist->SetLineWidth(2);
}

std::vector<double> parse_csv_to_doubles(const std::string& s) {
	std::vector<double> out;
	std::istringstream iss(s);
	std::string token;
	while (std::getline(iss, token, ',')) {
		size_t b = token.find_first_not_of(" \t");
		if (b == std::string::npos) {
			continue;
		}
		size_t e = token.find_last_not_of(" \t");
		token = token.substr(b, e - b + 1);
		try {
			out.push_back(std::stod(token));
		} catch (...) {
		}
	}
	return out;
}

std::string GetBaseName(const std::string& path) {
	std::filesystem::path p(path);
	return p.filename().string();
}

std::string StripExtension(const std::string& filename) {
	std::filesystem::path p(filename);
	return p.stem().string();
}

std::string BuildComboName(const std::string& mc_file, const std::string& data_file) {
	std::string mc_base = StripExtension(GetBaseName(mc_file));
	std::string data_base = StripExtension(GetBaseName(data_file));
	return "LossComp_" + mc_base + "_vs_" + data_base;
}

void EnsureDir(const std::string& path) {
	if (path.empty()) {
		return;
	}
	std::filesystem::create_directories(path);
}
