#include "Chi2_Fitter.h"
#include <iostream>
#include <typeinfo>
#include <string>
#include <stdexcept>

//constructor
Chi2_Fitter::Chi2_Fitter(const char* treename, const char* out_dir,
                const char* obs_s,
                const char* obs2_s,
                const char* in_file,
                std::vector<double> phibn_edges,
                std::vector<double> bn_edgs,
                std::vector<double> obs2bn,
                std::string fit_type)
    : TREENAME(treename),
    OUT_DIR(out_dir),
    OBS(obs_s),
    OBS2(obs2_s),
    IN_FILE(in_file),
    PHIBN_EDGES(phibn_edges),
    BN_EDGS(bn_edgs),
    OBS2BN(obs2bn),
    FIT_TYPE(fit_type)
    {
        // suppress RooFit messages
        RooMsgService* rms = &RooMsgService::instance();
        rms->setSilentMode(true);
        rms->setGlobalKillBelow(RooFit::ERROR);  

        //load the file and save the reference to the tree
        std::cout<<"Loading data from file: " << in_file << std::endl;
        TFile* file = TFile::Open(in_file);
        if (!file || file->IsZombie()) {
            throw std::runtime_error("Error: could not open ROOT file: " + std::string(in_file));
        }
        RAW_TREE = (TTree*)file->Get(treename);

        //populate the BN_CENTERS vector for later plotting
        BN_CENTERS.reserve(bn_edgs.size() - 1);
        for (size_t i = 0; i < bn_edgs.size() - 1; ++i) {
            double center = 0.5 * (bn_edgs[i] + bn_edgs[i + 1]);
            BN_CENTERS.push_back(center);
        }

        //populate the BN_CENTERS vector for later plotting
        PHIBN_CENTERS.reserve(phibn_edges.size() - 1);
        for (size_t i = 0; i < phibn_edges.size() - 1; ++i) {
            double center = 0.5 * (phibn_edges[i] + phibn_edges[i + 1]);
            PHIBN_CENTERS.push_back(center);
        }

        // Initialize 2D vectors for storing plot objects
        size_t n_obs_bins = BN_EDGS.size() - 1;
        size_t n_phi_bins = PHIBN_EDGES.size() - 1;
        
        N_sig_fitting_datathist.resize(n_obs_bins, std::vector<std::pair<TH1F*, TH1F*>>(n_phi_bins));
        N_sig_fitting_totalgraph.resize(n_obs_bins, std::vector<std::pair<TGraph*, TGraph*>>(n_phi_bins));
        N_sig_fitting_siggraph.resize(n_obs_bins, std::vector<std::pair<TGraph*, TGraph*>>(n_phi_bins));
        N_sig_fitting_bkggraph.resize(n_obs_bins, std::vector<std::pair<TGraph*, TGraph*>>(n_phi_bins));
        N_sig_fitting_legends.resize(n_obs_bins, std::vector<std::pair<TLegend*, TLegend*>>(n_phi_bins));
        N_sig_fitting_texts.resize(n_obs_bins, std::vector<std::pair<TLatex*, TLatex*>>(n_phi_bins));
        N_sig_fitting_paramboxes.resize(n_obs_bins, std::vector<std::pair<TPaveText*, TPaveText*>>(n_phi_bins));
        
        //initialize SinCanvas for sin fit plots
        SinCanvas = new TCanvas("SinCanvas", "Sin Fit Results", 1200, 800);


    }

//destructor
Chi2_Fitter::~Chi2_Fitter(){
        if (SinCanvas) {
        delete SinCanvas;
    }
}

void Chi2_Fitter::RunMhChi2Fit(int obs2bn_idx){
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning Mh chi2 phibinning Fit\033[0m\n"<<std::endl;
    
    //define Missing Mass cut
    double Mx_min = 0.85;
    double Mx_max = 1.05;

    //Pre-Cut the Tree for efficiency
    TCut Diphoton_cut = TCut("Mdiphoton<0.16 && 0.115<Mdiphoton");
    TCut Mx_cut = TCut(("Mx>" + std::to_string(Mx_min) + " && Mx<" + std::to_string(Mx_max)).c_str());
    TCut obs2_cut = TCut((OBS2 + ">" + std::to_string(OBS2BN[obs2bn_idx]) + " && " + OBS2 + "<" + std::to_string(OBS2BN[obs2bn_idx+1])).c_str());
    TCut pre_cut = Diphoton_cut && Mx_cut && obs2_cut;

    // Create temporary file for filtered trees to avoid memory issues
    std::cout << "  Creating pre-filtered tree..." << std::endl;
    std::string pid = std::to_string(getpid());
    TFile* tempFile = new TFile(("/tmp/filtered_tree_temp_MhFitMLM_" + pid + ".root").c_str(), "RECREATE");
    TTree* filteredTree = RAW_TREE->CopyTree(pre_cut.GetTitle());
    filteredTree->SetDirectory(tempFile); // Associate with temp file

    //plot binning scheme on filteredTree
    BinningSchemePlot(filteredTree);


    //Run Fitting Procedure  for neg and pos helicity data
    std::cout << "  Fitting Neg Helicity Data..." << std::endl;
    TCut* neg_hel =  new TCut(("hel==-1"));
    N_sig_neg = FitChi2(filteredTree,neg_hel,-1);

    std::cout << "  Fitting Pos Helicity Data..." << std::endl;
    TCut* pos_hel =  new TCut(("hel==1"));
    N_sig_pos = FitChi2(filteredTree,pos_hel,1);

    int n_plots = BN_EDGS.size() - 1;
    int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
    int rows = (int)TMath::Ceil((double)n_plots / cols);
    SinCanvas->Divide(cols, rows);

    //calculate alpha using polarization and N_sig values
    alpha.resize(BN_EDGS.size() - 1);
    A_sig.resize(BN_EDGS.size() - 1);
    
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        alpha[i].resize(PHIBN_EDGES.size() - 1);
        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            obs_bin_idx = i;
            phi_bin_idx = j;
            CalcAlpha(filteredTree, obs_bin_idx, phi_bin_idx);
        }
        //Fit alpha values to sin fit function to extract A_sig
        
        FitToSin(PHIBN_CENTERS,alpha[i],i);
    }
    
    
    SinCanvas->SaveAs((OUT_DIR + "SinFits.png").c_str());
    SinCanvas->Clear("D");
    tempFile->Close();
    delete tempFile;
    delete neg_hel;
    delete pos_hel;
}

void Chi2_Fitter::RunMxChi2Fit(int obs2bn_idx){
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning Mh chi2 phibinning Fit\033[0m\n"<<std::endl;
    
    //define Rho cut
    double Mh_min = 0.65;
    double Mh_max = 0.9;

    //Pre-Cut the Tree for efficiency
    TCut Diphoton_cut = TCut("Mdiphoton<0.16 && 0.115<Mdiphoton");
    TCut Mh_cut = TCut(("Mh>" + std::to_string(Mh_min) + " && Mh<" + std::to_string(Mh_max)).c_str());
    TCut obs2_cut = TCut((OBS2 + ">" + std::to_string(OBS2BN[obs2bn_idx]) + " && " + OBS2 + "<" + std::to_string(OBS2BN[obs2bn_idx+1])).c_str());
    TCut pre_cut = Diphoton_cut && Mh_cut && obs2_cut;

    // Create temporary file for filtered trees to avoid memory issues
    std::cout << "  Creating pre-filtered tree..." << std::endl;
    std::string pid = std::to_string(getpid());
    TFile* tempFile = new TFile(("/tmp/filtered_tree_temp_MxFitMLM_" + pid + ".root").c_str(), "RECREATE");
    TTree* filteredTree = RAW_TREE->CopyTree(pre_cut.GetTitle());
    filteredTree->SetDirectory(tempFile); // Associate with temp file

    //plot binning scheme on filteredTree
    BinningSchemePlot(filteredTree);

    //Run Fitting Procedure  for neg and pos helicity data
    std::cout << "  Fitting Neg Helicity Data..." << std::endl;
    TCut* neg_hel =  new TCut(("hel==-1"));
    N_sig_neg = FitChi2(filteredTree,neg_hel,-1);

    std::cout << "  Fitting Pos Helicity Data..." << std::endl;
    TCut* pos_hel =  new TCut(("hel==1"));
    N_sig_pos = FitChi2(filteredTree,pos_hel,1);

    int n_plots = BN_EDGS.size() - 1;
    int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
    int rows = (int)TMath::Ceil((double)n_plots / cols);
    SinCanvas->Divide(cols, rows);

    //calculate alpha using polarization and N_sig values
    alpha.resize(BN_EDGS.size() - 1);
    A_sig.resize(BN_EDGS.size() - 1);
    
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        alpha[i].resize(PHIBN_EDGES.size() - 1);
        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            obs_bin_idx = i;
            phi_bin_idx = j;
            CalcAlpha(filteredTree, obs_bin_idx, phi_bin_idx);
        }
        //Fit alpha values to sin fit function to extract A_sig
        
        FitToSin(PHIBN_CENTERS,alpha[i],i);
    }


    SinCanvas->SaveAs((OUT_DIR + "SinFits.png").c_str());
    SinCanvas->Clear("D");
    tempFile->Close();
    delete tempFile;
    delete neg_hel;
    delete pos_hel;
}

std::vector<std::vector<std::pair<double, double>>> Chi2_Fitter::FitChi2(TTree* filteredTree, TCut* hel_cut, int helicity) {

    //make helicity cut
    std::string pid = std::to_string(getpid());
    TFile* tempFile2 = new TFile(("/tmp/final_tree_temp_" + pid + ".root").c_str(), "RECREATE");
    TTree* final_tree = filteredTree->CopyTree(hel_cut->GetTitle());
    final_tree->SetDirectory(tempFile2); // Associate with temp file

    std::vector<std::vector<std::pair<double, double>>> N_sig_results;
    N_sig_results.resize(BN_EDGS.size() - 1); //resize to number of obs bins
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        N_sig_results[i].resize(PHIBN_EDGES.size() - 1); //resize to number of phi bins
        double obs_min = BN_EDGS[i];
        double obs_max = BN_EDGS[i + 1];
        std::cout << "    Fitting observable bin: (" << obs_min << "," << obs_max << ") -----------------------------------------------------------" << std::endl;

        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            obs_bin_idx = i;
            phi_bin_idx = j;
            
            std::cout<<"        Fitting "<<OBS<<" bin "<<i<<" phi bin "<<j<<std::endl;
            double phimin = PHIBN_EDGES[j];
            double phimax = PHIBN_EDGES[j + 1];

            //Call proper fitting function for this bin
            TCut bin_cut = TCut((OBS + ">" + std::to_string(obs_min) + " && " + OBS + "<" + std::to_string(obs_max) +
                                 " && phi>" + std::to_string(phimin) + " && phi<" + std::to_string(phimax)).c_str());
            if (FIT_TYPE.find("Mh") != std::string::npos){
                N_sig_results[i][j] = Mh_sig_fit(final_tree, bin_cut, helicity);
            }
            else if (FIT_TYPE.find("Mx") != std::string::npos){
                N_sig_results[i][j] = Mx_sig_fit(final_tree, bin_cut, helicity);
            }
            else{
                throw std::runtime_error("Error: Unknown FIT_TYPE specified: " + FIT_TYPE);}

            

        }

    }

  // Clean up final tree and temp file
  tempFile2->Close();
  delete tempFile2;

  return N_sig_results;

}

std::pair<double,double> Chi2_Fitter::Mh_sig_fit(TTree* binnedTree, TCut bin_cut,int helicity){
    //Define RooFit Variables
    RooRealVar Mh("Mh", "Mh", 0.4, 1.7); 
    RooRealVar mu("m_{0}", "mu", 0.78, 0.75, 0.9);
    RooRealVar sigma("#sigma", "sigma", 0.06, 0.00001, 0.1);
    RooRealVar p1("p1", "p1", 0, -5, 5);
    RooRealVar p2("p2", "p2", 0, -10, 10);
    RooRealVar p3("p3", "p3", 0, -5, 5);
    RooRealVar p4("p4", "p4", 0, -5, 5);
    
    int nEntries = binnedTree->GetEntries() /  BN_CENTERS.size(); //approximate number of entries per obs bin
    RooRealVar N_sig("N_{sig}", "N_sig", nEntries*0.7, 0, nEntries*1.2);
    RooRealVar N_bkg("N_{bkg}", "N_bkg", nEntries*0.3, 0, nEntries*1.2);
    
    RooRealVar obs(OBS.c_str(), OBS.c_str(),0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    
    //Define Roo Fitting Models
    RooGaussian sig("sig", "gaussian Fit", Mh, mu, sigma);
    RooArgList pars_pol(p1, p2, p3, p4);
    RooPolynomial background("background", "Background", Mh, pars_pol);
    RooAddPdf model_ext("model_ext", "Signal + Background", RooArgList(sig, background), RooArgList(N_sig, N_bkg));

    //Create Data Set for this bin
    RooDataSet binned_data("binned_data", "binned_data", RooArgSet(Mh,obs,phi),
                       RooFit::Import(*binnedTree),
                       RooFit::Cut(bin_cut));

    RooFitResult* fit_results = model_ext.fitTo(binned_data,
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1),
                                               RooFit::Extended(true));
    std::pair<double,double> N_sig_result(N_sig.getVal(), N_sig.getError());

    //Plot fit result. use helicity to assign graph to right entry in the pair
    
    PlotSigFitGraph(binned_data, Mh, sig, background, N_sig, N_bkg, model_ext, helicity);

    return N_sig_result;
}

std::pair<double,double> Chi2_Fitter::Mx_sig_fit(TTree* binnedTree, TCut bin_cut,int helicity){
    //Define RooFit Variables
    RooRealVar Mx("Mx", "Mx", 0.6, 1.7); 
    RooRealVar mu_sig("mu_{sig}", "mu", 0.94, 0.85, 1.2);
    RooRealVar sigma_sig("#sigma_{sig}", "sigma", 0.06, 0.01, 0.13);
    RooRealVar mu_bkg("mu_{bkg}", "mu", 2, 1.2, 3);
    RooRealVar sigma_bkg("#sigma_{bkg}", "sigma", 0.06, 0.01, 0.4);
    
    int nEntries = binnedTree->GetEntries() /  BN_CENTERS.size(); //approximate number of entries per obs bin
    RooRealVar N_sig("N_{sig}", "N_sig", nEntries*0.7, 0, nEntries*1.2);
    RooRealVar N_bkg("N_{bkg}", "N_bkg", nEntries*0.3, 0, nEntries*1.2);
    
    RooRealVar obs(OBS.c_str(), OBS.c_str(),0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    
    //Define Roo Fitting Models
    RooGaussian sig("sig", "gaussian Fit", Mx, mu_sig, sigma_sig);
    RooGaussian background("background", "Background", Mx, mu_bkg,sigma_bkg);
    RooAddPdf model_ext("model_ext", "Signal + Background", RooArgList(sig, background), RooArgList(N_sig, N_bkg));

    //Create Data Set for this bin
    RooDataSet binned_data("binned_data", "binned_data", RooArgSet(Mx,obs,phi),
                       RooFit::Import(*binnedTree),
                       RooFit::Cut(bin_cut));

    RooFitResult* fit_results = model_ext.fitTo(binned_data,
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1),
                                               RooFit::Extended(true));
    std::pair<double,double> N_sig_result(N_sig.getVal(), N_sig.getError());

    //Plot fit result. use helicity to assign graph to right entry in the pair
    
    PlotSigFitGraph(binned_data, Mx, sig, background, N_sig, N_bkg, model_ext, helicity);

    return N_sig_result;
}

void Chi2_Fitter::CalcAlpha(TTree* filteredTree, int obs_bin_idx, int phi_bin_idx){
    //Get N_sig values for this bin
    double N_sig_pos_val = N_sig_pos[obs_bin_idx][phi_bin_idx].first;
    double N_sig_neg_val = N_sig_neg[obs_bin_idx][phi_bin_idx].first;
    double N_sig_pos_err = N_sig_pos[obs_bin_idx][phi_bin_idx].second;
    double N_sig_neg_err = N_sig_neg[obs_bin_idx][phi_bin_idx].second;

    //find avg polarization for this bin
    TH1F hist1("hist1","hist1",200,0,1);
    TCut bin_cut = TCut((OBS + ">" + std::to_string(BN_EDGS[obs_bin_idx]) + " && " + OBS + "<" + std::to_string(BN_EDGS[obs_bin_idx+1]) +
                         " && phi>" + std::to_string(PHIBN_EDGES[phi_bin_idx]) + " && phi<" + std::to_string(PHIBN_EDGES[phi_bin_idx+1])).c_str());
    filteredTree->Draw("Pol>>hist1",bin_cut,"goff");
    double Pol_avg = hist1.GetMean();
    double Pol_stdv = hist1.GetStdDev();
    double Pol_N_hist = hist1.GetEntries();
    double Pol_avg_err = Pol_stdv / sqrt(Pol_N_hist);

    //Calculate alpha
    double denom = N_sig_pos_val + N_sig_neg_val;
    double alpha_val = (1/Pol_avg) * (N_sig_pos_val - N_sig_neg_val) / denom ;
    
    //Propagate error (assuming independent errors on N_sig and Pol)
    double dA_dNp = (1/Pol_avg)*2.0*N_sig_neg_val/(denom*denom);
    double dA_dNm = -(1/Pol_avg)*2.0*N_sig_pos_val/(denom*denom);
    double dA_dpol = -(1/Pol_avg)*alpha_val;
    double alpha_err = sqrt(pow(dA_dNp * N_sig_pos_err,2) + pow(dA_dNm * N_sig_neg_err,2)+pow(dA_dpol * Pol_avg_err,2));


    //Store alpha value with dummy error (to be filled later)
    std::pair<double,double> alpha_pair(alpha_val, alpha_err);
    alpha[obs_bin_idx][phi_bin_idx] = alpha_pair;
}

void Chi2_Fitter::FitToSin(std::vector<double>& x_vals,std::vector<std::pair<double,double>>& y,int obs_bin_idx){
    // Extract values and errors from pairs
    std::vector<double> y_vals;
    std::vector<double> y_errs;
    for(const auto& pair : y) {
        y_vals.push_back(pair.first);
        y_errs.push_back(pair.second);
    }
    std::vector<double> x_errs(x_vals.size(),0.0);
    SinCanvas->cd(obs_bin_idx+1);
    TGraphErrors* gr = new TGraphErrors(x_vals.size(),x_vals.data(),y_vals.data(),x_errs.data(),y_errs.data());
    gr->SetTitle("");

    TF1* fitFunc = new TF1(Form("fitFunc_%d",obs_bin_idx),"[0]*sin(x)",-3.14,3.14);
    fitFunc->SetParameter(0,0.2);

    gr->GetListOfFunctions()->Clear();
    gr->Fit(fitFunc,"RQ");
    double amplitude = fitFunc->GetParameter(0);
    double amplitude_err = fitFunc->GetParError(0);
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double chi2ndf = chi2 / ndf;

    
    gr->SetMarkerStyle(20);  // Common visible marker
    gr->SetMarkerSize(1);
    gr->GetXaxis()->SetTitle("#phi_{h}");
    gr->GetYaxis()->SetRangeUser(-0.3,0.3);
    gr->GetYaxis()->SetTitle("A_{LU}");
    gr->Draw("AP");
    fitFunc->Draw("same");
    
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextColor(kBlack);
    latex->DrawLatex(0.63, 0.85, Form("%.2f<%s<%.2f", BN_EDGS[obs_bin_idx], OBS.c_str(), BN_EDGS[obs_bin_idx+1]));
    latex->DrawLatex(0.15,0.85, Form("#chi^{2}/NDF = %.2f", chi2ndf));

    A_sig[obs_bin_idx] =  std::make_pair(amplitude,amplitude_err);

    SinCanvas->Update();

    // Write objects to ROOT file (UPDATE mode to avoid overwriting)
    std::string root_filename = OUT_DIR + "SinFits.root";
    TFile outFile(root_filename.c_str(), "UPDATE");
    gr->Write(Form("SinFit_graph_%d", obs_bin_idx));
    fitFunc->Write(Form("SinFit_func_%d", obs_bin_idx));
    latex->Write(Form("SinFit_latex_%d", obs_bin_idx));
    outFile.Close();
    
    
}
//Plotting Member Functions ---------------------------------------------------------------------------------------------------------------

void Chi2_Fitter::BinningSchemePlot(TTree* filteredTree) {
    std::cout << "  Generating Binning Scheme Plot..." << std::endl;
    TCanvas c2D;
    c2D.SetTickx();
    c2D.SetTicky();

    double bnmin = BN_EDGS.at(0) - std::abs(0.1*BN_EDGS.at(0));
    double bnmax = BN_EDGS.back() + 0.1*BN_EDGS.back();
    std::string hist_title = "binning scheme;" + OBS + " ;#phi_{h} [rad]";
    TH2F binning("binning",hist_title.c_str(),100,bnmin,bnmax,100,-3.2,3.2);

    filteredTree->Draw(("phi:" + OBS + ">>binning").c_str());
    
    
    binning.SetStats(0);
    binning.Draw();
    
       std::vector<TLine*> lines;
    for (double obsedg : BN_EDGS) {
        TLine* line = new TLine(obsedg, -3.14, obsedg, 3.14);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
    }
    for (double phiedg : PHIBN_EDGES){
        TLine* line = new TLine(bnmin, phiedg, bnmax, phiedg);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
    }
    
    c2D.Update();  // Force canvas to update with all drawn objects
    c2D.SaveAs((OUT_DIR + FIT_TYPE + "_2D" + OBS + "binningPlot.png").c_str());
    
    // Clean up
    for (TLine* line : lines) {
        delete line;
    }
}

void Chi2_Fitter::PlotToCanvas_N_sig_BarHist() {
    int n_plots = BN_EDGS.size() - 1;
    if (n_plots == 0) return;
  
    int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
    int rows = (int)TMath::Ceil((double)n_plots / cols);
    
    TCanvas* Nasym_c = new TCanvas("Nasym_c","Nasym_c",1200,800);
    Nasym_c->Divide(cols,rows);

    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        Nasym_c->cd(i+1);
        TH1D* h_pos = new TH1D(("h_pos_" + std::to_string(obs2_bin_idx) + "_" + std::to_string(i)).c_str(),"h_pos;#phi_{h};Signal Counts",PHIBN_EDGES.size()-1,PHIBN_EDGES.data());
        TH1D* h_neg = new TH1D(("h_neg_" + std::to_string(obs2_bin_idx) + "_" + std::to_string(i)).c_str(),"h_neg;#phi_{h};Signal Counts",PHIBN_EDGES.size()-1,PHIBN_EDGES.data());

        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            double N_sig_pos_val = N_sig_pos[i][j].first;
            double N_sig_neg_val = N_sig_neg[i][j].first;
            h_pos->SetBinContent(j+1, N_sig_pos_val);
            h_neg->SetBinContent(j+1, N_sig_neg_val);
        }

        h_pos->SetLineStyle(2);  
        h_pos->SetLineColor(kRed);
        h_pos->SetLineWidth(1);
        h_pos->SetStats(0); 
        h_neg->SetLineStyle(1);  
        h_neg->SetLineColor(kBlack);
        h_neg->SetLineWidth(1); 

        h_pos->Draw("HIST");
        h_neg->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.6,0.75,0.9,0.9);
        leg->AddEntry(h_pos,"Pos Helicity","l");
        leg->AddEntry(h_neg,"Neg Helicity","l");
        leg->Draw();

        std::string title = Form("%s(%.2f,%.2f) Bin", OBS.c_str(), BN_EDGS[i], BN_EDGS[i+1]);
        h_pos->SetTitle(title.c_str());
        
    }
    Nasym_c->SaveAs((OUT_DIR + OBS + "_Nsig_Asymmetry_BarHist.png").c_str());
    delete Nasym_c;
}

void Chi2_Fitter::PlotSigFitGraph(RooDataSet& binned_data, RooRealVar& x,
                                     RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext, int helicity){
  //Make ROOT objects for plotting and storing in member vectors -------------------------------------------------
  std::string idx_str = std::to_string(obs_bin_idx) + "_" + std::to_string(phi_bin_idx);
  
  // Create data histogram from RooDataSet
  std::string hist_name = "data_hist_" + idx_str;
  std::string hist_title = Form("Sig Extraction Fit %s(%.2f,%.2f)_phi(%.2f,%.2f)", 
                             OBS.c_str(), 
                             BN_EDGS[obs_bin_idx], 
                             BN_EDGS[obs_bin_idx+1],
                             PHIBN_EDGES[phi_bin_idx], 
                             PHIBN_EDGES[phi_bin_idx+1]);

  TH1F* temp_hist = (TH1F*)binned_data.createHistogram(hist_name.c_str(), x, RooFit::Binning(200, x.getMin(), x.getMax()));
  // Clone it to ensure complete independence from RooDataSet
  TH1F* data_hist = (TH1F*)temp_hist->Clone((hist_name + "_clone").c_str());
  data_hist->SetDirectory(0);  // Detach from any ROOT directory management
  delete temp_hist;  // Clean up temporary

  data_hist->SetTitle(hist_title.c_str());
  data_hist->SetMarkerStyle(20);
  data_hist->SetMarkerSize(0.5);
  data_hist->GetXaxis()->SetTitle((std::string(x.GetName()) + " (GeV)").c_str());
  data_hist->GetYaxis()->SetTitle("Events");
  data_hist->SetStats(0);
  
  // Create TGraph objects by evaluating PDFs
  const int nPoints = 200;
  double* xPoints = new double[nPoints];
  double* yTotal = new double[nPoints];
  double* ySig = new double[nPoints];
  double* yBkg = new double[nPoints];

  
  
  // Calculate bin width for proper normalization
  double binWidth = data_hist->GetBinWidth(1);
  
  RooArgSet args(x);
  for (int i = 0; i < nPoints; i++) {
    xPoints[i] = 0.4 + (1.7 - 0.4) * i / (nPoints - 1);
    x.setVal(xPoints[i]);
    
    // Evaluate PDFs and scale by number of events and bin width
    double sig_val = sig.getVal(args) * N_sig.getVal() * binWidth;
    double bkg_val = background.getVal(args) * N_bkg.getVal() * binWidth;
    
    ySig[i] = sig_val;
    yBkg[i] = bkg_val;
    yTotal[i] = sig_val + bkg_val;

  }
  
  // Create TGraph objects
  std::string total_graph_name = "totalFit_" + idx_str;
  std::string sig_graph_name = "sigFit_" + idx_str;
  std::string bkg_graph_name = "bkgFit_" + idx_str;
  
  TGraph* total_graph = new TGraph(nPoints, xPoints, yTotal);
  total_graph->SetName(total_graph_name.c_str());
  total_graph->SetLineStyle(kDashed);
  total_graph->SetLineColor(kBlack);
  total_graph->SetLineWidth(2);
  
  TGraph* sig_graph = new TGraph(nPoints, xPoints, ySig);
  sig_graph->SetName(sig_graph_name.c_str());
  sig_graph->SetLineColor(kBlue);
  sig_graph->SetLineWidth(2);
  
  TGraph* bkg_graph = new TGraph(nPoints, xPoints, yBkg);
  bkg_graph->SetName(bkg_graph_name.c_str());
  bkg_graph->SetLineColor(kRed);
  bkg_graph->SetLineWidth(2);
  
  // Clean up temporary arrays
  delete[] xPoints;
  delete[] yTotal;
  delete[] ySig;
  delete[] yBkg;

  //legend
  std::string leg_name = "legend_" + idx_str;
  TLegend* leg = new TLegend(0.45, 0.55, 0.75, 0.85);
  leg->SetName(leg_name.c_str());
  leg->SetBorderSize(0);
  leg->AddEntry(data_hist, "Data", "p");
  leg->AddEntry(total_graph, "sig + bkg", "l");
  leg->AddEntry(sig_graph, "Signal", "l");
  leg->AddEntry(bkg_graph, "Background", "l");

  //text annotation
  std::string text_name = "text_" + idx_str;
  std::string text = "";
  if(helicity == -1){
      text = "Helicity: -1";
  }
  else if(helicity == 1){
      text = "Helicity: +1";
  }
  TLatex* txt = new TLatex();
  txt->SetName(text_name.c_str());
  txt->SetNDC(true);
  txt->SetTextSize(0.06);
  txt->SetText(0.50, 0.45, text.c_str());
  //param box
  std::string param_box_name = "param_box_" + idx_str;
  TPaveText* param_box = new TPaveText(0.75, 0.15, 1, 0.85, "NDC"); //x.getMax()+0.05,0,x.getMax() + 0.55,data_hist->GetMaximum());
  param_box->SetName(param_box_name.c_str());
  param_box->SetFillColor(0);
  param_box->SetBorderSize(1);
  param_box->SetTextAlign(12);
  param_box->SetTextSize(0.06);
  RooArgSet* params = model_ext.getParameters(binned_data);
  for (auto* arg : *params) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
    if (!var) continue;
    
    TString name = var->GetName();
    double val = var->getVal();
    double err = var->getError();
    if (name.Contains("N_{")) {
      param_box->AddText(Form("%s: %.2e#pm%.2e", name.Data(), val, err));
    } else {
      param_box->AddText(Form("%s: %.2f#pm%.2f", name.Data(), val, err));
    }
  }

  // Store all objects for later use
  if (helicity == -1) {
        N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].second = data_hist;
        N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].second = total_graph;
        N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].second = sig_graph;
        N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].second = bkg_graph;
        N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].second = leg;
        N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].second = txt;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].second = param_box;
    }
    else if (helicity == 1) {
        N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].first = data_hist;
        N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].first = total_graph;
        N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].first = sig_graph;
        N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].first = bkg_graph;
        N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].first = leg;
        N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].first = txt;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].first = param_box;
    } else {
        throw std::invalid_argument("Helicity must be -1 or 1");
    } 
}

void Chi2_Fitter::PlotToCanvas_PostageStamp(std::vector<TH1F*>& data_hists,
                                            std::vector<TGraph*>& total_graphs,
                                            std::vector<TGraph*>& sig_graphs,
                                            std::vector<TGraph*>& bkg_graphs,
                                            std::vector<TLegend*>& legends,
                                            std::vector<TLatex*>& texts,
                                            std::vector<TPaveText*>& param_boxes,
                                            std::string title) {
  // Create postage stamp canvas with stored plots. 
  int n_plots = data_hists.size();
  if (n_plots == 0) return;
  
  int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
  int rows = (int)TMath::Ceil((double)n_plots / cols);
  
  TCanvas* c = new TCanvas("postage", title.c_str(), 1600, 800);
  c->Divide(cols, rows,0.01,0.01);
  for (int i = 0; i < n_plots; i++) {
    c->cd(i + 1);
    gPad->SetRightMargin(0.25); // Make room for external parameter box
    data_hists[i]->Draw();
    total_graphs[i]->Draw("L SAME");
    sig_graphs[i]->Draw("L SAME");
    bkg_graphs[i]->Draw("L SAME");
    legends[i]->Draw();
    texts[i]->Draw();

    gPad->Update(); 
    param_boxes[i]->Draw();
  }
  
  c->Update();
  
  // Save canvas as PNG
  std::string output_base = OUT_DIR + title;
  c->SaveAs((output_base + ".png").c_str());
  
  // Save canvas and all objects to ROOT file
  TFile f((output_base + ".root").c_str(), "RECREATE");
  c->Write();
  for (int i = 0; i < n_plots; i++) {
    data_hists[i]->Write(Form("hist_%d", i));
    total_graphs[i]->Write(Form("total_graph_%d", i));
    sig_graphs[i]->Write(Form("sig_graph_%d", i));
    bkg_graphs[i]->Write(Form("bkg_graph_%d", i));
    legends[i]->Write(Form("legend_%d", i));
    texts[i]->Write(Form("text_%d", i));
    param_boxes[i]->Write(Form("param_box_%d", i));
  }
  f.Close();
  
  delete c;
}

void Chi2_Fitter::MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
    const char* y_title, const char* x_title, 
    const char* title, Color_t color,  
    std::pair<double,double>* bounds=nullptr,std::vector<TGraph*>* out_graph=nullptr)
    {
    //Plot the given x and y data to a canvas and save it to OUT_DIR. THis method is meant to overlay
    //the plots over each other.

    TGraphErrors* gr = new TGraphErrors(x.size());
    for(size_t i=0;i<x.size();i++){
        gr->SetPoint(i, x[i], y[i].first);
        gr->SetPointError(i, 0, y[i].second);
    }
    gr->SetTitle(title);
    gr->SetName(title);
    gr->GetXaxis()->SetTitle(x_title);
    gr->GetYaxis()->SetTitle(y_title);
    if(bounds){
        gr->GetYaxis()->SetRangeUser(bounds->first,bounds->second);
    }
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(color);
    gr->SetLineColor(color);
    
    
    
    if(out_graph){
        out_graph->push_back(gr);
    } else {
        TCanvas* c = new TCanvas((std::string(title) + "_canvas").c_str(), (std::string(title) + "_canvas").c_str(),800,600);
        c->SetTickx();
        c->SetTicky();
        c->SetGridx();
        c->SetGridy();
        gr->Draw("AP");
        c->Update();
        c->SaveAs((OUT_DIR + "/"+ std::string(title) + ".png").c_str());
        
        TFile f((OUT_DIR + "/"+ std::string(title) + ".root").c_str(), "RECREATE");
        c->Write();   
        gr->Write();
        f.Close();
        delete c;
    }
}