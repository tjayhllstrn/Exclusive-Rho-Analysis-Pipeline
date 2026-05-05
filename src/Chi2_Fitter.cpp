#include "Chi2_Fitter.h"
#include <iostream>
#include <typeinfo>
#include <string>
#include <stdexcept>
#include <cstring>

//constructor
Chi2_Fitter::Chi2_Fitter(const char* treename, const char* out_dir,
                const char* obs_s,
                const char* obs2_s,
                const char* in_file,
                std::vector<double> phibn_edges,
                std::vector<double> bn_edgs,
                std::vector<double> obs2bn,
                std::string fit_type,
                bool rewrite_cache,
                double minPhoEnergy)
    : TREENAME(treename),
    OUT_DIR(out_dir),
    OBS(obs_s),
    OBS2(obs2_s),
    IN_FILE(in_file),
    PHIBN_EDGES(phibn_edges),
    BN_EDGS(bn_edgs),
    OBS2BN(obs2bn),
    FIT_TYPE(fit_type),
    REWRITE_CACHE(rewrite_cache),
    CACHE_DIR("/lustre24/expphy/volatile/clas12/users/tjhellst/cache"),
    MinPhoEnergy(minPhoEnergy) //GeV
    {
        USE_BIN_AVERAGE = true; //default to using data-weighted bin averages for fitting, can be turned off to use geometric bin centers instead (for testing systematics)
        USE_UNBINNED_SIG_FIT = false; //default uses unbinned signal fit using RooFit, set fault to fit to a histogram using Root.fit
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

        // Create cache directory if it doesn't exist
        gSystem->mkdir(CACHE_DIR.c_str(), true);
        std::cout << "Cache directory: " << CACHE_DIR << " (rewrite_cache=" << REWRITE_CACHE << ")" << std::endl;

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
        Chi2_values.resize(n_obs_bins, std::vector<std::pair<double, double>>(n_phi_bins));
        N_sig_fitting_paramboxes.resize(n_obs_bins, std::vector<std::pair<TPaveText*, TPaveText*>>(n_phi_bins));
        
        // initialize SinCanvas lazily per obs2 bin
        SinCanvas = nullptr;


    }

//destructor
Chi2_Fitter::~Chi2_Fitter(){
        if (SinCanvas) {
            delete SinCanvas;
            SinCanvas = nullptr;
        }
}

void Chi2_Fitter::CalculateBinAverages(TTree* filteredTree) {
    std::cout << "  Calculating data-weighted bin averages..." << std::endl;
    
    // Clear and resize the centers vectors
    BN_CENTERS.clear();
    BN_CENTERS.reserve(BN_EDGS.size() - 1);
    
    // Calculate OBS bin averages
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        TH1F hist("temp_obs_hist", "temp", 100, BN_EDGS[i], BN_EDGS[i+1]);
        TCut bin_cut = TCut((OBS + ">" + std::to_string(BN_EDGS[i]) + " && " + OBS + "<" + std::to_string(BN_EDGS[i+1])).c_str());
        filteredTree->Draw((OBS + ">>temp_obs_hist").c_str(), bin_cut, "goff");
        
        double bin_avg = hist.GetMean();
        int entries = hist.GetEntries();
        
        // Fall back to geometric center if bin is empty
        if (entries == 0 || bin_avg == 0) {
            bin_avg = 0.5 * (BN_EDGS[i] + BN_EDGS[i+1]);
            std::cout << "    Warning: Empty OBS bin [" << BN_EDGS[i] << ", " << BN_EDGS[i+1] 
                      << "], using geometric center: " << bin_avg << std::endl;
        } else {
            std::cout << "    OBS bin [" << BN_EDGS[i] << ", " << BN_EDGS[i+1] 
                      << "]: average = " << bin_avg << " (" << entries << " entries)" << std::endl;
        }
        
        BN_CENTERS.push_back(bin_avg);
    }
    
}


void Chi2_Fitter::RunMhChi2Fit(int obs2bn_idx){
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning Mh chi2 phibinning Fit\033[0m\n"<<std::endl;
    
    //define Missing Mass cut
    double Mx_min = 0.85;
    double Mx_max = 1.05;

    //Pre-Cut the Tree for efficiency
    TCut Diphoton_cut = TCut("Mdiphoton<0.164 && 0.104<Mdiphoton");
    TCut Mx_cut = TCut(("Mx>" + std::to_string(Mx_min) + " && Mx<" + std::to_string(Mx_max)).c_str());
    TCut obs2_cut = TCut((OBS2 + ">" + std::to_string(OBS2BN[obs2bn_idx]) + " && " + OBS2 + "<" + std::to_string(OBS2BN[obs2bn_idx+1])).c_str());
    TCut MinPhoCut = TCut(("pho1_E>" + std::to_string(MinPhoEnergy) + " && " + std::to_string(MinPhoEnergy) + "<pho2_E").c_str());
    TCut pre_cut = Diphoton_cut && Mx_cut && obs2_cut && MinPhoCut;

    // Create or open cached filtered tree
    std::cout << "  Creating/loading pre-filtered tree..." << std::endl;
    
    // Generate unique filename based on input file, fit_type, obs2 bin, and MinPhoEnergy
    std::string base_filename = GetBaseFilename(IN_FILE);
    std::string cache_filename = CACHE_DIR + "/filtered_tree_" + base_filename + "_" + FIT_TYPE + "_MhFit_" + OBS2 + "bin" + 
                                 std::to_string(obs2bn_idx) + "_" + 
                                 std::to_string(OBS2BN[obs2bn_idx]) + "_" + 
                                 std::to_string(OBS2BN[obs2bn_idx+1]) + "_" +
                                 "MinPhoE_" + std::to_string(MinPhoEnergy) +
                                 ".root";
    
    TFile* tempFile = nullptr;
    TTree* filteredTree = nullptr;
    
    if (!REWRITE_CACHE && !gSystem->AccessPathName(cache_filename.c_str())) {
        // File exists and we're not rewriting, open it
        std::cout << "    Using cached tree from: " << cache_filename << std::endl;
        tempFile = TFile::Open(cache_filename.c_str(), "READ");
        filteredTree = (TTree*)tempFile->Get("filtered_tree");
    } else {
        // Create new filtered tree
        std::cout << "    Creating new filtered tree at: " << cache_filename << std::endl;
        tempFile = new TFile(cache_filename.c_str(), "RECREATE");
        filteredTree = RAW_TREE->CopyTree(pre_cut.GetTitle());
        filteredTree->SetName("filtered_tree");
        filteredTree->SetDirectory(tempFile);
        tempFile->Write();
    }

    //plot binning scheme on filteredTree
    BinningSchemePlot(filteredTree);
    if (USE_BIN_AVERAGE) {
        CalculateBinAverages(filteredTree);
    }


    //Run Fitting Procedure  for neg and pos helicity data
    std::cout << "  \033[0;32mFitting Neg Helicity Data...\033[0m" << std::endl;
    TCut* neg_hel =  new TCut(("hel==-1"));
    N_sig_neg = FitChi2(filteredTree,neg_hel,-1);

    std::cout << "  \033[0;32mFitting Pos Helicity Data...\033[0m" << std::endl;
    TCut* pos_hel =  new TCut(("hel==1"));
    N_sig_pos = FitChi2(filteredTree,pos_hel,1);

    int n_plots = BN_EDGS.size() - 1;
    int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
    int rows = (int)TMath::Ceil((double)n_plots / cols);
    if (SinCanvas) {
        delete SinCanvas;
        SinCanvas = nullptr;
    }
    std::string canvas_name = "SinCanvas_" + std::to_string(obs2bn_idx);
    SinCanvas = new TCanvas(canvas_name.c_str(), "Sin Fit Results", 1200, 800);
    SinCanvas->Divide(cols, rows);

    //calculate A using polarization and N_sig values
    A.resize(BN_EDGS.size() - 1);
    A_sig.resize(BN_EDGS.size() - 1);
    
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        A[i].resize(PHIBN_EDGES.size() - 1);
        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            obs_bin_idx = i;
            phi_bin_idx = j;
            CalcA(filteredTree, obs_bin_idx, phi_bin_idx);
        }
        //Fit A values to sin fit function to extract A_sig
        
        FitToSin(PHIBN_CENTERS,A[i],i, filteredTree);
    }
    
    
    SinCanvas->Modified();
    SinCanvas->Update();
    SinCanvas->SaveAs((OUT_DIR + "SinFits.png").c_str());
    delete SinCanvas;
    SinCanvas = nullptr;
    tempFile->Close();
    delete tempFile;
    delete neg_hel;
    delete pos_hel;
}

void Chi2_Fitter::RunMxChi2Fit(int obs2bn_idx){
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning Mx chi2 phibinning Fit\033[0m\n"<<std::endl;
    
    //define Rho cut
    double Mh_min = 0.65;
    double Mh_max = 0.9;

    //Pre-Cut the Tree for efficiency
    TCut Diphoton_cut = TCut("Mdiphoton<0.164 && 0.104<Mdiphoton");
    TCut Mh_cut = TCut(("Mh>" + std::to_string(Mh_min) + " && Mh<" + std::to_string(Mh_max)).c_str());
    TCut obs2_cut = TCut((OBS2 + ">" + std::to_string(OBS2BN[obs2bn_idx]) + " && " + OBS2 + "<" + std::to_string(OBS2BN[obs2bn_idx+1])).c_str());
    TCut MinPhoCut = TCut(("pho1_E>" + std::to_string(MinPhoEnergy) + " && " + std::to_string(MinPhoEnergy) + "<pho2_E").c_str());
    TCut pre_cut = Diphoton_cut && Mh_cut && obs2_cut && MinPhoCut;

    // Create or open cached filtered tree
    std::cout << "  Creating/loading pre-filtered tree..." << std::endl;
    
    // Generate unique filename based on input file, fit_type and obs2 bin
    std::string base_filename = GetBaseFilename(IN_FILE);
    std::string cache_filename = CACHE_DIR + "/filtered_tree_" + base_filename + "_" + FIT_TYPE + "_MxFit_" + OBS2 + "bin" + 
                                 std::to_string(obs2bn_idx) + "_" + 
                                 std::to_string(OBS2BN[obs2bn_idx]) + "_" + 
                                 std::to_string(OBS2BN[obs2bn_idx+1]) + "_" +
                                 "MinPhoE_" + std::to_string(MinPhoEnergy) +
                                 ".root";
    
    TFile* tempFile = nullptr;
    TTree* filteredTree = nullptr;
    
    if (!REWRITE_CACHE && !gSystem->AccessPathName(cache_filename.c_str())) {
        // File exists and we're not rewriting, open it
        std::cout << "    Using cached tree from: " << cache_filename << std::endl;
        tempFile = TFile::Open(cache_filename.c_str(), "READ");
        filteredTree = (TTree*)tempFile->Get("filtered_tree");
    } else {
        // Create new filtered tree
        std::cout << "    Creating new filtered tree at: " << cache_filename << std::endl;
        tempFile = new TFile(cache_filename.c_str(), "RECREATE");
        filteredTree = RAW_TREE->CopyTree(pre_cut.GetTitle());
        filteredTree->SetName("filtered_tree");
        filteredTree->SetDirectory(tempFile);
        tempFile->Write();
    }

    // Calculate bin centers (either geometric or data-weighted average)
    if (USE_BIN_AVERAGE) {
        CalculateBinAverages(filteredTree);
    }
    
    //plot binning scheme on filteredTree
    BinningSchemePlot(filteredTree);

    //Run Fitting Procedure  for neg and pos helicity data
    std::cout << "  \033[0;32mFitting Neg Helicity Data...\033[0m" << std::endl;
    TCut* neg_hel =  new TCut(("hel==-1"));
    N_sig_neg = FitChi2(filteredTree,neg_hel,-1);

    std::cout << "  \033[0;32mFitting Pos Helicity Data...\033[0m" << std::endl;
    TCut* pos_hel =  new TCut(("hel==1"));
    N_sig_pos = FitChi2(filteredTree,pos_hel,1);

    int n_plots = BN_EDGS.size() - 1;
    int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
    int rows = (int)TMath::Ceil((double)n_plots / cols);
    if (SinCanvas) {
        delete SinCanvas;
        SinCanvas = nullptr;
    }
    std::string canvas_name = "SinCanvas_" + std::to_string(obs2bn_idx);
    SinCanvas = new TCanvas(canvas_name.c_str(), "Sin Fit Results", 1200, 800);
    SinCanvas->Divide(cols, rows);

    //calculate A using polarization and N_sig values
    A.resize(BN_EDGS.size() - 1);
    A_sig.resize(BN_EDGS.size() - 1);
    
    for (size_t i = 0; i < BN_EDGS.size() - 1; ++i) {
        A[i].resize(PHIBN_EDGES.size() - 1);
        for (size_t j = 0; j < PHIBN_EDGES.size() - 1; ++j) {
            obs_bin_idx = i;
            phi_bin_idx = j;
            CalcA(filteredTree, obs_bin_idx, phi_bin_idx);
        }
        //Fit A values to sin fit function to extract A_sig
        
        FitToSin(PHIBN_CENTERS,A[i],i, filteredTree);
    }


    SinCanvas->Modified();
    SinCanvas->Update();
    SinCanvas->SaveAs((OUT_DIR + "SinFits.png").c_str());
    delete SinCanvas;
    SinCanvas = nullptr;
    tempFile->Close();
    delete tempFile;
    delete neg_hel;
    delete pos_hel;
}

std::vector<std::vector<std::pair<double, double>>> Chi2_Fitter::FitChi2(TTree* filteredTree, TCut* hel_cut, int helicity) {

    //make helicity cut
    std::string pid = std::to_string(getpid());
    TFile* tempFile2 = new TFile(("/lustre24/expphy/volatile/clas12/users/tjhellst/final_tree_temp_" + pid + ".root").c_str(), "RECREATE");
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
                if (USE_UNBINNED_SIG_FIT) {
                    N_sig_results[i][j] = Mh_sig_fit(final_tree, bin_cut, helicity);
                } else {
                    // If not using unbinned fit, fall back to histogram-based fit
                    N_sig_results[i][j] = Mh_sig_fit_histogram(final_tree, bin_cut, helicity);
                }
            }
            else if (FIT_TYPE.find("Mx") != std::string::npos){
                if (USE_UNBINNED_SIG_FIT) {
                    N_sig_results[i][j] = Mx_sig_fit(final_tree, bin_cut, helicity);
                } else {
                    std::cerr << "Warning: Histogram signal fit not implemented for Mx. Falling back to unbinned fit." << std::endl;
                    N_sig_results[i][j] = Mx_sig_fit(final_tree, bin_cut, helicity);
                }
                N_sig_results[i][j] = Mx_sig_fit(final_tree, bin_cut, helicity);
            }
            else{
                throw std::runtime_error("Error: Unknown FIT_TYPE specified: " + FIT_TYPE);}

            

        }

    }

  // Clean up final tree and temp file
  tempFile2->Close();
  gSystem->Unlink(tempFile2->GetName());
  delete tempFile2;

  return N_sig_results;

}
std::pair<double,double> Chi2_Fitter::Mh_sig_fit_histogram(TTree* binnedTree, TCut bin_cut,int helicity){
    // Define fitting bounds
    double lb = 0.5; // GeV
    double ub = 1.25; // GeV

    // Create histogram for this bin
    TH1F hist("temp_hist", "temp", 100, lb, ub);
    binnedTree->Draw("Mh>>temp_hist", bin_cut, "goff");
    
    // Fit histogram with Voigtian signal + Chebychev background
    TF1 fit_func("fit_func", "[0]*TMath::Voigt((x-[1]), [2], [3]) + [4]*([5]*x + [6]*(2*x*x - 1))", lb, ub);
    
    // Set initial parameter values and limits
    fit_func.SetParameters(hist.GetMaximum()*0.7, 0.78, 0.06, 0.15, hist.GetMaximum()*0.3,0,0); // Initial guesses
    fit_func.SetParLimits(0, 1, hist.GetMaximum()); // N_sig
    fit_func.SetParLimits(1, 0.75, 0.9); // mu
    fit_func.SetParLimits(2, 0.01, 0.09); // sigma
    fit_func.SetParLimits(3, 0.145, 0.155); // gamma
    fit_func.SetParLimits(4, 1, hist.GetMaximum()); // N_bkg
    fit_func.SetParLimits(5, -1, 1); // Cheby coeff 1
    fit_func.SetParLimits(6, -1, 1); // Cheby coeff 2

    fit_func.SetParName(0, "N_sig");
    fit_func.SetParName(1, "mu");
    fit_func.SetParName(2, "sigma");
    fit_func.SetParName(3, "gamma");
    fit_func.SetParName(4, "N_bkg");
    fit_func.SetParName(5, "c1");
    fit_func.SetParName(6, "c2");

    hist.Fit(&fit_func, "RQ"); 

    double N_sig = fit_func.GetParameter(0);
    double N_sig_err = fit_func.GetParError(0);

    PlotSigFitGraph(hist, fit_func, helicity); 

    return std::make_pair(N_sig, N_sig_err);
}

std::pair<double,double> Chi2_Fitter::Mh_sig_fit(TTree* binnedTree, TCut bin_cut,int helicity){
    //Define RooFit Variables
    double lb = 0.5; //fitting bounds
    double ub = 1.25;


    // Create unique suffix for all RooFit objects to prevent parameter contamination between bins
    std::string idx = std::to_string(obs_bin_idx) + "_" + std::to_string(phi_bin_idx) + "_" + std::to_string(helicity);

    RooRealVar Mh("Mh","Mh",lb,ub); 
    RooRealVar obs(OBS.c_str(), OBS.c_str(),0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);

    //Create Data Set for this bin
    RooDataSet binned_data(("binned_data_" + idx).c_str(), "binned_data", RooArgSet(Mh,obs,phi),
                       RooFit::Import(*binnedTree),
                       RooFit::Cut(bin_cut));
    double N_total = binned_data.numEntries();

    //Define fit parameters for signal (Gaussian)
    RooRealVar mu(("m0_" + idx).c_str(), "#mu", 0.78, 0.75, 0.9);
    // RooRealVar sigma(("sigma_" + idx).c_str(), "#sigma", 0.06, 0.001, 0.3);
    RooRealVar sigma(("sigma_" + idx).c_str(), "#sigma", 0.06, 0.01, 0.09);
    RooRealVar gamma(("gamma_" + idx).c_str(), "#gamma", 0.15, 0.145, 0.155);

    RooRealVar mu2(("m02_" + idx).c_str(), "#mu_bkg", 0.5, 0.02, 0.66);
    RooRealVar sigma2(("sigma2_" + idx).c_str(), "#sigma_bkg", 0.05, 0.4, 0.59);

    //Define fit parameters for background (Chebychev polynomial)
    RooRealVar p1(("p1_" + idx).c_str(), "p1", 0, -1, 1);
    RooRealVar p2(("p2_" + idx).c_str(), "p2", 0, -1, 1);
    // RooRealVar p3(("p3_" + idx).c_str(), "p3", 0, -1, 1);
    // RooRealVar p4(("p4_" + idx).c_str(), "p4", 0, -2, 2);
    
    RooRealVar N_sig(("Nsig_" + idx).c_str(), "N_sig", N_total*0.7, 100, N_total);
    RooRealVar N_bkg_cheby(("Nbkg_cheby_" + idx).c_str(), "N_cheby", N_total*0.3, 50, N_total);
    RooRealVar N_bkg_gauss(("Nbkg_gauss_" + idx).c_str(), "N_gauss", N_total*0.05, 10, N_total);

    
    //Define Roo Fitting Models
    RooVoigtian sig(("sig_" + idx).c_str(), "Sig", Mh, mu, gamma, sigma);
    RooChebychev bkg_chebychev(("bkg_chebychev_" + idx).c_str(), "Bkg", Mh, RooArgList(p1, p2));
    RooGaussian bkg_gaus(("bkg_gaus_" + idx).c_str(), "Bkg Gauss", Mh, mu2, sigma2);
    // RooGenericPdf background(("background" + idx).c_str(),
    //                                    "@0 + @1",
    //                                    RooArgList(bkg_chebychev, bkg_gaus));

    RooAddPdf model_ext(("model_ext_" + idx).c_str(), "Sig + Bkg", RooArgList(sig, bkg_chebychev), RooArgList(N_sig, N_bkg_cheby));
    
    if (binned_data.numEntries() < 100) {
    std::cout << "WARNING: Low Statistics" << std::endl;
    }
    
    RooFitResult* fit_results = model_ext.fitTo(binned_data,
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1)); //,RooFit::Extended(true)

    //Plot fit result. use helicity to assign graph to right entry in the pair
    const int nPoints = 200;
    double* xPoints = new double[nPoints];
    double* yTotal = new double[nPoints];
    double* ySig = new double[nPoints];
    double* yBkg = new double[nPoints];
    
    
    // Calculate bin width for proper normalization
    int bin_number = 50;
    double binWidth = (Mh.getMax() - Mh.getMin()) / bin_number;
    RooArgSet args(Mh);
    for (int i = 0; i < nPoints; i++) {
        xPoints[i] = Mh.getMin() + (Mh.getMax() - Mh.getMin()) * i / (nPoints - 1);
        Mh.setVal(xPoints[i]);
        
        // Evaluate PDFs and scale by number of events and bin width
        double sig_val = sig.getVal(args) * N_sig.getVal() * binWidth;
        double bkg_cheby_val = bkg_chebychev.getVal(args) * N_bkg_cheby.getVal() * binWidth;
        double bkg_gauss_val = bkg_gaus.getVal(args) * N_bkg_gauss.getVal() * binWidth;
        
        ySig[i] = sig_val;
        yBkg[i] =  bkg_cheby_val; //bkg_gauss_val; //
        yTotal[i] = sig_val + bkg_cheby_val; // + bkg_gauss_val;
    }
    PlotSigFitGraph(binned_data, Mh, xPoints, ySig, yBkg,yTotal, model_ext, helicity, bin_number);
    
    std::pair<double,double> N_sig_result(N_sig.getVal(), N_sig.getError());
    return N_sig_result;
}

std::pair<double,double> Chi2_Fitter::Mx_sig_fit_histogram(TTree* binnedTree, TCut bin_cut,int helicity){
    //NOT FINISHED YET
    // Define fitting bounds
    // Scale Mx min/max values based on obs_bin_idx if OBS == "t_elec"
    double Mx_min_val = 0.55;
    double Mx_max_val = 1.3;
    
    if (OBS == "t_elec") {
        int n_obs_bins = BN_EDGS.size() - 1;
        double step_fraction = (n_obs_bins > 1) ? (double)obs_bin_idx / (n_obs_bins - 1) : 0.0;
        
        double min_start = 0.7;
        double min_end = 0.5;
        double max_start = 1.6;
        double max_end = 1.2;
        
        Mx_min_val = min_start - step_fraction * (min_start - min_end);
        Mx_max_val = max_start - step_fraction * (max_start - max_end);
    }
     std::cout<<"MX Fitting range: "<<Mx_min_val<<" to "<<Mx_max_val<<std::endl;
     double lb = Mx_min_val;
     double ub = Mx_max_val;

    // Create histogram for this bin
    TH1F hist("temp_hist", "temp", 100, lb, ub);
    binnedTree->Draw("Mx>>temp_hist", bin_cut, "goff");
    
    // Fit histogram with Voigtian signal + Chebychev background
    TF1 fit_func("fit_func", "gaus((x-[1]), [2]) + [4]*([5]*x + [6]*(2*x*x - 1))", lb, ub);
    
    // Set initial parameter values and limits
    fit_func.SetParameters(hist.GetMaximum()*0.7, 0.78, 0.06, 0.15, hist.GetMaximum()*0.3,0,0); // Initial guesses
    fit_func.SetParLimits(0, 1, hist.GetMaximum()); // N_sig
    fit_func.SetParLimits(1, 0.75, 0.9); // mu
    fit_func.SetParLimits(2, 0.01, 0.09); // sigma
    fit_func.SetParLimits(3, 0.145, 0.155); // gamma
    fit_func.SetParLimits(4, 1, hist.GetMaximum()); // N_bkg
    fit_func.SetParLimits(5, -1, 1); // Cheby coeff 1
    fit_func.SetParLimits(6, -1, 1); // Cheby coeff 2

    fit_func.SetParName(0, "N_sig");
    fit_func.SetParName(1, "mu");
    fit_func.SetParName(2, "sigma");
    fit_func.SetParName(3, "gamma");
    fit_func.SetParName(4, "N_bkg");
    fit_func.SetParName(5, "c1");
    fit_func.SetParName(6, "c2");

    hist.Fit(&fit_func, "RQ"); 

    double N_sig = fit_func.GetParameter(0);
    double N_sig_err = fit_func.GetParError(0);

    PlotSigFitGraph(hist, fit_func, helicity); 

    return std::make_pair(N_sig, N_sig_err);
}

std::pair<double,double> Chi2_Fitter::Mx_sig_fit(TTree* binnedTree, TCut bin_cut,int helicity){
    //Define RooFit Variables
    // Create unique suffix for all RooFit objects to prevent parameter contamination between bins
    std::string idx = std::to_string(obs_bin_idx) + "_" + std::to_string(phi_bin_idx) + "_" + std::to_string(helicity);

    // Scale Mx min/max values based on obs_bin_idx if OBS == "t_elec"
    double Mx_min_val = 0.55;
    double Mx_max_val = 1.3;
    
    if (OBS == "t_elec") {
        int n_obs_bins = BN_EDGS.size() - 1;
        double step_fraction = (n_obs_bins > 1) ? (double)obs_bin_idx / (n_obs_bins - 1) : 0.0;
        
        double min_start = 0.7;
        double min_end = 0.5;
        double max_start = 1.6;
        double max_end = 1.2;
        
        Mx_min_val = min_start - step_fraction * (min_start - min_end);
        Mx_max_val = max_start - step_fraction * (max_start - max_end);
    }
    std::cout<<"MX Fitting range: "<<Mx_min_val<<" to "<<Mx_max_val<<std::endl;
    RooRealVar Mx("Mx", "Mx", Mx_min_val, Mx_max_val); 
    RooRealVar mu_sig(("mu_sig_" + idx).c_str(), "#mu", 0.94, 0.92, 1.00);
    RooRealVar sigma_sig(("sigma_sig_" + idx).c_str(), "#sigma", 0.06, 0.02, 0.13);
    //Define fit parameters for background (Chebychev polynomial)
    RooRealVar p1(("p1_" + idx).c_str(), "p1", 0, -2, 2);
    RooRealVar p2(("p2_" + idx).c_str(), "p2", 0, -2, 2);
    RooRealVar p3(("p3_" + idx).c_str(), "p3", 0,-2,2);
    
    int nEntries = binnedTree->GetEntries() /  BN_CENTERS.size(); //approximate number of entries per obs bin
    RooRealVar N_sig(("Nsig_" + idx).c_str(), "N_sig", nEntries*0.7, 100, nEntries*1.2);
    RooRealVar N_bkg(("Nbkg_" + idx).c_str(), "N_bkg", nEntries*0.3, 50, nEntries*1.2);
    
    RooRealVar obs(OBS.c_str(), OBS.c_str(),0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    
    //Define Roo Fitting Models
    RooGaussian sig(("sig_" + idx).c_str(), "Sig", Mx, mu_sig, sigma_sig);
    RooChebychev background(("background_" + idx).c_str(), "Bkg", Mx, RooArgList(p1, p2, p3));
    RooAddPdf model_ext(("model_ext_" + idx).c_str(), "Sig + Bkg", RooArgList(sig, background), RooArgList(N_sig, N_bkg));

    //Create Data Set for this bin
    RooDataSet binned_data(("binned_data_" + idx).c_str(), "binned_data", RooArgSet(Mx,obs,phi),
                       RooFit::Import(*binnedTree),
                       RooFit::Cut(bin_cut));

    RooFitResult* fit_results = model_ext.fitTo(binned_data,
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1),
                                               RooFit::Extended(true));
    std::pair<double,double> N_sig_result(N_sig.getVal(), N_sig.getError());
    
    // Create TGraph objects by evaluating PDFs
    const int nPoints = 200;
    double* xPoints = new double[nPoints];
    double* yTotal = new double[nPoints];
    double* ySig = new double[nPoints];
    double* yBkg = new double[nPoints];

    // Calculate bin width for proper normalization
    int bin_number = 75;
    double binWidth = (Mx.getMax() - Mx.getMin()) / bin_number;
    RooArgSet args(Mx);
    for (int i = 0; i < nPoints; i++) {
        xPoints[i] = Mx.getMin() + (Mx.getMax() - Mx.getMin()) * i / (nPoints - 1);
        Mx.setVal(xPoints[i]);
        
        // Evaluate PDFs and scale by number of events and bin width
        double sig_val = sig.getVal(args) * N_sig.getVal() * binWidth;
        double bkg_val = background.getVal(args) * N_bkg.getVal() * binWidth;
        
        ySig[i] = sig_val;
        yBkg[i] = bkg_val;
        yTotal[i] = sig_val + bkg_val;
    }
    //Plot fit result. use helicity to assign graph to right entry in the pair
    PlotSigFitGraph(binned_data, Mx,xPoints, ySig, yBkg, yTotal, model_ext, helicity, bin_number);


    return N_sig_result;
}

void Chi2_Fitter::CalcA(TTree* filteredTree, int obs_bin_idx, int phi_bin_idx){
    //Get N_sig values for this bin
    double N_sig_pos_val = N_sig_pos[obs_bin_idx][phi_bin_idx].first;
    double N_sig_neg_val = N_sig_neg[obs_bin_idx][phi_bin_idx].first;
    double N_sig_pos_err = N_sig_pos[obs_bin_idx][phi_bin_idx].second;
    double N_sig_neg_err = N_sig_neg[obs_bin_idx][phi_bin_idx].second;

    //find avg polarization for this bin
    TH1F hist1("hist1","hist1",200,0,1);
    TH1F hist2("hist2","hist2",200,0,1);
    TCut bin_cut = TCut((OBS + ">" + std::to_string(BN_EDGS[obs_bin_idx]) + " && " + OBS + "<" + std::to_string(BN_EDGS[obs_bin_idx+1]) +
                         " && phi>" + std::to_string(PHIBN_EDGES[phi_bin_idx]) + " && phi<" + std::to_string(PHIBN_EDGES[phi_bin_idx+1])).c_str());
    filteredTree->Draw("Pol>>hist1",bin_cut,"goff");
    double Pol_avg = hist1.GetMean();
    double Pol_stdv = hist1.GetStdDev();
    double Pol_N_hist = hist1.GetEntries();
    double Pol_avg_err = Pol_stdv / sqrt(Pol_N_hist);
    filteredTree->Draw("eps>>hist2",bin_cut,"goff");
    double eps_avg = hist2.GetMean();
    double eps_stdv = hist2.GetStdDev();
    double eps_N_hist = hist2.GetEntries();
    double eps_avg_err = eps_stdv / sqrt(eps_N_hist);
    double depol = sqrt(2*eps_avg*(1-eps_avg));
    double depol_err = eps_avg_err*(1-2*eps_avg) / depol;

    //Calculate A
    double denom = N_sig_pos_val + N_sig_neg_val;
    
    double A_val = (1/Pol_avg) * (1/depol) * (N_sig_pos_val - N_sig_neg_val) / denom ;
    // double A_val = (N_sig_pos_val - N_sig_neg_val) / denom ;
    
    //Propagate error (assuming independent errors on N_sig, Pol, and eps)
    double dA_dNp = (1/Pol_avg)*(1/depol)*2.0*N_sig_neg_val/(denom*denom);
    double dA_dNm = -(1/Pol_avg)*(1/depol)*2.0*N_sig_pos_val/(denom*denom);
    double dA_dpol = -(1/Pol_avg)*A_val;
    double dA_depol = -(1/depol)*A_val;
    double A_err = sqrt(pow(dA_dNp * N_sig_pos_err,2) + pow(dA_dNm * N_sig_neg_err,2)+pow(dA_dpol * Pol_avg_err,2) + pow(dA_depol * depol_err,2));
    
    // double dA_dNp = 2.0*N_sig_neg_val/(denom*denom);
    // double dA_dNm = 2.0*N_sig_pos_val/(denom*denom);
    // double A_err = sqrt(pow(dA_dNp * N_sig_pos_err,2) + pow(dA_dNm * N_sig_neg_err,2));


    //Store A value 
    std::pair<double,double> A_pair(A_val, A_err);
    A[obs_bin_idx][phi_bin_idx] = A_pair;
}

void Chi2_Fitter::FitToSin(std::vector<double>& x_vals,std::vector<std::pair<double,double>>& y,int obs_bin_idx, TTree* filteredTree){
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

    //calculate A_sig = amplitude / depolarization*Pol factor
    // TH1F hist1 = TH1F("hist1","hist1",200,0,1);
    // TH1F hist2 = TH1F("hist2","hist2",200,0,1);
    // TCut bin_cut = TCut((OBS + ">" + std::to_string(BN_EDGS[obs_bin_idx]) + " && " + OBS + "<" + std::to_string(BN_EDGS[obs_bin_idx+1])).c_str());
    // filteredTree->Draw("Pol>>hist1",bin_cut,"goff");
    // double Pol_avg = hist1.GetMean();
    // double Pol_stdv = hist1.GetStdDev();
    // double Pol_N_hist = hist1.GetEntries();
    // double Pol_avg_err = Pol_stdv / sqrt(Pol_N_hist);
    // filteredTree->Draw("eps>>hist2",bin_cut,"goff");
    // double eps_avg = hist2.GetMean();
    // double eps_stdv = hist2.GetStdDev();
    // double eps_N_hist = hist2.GetEntries();
    // double eps_avg_err = eps_stdv / sqrt(eps_N_hist);
    // double depol = sqrt(2*eps_avg*(1-eps_avg));
    // double depol_err = eps_avg_err*(1-2*eps_avg) / depol;
    
    gr->SetMarkerStyle(20);  // Common visible marker
    gr->SetMarkerSize(1);
    gr->GetXaxis()->SetTitle("#phi_{h}");
    gr->GetYaxis()->SetRangeUser(-0.5,0.5);
    gr->GetYaxis()->SetTitle("A_{LU}");
    gr->Draw("AP");
    fitFunc->Draw("same");
    
    TLatex* latex = new TLatex();
    latex->SetNDC();
    latex->SetTextColor(kBlack);
    latex->DrawLatex(0.63, 0.85, Form("%.2f<%s<%.2f", BN_EDGS[obs_bin_idx], OBS.c_str(), BN_EDGS[obs_bin_idx+1]));
    latex->DrawLatex(0.15,0.85, Form("#chi^{2}/NDF = %.2f", chi2ndf));
    latex->DrawLatex(0.57,0.15, Form("amp = %.3f #pm %.3f", amplitude, amplitude_err));

    //calculate FLU and FLU error and store final results
    // double FLU = amplitude / (depol * Pol_avg);
    // double dFLU_damp = 1.0 / (depol * Pol_avg);
    // double dFLU_dpol = -amplitude / (depol * Pol_avg * Pol_avg);
    // double dFLU_depol = -amplitude / (depol * depol * Pol_avg);
    // double FLU_err = sqrt( pow(dFLU_damp * amplitude_err,2) + pow(dFLU_dpol * Pol_avg_err,2) + pow(dFLU_depol * depol_err,2) );
    
    double FLU = amplitude;
    double FLU_err = amplitude_err;

    A_sig[obs_bin_idx] =  std::make_pair(FLU, FLU_err);

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
void Chi2_Fitter::PlotSigFitGraph(TH1F& data_hist, TF1& fit_func, int helicity){
    // Create unique index string from current bin indices
    std::string idx_str = std::to_string(obs_bin_idx) + "_" + std::to_string(phi_bin_idx);

    const int nPoints = 200;
    double* xPoints = new double[nPoints];
    double* yTotal = new double[nPoints];
    double* ySig = new double[nPoints];
    double* yBkg = new double[nPoints];

    double xmin = data_hist.GetXaxis()->GetXmin();
    double xmax = data_hist.GetXaxis()->GetXmax();
    double sig_norm = fit_func.GetParameter(0);
    double mu = fit_func.GetParameter(1);
    double sigma = fit_func.GetParameter(2);
    double gamma = fit_func.GetParameter(3);
    double bkg_norm = fit_func.GetParameter(4);
    double p1 = fit_func.GetParameter(5);
    double p2 = fit_func.GetParameter(6);

    for (int i = 0; i < nPoints; ++i) {
        xPoints[i] = xmin + (xmax - xmin) * i / (nPoints - 1);
        yTotal[i] = fit_func.Eval(xPoints[i]);
        ySig[i] = sig_norm * TMath::Voigt(xPoints[i]-mu, sigma, gamma);
        yBkg[i] = (xPoints[i] * p1 + (2*xPoints[i]*xPoints[i] - 1) * p2) * bkg_norm; 
    }

    // Clone the provided data histogram to detach ownership
    std::string hist_name = std::string("data_hist_") + idx_str;
    TH1F* data_hist_clone = (TH1F*)data_hist.Clone(hist_name.c_str());
    data_hist_clone->SetDirectory(0);
    data_hist_clone->SetTitle(data_hist.GetTitle());
    data_hist_clone->SetMarkerStyle(20);
    data_hist_clone->SetMarkerSize(0.5);
    data_hist_clone->GetXaxis()->SetTitle(data_hist.GetXaxis()->GetTitle());
    data_hist_clone->GetYaxis()->SetTitle("Events");
    data_hist_clone->SetStats(0);

    // Create graphs
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

    double chi2NDF = CalculateChi2(data_hist_clone, total_graph);

    // Legend
    std::string leg_name = "legend_" + idx_str;
    TLegend* leg = new TLegend(0.55, 0.55, 0.85, 0.85);
    leg->SetName(leg_name.c_str());
    leg->SetBorderSize(0);
    leg->AddEntry(data_hist_clone, "Data", "p");
    leg->AddEntry(total_graph, "fit", "l");
    leg->AddEntry(sig_graph, "Signal", "l");
    leg->AddEntry(bkg_graph, "Background", "l");

    // Text annotation
    std::string text_name = "text_" + idx_str;
    std::string text = (helicity == -1) ? "Helicity: -1" : "Helicity: +1";
    TLatex* txt = new TLatex();
    txt->SetName(text_name.c_str());
    txt->SetNDC(true);
    txt->SetTextSize(0.06);
    txt->SetText(0.50, 0.47, Form("#splitline{%s}{#chi^{2}/NDF: %.2f}", text.c_str(), chi2NDF));

    // Parameter box from TF1 parameters
    std::string param_box_name = "param_box_" + idx_str;
    TPaveText* param_box = new TPaveText(0.75, 0.15, 1, 0.85, "NDC");
    param_box->SetName(param_box_name.c_str());
    param_box->SetFillColor(0);
    param_box->SetBorderSize(0);
    param_box->SetTextAlign(12);
    param_box->SetTextSize(0.055);
    int npar = fit_func.GetNpar();
    for (int ip = 0; ip < npar; ++ip) {
        const char* pname = fit_func.GetParName(ip);
        double pval = fit_func.GetParameter(ip);
        double perr = fit_func.GetParError(ip);
        if (pname && strlen(pname) > 0) {
            param_box->AddText(Form("%s: %.3g #pm %.3g", pname, pval, perr));
        } else {
            param_box->AddText(Form("p%d: %.3g #pm %.3g", ip, pval, perr));
        }
    }

    // Store objects in member vectors according to helicity
    if (helicity == -1) {
        N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].second = data_hist_clone;
        N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].second = total_graph;
        N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].second = sig_graph;
        N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].second = bkg_graph;
        N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].second = leg;
        N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].second = txt;
        Chi2_values[obs_bin_idx][phi_bin_idx].second = chi2NDF;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].second = param_box;
    }
    else if (helicity == 1) {
        N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].first = data_hist_clone;
        N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].first = total_graph;
        N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].first = sig_graph;
        N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].first = bkg_graph;
        N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].first = leg;
        N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].first = txt;
        Chi2_values[obs_bin_idx][phi_bin_idx].first = chi2NDF;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].first = param_box;
    } else {
        delete[] xPoints; delete[] yTotal; delete[] ySig; delete[] yBkg;
        throw std::invalid_argument("Helicity must be -1 or 1");
    }

    // Clean up temporary arrays
    delete[] xPoints;
    delete[] yTotal;
    delete[] ySig;
    delete[] yBkg;
}

void Chi2_Fitter::PlotSigFitGraph(RooDataSet& binned_data, RooRealVar& x, double xPoints[200],
                                     double ySig[200], double yBkg[200], double yTotal[200], RooAddPdf& model_ext, int helicity, int bin_number){
  //Make ROOT objects for plotting and storing in member vectors -------------------------------------------------
  std::string idx_str = std::to_string(obs_bin_idx) + "_" + std::to_string(phi_bin_idx);

  const int nPoints = 200;
  
  // Create data histogram from RooDataSet
  std::string hist_name = "data_hist_" + idx_str;
  std::string hist_title = Form("Sig Extraction Fit %s(%.2f,%.2f)_phi(%.2f,%.2f)", 
                             OBS.c_str(), 
                             BN_EDGS[obs_bin_idx], 
                             BN_EDGS[obs_bin_idx+1],
                             PHIBN_EDGES[phi_bin_idx], 
                             PHIBN_EDGES[phi_bin_idx+1]);

  TH1F* temp_hist = (TH1F*)binned_data.createHistogram(hist_name.c_str(), x, RooFit::Binning(bin_number, x.getMin(), x.getMax()));
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
  
//   // Create TGraph objects by evaluating PDFs
//   const int nPoints = 200;
//   double* xPoints = new double[nPoints];
//   double* yTotal = new double[nPoints];
//   double* ySig = new double[nPoints];
//   double* yBkg = new double[nPoints];

  
  
//   // Calculate bin width for proper normalization
//   double binWidth = data_hist->GetBinWidth(1);
//   RooArgSet args(x);
//   for (int i = 0; i < nPoints; i++) {
//     xPoints[i] = x.getMin() + (x.getMax() - x.getMin()) * i / (nPoints - 1);
//     x.setVal(xPoints[i]);
    
//     // Evaluate PDFs and scale by number of events and bin width
//     double sig_val = sig.getVal(args) * N_sig.getVal() * binWidth;
//     double bkg_cheby_val = background_cheby.getVal(args) * N_bkg_cheby.getVal() * binWidth;
//     double bkg_gauss_val = background_gauss.getVal(args) * N_bkg_gauss.getVal() * binWidth;
    
//     ySig[i] = sig_val;
//     yBkg[i] = bkg_cheby_val + bkg_gauss_val;
//     yTotal[i] = sig_val + bkg_cheby_val + bkg_gauss_val;
//   }
  
  // Create TGraph objects
  std::string total_graph_name = "totalFit_" + idx_str;
  std::string sig_graph_name = "sigFit_" + idx_str;
  std::string bkg_graph_name = "bkgFit_" + idx_str;
  
  TGraph* total_graph = new TGraph(nPoints, xPoints, yTotal);
  total_graph->SetName(total_graph_name.c_str());
  total_graph->SetLineStyle(kDashed);
  total_graph->SetLineColor(kBlack);
  total_graph->SetLineWidth(2);

  double chi2NDF = CalculateChi2(data_hist, total_graph);
  
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
  TLegend* leg = new TLegend(0.55, 0.55, 0.85, 0.85);
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
  txt->SetText(0.50, 0.47, Form("#splitline{%s}{#chi^{2}/NDF: %.2f}", text.c_str(), chi2NDF));


  //param box
  std::string param_box_name = "param_box_" + idx_str;
  TPaveText* param_box = new TPaveText(0.75, 0.15, 1, 0.85, "NDC"); //x.getMax()+0.05,0,x.getMax() + 0.55,data_hist->GetMaximum());
  param_box->SetName(param_box_name.c_str());
  param_box->SetFillColor(0);
  param_box->SetBorderSize(0);
  param_box->SetTextAlign(12);
  param_box->SetTextSize(0.055);
  RooArgSet* params = model_ext.getParameters(binned_data);
  for (auto* arg : *params) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
    if (!var) continue;
    
    TString name = var->GetTitle();
    double val = var->getVal();
    double err = var->getError();
    if (name.Contains("N_")) {
      param_box->AddText(Form("%s: %.1e#pm%.0f", name.Data(), val, err));
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
        Chi2_values[obs_bin_idx][phi_bin_idx].second = chi2NDF;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].second = param_box;
    }
    else if (helicity == 1) {
        N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].first = data_hist;
        N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].first = total_graph;
        N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].first = sig_graph;
        N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].first = bkg_graph;
        N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].first = leg;
        N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].first = txt;
        Chi2_values[obs_bin_idx][phi_bin_idx].first = chi2NDF;
        N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].first = param_box;
    } else {
        throw std::invalid_argument("Helicity must be -1 or 1");
    } 
}

double Chi2_Fitter::CalculateChi2(TH1F* data_hist, TGraph* total_graph) {
    double chi2 = 0.0;
    int nBins = data_hist->GetNbinsX();
    int ndf = 0;
    
    for (int i = 1; i <= nBins; i++) {
        double x = data_hist->GetBinCenter(i);
        double data = data_hist->GetBinContent(i);
        double data_err = sqrt(data); // Poisson error
        
        if (data == 0) continue; // Skip empty bins
        
        // Evaluate the fit at this x position
        double fit = total_graph->Eval(x);
        
        // Calculate chi2 contribution
        chi2 += pow((data - fit) / data_err, 2);
        ndf++;
    }
    
    std::cout << "      Chi2 = " << chi2 << ", NDF = " << ndf << ", Chi2/NDF = " << chi2/ndf << std::endl;
    
    return chi2/ndf;
}

void Chi2_Fitter::PlotToCanvas_PostageStamp(std::vector<TH1F*>& data_hists,
                                            std::vector<TGraph*>& total_graphs,
                                            std::vector<TGraph*>& sig_graphs,
                                            std::vector<TGraph*>& bkg_graphs,
                                            std::vector<TLegend*>& legends,
                                            std::vector<TLatex*>& texts,
                                            std::vector<double>& Chi2NDF,
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

  //avg Chi2NDF values and place on canvas
  double avg_chi2 = 0.0;
  for (double val : Chi2NDF) {
    avg_chi2 += val;
  }
  avg_chi2 /= Chi2NDF.size();
  
  c->cd(0); // Switch to main canvas (not individual pads)
  TLatex* avg_text = new TLatex();
  avg_text->SetNDC(true);
  avg_text->SetTextSize(0.02);
  avg_text->SetTextAlign(22); // Center alignment
  avg_text->DrawLatex(0.5, 0.02, Form("Average #chi^{2}/NDF: %.2f", avg_chi2));
  
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