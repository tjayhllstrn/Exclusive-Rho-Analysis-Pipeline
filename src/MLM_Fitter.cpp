#include "MLM_Fitter.h"
#include <iostream>
#include <typeinfo>


//constructor  
MLM_Fitter::MLM_Fitter(const char* treename, const char* out_dir,
                               const char* obs_s,
                               const char* obs2_s,
                              const char* in_file,
                            std::vector<double> bn_edgs,
                            std::vector<double> obs2bn,
                            std::string fit_type)
 : TREENAME(treename),
   OUT_DIR(out_dir),
   OBS(obs_s),
   OBS2(obs2_s),
   IN_FILE(in_file),
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
    
   }

  //destructor
  MLM_Fitter::~MLM_Fitter(){
    // std::cout<<"MLM_Fitter Destructor called"<<std::endl;
}
   
//---------------------------------------------------------------------------------------------------------------

void MLM_Fitter::RunMhFitMLM(int obs2bn_idx){
  //MLM fit to data in sigbkg and bkg regions, populating A_sigbkg and A_bkg data structures.
  std::cout << "Running Mh MLM Fit for obs2 bin: (" << OBS2BN[obs2bn_idx] << "," << OBS2BN[obs2bn_idx+1] << ")"<< std::endl;
  // Ensure vectors are sized to hold obs2 bins
  if(A_sigbkg.size() <= obs2bn_idx) {
    A_sigbkg.resize(obs2bn_idx + 1);
    A_bkg.resize(obs2bn_idx + 1);
  }

  //define sigbkg and bkg region edges
  sigbkg_min = 0.65;
  sigbkg_max = 0.9;
  bkg_min = 1.1;
  bkg_max = 1.7;

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
  
  //Run Fitting Procedure on each region
  std::cout << "  Fitting Signal+Background Region..." << std::endl;

  TCut* Mh_cut_sb =  new TCut(("Mh>" + std::to_string(sigbkg_min) + " && Mh<" + std::to_string(sigbkg_max)).c_str());
  A_sigbkg[obs2bn_idx] = FitMLM(filteredTree,Mh_cut_sb);

  std::cout << "  Fitting Background Region..." << std::endl;
  TCut* Mh_cut_b =  new TCut(("Mh>" + std::to_string(bkg_min) + " && Mh<" + std::to_string(bkg_max)).c_str());
  A_bkg[obs2bn_idx] = FitMLM(filteredTree,Mh_cut_b);

  PurityCalc(filteredTree);

  Calc_A_sig_from_A_sigbkg(obs2bn_idx);




  // Clean up filtered tree and temp file
  tempFile->Close();
  gSystem->Unlink(tempFile->GetName());
  delete tempFile;
  delete Mh_cut_sb;
  delete Mh_cut_b;
}

void MLM_Fitter::RunMxFitMLM(int obs2bn_idx){
  //MLM fit to data in sigbkg and bkg regions, populating A_sigbkg and A_bkg data structures.
  std::cout << "Running Mx MLM Fit for obs2 bin: (" << OBS2BN[obs2bn_idx] << "," << OBS2BN[obs2bn_idx+1] << ")"<< std::endl;
  // Ensure vectors are sized to hold obs2 bins
  if(A_sigbkg.size() <= obs2bn_idx) {
    A_sigbkg.resize(obs2bn_idx + 1);
    A_bkg.resize(obs2bn_idx + 1);
  }

  //define sigbkg and bkg region edges
  sigbkg_min = 0.85;
  sigbkg_max = 1.15;
  bkg_min = 1.4;
  bkg_max = 2.75;

  //define Rho Mh cut
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
  
  //Run Fitting Procedure on each region
  std::cout << "  Fitting Signal+Background Region..." << std::endl;
  
  TCut* Mx_cut_sb =  new TCut(("Mx>" + std::to_string(sigbkg_min) + " && Mx<" + std::to_string(sigbkg_max)).c_str());
  A_sigbkg[obs2bn_idx] = FitMLM(filteredTree,Mx_cut_sb);

  std::cout << "  Fitting Background Region..." << std::endl;
  TCut* Mx_cut_b =  new TCut(("Mx>" + std::to_string(bkg_min) + " && Mx<" + std::to_string(bkg_max)).c_str());
  A_bkg[obs2bn_idx] = FitMLM(filteredTree,Mx_cut_b);

  PurityCalc(filteredTree);

  Calc_A_sig_from_A_sigbkg(obs2bn_idx);

  // Clean up filtered tree and temp file
  tempFile->Close();
  gSystem->Unlink(tempFile->GetName());
  delete tempFile;
  delete Mx_cut_sb;
  delete Mx_cut_b;
}


//Private Member functions---------------------------------------------------------------------------------------------------------------
std::vector<std::pair<double,double>> MLM_Fitter::FitMLM(TTree* tree,TCut* bounds){
  //Perform the MLM fitting procedure on the given obs2 and Mh ranges
  //tree is already pre-filtered with Diphoton, Mx, and obs2 cuts
  
  //define RooFit Variables:
  RooRealVar A("A","A",-1,1);
  RooRealVar eps("eps","eps",0.6,0,1);
  RooRealVar phi("phi","phi",0,-3.25,3.25);
  RooRealVar hel("hel", "hel", -1, 1);
  RooRealVar Pol("Pol","Pol",0,1);
  RooRealVar obs(OBS.c_str(), OBS.c_str(), BN_EDGS.front(), BN_EDGS.back());
  RooRealVar Mx("Mx", "Mx", 0, 10);
  RooRealVar Mh("Mh", "Mh", 0, 3);
  

  //Define the Fitting Model
  std::string model_expr = "1+Pol*sqrt(2*eps*(1-eps))*A*hel*sin(phi)";
  RooGenericPdf model("model","model", model_expr.c_str(),
                          RooArgSet(A,eps,phi,Pol,hel));

  //now loop over observeable and calculate A for each bin.
  std::vector<std::pair<double,double>> results;
  std::string pid = std::to_string(getpid());
  TFile* tempFile2 = new TFile(("/tmp/final_tree_temp_" + pid + ".root").c_str(), "RECREATE");
  TTree* final_tree = tree->CopyTree(bounds->GetTitle());
  final_tree->SetDirectory(tempFile2); // Associate with temp file
  for(int i=0;i<BN_EDGS.size()-1;i++){
    double obs_min = BN_EDGS[i];
    double obs_max = BN_EDGS[i+1];
    std::cout << "    Fitting observable bin: (" << obs_min << "," << obs_max << ") -----------------------------------------------------------" << std::endl;

    // Create cut for this observable bin
    TCut bin_cut = TCut((OBS + ">" + std::to_string(obs_min) + " && " + OBS + "<" + std::to_string(obs_max)).c_str());
    
    // Create datasets for this bin GETS USED EVERYWHERE
    std::string dataset_name = "binned_data_bin" + std::to_string(i);
    RooDataSet binned_data(dataset_name.c_str(), dataset_name.c_str(), RooArgSet(eps,phi,obs,hel,Pol,Mh,Mx),
                       RooFit::Import(*final_tree),
                      RooFit::Cut(bin_cut));
  


    // Perform fit (suppress verbose output)
    RooFitResult* fitResult = (RooFitResult*)model.fitTo(binned_data,
                                                         RooFit::PrintLevel(-1),
                                                         RooFit::Save(true));

    // record fitted A and its error for this observable bin
    results.emplace_back(A.getVal(), A.getError());

  
    if(fitResult) delete fitResult;
    

  }

  // Clean up final tree and temp file
  tempFile2->Close();
  gSystem->Unlink(tempFile2->GetName());
  delete tempFile2;
  return results;

}

void MLM_Fitter::PurityCalc(TTree* tree){

  if(FIT_TYPE == "MhMLM"){
    PurityCalc_Mh(tree);
  }

  if (FIT_TYPE == "MxMLM"){
    PurityCalc_Mx(tree);
  }


}

void MLM_Fitter::PurityCalc_Mh(TTree* tree){
  //calculate purity by fitting to the Mh distribution
  //This implementation mirrors the Python purityCalc function
  
  std::cout << "Calculating Purity from Mh distribution..." << std::endl;
  double lb = 0.4; //fitting bounds
  double ub = 1.7;

  //Create RooRealVar for Mh fitting 
  RooRealVar Mh("Mh", "Mh", lb, ub);
  RooRealVar obs(OBS.c_str(), OBS.c_str(), BN_EDGS.front(), BN_EDGS.back());
  
  //Define fit parameters for signal (Gaussian)
  RooRealVar mu("m_{0}", "mu", 0.8, 0.6, 1);
  RooRealVar sigma("sigma_{sig}", "sigma", 0.06, 0.00001, 0.1);
  
  //Define fit parameters for background (Chebychev polynomial)
  RooRealVar p1("p1", "p1", 0, -1, 1);
  RooRealVar p2("p2", "p2", 0, -1, 1);
  RooRealVar p3("p3", "p3", 0,-1,1);
  RooRealVar p4("p4", "p4", 0,-1,1);
  
  //Create extended PDF parameters
  int nEntries = tree->GetEntries() /  BN_CENTERS.size(); //approximate number of entries per obs bin
  RooRealVar N_sig("N_{sig}", "N_sig", nEntries*0.7, 0, nEntries*1.2);
  RooRealVar N_bkg("N_{bkg}", "N_bkg", nEntries*0.3, 0, nEntries*1.2);
  
  //Create signal PDF (Gaussian)
  std::string sig_name = "sig" + std::to_string(bn_idx);
  RooGaussian sig(sig_name.c_str(), "gaussian Fit", Mh, mu, sigma);
  
  //Create background PDF (Chebychev) - reduced order
  RooArgList pars_pol(p1, p2, p3, p4);
  std::string bkg_name = "background" + std::to_string(bn_idx);
  RooChebychev background(bkg_name.c_str(), "ChebyChev", Mh, pars_pol);
  
  //Combine signal and background into extended model
  RooArgList components(sig, background);
  RooArgList yields(N_sig, N_bkg);
  std::string model_ext_name = "model_ext" + std::to_string(bn_idx);
  RooAddPdf model_ext(model_ext_name.c_str(), "Signal + Background", components, yields);

  //loop over each bin and perform fit
  for(int i=0;i<BN_EDGS.size()-1;i++){
    double obs_min = BN_EDGS[i];
    double obs_max = BN_EDGS[i+1];
    std::cout << "    Purity Fitting observable bin: (" << obs_min << "," << obs_max << ") -----------------------------------------------------------" << std::endl;
    bn_idx = i;
    // Create cut for this observable bin
    TCut bin_cut = TCut((OBS + ">" + std::to_string(obs_min) + " && " + OBS + "<" + std::to_string(obs_max)).c_str());
    
    std::string dataset_name = "binned_data_bin" + std::to_string(i);
    RooDataSet binned_data(dataset_name.c_str(), dataset_name.c_str(), RooArgSet(Mh,obs),
                       RooFit::Import(*tree),
                      RooFit::Cut(bin_cut));
  


    // Perform fit (suppress verbose output)
    RooFitResult* fit_results = model_ext.fitTo(binned_data, 
                                               RooFit::Range("fullRange"),
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1),
                                               RooFit::Extended(true));


    PurityFromFit(sig,background,N_sig,N_bkg,Mh, fit_results);
    PlotPurityGraph(binned_data,Mh,sig,background,N_sig,N_bkg,model_ext);
    
    if (fit_results) delete fit_results;
  }
}

void MLM_Fitter::PurityCalc_Mx(TTree* tree){
  //calculate purity by fitting to the Mx distribution
  //This implementation mirrors the Python purityCalc function
  
  std::cout << "Calculating Purity from Mx distribution..." << std::endl;
  double lb = 0.6; //fitting bounds
  double ub = 1.7;

  //Create RooRealVar for Mx fitting 
  RooRealVar Mx("Mx", "Mx", lb, ub);
  RooRealVar obs(OBS.c_str(), OBS.c_str(), BN_EDGS.front(), BN_EDGS.back());
  
  //Define fit parameters for signal (Gaussian)
  RooRealVar mu_sig("m_{sig}", "mu", 0.94, 0.85, 1.2);
  RooRealVar sigma_sig("sigma_{sig}", "sigma", 0.06, 0.01, 0.13);
  
  //Define fit parameters for background (Gaussian)
  // Reduced to 2nd order and tighter constraints to prevent negative PDF values
   RooRealVar mu_bkg("mu_{bkg}", "mu", 2, 1.2, 3);
   RooRealVar sigma_bkg("#sigma_{bkg}", "sigma", 0.06, 0.01, 0.4);

  //Create extended PDF parameters
  int nEntries = tree->GetEntries() /  BN_CENTERS.size(); //approximate number of entries per obs bin
  RooRealVar N_sig("N_{sig}", "N_sig", nEntries*0.7, 0, nEntries*1.2);
  RooRealVar N_bkg("N_{bkg}", "N_bkg", nEntries*0.3, 0, nEntries*1.2);
  
  //Create signal PDF (Gaussian)
  std::string sig_name = "sig" + std::to_string(bn_idx);
  RooGaussian sig(sig_name.c_str(), "gaussian Fit", Mx, mu_sig, sigma_sig);
  
  //Create background PDF (Gaussian)
  std::string bkg_name = "background" + std::to_string(bn_idx);
  RooGaussian background(bkg_name.c_str(), "Gaussian Background", Mx, mu_bkg, sigma_bkg);
  
  //Combine signal and background into extended model
  RooArgList components(sig, background);
  RooArgList yields(N_sig, N_bkg);
  std::string model_ext_name = "model_ext" + std::to_string(bn_idx);
  RooAddPdf model_ext(model_ext_name.c_str(), "Signal + Background", components, yields);

  //loop over each bin and perform fit
  for(int i=0;i<BN_EDGS.size()-1;i++){
    double obs_min = BN_EDGS[i];
    double obs_max = BN_EDGS[i+1];
    std::cout << "    Purity Fitting observable bin: (" << obs_min << "," << obs_max << ") -----------------------------------------------------------" << std::endl;
    bn_idx = i;
    // Create cut for this observable bin
    TCut bin_cut = TCut((OBS + ">" + std::to_string(obs_min) + " && " + OBS + "<" + std::to_string(obs_max)).c_str());
    
    std::string dataset_name = "binned_data_bin" + std::to_string(i);
    RooDataSet binned_data(dataset_name.c_str(), dataset_name.c_str(), RooArgSet(Mx,obs),
                       RooFit::Import(*tree),
                      RooFit::Cut(bin_cut));
  


    // Perform fit (suppress verbose output)
    RooFitResult* fit_results = model_ext.fitTo(binned_data, 
                                               RooFit::Range("fullRange"),
                                               RooFit::Save(true),
                                               RooFit::PrintLevel(-1),
                                               RooFit::Extended(true));


    PurityFromFit(sig,background,N_sig,N_bkg,Mx, fit_results);
    PlotPurityGraph(binned_data,Mx,sig,background,N_sig,N_bkg,model_ext);
    
    if (fit_results) delete fit_results;
  }
}

void MLM_Fitter::PurityFromFit(RooAbsPdf& sig, RooAbsPdf& background,
                                   RooRealVar& N_sig, RooRealVar& N_bkg,RooRealVar& x,
                                  RooFitResult* fit_results){
  //Integrate to find purity (mirroring integrate_u function from Python)
  std::cout << "  Calculating purity from integrals..." << std::endl;
                                    
  x.setRange("sigbkgRange",sigbkg_min, sigbkg_max);

  //Create integrals over the signal+background region
  RooAbsReal* bkg_int = background.createIntegral(RooArgSet(x), 
                                                   RooArgSet(x), 
                                                   "sigbkgRange");
  double bkg_perc = bkg_int->getVal();
  
  RooAbsReal* sig_int = sig.createIntegral(RooArgSet(x), 
                                            RooArgSet(x), 
                                            "sigbkgRange");
  double sig_perc = sig_int->getVal();
  
  //Calculate local number of signal and background events in the region
  double sig_N_local = sig_perc * N_sig.getVal();
  double bkg_N_local = bkg_perc * N_bkg.getVal();
  
  //Calculate purity
  double denom = bkg_N_local + sig_N_local;
  double num = sig_N_local;


  u = num / denom;
  
  //Calculate error
  double bkg_perc_err = bkg_int->getPropagatedError(*fit_results);
  double sig_perc_err = sig_int->getPropagatedError(*fit_results);
  
  double num_err = (sig_perc_err/sig_perc + N_sig.getError()/N_sig.getVal()) * num;
  double denom_err = num_err + (bkg_perc_err/bkg_perc + N_bkg.getError()/N_bkg.getVal()) * bkg_N_local;
  
  u_err = u * TMath::Sqrt(TMath::Power(num_err/num, 2) + TMath::Power(denom_err/denom, 2));
  
  std::cout << "  Purity: u = " << u << " +/- " << u_err << std::endl;
  Purities.emplace_back(u, u_err);
  

}

void MLM_Fitter::Calc_A_sig_from_A_sigbkg(int obs2bn_idx){
  // Calculate signal asymmetry from A_sigbkg = u*A_sig + (1-u)*A_bkg
  // A_sig = (1/u) * (A_sigbkg - (1-u)*A_bkg)
  
  std::cout << "  Calculating A_sig from A_sigbkg, A_bkg, and purity..." << std::endl;
  
  // Ensure A_sig vector is sized to hold obs2 bins
  if(A_sig.size() <= obs2bn_idx) {
    A_sig.resize(obs2bn_idx + 1);
  }
  
  // Loop through each observable bin
  for(size_t i = 0; i < A_sigbkg[obs2bn_idx].size(); i++){
    double A_sigbkg_val = A_sigbkg[obs2bn_idx][i].first;
    double A_sigbkg_err = A_sigbkg[obs2bn_idx][i].second;
    double A_bkg_val = A_bkg[obs2bn_idx][i].first;
    double A_bkg_err = A_bkg[obs2bn_idx][i].second;
    double u_val = Purities[i].first;
    double u_err = Purities[i].second;
    
    // Calculate A_sig value
    double A_sig_val = (1.0/u_val) * (A_sigbkg_val - (1.0 - u_val) * A_bkg_val);
    
    // Error propagation based on general formula
    double a = (1.0/u_val) * A_sigbkg_err;
    double b = ((1.0 - u_val)/u_val) * A_bkg_err;
    double c = TMath::Power(u_val, -2) * u_err * (A_bkg_val - A_sigbkg_val);
    double A_sig_err = TMath::Sqrt(TMath::Power(a, 2) + TMath::Power(b, 2) + TMath::Power(c, 2));
    
    A_sig[obs2bn_idx].emplace_back(A_sig_val, A_sig_err);
    
    std::cout << "    Bin " << i << ": A_sig = " << A_sig_val << " +/- " << A_sig_err << std::endl;
  }
}

//PLotting FUnctions---------------------------------------------------------------------------------------------------------------
void MLM_Fitter::PlotPurityGraph(RooDataSet& binned_data, RooRealVar& x,
                                     RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext){
  //Make ROOT objects for plotting and storing in member vectors -------------------------------------------------
  std::string idx_str = std::to_string(bn_idx);
  
  // Create data histogram from RooDataSet
  std::string hist_name = "data_hist_" + idx_str;
  std::string hist_title = "Purity Calculation Fit " + OBS + "(" + 
                        Form("%.2f", BN_EDGS[bn_idx]) + "," + 
                        Form("%.2f", BN_EDGS[bn_idx+1]) + ")";
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
  TLegend* leg = new TLegend(0.45, 0.55, 0.75, 0.85);
  leg->SetName(leg_name.c_str());
  leg->SetBorderSize(0);
  leg->AddEntry(data_hist, "Data", "p");
  leg->AddEntry(total_graph, "Gauss+Cheby", "l");
  leg->AddEntry(sig_graph, "Signal", "l");
  leg->AddEntry(bkg_graph, "Background", "l");

  //text annotation
  std::string text_name = "text_" + idx_str;
  TLatex* text = new TLatex();
  text->SetName(text_name.c_str());
  text->SetNDC(true);
  text->SetTextSize(0.06);
  text->SetText(0.50, 0.45, Form("#splitline{u = %.4f #pm %.4f}{#chi^{2}/NDF: %.2f}", u, u_err, chi2NDF));

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
  purity_data_hists.push_back(data_hist);
  purity_total_graphs.push_back(total_graph);
  purity_sig_graphs.push_back(sig_graph);
  purity_bkg_graphs.push_back(bkg_graph);
  purity_legends.push_back(leg);
  purity_texts.push_back(text);
  purity_param_boxes.push_back(param_box);


    
  
}

double MLM_Fitter::CalculateChi2(TH1F* data_hist, TGraph* total_graph) {
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


void MLM_Fitter::PlotToCanvas_PostageStamp(){
  // Create postage stamp canvas with stored plots
  int n_plots = purity_data_hists.size();
  if (n_plots == 0) return;
  
  int cols = (int)TMath::Ceil(TMath::Sqrt(n_plots));
  int rows = (int)TMath::Ceil((double)n_plots / cols);
  
  TCanvas* c = new TCanvas("purity_postage", "Purity Fits", 1600, 800);
  c->Divide(cols, rows,0.01,0.01);
  for (int i = 0; i < n_plots; i++) {
    c->cd(i + 1);
    gPad->SetRightMargin(0.25); // Make room for external parameter box
    purity_data_hists[i]->Draw();
    purity_total_graphs[i]->Draw("L SAME");
    purity_sig_graphs[i]->Draw("L SAME");
    purity_bkg_graphs[i]->Draw("L SAME");
    purity_legends[i]->Draw();
    purity_texts[i]->Draw();

    gPad->Update(); 
    purity_param_boxes[i]->Draw();
  }
  
  c->Update();
  
  // Save canvas as PNG
  std::string output_base = OUT_DIR + "/purity_postage_stamp";
  c->SaveAs((output_base + ".png").c_str());
  
  // Save canvas and all objects to ROOT file
  TFile f((output_base + ".root").c_str(), "RECREATE");
  c->Write();
  for (int i = 0; i < n_plots; i++) {
    purity_data_hists[i]->Write(Form("hist_%d", i));
    purity_total_graphs[i]->Write(Form("total_graph_%d", i));
    purity_sig_graphs[i]->Write(Form("sig_graph_%d", i));
    purity_bkg_graphs[i]->Write(Form("bkg_graph_%d", i));
    purity_legends[i]->Write(Form("legend_%d", i));
    purity_texts[i]->Write(Form("text_%d", i));
    purity_param_boxes[i]->Write(Form("param_box_%d", i));
  }
  f.Close();
  
  delete c;
}

void MLM_Fitter::MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
    const char* y_title, const char* x_title, 
    const char* title, Color_t color,  
    std::pair<double,double>* bounds/*=nullptr*/,std::vector<TGraph*>* out_graph/*=nullptr*/)
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
void MLM_Fitter::PlotToCanvas_overlayed(std::vector<TGraph*> graphs, const char* canvas_name, const char* canvas_title){
    //Plot the given x and y data to a canvas and save it to OUT_DIR. THis method is meant to overlay
    //the plots over each other.
    TCanvas* c = new TCanvas((std::string(canvas_name) + "_canvas").c_str(), (std::string(canvas_name) + "_canvas").c_str(),800,600); 
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    TLegend* legend = new TLegend(0.65,0.65,0.85,0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    bool first = true;
    for(auto gr : graphs){
        if(first){
            gr->Draw("AP");
            first = false;
        } else {
            gr->Draw("P SAME");
        }
        legend->AddEntry(gr, gr->GetTitle(), "p");
    }
    
    c->cd();
    legend->Draw();
    
    c->Update();
    c->SaveAs((OUT_DIR + "/"+ std::string(canvas_name) + std::string(canvas_title) + ".png").c_str());
    
    TFile f((OUT_DIR + "/"+ std::string(canvas_name) + std::string(canvas_title) + ".root").c_str(), "RECREATE");
    c->Write();   
    for(auto gr : graphs){
        gr->Write();
    }
    legend->Write();   
    f.Close();
    delete c;
}