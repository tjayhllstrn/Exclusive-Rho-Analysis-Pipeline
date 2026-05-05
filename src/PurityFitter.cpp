#include "MLM_Fitter.h"
#include <iostream>

PurityFitter::PurityFitter(TTree* input_tree,std::string out_dir){
    //constructor
    input_data = input_tree;
    OUT_DIR = out_dir;
    RooMsgService* rms = &RooMsgService::instance();
    rms->setSilentMode(true);
    rms->setGlobalKillBelow(RooFit::ERROR);
}

PurityFitter::~PurityFitter(){
    //destructor
}

void PurityFitter::FitPurity(){
  //This implementation mirrors the Python purityCalc function

  double lb = 0.08; //fitting bounds
  double ub = 0.21;

  //Create RooRealVar for Mdiphoton fitting 
  RooRealVar Mdiphoton("Mdiphoton", "Mdiphoton", lb, ub);
  RooRealVar pho1_E("pho1_E", "pho1_E", minPhoE, 100);
  RooRealVar pho2_E("pho2_E", "pho2_E", minPhoE, 100);
  
  //Define fit parameters for signal (Voigtian)
  RooRealVar mu("m_{0}", "#mu", 0.13, 0.1, 0.16);
  RooRealVar sigma("sigma_{sig}", "#sigma", 0.01, 0.001, 0.02);
//   RooRealVar gamma("#gamma", "#gamma", 0.008, 0.007, 0.009); FWHM for pi0 is near 8 eV, which is basically 0 so don't need voigtian
  
  //Define fit parameters for background (Chebychev polynomial)
  RooRealVar p1("p1", "p1", 0, -1, 1);
  RooRealVar p2("p2", "p2", 0, -1, 1);
  // RooRealVar p3("p3", "p3", 0,-1,1);
  
  //Create extended PDF parameters
  int nEntries = input_data->GetEntries();
  RooRealVar N_sig("N_{sig}", "N_{sig}", nEntries*0.8, 0, nEntries*1.2);
  RooRealVar N_bkg("N_{bkg}", "N_{bkg}", nEntries*0.2, 0, nEntries*1.2);
  
  //Create signal PDF (Gaussian)
  RooGaussian sig("sig", "Signal", Mdiphoton, mu, sigma);
  
  //Create background PDF (Chebychev) 
  RooChebychev background("background", "Bkg Cheby", Mdiphoton, RooArgList(p1, p2));
  
  //Combine signal and background into extended model
  RooArgList components(sig, background);
  RooArgList yields(N_sig, N_bkg);
  RooAddPdf model_ext("model_ext", "Signal + Background", components, yields);

  //Create data set and perform fit
  RooDataSet cut_data("cut_data", "cut_data", RooArgSet(Mdiphoton,pho1_E,pho2_E),
                        RooFit::Import(*input_data));
                        
  double N_total = cut_data.numEntries();
  N_sig.setMax(N_total);
  N_sig.setVal(N_total*0.8);
  N_sig.setMax(N_total);
  N_bkg.setVal(N_total*0.2);

// Perform fit (suppress verbose output)
  RooFitResult* fit_results = model_ext.fitTo(cut_data,
                                            RooFit::Save(true),
                                            RooFit::PrintLevel(-1),
                                            RooFit::Extended(true));


  PurityFromFit(sig,background,N_sig,N_bkg,Mdiphoton, fit_results);
  PlotPurityGraph(cut_data,Mdiphoton,sig,background,N_sig,N_bkg,model_ext);
    
  if (fit_results) delete fit_results;
}

void PurityFitter::PurityFromFit(RooAbsPdf& sig, RooAbsPdf& background,
                                   RooRealVar& N_sig, RooRealVar& N_bkg,RooRealVar& x,
                                  RooFitResult* fit_results){
  //Integrate to find purity (mirroring integrate_u function from Python)
  std::cout << "  Calculating purity from integrals..." << std::endl;

  //Create integrals over the signal+background region
  RooAbsReal* bkg_int = background.createIntegral(RooArgSet(x), 
                                                   RooArgSet(x));
  double bkg_perc = bkg_int->getVal();
  
  RooAbsReal* sig_int = sig.createIntegral(RooArgSet(x), 
                                            RooArgSet(x));
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

void PurityFitter::PlotPurityGraph(RooDataSet& binned_data, RooRealVar& x,
                                     RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext){
  std::string idx_str = std::to_string(minPhoE).substr(0,4);
  // Create data histogram from RooDataSet
  std::string hist_name = std::string("data_hist_") + idx_str;
  std::string hist_title = std::string("Purity Calculation Fit E_{#gamma} > ") + idx_str;
  TH1F* temp_hist = (TH1F*)binned_data.createHistogram(hist_name.c_str(), x, RooFit::Binning(75, x.getMin(), x.getMax()));

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
  data_hist->GetYaxis()->SetRangeUser(0, data_hist->GetMaximum() * 1.2);  // Set y-axis min to 0
  
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
    xPoints[i] = x.getMin() + (x.getMax() - x.getMin()) * i / (nPoints - 1);
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
  TLegend* leg = new TLegend(0.5, 0.55, 0.85, 0.85);
  leg->SetName(leg_name.c_str());
  leg->SetBorderSize(0);
  leg->AddEntry(data_hist, "Data", "p");
  leg->AddEntry(total_graph, "Sig+Bkg", "l");
  leg->AddEntry(sig_graph, "Sig", "l");
  leg->AddEntry(bkg_graph, "Bkg", "l");

  //text annotation
  std::string text_name = "text_" + idx_str;
  TLatex* text = new TLatex();
  text->SetName(text_name.c_str());
  text->SetNDC(true);
  text->SetTextSize(0.06);
  text->SetText(0.5, 0.45, Form("#splitline{u = %.2f #pm %.3f}{#chi^{2}/NDF: %.2f}", u, u_err, chi2NDF));

  //param box
  std::string param_box_name = "param_box_" + idx_str;
  TPaveText* param_box = new TPaveText(0.75, 0.15, 1, 0.85, "NDC"); //x.getMax()+0.05,0,x.getMax() + 0.55,data_hist->GetMaximum());
  param_box->SetName(param_box_name.c_str());
  param_box->SetFillColor(0);
  param_box->SetBorderSize(1);
  param_box->SetTextAlign(12);
  param_box->SetTextSize(0.05);
  RooArgSet* params = model_ext.getParameters(binned_data);
  for (auto* arg : *params) {
    RooRealVar* var = dynamic_cast<RooRealVar*>(arg);
    if (!var) continue;
    
    TString name = var->GetTitle();
    double val = var->getVal();
    double err = var->getError();
    if (name.Contains("N_{")) {
      param_box->AddText(Form("%s: %.1e#pm%.0f", name.Data(), val, err));
    } else {
      param_box->AddText(Form("%s: %.3f#pm%.3f", name.Data(), val, err));
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

double PurityFitter::CalculateChi2(TH1F* data_hist, TGraph* total_graph) {
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


void PurityFitter::PlotToCanvas_PostageStamp(){
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
    purity_data_hists[i]->GetYaxis()->SetRangeUser(0, purity_data_hists[i]->GetMaximum() * 1.2);  // Force y-axis to start at 0
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
  std::string output_base = OUT_DIR + "/purity_postage_stamp" + pre_cut_name;
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

void PurityFitter::MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
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