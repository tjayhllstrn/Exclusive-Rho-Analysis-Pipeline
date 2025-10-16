#include "PhiBinnedFitRunner.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
//something to look into : RooWorkspace that stores all the objects in one place and they get called by their name only

std::pair<double,double> FitToSin(std::vector<double> x_vals,std::vector<double> y_vals,std::vector<double> y_errs,int i,std::string outputDir,const char* obs_str,double obsmin,double obsmax,TCanvas* SinCanvas){
    std::vector<double> x_errs(x_vals.size(),0.0);
    SinCanvas->cd(i+1);
    TGraphErrors* gr = new TGraphErrors(x_vals.size(),x_vals.data(),y_vals.data(),x_errs.data(),y_errs.data());
    gr->SetTitle("");

    TF1* fitFunc = new TF1(Form("fitFunc_%d",i),"[0]*sin(x)",-3.14,3.14);
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
    latex->DrawLatex(0.63, 0.85, Form("%.2f<%s<%.2f", obsmin,obs_str,obsmax));
    latex->DrawLatex(0.15,0.85, Form("#chi^{2}/NDF = %.2f", chi2ndf));

    return std::make_pair(amplitude,amplitude_err);
    
}
void PlotSigFit(RooAddPdf model_ext, RooDataSet DS,RooPlot* frame,std::string hel,int i, int j,std::string outputDir,const char* obs_str,double obsmin,double obsmax,double phimin,double phimax,const char* flag,const char* obs2_str,double obs2min,double obs2max){
    using namespace RooFit;
    DS.plotOn(frame, MarkerSize(0.5), Name("data"));
    model_ext.plotOn(frame, LineStyle(kDashed), LineColor(kBlack), Name("totalFit"));
    model_ext.plotOn(frame, Components("background"), LineColor(kRed), Name("bkgFit"));
    model_ext.plotOn(frame, Components("sig"), LineColor(kBlue), Name("sigFit"));
    frame->SetTitle(Form("%s hel %sbin %d phibin %d",hel.c_str(),obs_str,(int)i,(int)j));

    TCanvas* c = new TCanvas();
    c->Divide(2,1);
    TPad* pad1 = (TPad*)c->cd(1);
    pad1->SetPad(0.0,0.1,0.9,1.0);
    frame->Draw();
    TLatex* obsbinLabel = new TLatex(0.65,0.48,Form("%.2f < %s < %.2f",obsmin,obs_str,obsmax));
    obsbinLabel->SetNDC(true);
    obsbinLabel->SetTextSize(0.04);
    obsbinLabel->SetTextColor(kBlack);
    
    TLatex* phibinLabel = new TLatex(0.65,0.43,Form("%.2f < #phi_{h} < %.2f",phimin,phimax));
    phibinLabel->SetNDC(true);
    phibinLabel->SetTextSize(0.04);
    phibinLabel->SetTextColor(kBlack);

    TLegend* leg = new TLegend(0.65,0.55,0.85,0.85);
    leg->AddEntry(frame->findObject("data"),"Data","p");
    leg->AddEntry(frame->findObject("totalFit"),"sig+bkg Fit","l");
    leg->AddEntry(frame->findObject("sigFit"),"signal","l");
    leg->AddEntry(frame->findObject("bkgFit"),"Background","l");
    leg->SetBorderSize(0);

    leg->Draw();
    obsbinLabel->Draw();
    phibinLabel->Draw();
    
    TPad* pad2 = (TPad*)c->cd(2);
    pad2->SetPad(0.815,0.0,1.0,1.0);
    TPaveText* param_box = new TPaveText(0.05, 0.2, 0.95, 0.8, "NDC");
    param_box->SetFillColor(0);
    param_box->SetTextAlign(22);
    param_box->SetTextSize(0.1);

    RooArgSet* params = model_ext.getParameters(DS);
    TIterator* it = params->createIterator();
    RooRealVar* p = nullptr;

    while ((p = (RooRealVar*)it->Next())) {
        TString name = p->GetName();
        double val = p->getVal();
        double val_err = p->getError();

        TString text;
        if (name == "N_{bkg}" || name == "N_{sig}") {
             text = Form("%s: %.0e#pm%.0e", name.Data(), val, val_err);
        } else {
            text = Form("%s: %.2f#pm%.2f", name.Data(), val, val_err);
        }

        param_box->AddText(text);
    }
    
    param_box->Draw();
    c->SaveAs(Form("%s/phiBinningFits_%s/%s_fit_%s%d_phi%d_%s%.2f-%.2f.png", outputDir.c_str(),flag,hel.c_str(),obs_str,(int)i, (int)j,obs2_str,obs2min,obs2max));
}


std::vector<double> PhiBinnedFitRunner::Run(const std::vector<double>& bn_edges,const std::vector<double>& phibn_edges, const std::vector<std::string>& inputFiles, const std::string& outputDir,const char* obs_str,const char* treeName, const char* obs2_str, double obs2min, double obs2max) {
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning chi2 phibinning Fit\033[0m\n"<<std::endl;
    // Load TChain
    TChain uncut(treeName);
    for (const auto& str : inputFiles) {
        uncut.Add(str.c_str());
    }
    
    const char* cut_str;
    if (strcmp(treeName,"pippi0")==0){
        cut_str = Form("Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx<1.05 && %f<%s&&%s<%f",obs2min,obs2_str,obs2_str,obs2max);
    }

    else{
        cut_str = Form("0.85<Mx && Mx<1.05&& %f<%s&&%s<%f",obs2min,obs2_str,obs2_str,obs2max);
    }
    TTree* chain = uncut.CopyTree(cut_str);
    
    RooRealVar Mh("Mh", "Mh", 0.4, 1.7); 
    Mh.setRange("fullRange", 0.4, 1.4); //0.4-1.7
    
    RooRealVar mu("m_{0}", "mu", 0.78, 0.75, 0.9);
    RooRealVar sigma("#sigma", "sigma", 0.06, 0.00001, 0.1);
    RooRealVar p1("p1", "p1", 0, -5, 5);
    RooRealVar p2("p2", "p2", 0, -10, 10);
    RooRealVar p3("p3", "p3", 0, -5, 5);
    RooRealVar p4("p4", "p4", 0, -5, 5);
    
    RooRealVar N_sig("N_{sig}", "N_sig", 10000, 0, chain->GetEntries());
    RooRealVar N_bkg("N_{bkg}", "N_bkg", 10000, 0, chain->GetEntries());
    
    RooRealVar obs(obs_str, obs_str,0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    RooRealVar hel("hel", "hel", -1, 1);
    
    RooGaussian sig("sig", "gaussian Fit", Mh, mu, sigma);
    RooArgList pars_pol(p1, p2, p3, p4);
    RooPolynomial background("background", "Background", Mh, pars_pol);
    RooAddPdf model_ext("model_ext", "Signal + Background", RooArgList(sig, background), RooArgList(N_sig, N_bkg));
    
    TCanvas c2D;
    c2D.SetTickx();
    c2D.SetTicky();

    double bnmin = bn_edges.at(0) - std::abs(0.1*bn_edges.at(0));
    double bnmax = bn_edges.back() + 0.1*bn_edges.back();
    TH2F binning("binning",Form("binning scheme;%s ;#phi_{h} [rad]",obs_str),100,bnmin,bnmax,100,-3.2,3.2);

    chain->Draw(Form("phi:%s>>binning",obs_str));
    
    
    binning.SetStats(0);
    binning.Draw();
    
    std::vector<TLine*> lines;
    for (double obsedg : bn_edges) {
        TLine* line = new TLine(obsedg, -3.14, obsedg, 3.14);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    for (double phiedg : phibn_edges){
        TLine* line = new TLine(bnmin, phiedg, bnmax, phiedg);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    c2D.SaveAs(Form("%s/MhChi2_2D%sbinningPlot_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));

    //Loop over bins
    std::vector<double> results;
    std::vector<double> phibn_avgVals;
    std::vector<double> depolarizations;
    std::vector<double> depolarizations_err;
    TCanvas* SinCanvas = new TCanvas("SinCanvas","Sin fitting",800,800);
    int grid_n = (int) std::round(sqrt(bn_edges.size()-1));
    SinCanvas->Divide(grid_n,grid_n);
    TCanvas* Nasym_c = new TCanvas("Nasym_c","Nasym_c",800,800);
    Nasym_c->Divide(grid_n,grid_n);
    
    for (size_t i = 0; i < bn_edges.size() - 1; ++i) {
        std::cout<<"--------------------------------------------------------------------------------"<<std::endl;
        std::vector<double> Aphi;
        std::vector<double> Aphi_errs;
        std::vector<double> N_pos;
        std::vector<double> N_neg;
        phibn_avgVals.clear();
        double obsmin = bn_edges[i];
        double obsmax = bn_edges[i + 1];
        for (size_t j = 0; j < phibn_edges.size() - 1; ++j) {
            std::cout<<"\033[0;32mFitting "<<obs_str<<"bin "<<i<<" phi bin "<<j<<"\033[0m\n"<<std::endl;
            double phimin = phibn_edges[j];
            double phimax = phibn_edges[j + 1];

            TString cut = Form("%s > %f && %s < %f && phi > %f && phi < %f && hel == -1", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            RooDataSet neg_DS("neg_DS", "neg_DS", chain, RooArgSet(Mh, obs, phi, hel), cut);

            cut = Form("%s > %f && %s < %f && phi > %f && phi < %f && hel == 1", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            RooDataSet pos_DS("pos_DS", "pos_DS", chain, RooArgSet(Mh, obs, phi, hel), cut);

            //calculate phi location of the average number of events
            TH1F* hphi = new TH1F(Form("h_%d_%d",i,j),"histogram of phi in this bin",200,-3.14,3.14);
            cut = Form("%s > %f && %s < %f && phi > %f && phi < %f", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            chain->Draw(Form("phi>>h_%d_%d",i,j),cut,"goff");
            double mean = hphi->GetMean();
            phibn_avgVals.push_back(mean);

            //fit Mh to separate sig and bkg and calculate N_sig for pos and neg helicities
            model_ext.fitTo(neg_DS, Range("fullRange"), Save(), PrintLevel(-1), Verbose(false),Warnings(false),Extended(kTRUE));
            double N_sig_neg = N_sig.getVal();
            double N_sig_neg_err = N_sig.getError();
            RooPlot* frame_neg = Mh.frame(0.4, 1.7);
            PlotSigFit(model_ext, neg_DS,frame_neg,"neg",i,j,outputDir,obs_str,obsmin,obsmax,phimin,phimax,"Mh",obs2_str,obs2min,obs2max);
            delete frame_neg;
            
            //long N_bkg_neg = N_bkg.getVal();
            model_ext.fitTo(pos_DS, Range("fullRange"), Save(), PrintLevel(-1), Verbose(false),Warnings(false),Extended(kTRUE));
            double N_sig_pos = N_sig.getVal();
            double N_sig_pos_err = N_sig.getError();
            RooPlot* frame_pos = Mh.frame(0.4, 1.7);
            PlotSigFit(model_ext, pos_DS,frame_pos,"pos",i,j,outputDir,obs_str,obsmin,obsmax,phimin,phimax,"Mh",obs2_str,obs2min,obs2max);
            delete frame_pos;

            //beam polarization taken to be the avg of polarizations in the sample
            TH1F hist1("hist1","hist1",200,0,1);
            chain->Draw("Pol>>hist1",Form("%f<%s && %s<%f",obsmin,obs_str,obs_str,obsmax),"goff");
            double Pol_avg = hist1.GetMean();
            double Pol_stdv = hist1.GetStdDev();
            double Pol_N_hist = hist1.GetEntries();
            double Pol_avg_err = Pol_stdv / sqrt(Pol_N_hist);
            
            //calculate A for this phi bin 
            double denom = N_sig_pos + N_sig_neg;
            double A = (1/Pol_avg)*(N_sig_pos - N_sig_neg) / denom;
            double dA_dNp = (1/Pol_avg)*2.0*N_sig_neg/(denom*denom);
            double dA_dNm = -(1/Pol_avg)*2.0*N_sig_pos/(denom*denom);
            double dA_dpol = -(1/Pol_avg)*A;
            double A_err = sqrt(pow(dA_dNp * N_sig_pos_err,2) + pow(dA_dNm * N_sig_neg_err,2)+pow(dA_dpol * Pol_avg_err,2));
            Aphi.push_back(A);
            Aphi_errs.push_back(A_err);
            N_pos.push_back(N_sig_pos);
            N_neg.push_back(N_sig_neg);
        
        }
        auto res = FitToSin(phibn_avgVals,Aphi,Aphi_errs,i,outputDir,obs_str,obsmin,obsmax,SinCanvas);
        double alpha = res.first;
        double alpha_err = res.second;
        
        //define depolarization and error
        TH1F hist("hist","hist",200,0,1);
        chain->Draw("eps>>hist",Form("%f<%s && %s<%f",obsmin,obs_str,obs_str,obsmax),"goff");
        double eps_avg = hist.GetMean();
        double eps_stdv = hist.GetStdDev();
        double N_hist = hist.GetEntries();
        double eps_avg_err = eps_stdv / sqrt(N_hist);

        double avg_depolarization = sqrt(2*eps_avg*(1-eps_avg));
        double avg_depolarization_err = eps_avg_err * (abs(1-2*eps_avg)/sqrt(2*eps_avg*(1-eps_avg)));
        
        double FLU_UU = alpha / avg_depolarization;
        double FLU_UU_err = FLU_UU*sqrt(pow(alpha_err/alpha,2)+pow(avg_depolarization_err/avg_depolarization,2));
        results.push_back(FLU_UU);
        results.push_back(FLU_UU_err);
        
        /*
        TGraph* gr3 = new TGraph(phibn_avgVals.size(),phibn_avgVals.data(),N_pos.data());
        TGraph* gr4 = new TGraph(phibn_avgVals.size(),phibn_avgVals.data(),N_neg.data());
        gr4->SetLineStyle(2);  // Common visible marker
        gr4->SetLineColor(kRed);
        gr4->SetLineWidth(1.2); 
        gr3->SetLineStyle(1);  // Common visible marker
        gr3->SetLineColor(kBlack);
        gr3->SetLineWidth(1.2); 
        */
        Nasym_c->cd(i+1);
        TH1D* h_pos = new TH1D("h_pos","h_pos;#phi_{h};Signal Counts",phibn_edges.size()-1,phibn_edges.data());
        TH1D* h_neg = new TH1D("h_neg","h_neg;#phi_{h};Signal Counts",phibn_edges.size()-1,phibn_edges.data());
        for (int k = 0; k<phibn_edges.size()-1; k++){
            h_pos->SetBinContent(k+1,N_pos[k]);
            h_neg->SetBinContent(k+1,N_neg[k]);
        }
        h_pos->SetLineStyle(2);  
        h_pos->SetLineColor(kRed);
        h_pos->SetLineWidth(1.2); 
        h_neg->SetLineStyle(1);  
        h_neg->SetLineColor(kBlack);
        h_neg->SetLineWidth(1.2); 
        
        TLegend* leg = new TLegend(0.60,0.75,0.9,0.9);
        leg->AddEntry(h_pos,"Pos Helicity");
        leg->AddEntry(h_neg,"Neg Helicity");

        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextColor(kBlack);
        

        h_neg->SetStats(0);
        h_neg->SetTitle("");
        h_neg->Draw();
        h_pos->Draw("SAME");
        latex->DrawLatex(0.35, 0.85, Form("%.2f<%s<%.2f", obsmin,obs_str,obsmax));
        leg->Draw();
        
        
        
    }
    Nasym_c->SaveAs(Form("%s/Mh%sNasym_grid_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));
    SinCanvas->SaveAs(Form("%s/MhPhiBinningSinFits_%s_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));
    delete SinCanvas;
    delete Nasym_c;
    return results;
    
}
//------------------------------------------------------------------------------------------------------------------------------------------------------


std::vector<double> PhiBinnedFitRunner::Run_mxFit(const std::vector<double>& bn_edges,const std::vector<double>& phibn_edges, const std::vector<std::string>& inputFiles, const std::string& outputDir,const char* obs_str,const char* treeName, const char* obs2_str, double obs2min, double obs2max) {
    using namespace RooFit;
    std::cout<<"\033[0;32mRunning chi2 phibinning Fit\033[0m\n"<<std::endl;
    // Load TChain
    TChain* uncut = new TChain(treeName);
    for (const auto& str : inputFiles) {
        uncut->Add(str.c_str());
    }
    
    std::string cut_str;
    cut_str = Form("Mdiphoton<0.16 && 0.115<Mdiphoton && 0.65<Mh&&Mh<0.9&& %f<%s&&%s<%f",obs2min,obs2_str,obs2_str,obs2max);
    // RooRealVar Mh("Mh", "Mh", 0.65, 0.9); 
    // RooRealVar Mdiphoton("Mdiphoton", "Mdiphoton", 0.115, 0.65); 
    // RooRealVar obs2(obs2_str, obs2_str,0); 
    TFile temp_file(Form("%s/Mxtempcut_tree.root",outputDir.c_str()),"RECREATE");
    TTree* cut_tree = uncut->CopyTree(cut_str.c_str());
    double N_sig_max = cut_tree->GetEntries();
    double N_bkg_max = cut_tree->GetEntries();
    TCanvas c2D;
    c2D.SetTickx();
    c2D.SetTicky();

    double bnmin = bn_edges.at(0) - std::abs(0.1*bn_edges.at(0));
    double bnmax = bn_edges.back() + 0.1*bn_edges.back();
    
    TH2F binning("binning",Form("binning scheme;%s ;#phi_{h} [rad]",obs_str),100,bnmin,bnmax,100,-3.2,3.2);
    cut_tree->Draw(Form("phi:%s>>binning",obs_str));
    
    binning.SetStats(0);
    binning.Draw();
    
    std::vector<TLine*> lines;
    for (double obsedg : bn_edges) {
        TLine* line = new TLine(obsedg, -3.14, obsedg, 3.14);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    for (double phiedg : phibn_edges){
        TLine* line = new TLine(bnmin, phiedg, bnmax, phiedg);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    c2D.SaveAs(Form("%s/MxChi2_2D%sbinningPlot_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));
    cut_tree->Write();
    temp_file.Close();

    
    RooRealVar Mx("Mx", "Mx", 0.6, 1.7); 
    Mx.setRange("fullRange", 0.6, 1.7); 
    
    RooRealVar mu_sig("mu_{sig}", "mu", 0.94, 0.85, 1.2);
    RooRealVar sigma_sig("#sigma_{sig}", "sigma", 0.06, 0.01, 0.13);
    RooRealVar mu_bkg("mu_{bkg}", "mu", 2, 1.2, 3);
    RooRealVar sigma_bkg("#sigma_{bkg}", "sigma", 0.06, 0.01, 0.4);
    
    RooRealVar N_sig("N_{sig}", "N_sig", 10000, 0, N_sig_max);
    RooRealVar N_bkg("N_{bkg}", "N_bkg", 10000, 900, N_bkg_max); //high min bc when N_bkg=0, the fit fails
    
    RooRealVar obs(obs_str, obs_str,0);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    RooRealVar hel("hel", "hel", -1, 1);
    
    RooGaussian sig("sig", "gaussian Fit", Mx, mu_sig, sigma_sig);
    RooGaussian background("background", "Background", Mx, mu_bkg,sigma_bkg);
    RooAddPdf model_ext("model_ext", "Signal + Background", RooArgList(sig, background), RooArgList(N_sig, N_bkg));
    
    

    //Loop over bins
    std::vector<double> results;
    std::vector<double> phibn_avgVals;
    std::vector<double> depolarizations;
    std::vector<double> depolarizations_err;
    TCanvas* SinCanvas = new TCanvas("SinCanvas","Sin fitting",800,800);
    int grid_n = (int) std::round(sqrt(bn_edges.size()-1));
    SinCanvas->Divide(grid_n,grid_n);
    TCanvas* Nasym_c = new TCanvas("Nasym_c","Nasym_c",800,800);
    Nasym_c->Divide(grid_n,grid_n);

    TFile* infile = TFile::Open(Form("%s/Mxtempcut_tree.root",outputDir.c_str()),"READ");
    TTree* chain = (TTree*)infile->Get(treeName);
    
    for (size_t i = 0; i < bn_edges.size() - 1; ++i) {
        std::cout<<"--------------------------------------------------------------------------------"<<std::endl;
        std::vector<double> Aphi;
        std::vector<double> Aphi_errs;
        std::vector<double> N_pos;
        std::vector<double> N_neg;
        phibn_avgVals.clear();
        double obsmin = bn_edges[i];
        double obsmax = bn_edges[i + 1];
        for (size_t j = 0; j < phibn_edges.size() - 1; ++j) {
            std::cout<<"\033[0;32mFitting "<<obs_str<<"bin "<<i<<" phi bin "<<j<<"\033[0m\n"<<std::endl;
            
            double phimin = phibn_edges[j];
            double phimax = phibn_edges[j + 1];
            
            TString cut = Form("%s > %f && %s < %f && phi > %f && phi < %f && hel == -1", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            RooDataSet neg_DS("neg_DS", "neg_DS", chain, RooArgSet(Mx, obs, phi, hel), cut);
            
            cut = Form("%s > %f && %s < %f && phi > %f && phi < %f && hel == 1", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            RooDataSet pos_DS("pos_DS", "pos_DS", chain, RooArgSet(Mx, obs, phi, hel), cut);
            
            
            //calculate phi location of the average number of events
            TH1F* hphi = new TH1F(Form("h_%d_%d",i,j),"histogram of phi in this bin",200,-3.14,3.14);
            cut = Form("%s > %f && %s < %f && phi > %f && phi < %f", obs_str,obsmin, obs_str,obsmax, phimin, phimax);
            chain->Draw(Form("phi>>h_%d_%d",i,j),cut,"goff");
            double mean = hphi->GetMean();
            phibn_avgVals.push_back(mean);

            //fit Mx to separate sig and bkg and calculate N_sig for pos and neg helicities
            model_ext.fitTo(neg_DS, Range("fullRange"), Save(), PrintLevel(-1), Verbose(false),Warnings(false),Extended(kTRUE));
            double N_sig_neg = N_sig.getVal();
            double N_sig_neg_err = N_sig.getError();
            RooPlot* frame_neg = Mx.frame(0.6, 2.3);
            PlotSigFit(model_ext, neg_DS,frame_neg,"neg",i,j,outputDir,obs_str,obsmin,obsmax,phimin,phimax,"Mx",obs2_str,obs2min,obs2max);
            delete frame_neg;
            
            //long N_bkg_neg = N_bkg.getVal();
            model_ext.fitTo(pos_DS, Range("fullRange"), Save(), PrintLevel(-1), Verbose(false),Warnings(false),Extended(kTRUE));
            double N_sig_pos = N_sig.getVal();
            double N_sig_pos_err = N_sig.getError();
            RooPlot* frame_pos = Mx.frame(0.6, 2.3);
            PlotSigFit(model_ext, pos_DS,frame_pos,"pos",i,j,outputDir,obs_str,obsmin,obsmax,phimin,phimax,"Mx",obs2_str,obs2min,obs2max);
            delete frame_pos;

            //beam polarization taken to be the avg of polarizations in the sample
            TH1F hist1("hist1","hist1",200,0,1);
            chain->Draw("Pol>>hist1",Form("%f<%s && %s<%f",obsmin,obs_str,obs_str,obsmax),"goff");
            double Pol_avg = hist1.GetMean();
            double Pol_stdv = hist1.GetStdDev();
            double Pol_N_hist = hist1.GetEntries();
            double Pol_avg_err = Pol_stdv / sqrt(Pol_N_hist);
            
            //calculate A for this phi bin 
            double denom = N_sig_pos + N_sig_neg;
            double A = (1/Pol_avg)*(N_sig_pos - N_sig_neg) / denom;
            double dA_dNp = (1/Pol_avg)*2.0*N_sig_neg/(denom*denom);
            double dA_dNm = -(1/Pol_avg)*2.0*N_sig_pos/(denom*denom);
            double dA_dpol = -(1/Pol_avg)*A;
            double A_err = sqrt(pow(dA_dNp * N_sig_pos_err,2) + pow(dA_dNm * N_sig_neg_err,2)+pow(dA_dpol * Pol_avg_err,2));
            Aphi.push_back(A);
            Aphi_errs.push_back(A_err);
            N_pos.push_back(N_sig_pos);
            N_neg.push_back(N_sig_neg);
        
        }
        auto res = FitToSin(phibn_avgVals,Aphi,Aphi_errs,i,outputDir,obs_str,obsmin,obsmax,SinCanvas);
        double alpha = res.first;
        double alpha_err = res.second;
        
        //define depolarization and error
        TH1F hist("hist","hist",200,0,1);
        chain->Draw("eps>>hist",Form("%f<%s && %s<%f",obsmin,obs_str,obs_str,obsmax),"goff");
        double eps_avg = hist.GetMean();
        double eps_stdv = hist.GetStdDev();
        double N_hist = hist.GetEntries();
        double eps_avg_err = eps_stdv / sqrt(N_hist);

        double avg_depolarization = sqrt(2*eps_avg*(1-eps_avg));
        double avg_depolarization_err = eps_avg_err * (abs(1-2*eps_avg)/sqrt(2*eps_avg*(1-eps_avg)));
        
        double FLU_UU = alpha / avg_depolarization;
        double FLU_UU_err = FLU_UU*sqrt(pow(alpha_err/alpha,2)+pow(avg_depolarization_err/avg_depolarization,2));
        results.push_back(FLU_UU);
        results.push_back(FLU_UU_err);
        

        Nasym_c->cd(i+1);
        TH1D* h_pos = new TH1D("h_pos","h_pos;#phi_{h};Signal Counts",phibn_edges.size()-1,phibn_edges.data());
        TH1D* h_neg = new TH1D("h_neg","h_neg;#phi_{h};Signal Counts",phibn_edges.size()-1,phibn_edges.data());
        for (int k = 0; k<phibn_edges.size()-1; k++){
            h_pos->SetBinContent(k+1,N_pos[k]);
            h_neg->SetBinContent(k+1,N_neg[k]);
        }
        h_pos->SetLineStyle(2);  
        h_pos->SetLineColor(kRed);
        h_pos->SetLineWidth(1.2); 
        h_neg->SetLineStyle(1);  
        h_neg->SetLineColor(kBlack);
        h_neg->SetLineWidth(1.2); 
        
        TLegend* leg = new TLegend(0.60,0.75,0.9,0.9);
        leg->AddEntry(h_pos,"Pos Helicity");
        leg->AddEntry(h_neg,"Neg Helicity");

        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextColor(kBlack);
        

        h_neg->SetStats(0);
        h_neg->SetTitle("");
        h_neg->Draw();
        h_pos->Draw("SAME");
        latex->DrawLatex(0.35, 0.85, Form("%.2f<%s<%.2f", obsmin,obs_str,obsmax));
        leg->Draw();
        
        
        
    }
    std::cout<<"\033[0;32mSuccessfully completed phiBinning\033[0m\n"<<std::endl;
    Nasym_c->SaveAs(Form("%s/Mx%sNasym_grid_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));
    SinCanvas->SaveAs(Form("%s/MxPhiBinningSinFits_%s_%s%.2f-%.2f.png",outputDir.c_str(),obs_str,obs2_str,obs2min,obs2max));
    delete SinCanvas;
    delete Nasym_c;
    return results;
    
}



