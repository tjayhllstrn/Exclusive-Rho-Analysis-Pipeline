#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFit.h"
#include <math.h>
#include <vector>
#include <iostream>

double FitToSin(vector<double> x_vals,vector<double> y_vals,int i){
    TGraph* gr = new TGraph(x_vals.size(),x_vals.data(),y_vals.data());
    gr->SetTitle(Form("Sin Fit z%d",(int)i));

    TF1* fitFunc = new TF1("fitFunc","[0]*sin(x)",-3.14,3.14);
    fitFunc->SetParameter(0,0.1);

    gr->Fit(fitFunc,"R");

    TCanvas* c1 = new TCanvas();
    
    gr->SetMarkerStyle(20);  // Common visible marker
    gr->SetMarkerSize(1.2);
    gr->Draw("AP");
    fitFunc->Draw("same");
    c1->SaveAs(Form("out/test/SinFit_z%d.png",(int)i));

    double amplitude = fitFunc->GetParameter(0);

    return amplitude;
    
}

void BSA_phiBinning_RooFit() {
    using namespace RooFit;
    gROOT->SetBatch(kFALSE);

    // Load TChain
    TChain uncut("pippi0");
    uncut.Add("out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root");
    uncut.Add("out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root");
    
    TTree* chain = uncut.CopyTree("Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx<1.05");
    
    RooRealVar Mh("Mh", "Mh", 0.4, 1.7);
    Mh.setRange("fullRange", 0.4, 1.4);
    Mh.setRange("sigbkgRange", 0.65, 0.9);
    
    RooRealVar mu("m_{0}", "mu", 0.8, 0.6, 1);
    RooRealVar sigma("#sigma", "sigma", 0.06, 0.00001, 0.4);
    RooRealVar p1("p1", "p1", 0, -5, 5);
    RooRealVar p2("p2", "p2", 0, -5, 5);
    RooRealVar p3("p3", "p3", 0, -5, 5);
    RooRealVar p4("p4", "p4", 0, -5, 5);
    
    RooRealVar N_sig("N_{sig}", "N_sig", 10000, 0, chain->GetEntries());
    RooRealVar N_bkg("N_{bkg}", "N_bkg", 10000, 0, chain->GetEntries());
    
    RooRealVar z("z", "z", 0, 1);
    RooRealVar phi("phi", "phi", -3.14, 3.14);
    RooRealVar hel("hel", "hel", -1, 1);
    
    RooGaussian sig("sig", "gaussian Fit", Mh, mu, sigma);
    RooArgList pars_pol(p1, p2, p3, p4);
    RooChebychev background("background", "Background", Mh, pars_pol);
    RooAddPdf model_ext("model_ext", "Signal + Background", RooArgList(sig, background), RooArgList(N_sig, N_bkg));
    
    std::vector<double> zbn_edges = {0.96, 1.0}; //{0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 
    std::vector<double> phibn_edges = {-3.14, -2.57, -2,-0.75, 0.75,2,2.57, 3.14};
    std::vector<double> phibn_centers;
    TCanvas c2D;
    c2D.SetTickx();
    c2D.SetTicky();
    
    TH2F binning("binning","binning scheme;z [GeV];#phi_{h} [rad]",100,0.6,1,100,-3.2,3.2);
    
    chain->Draw("phi:z>>binning");
    
    binning.SetStats(0);
    binning.Draw();
    
    std::vector<TLine*> lines;
    for (double z : zbn_edges) {
        TLine* line = new TLine(z, -3.14, z, 3.14);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    for (double phi : phibn_edges){
        TLine* line = new TLine(0.6, phi, 1, phi);
        line->SetLineColor(kBlack);
        line->SetLineWidth(3);
        line->Draw();
        lines.push_back(line);
            }
    c2D.Update();
    c2D.Draw();

    //Loop over bins
    vector<double> alphaz;
    vector<double> zbin_centers;
    
    for (size_t i = 0; i < zbn_edges.size() - 1; ++i) {
        std::vector<double> Aphi;
        for (size_t j = 0; j < phibn_edges.size() - 1; ++j) {
            double zmin = zbn_edges[i];
            double zmax = zbn_edges[i + 1];
            double phimin = phibn_edges[j];
            double phimax = phibn_edges[j + 1];
            double phi_absval = fabs(phimax-phimin);
            phibn_centers.push_back(phi_absval/2 + phimin);
            double z_absval = fabs(zmax-zmin);
            zbn_centers.push_back(z_absval/2 + zmin);

            TString cut = Form("z > %f && z < %f && phi > %f && phi < %f && hel == -1", zmin, zmax, phimin, phimax);
            RooDataSet neg_DS("neg_DS", "neg_DS", chain, RooArgSet(Mh, z, phi, hel), cut);

            cut = Form("z > %f && z < %f && phi > %f && phi < %f && hel == 1", zmin, zmax, phimin, phimax);
            RooDataSet pos_DS("pos_DS", "pos_DS", chain, RooArgSet(Mh, z, phi, hel), cut);

            //fit Mh to separate sig and bkg and calculate N_sig for pos and neg helicities
            model_ext.fitTo(neg_DS, Range("fullRange"), Save(), PrintLevel(-1), Extended(kTRUE));
            RooAbsReal* sig_int = sig.createIntegral(RooArgSet(Mh), RooArgSet(Mh), "sigbkgRange");
            double sig_perc = sig_int->getVal();
            double N_sig_neg = sig_perc * N_sig.getVal();

            //long N_bkg_neg = N_bkg.getVal();
            model_ext.fitTo(pos_DS, Range("fullRange"), Save(), PrintLevel(-1), Extended(kTRUE));
            sig_int = sig.createIntegral(RooArgSet(Mh), RooArgSet(Mh), "sigbkgRange");
            sig_perc = sig_int->getVal();
            double N_sig_pos = sig_perc * N_sig.getVal();

            // Optional: plot
            RooPlot* frame = Mh.frame(0.4, 1.7);
            neg_DS.plotOn(frame, MarkerSize(0.5), Name("data"));
            model_ext.plotOn(frame, LineStyle(kDashed), LineColor(kBlack), Name("totalFit"));
            model_ext.plotOn(frame, Components("background"), LineColor(kRed), Name("chebychevFit"));
            model_ext.plotOn(frame, Components("sig"), LineColor(kBlue), Name("sigFit"));

            TCanvas* c = new TCanvas();
            c->Divide(2,1);
            TPad* pad1 = (TPad*)c->cd(1);
            pad1->SetPad(0.0,0.1,0.9,1.0);
            frame->Draw();
            
            TPad* pad2 = (TPad*)c->cd(2);
            pad2->SetPad(0.815,0.0,1.0,1.0);
            TPaveText* param_box = new TPaveText(0.05, 0.2, 0.95, 0.8, "NDC");
            param_box->SetFillColor(0);
            param_box->SetTextAlign(22);
            param_box->SetTextSize(0.1);

            RooArgSet* params = model_ext.getParameters(pos_DS);
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

            c->SaveAs(Form("out/test/fit_z%d_phi%d.png", (int)i, (int)j));

            //calculate A for this phi bin
            double A = (N_sig_pos - N_sig_neg) / (N_sig_pos + N_sig_neg);
            Aphi.push_back(A);
        
        }
        double alpha = FitToSin(phibn_centers,Aphi,i);
        alphaz.push_back(alpha);
        //printf("[%.3f, %.3f, %.3f,%.3f]\n", phibn_centers[0], phibn_centers[1], phibn_centers[2],phibn_centers[3]);
        //std::cout << "Fitted amplitude: " << alpha << std::endl;
        
    }

    
}

