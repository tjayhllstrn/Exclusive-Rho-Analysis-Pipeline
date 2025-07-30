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
#include <vector>
#include <iostream>

void FitToSin(){
    int i = 100;
    std::vector<double> x_vals = {-3.14,-1.01,0,1.01,3.14};
    std::vector<double> y_vals = {-1,-0.5,0,0.5,1};
    TGraph* gr = new TGraph(x_vals.size(),x_vals.data(),y_vals.data());
    gr->SetTitle(Form("Sin Fit z%d",(int)i));

    TF1* fitFunc = new TF1("fitFunc","[0]*sin(x)",-3.14,3.14);
    fitFunc->SetParameter(0,0.1);

    gr->Fit(fitFunc,"R");

    TCanvas* c1 = new TCanvas();

    gr->SetMarkerStyle(20);  // Common visible marker
    gr->SetMarkerSize(1.2);  // Optional: make it more visible
    gr->Draw("AP");
    fitFunc->Draw("same");
    c1->SaveAs(Form("out/test/SinFit_z%d.png",(int)i));

    double amplitude = fitFunc->GetParameter(0);
    
    std::cout<<"fitted amplitude"<<amplitude<<std::endl;
    
}