#!/usr/bin/python3
#This macro inputs a group of root files with TTrees, an output directory, a variable to bin in, and the bin edges. It outputs Asymmetry fitting values and errors based on the single maximum liklihood fit to the probability fxn deroved from the cross section for exclusive rho production.
#macros/AsymmetryFitting.py RGAinbending_thetabinning_MLM
import sys
import os
from importlib import import_module
sys.path.append('./src') 
sys.path.append('macros/config') 
from AsymmetryUtils import *
module_name = sys.argv[1]
inputs = import_module(module_name).inputs
import numpy as np
import ROOT

def main(inputFiles,outputDir,
         fit_type,obs,bin_edges,
        graph_title,
        treeName,
        obs2,obs2bin):
    
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libRooFit")
    ROOT.gInterpreter.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);")
    ROOT.gInterpreter.ProcessLine("RooMsgService::instance().setStreamStatus(0, false);")
    
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles] 

    if not os.path.isdir(outputDir):
        print(f"The directory '{outputDir}' does not exist.")
        sys.exit(1)

    phibn_edges = np.linspace(-3.14,3.14,7) #[-3.14,-2.38,-1.57, 0.0,1.57,2.38, 3.14]
    
    results = []
    legend_names = []
    for bn_idx in range(len(obs2bin)-1): 
        if fit_type[0]=="1":
            results.append(RunMhFitMLM(obs,bin_edges,inputFiles,outputDir,treeName,obs2,obs2bin[bn_idx],obs2bin[bn_idx+1]))
            legend_names.append(f"A_{{sig}} MLM Fit to M_{{h}} {obs2} ({obs2bin[bn_idx]},{obs2bin[bn_idx+1]})")
        if fit_type[1]=="1":
            results.append(RunMhFitChi2(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2bin[bn_idx],obs2bin[bn_idx+1]))
            legend_names.append(f"A_{{sig}} Chi2 Fit to M_{{h}} {obs2} ({obs2bin[bn_idx]},{obs2bin[bn_idx+1]})")
        if fit_type[2]=="1":
            results.append(RunMxFitMLM(obs,bin_edges,inputFiles,outputDir,treeName,obs2,obs2bin[bn_idx],obs2bin[bn_idx+1]))
            legend_names.append(f"A_{{sig}} MLM Fit to M_{{x}} {obs2} ({obs2bin[bn_idx]},{obs2bin[bn_idx+1]})")
        if fit_type[3]=="1":
            results.append(RunMxFitChi2(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2bin[bn_idx],obs2bin[bn_idx+1]))
            legend_names.append(f"A_{{sig}} Chi2 Fit to M_{{x}} {obs2} ({obs2bin[bn_idx]},{obs2bin[bn_idx+1]})")

        
    #--------------------------------------------------------------------
    #Plot the results
    bin_centers = []
    for i in range(len(bin_edges)-1):
        absval = abs(bin_edges[i]-bin_edges[i+1])
        bin_centers.append(bin_edges[i]+absval/2)
    outputFile = os.path.join(outputDir,f"RESULTS_{obs}binning_{fit_type}")
    PlotGraphs(results,
                legend_names, graph_title,
                bin_centers,outputFile,obs)
    return 0

def RunMxFitMLM(obs,bin_edges,inputFiles,outputDir,treeName,obs2,obs2min,obs2max):
    sigbkg_ub = 1.15
    sigbkg_lb = 0.85
    bkg_ub = 2.75
    bkg_lb = 1.4

    print("\033[92mfitting sig + bkg region\033[0m")
    A_sigbkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                                 inputFiles,treeName,
                                                 sigbkg_lb,sigbkg_ub,
                                                 obs2,obs2min,obs2max,
                                                 "Mx","sim") #runs "sim" or "single"
    print("\033[92mfitting bkg region\033[0m")
    A_bkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                             inputFiles,treeName,
                                             bkg_lb,bkg_ub,
                                             obs2,obs2min,obs2max,
                                             "Mx","sim") #runs "sim" or "single"
    
    #--------------------------------------------------------------------
    #calculate purity based on fitting Mh with double crystal ball + poly4 and integrating:
    print("\033[92mCalculating Exclusive Event Purity\033[0m")
    u_results = []
    for idx in range(len(bin_edges)-1):
        obsmin = bin_edges[idx]
        obsmax = bin_edges[idx+1]
        u_results.append(purityCalc_mxFit(inputFiles,outputDir,"MxMLM",sigbkg_lb,sigbkg_ub,obs,obsmin,obsmax,obs2,obs2min,obs2max,treeName,True))
    
    #plot purity by bin:
    u_vals = [entry[0] for entry in u_results]
    u_errs = [entry[1] for entry in u_results]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(u_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(u_errs,dtype='float64')
        
    gr = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(ROOT.kBlue)
    gr.SetLineColor(ROOT.kBlue)
    gr.SetLineWidth(2)
    gr.SetTitle("")
    if obs =="cth":
        gr.GetXaxis().SetTitle(f"cos#theta")
    else:
        gr.GetXaxis().SetTitle(f"{obs}")

    leg = ROOT.TPaveText(0.15, 0.75, 0.55, 0.85, "NDC") #ROOT.TLegend(0.15,0.75,0.55,0.85) #xmin,ymin,xmax.ymax
    #leg.SetBorderSize(1)
    leg.AddText("Exclusive Signal Purity")#leg.AddEntry(gr,"Exclusive Signal Purity","p")
    leg.SetTextSize(0.04)

    PurityOutputFile = os.path.join(outputDir,f"PurityFits/Mx_purityPerBin_{obs2}{obs2min}-{obs2max}")
    outroot = ROOT.TFile(f"{PurityOutputFile}.root","RECREATE")
    c = ROOT.TCanvas()
    gr.Draw("AP")
    leg.Draw("SAME")
    gr.Write()
    leg.Write()
    outroot.Close()
    c.SaveAs(f"{PurityOutputFile}.png")

    #--------------------------------------------------------------------
    #find the signal Assymetry based on A_sigbkg = u*A_sig + (1-u)*A_bkg
    A_sig_results = []
    for i in range(len(A_sigbkg_results)):
        A_sigbkg_val = A_sigbkg_results[i][0]
        A_sigbkg_err = A_sigbkg_results[i][1]
        A_bkg_val = A_bkg_results[i][0]
        A_bkg_err = A_bkg_results[i][1]
        u_val = u_results[i][0]
        u_err = u_results[i][1]
            
        A_sig_val = 1/u_val*(A_sigbkg_val-(1-u_val)*A_bkg_val)
        #error propogation based on general formula
        a = 1/u_val*A_sigbkg_err
        b = (1-u_val)/u_val * A_bkg_err
        c = ROOT.TMath.Power(u_val,-2)*u_err*(A_bkg_val-A_sigbkg_val)
        A_sig_err = ROOT.TMath.Sqrt(ROOT.TMath.Power(a,2)+ROOT.TMath.Power(b,2)+ROOT.TMath.Power(c,2))
    
        A_sig_results.append((A_sig_val,A_sig_err))

    return A_sig_results

def RunMxFitChi2(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2min,obs2max):
    ROOT.gSystem.Load("/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/macros/libPhiBinnedFitRunner.so")
    runner = ROOT.PhiBinnedFitRunner()

    alpha = runner.Run_mxFit(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2min,obs2max)
            
    results = list(alpha)
    A_sig_results = []
    for i in range(0,len(results),2):
        A_sig_results.append((results[i],results[i+1]))

    return A_sig_results
        


def RunMhFitMLM(obs,bin_edges,inputFiles,outputDir,treeName,obs2,obs2min,obs2max):
    sigbkg_ub = 0.9
    sigbkg_lb = 0.65
    bkg_ub = 1.7
    bkg_lb = 1.1

    print("\033[92mfitting sig + bkg region\033[0m")
    A_sigbkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                                 inputFiles,treeName,
                                                 sigbkg_lb,sigbkg_ub,
                                                 obs2,obs2min,obs2max,
                                                 "Mh","sim") #runs "sim" or "single"
    print("\033[92mfitting bkg region\033[0m")
    A_bkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                             inputFiles,treeName,
                                             bkg_lb,bkg_ub,
                                              obs2,obs2min,obs2max,
                                             "Mh","sim") #runs "sim" or "single"

    #--------------------------------------------------------------------
    #calculate purity based on fitting Mh with double crystal ball + poly4 and integrating:
    print("\033[92mCalculating Exclusive Event Purity\033[0m")
    u_results = []
    for idx in range(len(bin_edges)-1):
        obsmin = bin_edges[idx]
        obsmax = bin_edges[idx+1]
        u_results.append(purityCalc(inputFiles,outputDir,"MhMLM",sigbkg_lb,sigbkg_ub,obs,obsmin,obsmax,obs2,obs2min,obs2max,treeName,True))
    
    #plot purity by bin:
    u_vals = [entry[0] for entry in u_results]
    u_errs = [entry[1] for entry in u_results]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(u_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(u_errs,dtype='float64')
        
    gr = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(ROOT.kBlue)
    gr.SetLineColor(ROOT.kBlue)
    gr.SetLineWidth(2)
    gr.SetTitle("")
    if obs =="cth":
        gr.GetXaxis().SetTitle(f"cos#theta")
    else:
        gr.GetXaxis().SetTitle(f"{obs}")

    leg = ROOT.TPaveText(0.15, 0.75, 0.55, 0.85, "NDC") #ROOT.TLegend(0.15,0.75,0.55,0.85) #xmin,ymin,xmax.ymax
    #leg.SetBorderSize(1)
    leg.AddText("Exclusive Signal Purity")#leg.AddEntry(gr,"Exclusive Signal Purity","p")
    leg.SetTextSize(0.04)
    
    
    PurityOutputFile = os.path.join(outputDir,f"PurityFits/Mh_purityPerBin_{obs2}{obs2min}-{obs2max}")
    outroot = ROOT.TFile(f"{PurityOutputFile}.root","RECREATE")
    c = ROOT.TCanvas()
    gr.Draw("AP")
    leg.Draw("SAME")
    gr.Write()
    leg.Write()
    outroot.Close()
    c.SaveAs(f"{PurityOutputFile}.png")

    #--------------------------------------------------------------------
    #find the signal Assymetry based on A_sigbkg = u*A_sig + (1-u)*A_bkg
    A_sig_results = []
    for i in range(len(A_sigbkg_results)):
        A_sigbkg_val = A_sigbkg_results[i][0]
        A_sigbkg_err = A_sigbkg_results[i][1]
        A_bkg_val = A_bkg_results[i][0]
        A_bkg_err = A_bkg_results[i][1]
        u_val = u_results[i][0]
        u_err = u_results[i][1]
            
        A_sig_val = 1/u_val*(A_sigbkg_val-(1-u_val)*A_bkg_val)
        #error propogation based on general formula
        a = 1/u_val*A_sigbkg_err
        b = (1-u_val)/u_val * A_bkg_err
        c = ROOT.TMath.Power(u_val,-2)*u_err*(A_bkg_val-A_sigbkg_val)
        A_sig_err = ROOT.TMath.Sqrt(ROOT.TMath.Power(a,2)+ROOT.TMath.Power(b,2)+ROOT.TMath.Power(c,2))
    
        A_sig_results.append((A_sig_val,A_sig_err))

    return A_sig_results


def RunMhFitChi2(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2min,obs2max):
    ROOT.gSystem.Load("/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/macros/libPhiBinnedFitRunner.so")
    runner = ROOT.PhiBinnedFitRunner()
    alpha = runner.Run(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName,obs2,obs2min,obs2max)            
    results = list(alpha) 
    A_sig_results = []
    for i in range(0,len(results),2):
        A_sig_results.append((results[i],results[i+1]))

    return A_sig_results

                   


if __name__ == "__main__":
    input_params = inputs()
    main(*input_params)
