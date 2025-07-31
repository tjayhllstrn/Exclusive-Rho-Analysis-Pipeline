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

def main(inputFiles,outputDir,filename,
         fit_type,obs,bin_edges,
        graph_title,plot_type,save_purity_plot,
        treeName,exclusiveFit):
    
    ROOT.gROOT.SetBatch(True)
    ROOT.gSystem.Load("libRooFit")
    ROOT.gInterpreter.ProcessLine("RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);")
    ROOT.gInterpreter.ProcessLine("RooMsgService::instance().setStreamStatus(0, false);")
    
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles] 

    if not os.path.isdir(outputDir):
        print(f"The directory '{outputDir}' does not exist.")
        sys.exit(1)

    if exclusiveFit == "Mh":
        sigbkg_ub = 0.9
        sigbkg_lb = 0.65
        bkg_ub = 1.7
        bkg_lb = 1.1
    elif exclusiveFit == "Mx": #when fitting with Mx as main independed var
        sigbkg_ub = 1.15
        sigbkg_lb = 0.85
        bkg_ub = 2.75
        bkg_lb = 1.4

    #--------------------------------------------------------------------
    #run the fit depending on fit_type 
    A_sig_results = []
    if fit_type == "MLM" or fit_type == "MLMchi2":
        print("\033[92mfitting sig + bkg region\033[0m")
        A_sigbkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                                     inputFiles,treeName,
                                                     sigbkg_lb,sigbkg_ub,
                                                     exclusiveFit,"sim") #runs "sim" or "single"
        print("\033[92mfitting bkg region\033[0m")
        A_bkg_results,bin_centers = RunFittingMLM(obs,bin_edges,
                                                  inputFiles,treeName,
                                                  bkg_lb,bkg_ub,
                                                  exclusiveFit,"sim") #runs "sim" or "single"

        #--------------------------------------------------------------------
        #calculate purity based on fitting Mh with double crystal ball + poly4 and integrating:
        print("\033[92mCalculating Exclusive Event Purity\033[0m")
        u_results = []
        if exclusiveFit == "Mh":
            for idx in range(len(bin_edges)-1):
                obsmin = bin_edges[idx]
                obsmax = bin_edges[idx+1]
                u_results.append(purityCalc(inputFiles,outputDir,filename,sigbkg_lb,sigbkg_ub,obs,obsmin,obsmax,treeName,save_purity_plot))
        elif exclusiveFit == "Mx":
            for idx in range(len(bin_edges)-1):
                obsmin = bin_edges[idx]
                obsmax = bin_edges[idx+1]
                u_results.append(purityCalc_mxFit(inputFiles,outputDir,filename,sigbkg_lb,sigbkg_ub,obs,obsmin,obsmax,treeName,save_purity_plot))
    
        #plot purity by bin:
        u_vals = [entry[0] for entry in u_results]
        u_errs = [entry[1] for entry in u_results]
        #shift points to avoid overlapping error bars:
        jittered_bin_centers = [entry*(1+0.005*idx) for entry in bin_centers]
        x = np.array(jittered_bin_centers, dtype='float64')
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
        
        PurityOutputFile = os.path.join(outputDir,f"PurityFits/{filename}_purityPerBin.png")
        c = ROOT.TCanvas()
        gr.Draw("AP")
        leg.Draw("SAME")
        c.SaveAs(PurityOutputFile)
        #--------------------------------------------------------------------
        #find the signal Assymetry based on A_sigbkg = u*A_sig + (1-u)*A_bkg
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
    
    phibn_edges = np.linspace(-3.14,3.14,7) #[-3.14,-2.38,-1.57, 0.0,1.57,2.38, 3.14]
 
    if fit_type == "chi2":
        plot_type = "sig_only"
        ROOT.gSystem.Load("/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/macros/libPhiBinnedFitRunner.so")
        runner = ROOT.PhiBinnedFitRunner()
        if exclusiveFit =="Mh":
            alpha = runner.Run(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName)
        elif exclusiveFit=="Mx":
            alpha = runner.Run_mxFit(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName)
            
        results = list(alpha) 
        for i in range(0,len(results),2):
            A_sig_results.append((results[i],results[i+1]))
        bin_centers = []
        for i in range(len(bin_edges)-1):
            bin_centers.append(abs(bin_edges[i]-bin_edges[i+1])/2+bin_edges[i])

    if fit_type == "MLMchi2":
        ROOT.gSystem.Load("/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/macros/libPhiBinnedFitRunner.so")
        runner = ROOT.PhiBinnedFitRunner()
        if exclusiveFit =="Mh":
            alpha = runner.Run(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName)
        elif exclusiveFit=="Mx":
            alpha = runner.Run_mxFit(bin_edges,phibn_edges,inputFiles,outputDir,obs,treeName)
        results = list(alpha)
        A_sig_results_chi2 = []
        for i in range(0,len(results),2):
            A_sig_results_chi2.append((results[i],results[i+1]))

    #--------------------------------------------------------------------
    #Plot the results
    outputFile = os.path.join(outputDir,f"{filename}.png")

    if plot_type == "sideband_comp":
        results = [A_sig_results,A_sigbkg_results,A_bkg_results]
        legend_names = [f"A_{{sig}} {fit_type} Fit",f"A_{{sig+bkg}} {fit_type} Fit",f"A_{{bkg}} {fit_type} Fit"]

    if plot_type == "sig_only":
        results = [A_sig_results]
        legend_names = [f"A_{{sig}} {fit_type} Fit"]

    if plot_type == "MLMchi2_comp":
        results = [A_sig_results_chi2,A_sig_results]
        legend_names = [f"A_{{sig}} MLM Fit",f"A_{{sig}} chi2 Fit"]
    if plot_type == "MLMchi2_sideband_comp":
        results = [A_sig_results_chi2,A_sig_results,A_sigbkg_results,A_bkg_results]
        legend_names = [f"A_{{sig}} chi2 Fit",f"A_{{sig}} MLM Fit",f"A_{{sig+bkg}} {fit_type} Fit",f"A_{{bkg}} {fit_type} Fit"]
        
    
    PlotGraphs(results,
                legend_names, graph_title,
                bin_centers,outputFile,obs)
    return 0
                   


if __name__ == "__main__":
    input_params = inputs()
    main(*input_params)
