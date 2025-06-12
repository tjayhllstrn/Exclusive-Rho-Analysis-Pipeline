#!/usr/bin/python3
#This macro inputs a group of root files with TTrees, an output directory, a variable to bin in, and the bin edges. It outputs Asymmetry fitting values and errors based on the single maximum liklihood fit to the probability fxn deroved from the cross section for exclusive rho production.
import sys
import os
sys.path.append('./src') 
from AsymmetryUtils import *
import numpy as np
import ast

#This function defines custom options that you might want to change per fitting run and then runs them using functions defined in AsymmetyryUtils. Fit Type options are "MLM" for simultaneous MLM itting method, and a "chi2" type that does a chi2 fit in bins of phi rather than uses MLM.
#macros/AsymmetryFitting.py '["out/pippi0_fall2018_in_pass1/pippi0_fall2018_in_pass1.root","out/pippi0_spring2019_in_pass1/pippi0_spring2019_in_pass1.root"]' "./out/BSA_Plots/" "RGA_in_RegionComparison" "MLM"

def main(inputFiles,outputDir,filename,
         fit_type = "chi2",
         obs = "z",bin_edges=np.linspace(0.5,1,11)):
         
    if not isinstance(inputFiles, list):
        inputFiles = [inputFiles] 

    if not os.path.isdir(outputDir):
        print(f"The directory '{outputDir}' does not exist.")
        sys.exit(1)

    sigbkg_ub = 0.9
    sigbkg_lb = 0.65
    bkg_ub = 1.7
    bkg_lb = 1.1

    #run the fit depending on fit_type
    if fit_type == "MLM":
        print("\033[92mfitting sig + bkg region\033[0m")
        A_sigbkg_results,bin_centers = RunFittingMLM(obs,bin_edges,inputFiles,sigbkg_lb,sigbkg_ub,"sim") #runs "sim" or "single"
        print("\033[92mfitting bkg region\033[0m")
        A_bkg_results,bin_centers = RunFittingMLM(obs,bin_edges,inputFiles,bkg_lb,bkg_ub,"sim") #runs "sim" or "single"
        
    if fit_type == "chi2":
        A_results,bin_centers = RunFittingChi2(obs,bin_edges,inputFiles)

    #calculate purity based on fitting Mh with double crystal ball + poly4 and integrating:
    print("\033[92mCalculating Exclusive Event Purity\033[0m")
    u = purityCalc(inputFiles,outputDir,filename,sigbkg_lb,sigbkg_ub, "Q")

    A_sig_results = []
    for i in range(len(A_sigbkg_results)):
        A_sigbkg_val = A_sigbkg_results[i][0]
        A_sigbkg_err = A_sigbkg_results[i][1]
        A_bkg_val = A_bkg_results[i][0]
        A_bkg_err = A_bkg_results[i][1]
        
        A_sig_val = 1/u[0]*(A_sigbkg_val-(1-u[0])*A_bkg_val)
        A_sig_err = 0 #Change later ---------------------------------???????????????????????????????????????????

        A_sig_results.append((A_sig_val,A_sig_err))
    
    #Plot the results
    graph_title = "RGA Inbending A"
    outputFile = os.path.join(outputDir,f"{filename}.png")
    #RunPlotting(A_results,bin_centers,
    #            leg_title,graph_title,outputFile,
    #            obs)

    ROOT.gROOT.SetBatch(True)
    c = ROOT.TCanvas()
    c.SetTickx()
    c.SetTicky()
    c.SetGridx()
    c.SetGridy()

    #prepare arrays for graphing sig+bkg
    A_sigbkg_vals = [entry[0] for entry in A_sigbkg_results]
    A_sigbkg_errs = [entry[1] for entry in A_sigbkg_results]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(A_sigbkg_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(A_sigbkg_errs,dtype='float64')
    
    gr_sigbkg = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr_sigbkg.SetMarkerColor(ROOT.kGreen)
    gr_sigbkg.SetMarkerStyle(20)
    gr_sigbkg.SetTitle(graph_title)
    gr_sigbkg.GetXaxis().SetTitle(obs)
    gr_sigbkg.GetYaxis().SetTitle("F_{LU}^{sin#phi}/F_{UU}")
    

    #prepare arrays for graphing bkg
    A_bkg_vals = [entry[0] for entry in A_bkg_results]
    A_bkg_errs = [entry[1] for entry in A_bkg_results]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(A_bkg_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(A_bkg_errs,dtype='float64')

    gr_bkg = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr_bkg.SetMarkerColor(ROOT.kRed)
    gr_bkg.SetMarkerStyle(20)


    #prepare arrays for graphing sig
    A_sig_vals = [entry[0] for entry in A_sig_results]
    A_sig_errs = [entry[1] for entry in A_sig_results]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(A_sig_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(A_sig_errs,dtype='float64')

    gr_sig = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr_sig.SetMarkerColor(ROOT.kBlue)
    gr_sig.SetMarkerStyle(20)

    leg = ROOT.TLegend(0.55,0.15,0.85,0.35) #xmin,ymin,xmax.ymax
    leg.AddEntry(gr_sigbkg,"A_{sig+bkg} MLM Fit","p") #name of object being referred, "title shown", (fit line or point)
    leg.AddEntry(gr_bkg,"A_{bkg} MLM Fit","p") #name of object being referred, "title shown", (fit line or point)
    leg.AddEntry(gr_sig,"A_{sig} MLM Fit","p") #name of object being referred, "title shown", (fit line or point)
    leg.SetBorderSize(0)
    
    gr_sigbkg.Draw("AP")
    gr_bkg.Draw("P SAME")
    gr_sig.Draw("P SAME")
    leg.Draw()
    c.SaveAs(outputFile)
    
if __name__ == "__main__":
    args = [arg for arg in sys.argv[2:]]
    inputFiles=ast.literal_eval(sys.argv[1])
    main(inputFiles,*args)
