#!/usr/bin/python3
import numpy as np

#configuration file for BSA Fitting using MLM on inbending data with theta binning scheme 

def inputs():
    inputFiles = ["out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root","out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippi0_RGA_inbending_pass2_merged"
    filename = "chi2_cthbinning"
    fit_type = "chi2"    #this can be "MLM" or "Chi2"
    obs = "cth"
    bin_edges = [-0.8,-0.53,-0.3,-0.15,0.0,0.15,0.3,0.53,0.8]
    graph_title = "RGA #pi^{+}#pi^{0} inbending Merged"
    plot_type = "sig_only"  #this can be "sideband_comp" or "sig_only"
    save_purity_plot = True
    treeName = "pippi0"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName