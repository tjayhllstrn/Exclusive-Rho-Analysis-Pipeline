#!/usr/bin/python3


#configuration file for BSA Fitting using MLM on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippi0_fall2018_out_pass1/pippi0_fall2018_out_pass1.root"]
    outputDir = "./out/BSA_Plots/RGA_outbending_pass1"
    filename = "MLM_zbinning"
    fit_type = "MLM"    #this can be "MLM" or "Chi2"
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "RGA_outbending"
    plot_type = "sig_only"  #this can be ""sideband_comp"" or "sig_only"
    save_purity_plot = True
    treeName = "pippi0"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName