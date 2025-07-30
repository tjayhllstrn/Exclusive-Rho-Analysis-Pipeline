#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with cth binning scheme 

def inputs():
    inputFiles = ["out/pippim_spring2019_in_pass2/pippim_spring2019_in_pass2.root","out/pippim_fall2018_in_pass2/pippim_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippim_RGA_inbending_pass2_merged"
    filename = "MLMchi2_cthbinning"
    fit_type = "MLMchi2"    #this can be "MLM" "chi2" or "MLMchi2
    obs = "cth"
    bin_edges = [-0.8,-0.53,-0.3,-0.15,0.0,0.15,0.3,0.53,0.8]
    graph_title = "RGA #pi^{+}#pi^{-} inbending Merged"
    plot_type = "MLMchi2_comp"  #this can be "sideband_comp" "sig_only" or "MLMchi2_comp"
    save_purity_plot = True
    treeName = "pippim"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName