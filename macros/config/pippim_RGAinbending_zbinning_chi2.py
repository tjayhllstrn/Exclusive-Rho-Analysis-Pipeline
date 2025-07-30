#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippim_spring2019_in_pass2/pippim_spring2019_in_pass2.root","out/pippim_fall2018_in_pass2/pippim_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippim_RGA_inbending_pass2_merged"
    filename = "chi2_zbinning"
    fit_type = "chi2"    #this can be "MLM" or "Chi2"
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "RGA #pi^{+}#pi^{-} inbending Merged"
    plot_type = "sig_only"  #this can be "sideband_comp" "sig_only", "MLMchi2_comp", or "MLMchi2_sideband_comp"
    save_purity_plot = False
    treeName = "pippim"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName