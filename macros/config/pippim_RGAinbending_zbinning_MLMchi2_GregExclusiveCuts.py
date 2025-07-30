#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/Greg_exclusive_pippim/Fall2018Spring2019_RGA_inbending_merged_exclusive_cuts.root"]
    outputDir = "./out/BSA_Plots/GregsCuts_merged_exclusive"
    filename = "MLMchi2_zbinning"
    fit_type = "MLMchi2"    #this can be "MLM" "chi2" or "MLMchi2
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "RGA #pi^{+}#pi^{-} inbending Merged"
    plot_type = "MLMchi2_comp"  #this can be "sideband_comp" "sig_only" or "MLMchi2_comp"
    save_purity_plot = True
    treeName = "dihadron_exclusive_cuts"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName