#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root","out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippi0_RGA_inbending_pass2_merged"
    filename = "MLMchi2_zbinning"
    fit_type = "MLMchi2"    #this can be "MLM" "chi2" or "MLMchi2
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "RGA #pi^{+}#pi^{0} inbending Merged"
    plot_type = "MLMchi2_comp"  #this can be "sideband_comp" "sig_only", "MLMchi2_comp", or "MLMchi2_sideband_comp"
    save_purity_plot = True
    treeName = "pippi0"
    return inputFiles,outputDir,filename,fit_type,obs,bin_edges,graph_title,plot_type,save_purity_plot,treeName