#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root","out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippi0_RGA_inbending_pass2_merged"
    fit_type = "1111"    #indeces turn on or off the different fits (MhMLM,MhChi2,MxMLM,MxChi2)
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "RGA #pi^{+}#pi^{0} inbending Merged"
    treeName = "pippi0"
    return inputFiles,outputDir,fit_type,obs,bin_edges,graph_title,treeName
