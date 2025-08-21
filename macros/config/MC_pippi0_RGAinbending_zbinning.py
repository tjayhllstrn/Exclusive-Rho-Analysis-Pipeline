#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/MC_pippi0_fall2018_in_pass2/MC_pippi0_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/MC_pippi0_RGA_inbending_pass2_fall"
    fit_type = "1111"    #indeces turn on or off the different fits (MhMLM,MhChi2,MxMLM,MxChi2)
    obs = "z"
    bin_edges = [0.66, 0.74, 0.79, 0.82, 0.85, 0.88, 0.91, 0.93, 0.96, 1]
    graph_title = "MonteCarlo Fall Inbending"
    treeName = "pippi0"
    return inputFiles,outputDir,fit_type,obs,bin_edges,graph_title,treeName