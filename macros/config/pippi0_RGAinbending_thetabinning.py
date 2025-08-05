#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with cth binning scheme 

def inputs():
def inputs():
    inputFiles = ["out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root","out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippi0_RGA_inbending_pass2_merged"
    fit_type = "1111"    #indeces turn on or off the different fits (MhMLM,MhChi2,MxMLM,MxChi2)
    obs = "cth"
    bin_edges = [-0.8,-0.53,-0.3,-0.15,0.0,0.15,0.3,0.53,0.8]
    graph_title = "RGA #pi^{+}#pi^{0} inbending Merged"
    treeName = "pippi0"
    return inputFiles,outputDir,fit_type,obs,bin_edges,graph_title,treeName