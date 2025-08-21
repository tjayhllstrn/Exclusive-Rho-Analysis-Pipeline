#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root","out/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root"]
    outputDir = "./out/BSA_Plots/pippi0_RGA_inbending_pass2_merged/tBinning_cthDep"
    fit_type = "0001"    #indeces turn on or off the different fits (MhMLM,MhChi2,MxMLM,MxChi2)
    obs = "t_elec"
    bin_edges = [-3.21, -2.53, -1.99, -1.64, -1.35, -1.08, -0.83, -0.62, -0.4, -0.17]
    graph_title = "RGA #pi^{+}#pi^{0} inbending Merged"
    treeName = "pippi0"
    obs2 = "cth"     #obs2 is the one that bins the obs dependence
    obs2bin = [-1,1]#[-1,-0.17,0.12,1]
    return inputFiles,outputDir,fit_type,obs,bin_edges,graph_title,treeName,obs2,obs2bin