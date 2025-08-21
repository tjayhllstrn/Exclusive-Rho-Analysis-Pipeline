#!/usr/bin/python3


#configuration file for BSA Fitting using chi2 on RGA inbending data with z binning scheme used in Greg's Analysis Note

def inputs():
    inputFiles = ["out/pippi0_fall2018_in_pass2/nSidis_005032.root"]#,"out/pippi0_spring2019_in_pass2/pippi0_spring2019_in_pass2.root"]
    outputDir = "./out/test"
    fit_type = "0001"    #indeces turn on or off the different fits (MhMLM,MhChi2,MxMLM,MxChi2)
    obs = "z"
    bin_edges = [0.66,0.72]
    graph_title = "RGA #pi^{+}#pi^{0} inbending"
    treeName = "pippi0"
    obs2 = "cth"     #obs2 is the one that bins the obs dependence
    obs2bin = [-1,1]#[-1,-0.17,0.12,1]
    return inputFiles,outputDir,fit_type,obs,bin_edges,graph_title,treeName,obs2,obs2bin
