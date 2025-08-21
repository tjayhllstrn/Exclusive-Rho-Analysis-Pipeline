#ifndef PHIBINNEDFIT_H
#define PHIBINNEDFIT_H

#include <vector>
#include <string>
#include <RooFit.h>

class PhiBinnedFitRunner {
public:
    PhiBinnedFitRunner(){}
    std::vector<double> Run(const std::vector<double>& bn_edges,const std::vector<double>& phibn_edges, const std::vector<std::string>& inputFiles, const std::string& outputDir,const char* obs_str,const char* treeName, const char* obs2_str, double obs2min,double obs2max);
    std::vector<double> Run_mxFit(const std::vector<double>& bn_edges,const std::vector<double>& phibn_edges, const std::vector<std::string>& inputFiles, const std::string& outputDir,const char* obs_str,const char* treeName, const char* obs2_str, double obs2min,double obs2max);
};

#endif
