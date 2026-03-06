#include <iostream>
#include <string>
#include <glob.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "../src/PurityFitter.h"
#include "../src/PurityFitter.cpp"

//clas12root -l -b -q 'macros/pi0Purity.cpp("out/pippi0_fall2018_pi0NoiseThreshold","pippi0_fall2018_pi0NoiseThreshold")'

void CalculatePurity(const char* input_dir, const char* inFile_pattern,bool REWRITE_CACHE = false){
    //set variables:
    std::string CACHE_DIR = "/lustre24/expphy/volatile/clas12/users/tjhellst/cache/pi0Purities/";
    std::string OBS2 = "cth";
    std::string OBS2_min = "-1";
    std::string OBS2_max = "1";
    std::vector<double> MinPhoE_Cuts = {0.0,0.05,0.1,0.15, 0.2,0.25, 0.3, 0.35, 0.4}; // Min Pho Cuts to fit (in GeV)



    // Output results (you can modify this part as needed)
    std::cout << "input: " << input_dir << ", file pattern: " << inFile_pattern << std::endl;

    // Create a TChain to hold the tree
    TChain* chain = new TChain("pippi0");
    // Use glob to find matching files
    std::string pattern = std::string(input_dir) + "/" + std::string(inFile_pattern) + ".root";
    glob_t glob_result;
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    // Add matched files to the chain
    for(unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        std::cout << "Adding file: " << glob_result.gl_pathv[i] << std::endl;
        chain->Add(glob_result.gl_pathv[i]);
    }
    const size_t nMatchedFiles = glob_result.gl_pathc;
    globfree(&glob_result);
    std::cout<<"Input TChain has " << chain->GetEntries() << " entries from " << nMatchedFiles << " files." << std::endl;

    //pre-cut the tree 
    TCut Mx_cut = TCut("0.85<Mx&&Mx<1.05");
    TCut RhoCut = TCut("0.65<Mh&&Mh<0.9");
    TCut obs2_cut = TCut((OBS2 + ">" + OBS2_min + " && " + OBS2 + "<" + OBS2_max).c_str());
    TCut pre_cut = RhoCut && obs2_cut;
    std::string pre_cut_name = "Mx_Rho_cth"; //"Rho_cth" "MxCut" "Mx_Rho_cth" "Mx_Rho_cth_noGBT"
    // Create or open cached filtered tree
    std::cout << "  Creating/loading pre-filtered tree..." << std::endl;
    std::string safe_cache_ID = std::string(input_dir) + "_" + std::string(inFile_pattern);
    for (char& c : safe_cache_ID) {
            // Replace glob + path-unfriendly chars
            if (c == '*' || c == '?' || c == '[' || c == ']' ||
                c == '{' || c == '}' || c == ':' || c=='/' || c == ' ' ) {
                c = '-';
            }
        }
    // Generate unique filename based on input file, fit_type and obs2 bin
    std::string base_filename = safe_cache_ID + "_" + std::to_string(nMatchedFiles);
    std::string cache_filename = CACHE_DIR + "filtered_tree_" + base_filename + "_" + pre_cut_name + "_pi0Purity_" + 
                                OBS2 + "_" + OBS2_min + "_" + OBS2_max + "_" +
                                ".root";
    TFile* tempFile = nullptr;
    TTree* filteredTree = nullptr;
    
    if (!REWRITE_CACHE && !gSystem->AccessPathName(cache_filename.c_str())) {
        // File exists and we're not rewriting, open it
        std::cout << "    Using cached tree from: " << cache_filename << std::endl;
        tempFile = TFile::Open(cache_filename.c_str(), "READ");
        filteredTree = (TTree*)tempFile->Get("filtered_tree");
    } else {
        // Create new filtered tree
        std::cout << "    Creating new filtered tree at: " << cache_filename << std::endl;
        tempFile = new TFile(cache_filename.c_str(), "RECREATE");
        filteredTree = chain->CopyTree(pre_cut.GetTitle());
        filteredTree->SetName("filtered_tree");
        filteredTree->SetDirectory(tempFile); // Associate with temp file
        tempFile->Write();  // Save to cache
    }

    // Configure PurityFitter Object with TChain as the input tree
    std::string outdir = std::string(input_dir) + "/pi0_purity_results"; // You can modify this to specify a different output directory
    PurityFitter purityFitter(filteredTree, outdir);

    // loop through MinPhoE cuts and calculate purity for each
    for (double minPhoE : MinPhoE_Cuts) {
        std::cout << "Fitting purity with MinPhoE cut: " << minPhoE << " GeV" << std::endl;
        purityFitter.minPhoE = minPhoE;
        purityFitter.FitPurity();
    }
    purityFitter.pre_cut_name = pre_cut_name; // Store pre-cut name for labeling plots
    purityFitter.PlotToCanvas_PostageStamp();
    
    //Plot purity results
    purityFitter.MakeGraphLinePlot(purityFitter.Purities, MinPhoE_Cuts,
            "PurityPer Pho E Cut",
            "Min Pho E Cut (GeV)",
            ("Purity_perPhoCut_" + pre_cut_name).c_str(),
            kBlack);
}



int pi0Purity(const char* input_dir, const char* inFile_pattern,bool REWRITE_CACHE = false){
    CalculatePurity(input_dir, inFile_pattern, REWRITE_CACHE);
    return 0;
}