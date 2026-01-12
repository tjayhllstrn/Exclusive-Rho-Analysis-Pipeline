#include "../src/MLM_Fitter.h"
#include "../src/MLM_Fitter.cpp"
#include "../src/Chi2_Fitter.h"
#include "../src/Chi2_Fitter.cpp"
#include <TEnv.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>

// helper: parse a comma-separated string into doubles
static std::vector<double> parse_csv_to_doubles(const std::string &s){
    std::vector<double> out;
    std::istringstream iss(s);
    std::string token;
    while(std::getline(iss, token, ',')){
        size_t b = token.find_first_not_of(" \t");
        if(b == std::string::npos) continue;
        size_t e = token.find_last_not_of(" \t");
        token = token.substr(b, e - b + 1);
        try{ out.push_back(std::stod(token)); } catch(...){}
    }
    return out;
}

int AsymmetryMLM(const char* file_name, 
    std::string fit_type, const char* config_file
    ){
        // load configuration
    TEnv env; 
    env.ReadFile(config_file, kEnvLocal);
    
    //Initialize and parse fitter
    std::vector<double> bn_edgs = parse_csv_to_doubles(env.GetValue("bn_edgs", ""));
    std::vector<double> obs2bn = parse_csv_to_doubles(env.GetValue("obs2bn", ""));
    std::string file_base = env.GetValue("file_base", "./out/pippi0_fall2018_in_pass2/");
    MLM_Fitter fitter(
        env.GetValue("treename", "pippi0"),
        (file_base + "BSA_Plots/" + fit_type + "/").c_str(),
        env.GetValue("obs", "z"),
        env.GetValue("obs2", "cth"),
        (file_base + file_name).c_str(),
        bn_edgs,
        obs2bn,
        fit_type
    );
    
    //Loop through obs2 bins and run MLM fitting based on the fit_type
    for(size_t bn_idx = 0; bn_idx < fitter.OBS2BN.size() - 1; bn_idx++){

        fitter.OUT_DIR = file_base + "BSA_Plots/" + fit_type + "/" + fitter.OBS2 + "_" + std::to_string(fitter.OBS2BN[bn_idx]).substr(0,5) + "_" + std::to_string(fitter.OBS2BN[bn_idx+1]).substr(0,5) + "/";
        // Check if directory exists and create if needed
        if (!std::filesystem::exists(fitter.OUT_DIR)) {
            std::filesystem::create_directories(fitter.OUT_DIR);
        }
        
        // MhMLM fit
        if(fit_type == "MhMLM"){
            fitter.RunMhFitMLM(bn_idx);
        }
        
        // MxMLM fit
        if(fit_type == "MxMLM"){
            fitter.RunMxFitMLM(bn_idx);
        }
        
        //Plot Asymmetry Results
        std::vector<TGraph*> asymmetry_graphs;
        std::pair<double,double> ybounds(-0.2, 0.5);
        
        fitter.MakeGraphLinePlot(fitter.A_sigbkg[bn_idx], fitter.BN_CENTERS,
            "F_{LU}^{sin#phi}/F_{UU} [arb. units]",
            fitter.OBS.c_str(),
            "A_sig_bkg",
            kRed-6,
            &ybounds,
            &asymmetry_graphs);

        fitter.MakeGraphLinePlot(fitter.A_bkg[bn_idx], fitter.BN_CENTERS,
            "F_{LU}^{sin#phi}/F_{UU} [arb. units]",
            fitter.OBS.c_str(),
            "A_bkg",
            kCyan-6,
            &ybounds,
            &asymmetry_graphs);

        fitter.MakeGraphLinePlot(fitter.A_sig[bn_idx], fitter.BN_CENTERS,
            "F_{LU}^{sin#phi}/F_{UU} [arb. units]",
            fitter.OBS.c_str(),
            "A_sig",
            kBlack,
            &ybounds,
            &asymmetry_graphs);

        fitter.PlotToCanvas_overlayed(asymmetry_graphs, "Asymmetries", ("Asymmetries_per_" + fitter.OBS).c_str());

        //Plot purity results
        fitter.MakeGraphLinePlot(fitter.Purities, fitter.BN_CENTERS,
            "PurityPerBin",
            fitter.OBS.c_str(),
            "Purity_perBin",
            kBlack);

        fitter.PlotToCanvas_PostageStamp();

        //clear purity results for next obs2 bin
        fitter.Purities.clear(); 
        fitter.purity_data_hists.clear();      
        fitter.purity_total_graphs.clear();  
        fitter.purity_sig_graphs.clear();    
        fitter.purity_bkg_graphs.clear();    
        fitter.purity_legends.clear();      
        fitter.purity_texts.clear();         
        fitter.purity_param_boxes.clear(); 

    }
    return 0;
    }

int AsymmetryChi2(const char* file_name, 
    std::string fit_type, const char* config_file
    ){
    // Placeholder for Chi2 fitting implementation
    TEnv env; 
    env.ReadFile(config_file, kEnvLocal);
    
    //Initialize and parse fitter
    std::vector<double> bn_edgs = parse_csv_to_doubles(env.GetValue("bn_edgs", ""));
    std::vector<double> obs2bn = parse_csv_to_doubles(env.GetValue("obs2bn", ""));
    std::string file_base = env.GetValue("file_base", "./out/pippi0_fall2018_in_pass2/");
    Chi2_Fitter fitter(
        env.GetValue("treename", "pippi0"),
        (file_base + "BSA_Plots/" + fit_type + "/").c_str(),
        env.GetValue("obs", "z"),
        env.GetValue("obs2", "cth"),
        (file_base + file_name).c_str(),
        bn_edgs,
        obs2bn,
        fit_type
    );


    return 0;
    }


int AsymmetryFitting(const char* file_name = "nSidis_005032.root", 
    std::string fit_type = "MhMLM", const char* config_file = "config/pippi0_RGAinbending_zbinning.txt"
    ){
    
    if (fit_type.find("MLM") != std::string::npos) {
    // MLM fit type
    return AsymmetryMLM(file_name, fit_type, config_file);
    }

    if (fit_type.find("Chi2") != std::string::npos) {
    return AsymmetryChi2(file_name, fit_type, config_file);
    }
    

    return 0;
}
