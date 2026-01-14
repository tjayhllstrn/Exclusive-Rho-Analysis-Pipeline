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
    // Read config file
    TEnv env; 
    env.ReadFile(config_file, kEnvLocal);
    
    //Initialize and parse fitter
    std::vector<double> phibn_edges = parse_csv_to_doubles(env.GetValue("phibn_edges", ""));
    std::vector<double> bn_edgs = parse_csv_to_doubles(env.GetValue("bn_edgs", ""));
    std::vector<double> obs2bn = parse_csv_to_doubles(env.GetValue("obs2bn", ""));
    std::string file_base = env.GetValue("file_base", "./out/pippi0_fall2018_in_pass2/");
    Chi2_Fitter fitter(
        env.GetValue("treename", "pippi0"),
        (file_base + "BSA_Plots/" + fit_type + "/").c_str(),
        env.GetValue("obs", "z"),
        env.GetValue("obs2", "cth"),
        (file_base + file_name).c_str(),
        phibn_edges,
        bn_edgs,
        obs2bn,
        fit_type
    );

    //Loop through obs2 bins and run Chi2 fitting based on the fit_type
    for(size_t bn_idx = 0; bn_idx < fitter.OBS2BN.size() - 1; bn_idx++){
        fitter.OUT_DIR = file_base + "BSA_Plots/" + fit_type + "/" + fitter.OBS2 + "_" + std::to_string(fitter.OBS2BN[bn_idx]).substr(0,5) + "_" + std::to_string(fitter.OBS2BN[bn_idx+1]).substr(0,5) + "/";
        // Check if directory exists and create if needed
        if (!std::filesystem::exists(fitter.OUT_DIR)) {
            std::filesystem::create_directories(fitter.OUT_DIR);
        }
        //Make sure there is no SinFits.root file already existing in the output directory
        std::string root_filename = fitter.OUT_DIR + "SinFits.root";
        if (std::filesystem::exists(root_filename)) {
            std::filesystem::remove(root_filename);
        }
        fitter.obs2_bin_idx = bn_idx;
        // MhChi2 fit
        if(fit_type == "MhChi2"){
            fitter.RunMhChi2Fit(bn_idx);
        }
        
        // MxChi2 fit
        if(fit_type == "MxChi2"){
            fitter.RunMxChi2Fit(bn_idx);
        }

        //Plot N_sig values to compare neg and pos helicity
        fitter.PlotToCanvas_N_sig_BarHist();

        //plot A_sig!!
        fitter.MakeGraphLinePlot(fitter.A_sig, fitter.BN_CENTERS,
            "F_{LU}^{sin#phi}/F_{UU} [arb. units]",
            fitter.OBS.c_str(),
            "A_sig",
            kBlack);

        //Plot N_sig extraction fits. One postage stamp plot per obs bin - first need to flatten arrays for easier plotting
        for(size_t obs_bin_idx = 0; obs_bin_idx < fitter.BN_EDGS.size() -1; obs_bin_idx++){
            std::vector<TH1F*> data_hists;
            std::vector<TGraph*> total_graphs;
            std::vector<TGraph*> sig_graphs;
            std::vector<TGraph*> bkg_graphs;
            std::vector<TLegend*> legends;
            std::vector<TLatex*> texts;
            std::vector<TPaveText*> param_boxes;
            std::string title = "ExtractionFits_" + fitter.OBS + "(" + std::to_string(fitter.BN_EDGS[obs_bin_idx]).substr(0,5) + ", " + std::to_string(fitter.BN_EDGS[obs_bin_idx+1]).substr(0,5) + ")";

            for(size_t phi_bin_idx = 0; phi_bin_idx < fitter.PHIBN_EDGES.size() -1; phi_bin_idx++){
                data_hists.push_back(fitter.N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].first);
                data_hists.push_back(fitter.N_sig_fitting_datathist[obs_bin_idx][phi_bin_idx].second);
                total_graphs.push_back(fitter.N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].first);
                total_graphs.push_back(fitter.N_sig_fitting_totalgraph[obs_bin_idx][phi_bin_idx].second);
                sig_graphs.push_back(fitter.N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].first);
                sig_graphs.push_back(fitter.N_sig_fitting_siggraph[obs_bin_idx][phi_bin_idx].second);
                bkg_graphs.push_back(fitter.N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].first);
                bkg_graphs.push_back(fitter.N_sig_fitting_bkggraph[obs_bin_idx][phi_bin_idx].second);
                legends.push_back(fitter.N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].first);
                legends.push_back(fitter.N_sig_fitting_legends[obs_bin_idx][phi_bin_idx].second);
                texts.push_back(fitter.N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].first);
                texts.push_back(fitter.N_sig_fitting_texts[obs_bin_idx][phi_bin_idx].second);
                param_boxes.push_back(fitter.N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].first);
                param_boxes.push_back(fitter.N_sig_fitting_paramboxes[obs_bin_idx][phi_bin_idx].second);
            }

            fitter.PlotToCanvas_PostageStamp(data_hists,
                                            total_graphs,
                                            sig_graphs,
                                            bkg_graphs,
                                            legends,
                                            texts,
                                            param_boxes,
                                            title);
        
        //clear results containers
        fitter.N_sig_pos.clear(); 
        fitter.N_sig_neg.clear(); 
        fitter.alpha.clear();      
        fitter.A_sig.clear();  
        }
    }
    return 0;
    }


int AsymmetryFitting(const char* file_name = "nSidis_005032.root", 
    std::string fit_type = "MhChi2", const char* config_file = "config/pippi0_RGAinbending_zbinning.txt")
    {
    
    if (fit_type.find("MLM") != std::string::npos) {
    // MLM fit type
    return AsymmetryMLM(file_name, fit_type, config_file);
    }

    if (fit_type.find("Chi2") != std::string::npos) {
    return AsymmetryChi2(file_name, fit_type, config_file);
    }
    

    return 0;
}
