#include "src/InjectionStudy_utils.C"
#include "src/InjectionStudy_utils.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TH1.h>
#include <TKey.h>
#include <TObject.h>

#include <algorithm>
#include <glob.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

std::vector<std::string> get_root_file_names(const char* Asym){
    std::vector<std::string> root_file_names;
    const std::string pattern = "out/pippi0_rgaMC_in_fa18_pass2_injAsym_" + std::string(Asym) + "*";

    glob_t glob_result;
    if (glob(pattern.c_str(), GLOB_TILDE, nullptr, &glob_result) != 0) {
        return root_file_names;
    }

    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        for (const auto& entry : fs::recursive_directory_iterator(glob_result.gl_pathv[i])) {
            if (!entry.is_regular_file()) {
                continue;
            }

            const std::string file_name = entry.path().filename().string();
            const bool is_a_sig_file = file_name == "A_sig.root";
            const bool is_asym_file = file_name.rfind("AsymmetriesAsymmetries_per", 0) == 0 && entry.path().extension() == ".root";

            if (is_a_sig_file || is_asym_file) {
                root_file_names.push_back(entry.path().string());
            }
        }
    }

    globfree(&glob_result);

    std::sort(root_file_names.begin(), root_file_names.end());
    return root_file_names;
}



static TGraph* build_asymmetry_line(const TGraph* original, const char* asym_type, int n_points = 200){
    if (original == nullptr || original->GetN() == 0) {
        return nullptr;
    }

    double xmin = 0.0;
    double xmax = 0.0;
    double x = 0.0;
    double y = 0.0;
    original->GetPoint(0, x, y);
    xmin = x;
    xmax = x;

    for (int i = 1; i < original->GetN(); ++i) {
        original->GetPoint(i, x, y);
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
    }

    if (xmin == xmax) {
        xmin -= 1.0;
        xmax += 1.0;
    }

    TGraph* asymmetry_line = new TGraph(n_points);
    asymmetry_line->SetName("A_sig_injected_asymmetry_line");
    asymmetry_line->SetTitle("Injected asymmetry function");

    cross_section cc(asym_type);
    const double step = (xmax - xmin) / static_cast<double>(n_points - 1);
    for (int i = 0; i < n_points; ++i) {
        const double x_value = xmin + step * static_cast<double>(i);
        asymmetry_line->SetPoint(i, x_value, cc.Asymmetry(x_value));
    }

    return asymmetry_line;
}
static void draw_root_object(TFile& input_file, const char* asym_type){
    TGraph* original = dynamic_cast<TGraph*>(input_file.Get("A_sig"));
    TGraph* A_bkg = dynamic_cast<TGraph*>(input_file.Get("A_bkg"));
    TGraph* A_sig_bkg = dynamic_cast<TGraph*>(input_file.Get("A_sig_bkg"));
    TLegend* leg = dynamic_cast<TLegend*>(input_file.Get("TPave"));

    std::string base_name = fs::path(input_file.GetName()).stem().string();
    std::string parent_dir = fs::path(input_file.GetName()).parent_path().string();
    std::string out_png = parent_dir + "/" + base_name + std::string("_withFunction.png");

    if (original == nullptr) {
        std::cerr << "No drawable object found in " << input_file.GetName() << std::endl;
        return;
    }

    TGraph* asymmetry_line = build_asymmetry_line(original, asym_type);
    if (asymmetry_line == nullptr) {
        std::cerr << "Failed to build asymmetry line for " << input_file.GetName() << std::endl;
        return;
    }

    // original->SetMarkerStyle(kFullCircle);
    // original->SetMarkerColor(kBlue + 1);
    // original->SetLineColor(kBlue + 1);

    asymmetry_line->SetLineColor(kRed + 1);
    asymmetry_line->SetLineWidth(2);
    asymmetry_line->SetLineStyle(2);

    TCanvas canvas("c_view_Asym_inj_results", "c_view_Asym_inj_results", 800, 600);
    canvas.SetGrid();

    const char* label = "Chi2 Fit";
    if (A_bkg) {
        label = "MLM Fit";
    }

    original->SetTitle((std::string("Extracted ") + label + std::string(" vs Injected Signal Asymmetry")).c_str());
    original->Draw("AP");
    asymmetry_line->Draw("L SAME");

    
    if (leg) {
        leg->AddEntry(asymmetry_line, "Injected Signal Asymmetry", "l");
        leg->SetX1NDC(0.45);

        TLegendEntry* entry = (TLegendEntry*)leg->GetListOfPrimitives()->At(0);
        entry->SetLabel("A sig + bkg");
        TLegendEntry* entry1 = (TLegendEntry*)leg->GetListOfPrimitives()->At(1);
        entry1->SetLabel("A Background");
        TLegendEntry* entry2 = (TLegendEntry*)leg->GetListOfPrimitives()->At(2);
        entry2->SetLabel("A Signal");

        leg->Draw();
    } else {
        TLegend* legend = new TLegend(0.45, 0.75, 0.88, 0.88);
        legend->AddEntry(original, "Extracted Signal Asymmetry", "p");
        legend->AddEntry(asymmetry_line, "Injected Signal Asymmetry", "l");
        legend->SetBorderSize(0);
        legend->Draw();
    }
    if (A_bkg) {
        A_bkg->Draw("P SAME");
    }
    if (A_sig_bkg) {
        A_sig_bkg->Draw("P SAME");
    }

    canvas.SaveAs(out_png.c_str());
    input_file.cd();
    asymmetry_line->Write("A_injected");
    input_file.Close();


    delete asymmetry_line;
}

void view_Asym_inj_results(){
    const char* Asym_types[3] = {"A10", "Alin_t","A04"};

    std::cout << "plotting Asymmetry results for: " << std::endl;
    for (const char* asym_type : Asym_types) {
        std::cout << asym_type << std::endl;
        std::vector<std::string> files_for_asym = get_root_file_names(asym_type);

        std::cout << "Found " << files_for_asym.size() << " matching ROOT files." << std::endl;
        for (const std::string& root_file_name : files_for_asym) {
            fs::path input_path(root_file_name);
            TFile input_file(root_file_name.c_str(), "UPDATE");
            if (input_file.IsZombie()) {
                std::cerr << "Failed to open: " << root_file_name << std::endl;
                continue;
            }

            draw_root_object(input_file, asym_type);
            std::cout << "Saved plot: " << input_file.GetName() << std::endl;
        }
    }

}