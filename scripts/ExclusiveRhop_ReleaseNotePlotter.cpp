//full data: root -l -b -q 'scripts/ExclusiveRhop_ReleaseNotePlotter.cpp("/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/pippi0_merged_in_pass2/pippi0_merged_in_pass2.root")'
//Mh cut short: root -l -b -q 'scripts/ExclusiveRhop_ReleaseNotePlotter.cpp("/volatile/clas12/users/tjhellst/cache/filtered_tree_pippi0_merged_in_pass2_MhChi2_MhFit_cthbin0_-0.190000_0.270000_MinPhoE_0.200000.root")'
//just fall: root -l -b -q 'scripts/ExclusiveRhop_ReleaseNotePlotter.cpp("/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/pippi0_fall2018_in_pass2/pippi0_fall2018_in_pass2.root")'
#include "../out/ExclusiveRhop_ReleaseNotePlots/style.C"
//Global varialbes for cuts
TCut MX_CUT = TCut("MX_CUT", "0.85 < Mx&&Mx < 1.05");
TCut MH_CUT = TCut("MH_CUT", "0.65 < Mh&&Mh < 0.9");
TCut MDIPHOTON_CUT = TCut("MDIPHOTON_CUT", "0.104 < Mdiphoton&&Mdiphoton < 0.164");
TCut MINPHOE = TCut("MINPHOE", "pho1_E > 0.2 && pho2_E > 0.2");


TCut SIDEBAND_CUT = TCut("SIDEBAND_CUT", "(1.1<Mh&&Mh<1.7) | (0.25<Mh&&Mh<0.5)");
TCut SIGNAL_CUT = TCut("SIGNAL_CUT", "0.65<Mh&&Mh<0.9");
TCut SIDEBAND_CUT_Mx = TCut("SIDEBAND_CUT_Mx", "1.4<Mx&&Mx<2.75");
TCut SIGNAL_CUT_Mx = TCut("SIGNAL_CUT_Mx", "0.85<Mx&&Mx<1.15");


int FillHistogram(TH1F* h, const std::string& var, TTree* tree, const std::string& selection = "") {
    // Fill the histogram with data (this is just a placeholder)
    if (selection.empty()) {
        tree->Draw((var + ">>" + h->GetName()).c_str());
    } else {
        tree->Draw((var + ">>" + h->GetName()).c_str(), selection.c_str());
    }
    if (tree->GetEntries() == 0) {
        std::cerr << "Warning: No entries found for variable " << var << " with selection: " << selection << std::endl;
        return 0;
    }
    return 1;
}

int DrawRegionsHistogram(TTree* tree){
    //Plot Mh regions
    TH1F Mh = TH1F("Mh", "Neutron Mass Cut;M_{h} [GeV];Events", 100, 0.2, 1.8);
    FillHistogram(&Mh, "Mh", tree, (MX_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());

    TH1F Mh_sideband = TH1F("Mh_sideband", "", 100, 0.2, 1.8);
    FillHistogram(&Mh_sideband, "Mh", tree, (SIDEBAND_CUT+MX_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());
    Mh_sideband.SetLineColor(MAGENTA);
    Mh_sideband.SetFillColor(MAGENTA);

    TH1F Mh_signal = TH1F("Mh_signal", "", 100, 0.2, 1.8);
    FillHistogram(&Mh_signal, "Mh", tree, (SIGNAL_CUT+MX_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());
    Mh_signal.SetLineColor(BLUE);
    Mh_signal.SetFillColor(BLUE);

    TLegend legend(0.6, 0.7, 1, 0.85);
    legend.AddEntry(&Mh_signal, "Signal Region", "f");
    legend.AddEntry(&Mh_sideband, "Sideband Region", "f");

    //Plot Mx Regions
    TH1F Mx = TH1F("Mx", "#rho Mass Cut;M_{x} [GeV];Events", 100, 0.2, 2.8);
    FillHistogram(&Mx, "Mx", tree, (MH_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());

    TH1F Mx_sideband = TH1F("Mx_sideband", "", 100, 0.2, 2.8);
    FillHistogram(&Mx_sideband, "Mx", tree, (SIDEBAND_CUT_Mx+MH_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());
    Mx_sideband.SetLineColor(MAGENTA);
    Mx_sideband.SetFillColor(MAGENTA);

    TH1F Mx_signal = TH1F("Mx_signal", "", 100, 0.2, 2.8);
    FillHistogram(&Mx_signal, "Mx", tree, (SIGNAL_CUT_Mx+MH_CUT+MDIPHOTON_CUT+MINPHOE).GetTitle());
    Mx_signal.SetLineColor(BLUE);
    Mx_signal.SetFillColor(BLUE);

    TCanvas RegionsCanvas("RegionsCanvas", "Regions Canvas", 1400, 800);
    RegionsCanvas.Divide(2,1);
    RegionsCanvas.cd(1);
    Mh.Draw();
    Mh_sideband.Draw("same");
    Mh_signal.Draw("same");
    Mh.Draw("same");
    legend.Draw();
    RegionsCanvas.cd(2);
    Mx.Draw();
    Mx_sideband.Draw("same");
    Mx_signal.Draw("same");
    Mx.Draw("same");
    RegionsCanvas.SaveAs("out/ExclusiveRhop_ReleaseNotePlots/Regions.png");
    RegionsCanvas.SaveAs("out/ExclusiveRhop_ReleaseNotePlots/Regions.root");
    RegionsCanvas.SaveAs("../6a0369d69f5534ce43bd2e59/figures/Regions.pdf");
    return 1;
}

int DrawMdiphotonHistogram(TTree* tree){
    TH1F Mdiphoton_Mx = TH1F("Mdiphoton_Mx", ";M_{#gamma#gamma} [GeV];Events", 100, 0.0, 0.5);
    TH1F Mdiphoton_Mh = TH1F("Mdiphoton_Mh", ";M_{#gamma#gamma} [GeV];Events", 100, 0.0, 0.5);
    FillHistogram(&Mdiphoton_Mx, "Mdiphoton", tree, (MH_CUT+MINPHOE).GetTitle());
    FillHistogram(&Mdiphoton_Mh, "Mdiphoton", tree, (MX_CUT+MINPHOE).GetTitle());


    TCanvas MdiphotonCanvas("MdiphotonCanvas", "Mdiphoton Canvas", 1400, 800);
    MdiphotonCanvas.Divide(2,1);
    MdiphotonCanvas.cd(1);
    Mdiphoton_Mh.Draw();
    MdiphotonCanvas.cd(2);
    Mdiphoton_Mx.Draw();
    MdiphotonCanvas.SaveAs("out/ExclusiveRhop_ReleaseNotePlots/Mdiphoton.png");
    MdiphotonCanvas.SaveAs("out/ExclusiveRhop_ReleaseNotePlots/Mdiphoton.root");
    MdiphotonCanvas.SaveAs("../6a0369d69f5534ce43bd2e59/figures/Mdiphoton.pdf");
    return 1;
}

int DrawStandardHistogram(TTree* tree, const std::string& var, const std::string& xlabel, double min, double max) {
    TH1F Mx = TH1F("Mx", (" ;" + xlabel + ";Events").c_str(), 100, min, max);
    TH1F Mh = TH1F("Mh", (" ;" + xlabel + ";Events").c_str(), 100, min, max);
    FillHistogram(&Mx, var, tree, (MH_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());
    FillHistogram(&Mh, var, tree, (MX_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());

    TH1F Mx_sideband = TH1F("Mx_sideband", "", 100, min, max);
    TH1F Mh_sideband = TH1F("Mh_sideband", "", 100, min, max);
    FillHistogram(&Mx_sideband, var, tree, (SIDEBAND_CUT_Mx+MH_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());
    FillHistogram(&Mh_sideband, var, tree, (SIDEBAND_CUT+MX_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());
    Mx_sideband.SetLineColor(MAGENTA);
    Mh_sideband.SetLineColor(MAGENTA);

    TH1F Mx_signal = TH1F("Mx_signal", "", 100, min, max);
    TH1F Mh_signal = TH1F("Mh_signal", "", 100, min, max);
    FillHistogram(&Mx_signal, var, tree, (SIGNAL_CUT_Mx+MH_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());
    FillHistogram(&Mh_signal, var, tree, (SIGNAL_CUT+MX_CUT+MINPHOE+MDIPHOTON_CUT).GetTitle());
    Mx_signal.SetLineColor(BLUE);
    Mh_signal.SetLineColor(BLUE);

    TLegend leg = TLegend(0.2, 0.9, 1, 0.95);
    //dummy objects for legend
    TH1F sideband_box = TH1F("sideband_box", "", 1, 0, 1);
    sideband_box.SetFillColor(Mh_sideband.GetLineColor());
    TH1F signal_box = TH1F("signal_box", "", 1, 0, 1);
    signal_box.SetFillColor(Mh_signal.GetLineColor());
    TH1F total_box = TH1F("total_box", "", 1, 0, 1);
    total_box.SetFillColor(Mh.GetLineColor());
    leg.AddEntry(&total_box, "Total", "f");
    leg.AddEntry(&signal_box, "Signal Region", "f");
    leg.AddEntry(&sideband_box, "Sideband Region", "f");
    leg.SetNColumns(3);

    Mh.SetMinimum(0);
    Mx.SetMinimum(0);
    
    TCanvas Canvas("Canvas", "Canvas", 1400, 800);
    Canvas.Divide(2,1);
    Canvas.cd(1);
    Mh.Draw();
    Mh_sideband.Draw("same");
    Mh_signal.Draw("same");
    leg.Draw();
    Canvas.cd(2);
    Mx.Draw();
    Mx_sideband.Draw("same");
    Mx_signal.Draw("same");
    leg.Draw();
    Canvas.SaveAs(("out/ExclusiveRhop_ReleaseNotePlots/" + var + ".png").c_str());
    Canvas.SaveAs(("out/ExclusiveRhop_ReleaseNotePlots/" + var + ".root").c_str());
    Canvas.SaveAs(("../6a0369d69f5534ce43bd2e59/figures/" + var + ".pdf").c_str());

    return 1;
}

void ExclusiveRhop_ReleaseNotePlotter(const char* root_file, bool update_overleaf = false) {
    SetReleaseNoteStyle();
    DefineColors();
    std::cout << "Plotting from file: " << root_file << std::endl;
    TFile* file = TFile::Open(root_file, "READ");
    TTree* pippi0 = (TTree*)file->Get("pippi0");
    TTree* EventTree = (TTree*)file->Get("EventTree");

    if (update_overleaf) {
        std::cout << "Pulling latest changes from Overleaf..." << std::endl;
        system("cd ../6a0369d69f5534ce43bd2e59 && git pull");
    }


    //Draw Kinematic Plots
    // DrawRegionsHistogram(pippi0);
    // DrawMdiphotonHistogram(pippi0);

    // std::vector<std::string> vars       = {"Q2", "W", "x", "y", "phi", "t_elec","z","varphi", "cth"};
    // std::vector<std::string> var_labels = {"Q^{2} [GeV^{2}]", "W [GeV]", "x", "y", "#phi [rad]", "t [GeV^{2}]", "z", "#varphi [rad]", "cos#theta"};
    // std::vector<double> min_vals        = {0, 2.0, 0, 0.2, -3.14,-6,0,-3.14,-1};
    // std::vector<double> max_vals        = {10, 4, 1, 0.8, 3.14,0,1,3.14,1};
    std::vector<std::string> vars       = {"phi", "t_elec","z","varphi", "cth"};
    std::vector<std::string> var_labels = {"#phi [rad]", "t [GeV^{2}]", "z", "#varphi [rad]", "cos#theta"};
    std::vector<double> min_vals        = {-3.14,-6,0,-3.14,-1};
    std::vector<double> max_vals        = {3.14,0,1,3.14,1};
    for (int i = 0; i < vars.size(); ++i) {
        DrawStandardHistogram(pippi0, vars[i], var_labels[i], min_vals[i], max_vals[i]);
    }
    if (update_overleaf) {
        std::cout << "Updating Overleaf..." << std::endl;
        system("cd ../6a0369d69f5534ce43bd2e59 && git add figures/*.pdf");
        system("cd ../6a0369d69f5534ce43bd2e59 && git commit -m \"update plots\"");
        system("cd ../6a0369d69f5534ce43bd2e59 && git push");
    }
}
