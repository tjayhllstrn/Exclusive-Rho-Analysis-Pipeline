#include "../src/InjectionStudy_utils.C"
#include "../src/InjectionStudy_utils.h"
#include <iostream>
#include <cmath>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TRandom.h>
//clas12root -l -b -q 'macros/InjectionStudy.cpp("/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/pippi0_rgaMC_in_fa18_pass2/clasdis_rga_fa18_inb_45nA_10604MeV-0001.root","A10","A01")'

bool is_RhoSig(inj_event event){
    bool pipParentIsRho = event.truepipparent_pid==213;
    bool phoParentsParentIsRho = event.truepho1parentparent_pid==213 && event.truepho2parentparent_pid==213;
    bool phoParentsSame = event.truepho1parent_id==event.truepho2parent_id;
    bool phoParentsPi0 = event.truepho1parent_pid == 111 && event.truepho2parent_pid == 111;
    bool AllRhoCuts = pipParentIsRho && phoParentsParentIsRho && phoParentsSame && phoParentsPi0;

    return AllRhoCuts;
}

void InjectionStudy(const char* inputFile, const char* InjAsym, const char* BkgInjAsym){
    const char* GREEN="\033[0;32m";
    const char* RESET="\033[0m";
    std::string inputPath(inputFile);
    std::string::size_type slashPos = inputPath.find_last_of('/');
    std::string inputDir = (slashPos == std::string::npos) ? std::string(".") : inputPath.substr(0, slashPos);
    std::string inputBase = (slashPos == std::string::npos) ? inputPath : inputPath.substr(slashPos + 1);
    std::string outputDir = inputDir + "_injAsym_" + InjAsym + "_bkg" + BkgInjAsym;
    std::string outputPath = outputDir + "/" + inputBase;

    std::cout << GREEN << "Starting injection study with input file: " << RESET << inputFile
              << GREEN << ", output file: " << RESET << outputPath.c_str()
              << GREEN << ", and injected asymmetry: " << RESET << InjAsym << std::endl;

    if (gSystem) {
        gSystem->mkdir(outputDir.c_str(), true);
    }

    TFile* inFile = TFile::Open(inputFile, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Failed to open input file: " << inputFile << std::endl;
        return;
    }

    TTree* inTree = dynamic_cast<TTree*>(inFile->Get("pippi0"));
    if (!inTree) {
        std::cerr << "Could not find tree 'pippi0' in input file." << std::endl;
        inFile->Close();
        return;
    }

    double phi = 0.0;
    double Pol = 0.0;
    double eps = 0.0;
    double t = 0.0;
    double truepipparent_pid = 0.0;
    double truepho1parentparent_pid = 0.0;
    double truepho2parentparent_pid = 0.0;
    double truepho1parent_id = 0.0;
    double truepho2parent_id = 0.0;
    double truepho1parent_pid = 0.0;
    double truepho2parent_pid = 0.0;
    inTree->SetBranchAddress("phi_true", &phi);
    inTree->SetBranchAddress("Pol", &Pol);
    inTree->SetBranchAddress("eps", &eps);
    inTree->SetBranchAddress("t_elec_true", &t);
    inTree->SetBranchAddress("truepipparent_pid", &truepipparent_pid);
    inTree->SetBranchAddress("truepho1parentparent_pid", &truepho1parentparent_pid);
    inTree->SetBranchAddress("truepho2parentparent_pid", &truepho2parentparent_pid);
    inTree->SetBranchAddress("truepho1parent_id", &truepho1parent_id);
    inTree->SetBranchAddress("truepho2parent_id", &truepho2parent_id);
    inTree->SetBranchAddress("truepho1parent_pid", &truepho1parent_pid);
    inTree->SetBranchAddress("truepho2parent_pid", &truepho2parent_pid);


    TFile* outFile = TFile::Open(outputPath.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Failed to create output file: " << outputPath << std::endl;
        inFile->Close();
        return;
    }

    outFile->cd();
    TTree* outTree = inTree->CloneTree(0);
    if (!outTree) {
        std::cerr << "Failed to clone tree structure." << std::endl;
        outFile->Close();
        inFile->Close();
        return;
    }

    int hel = 0;
    TBranch* helBranch = outTree->GetBranch("hel");
    if (!helBranch) {
        std::cerr << "Output tree does not contain branch 'hel'." << std::endl;
        outFile->Close();
        inFile->Close();
        return;
    }
    helBranch->SetAddress(&hel);

    TRandom3 seededRng(12345);
    gRandom = &seededRng;

    cross_section Rho_sigma(InjAsym);
    cross_section bkg_sigma(BkgInjAsym);

    const Long64_t nEntries = inTree->GetEntries();
    int sig_events = 0;
    for (Long64_t ev = 0; ev < nEntries; ++ev){
        if (nEntries > 0 && ev % (nEntries >= 100 ? nEntries / 100 : 1) == 0) {
            int progress = static_cast<int>((100 * ev) / nEntries);
            std::cout << progress << "% complete.\r" << std::flush;
        }

        inTree->GetEntry(ev);
        inj_event event(phi, Pol, eps, t, truepipparent_pid, truepho1parentparent_pid, truepho2parentparent_pid, truepho1parent_id, truepho2parent_id, truepho1parent_pid, truepho2parent_pid);

        if (is_RhoSig(event)) {
            sig_events++;
            hel = static_cast<int>(std::round(Rho_sigma.AssignmentRule(event)));
        }
        else{
            hel = static_cast<int>(std::round(bkg_sigma.AssignmentRule(event)));
        }
        
        outTree->Fill();
    }

    std::cout << "100% complete." << std::endl;
    std::cout << GREEN << "Writing output tree with " << RESET << outTree->GetEntries() << GREEN << " entries and " << RESET << sig_events << GREEN << " signal events to: " << RESET << outputPath << std::endl;
    outTree->Write();

    outFile->Close();
    inFile->Close();
    std::cout << GREEN << "Injection study complete." << RESET << std::endl;
}