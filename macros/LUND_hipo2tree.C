#include "../src/Structs.h"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "../src/Kinematics.C"
Kinematics kin;

//this macro searches the LUND hipo bank for Rho Plus events and calculates relvenant variables related to these events
int LUND_hipo2tree(const char* hipoFile = "",
    const char* outputFile = "LUND_output.root",
    const int maxEvents = 100){

    if(strlen(hipoFile)==0){
        std::cerr << "No input hipo file provided. Usage: LUND_hipo2tree(\"in.hipo\", \"out.root\", maxEvents)" << std::endl;
        return 1;
    }

    // Create a TFile to save the data
    TFile* fOut = new TFile(outputFile, "RECREATE");
    // Create a TTree to store the data
    TTree* outTree = new TTree("EventTree", "LUND-derived event tree");
    double Mh = -1;
    double Mdiphoton = -1;
    double rho_M = -1;
    double rho_phi = -1;
    double pion_M = -1;
    double t = -1;
    double Q2 = -1;
    double Mx = -1;
    double W = -1;
    double y = -1;
    double xF = -1;
    double z = -1;
    int neutron_parent_pindex = -1;
    int neutron_parent_pid = -1;
    outTree->Branch("Mh", &Mh, "Mh/D");
    outTree->Branch("Mdiphoton", &Mdiphoton, "Mdiphoton/D");
    outTree->Branch("rho_M", &rho_M, "rho_M/D");
    outTree->Branch("rho_phi", &rho_phi, "rho_phi/D");
    outTree->Branch("pion_M", &pion_M, "pion_M/D");
    outTree->Branch("t", &t, "t/D");
    outTree->Branch("Q2", &Q2, "Q2/D");
    outTree->Branch("Mx", &Mx, "Mx/D");
    outTree->Branch("W", &W, "W/D");
    outTree->Branch("y", &y, "y/D");
    outTree->Branch("neutron_parent_pindex", &neutron_parent_pindex, "neutron_parent_pindex/I");
    outTree->Branch("neutron_parent_pid", &neutron_parent_pid, "neutron_parent_pid/I");
    outTree->Branch("xF", &xF, "xF/D");
    outTree->Branch("z",&z,"z/D");



    // Configure CLAS12 Reader and HipoChain
    // -------------------------------------
    clas12root::HipoChain chain;
    chain.Add(hipoFile);
    auto& c12 = chain.C12ref();

    // Add the LUND bank and get field orders
    int idx_LUNDBank = chain.GetC12Reader()->addBank("MC::Lund");
    int ip_pid = c12->getBankOrder(idx_LUNDBank, "pid");
    int ip_type = c12->getBankOrder(idx_LUNDBank, "type");
    int ip_parent = c12->getBankOrder(idx_LUNDBank, "parent");
    int ip_daughter = c12->getBankOrder(idx_LUNDBank, "daughter");
    int ip_px = c12->getBankOrder(idx_LUNDBank, "px");
    int ip_py = c12->getBankOrder(idx_LUNDBank, "py");
    int ip_pz = c12->getBankOrder(idx_LUNDBank, "pz");
    int ip_E = c12->getBankOrder(idx_LUNDBank, "energy");
    int ip_mass = c12->getBankOrder(idx_LUNDBank, "mass");
    int ip_pindex = c12->getBankOrder(idx_LUNDBank, "index");

    int eventIdx = 0;
    int writtenEvents = 0;
    // Loop over events in the chain
    while(chain.Next() == true && (eventIdx < maxEvents || maxEvents < 0)){
        auto lund = c12->getBank(idx_LUNDBank);
        ++eventIdx;
        if (eventIdx % 100000 == 0 && eventIdx != 0) {
            std::cout << eventIdx << " events read | "
                      << writtenEvents * 100.0 / eventIdx << "% passed event selection"
                      << std::endl;
        }

        const int nrows = lund->getRows();

        std::vector<TLorentzVector> photons;
        TLorentzVector pion(0,0,0,0);
        TLorentzVector neutron(0,0,0,0);
        bool hasPion=false;
        bool hasNeutron=false;

        //check if this event has a rho+ before continuing - capture the pindex
        int Rhoindx = -1;
        for(int i=0; i<nrows;i++){
            int pid = -999;
            if(ip_pid>=0) pid = lund->getInt(ip_pid,i);
            if(pid==213){ // rho+ in final state
                
                Rhoindx = lund->getInt(ip_pindex,i);
                rho_M = lund->getFloat(ip_mass,i);
                rho_phi = std::atan2(lund->getFloat(ip_py,i), lund->getFloat(ip_px,i));
                break;
            }
        }
        if(Rhoindx==-1) continue;

        TLorentzVector target(0,0,0,0.938272); // proton target at rest
        for(int i=0; i<nrows;i++){
            int pid = -999;
            int type = -999;
            if(ip_pid>=0) pid = lund->getInt(ip_pid,i);
            if(ip_type>=0) type = lund->getInt(ip_type,i);
            if(pid==2212 && type==21){ // find the target proton
                target.SetPxPyPzE(lund->getFloat(ip_px,i), lund->getFloat(ip_py,i), lund->getFloat(ip_pz,i), lund->getFloat(ip_E,i));
                break;
            }
        }

        //find q:
        TLorentzVector q(0,0,0,0);
        TLorentzVector ebeam(0,0,0,0);
        TLorentzVector escatter(0,0,0,0);
        int e_pindex = -999;
        bool foundElectron = false;
        bool foundScatteredElectron = false;
        for(int i=0; i<nrows; i++){
            int pid = -999;
            int type = -999;
            pid = lund->getInt(ip_pid,i);
            type = lund->getInt(ip_type,i);
            if(type==21 && pid==11 && !foundElectron){ //beam electron
                ebeam; ebeam.SetPxPyPzE(lund->getFloat(ip_px,i), lund->getFloat(ip_py,i), lund->getFloat(ip_pz,i), lund->getFloat(ip_E,i));
                e_pindex = lund->getInt(ip_pindex,i);
                foundElectron = true;
            }
            if(type==1 && pid==11 && lund->getInt(ip_parent,i)==e_pindex && !foundScatteredElectron){ // scattered electron
                escatter; escatter.SetPxPyPzE(lund->getFloat(ip_px,i), lund->getFloat(ip_py,i), lund->getFloat(ip_pz,i), lund->getFloat(ip_E,i));
                foundScatteredElectron = true;
            }
            if(foundElectron && foundScatteredElectron){
                q = ebeam - escatter; //l - l'
                Q2 = -q.M2();
                break;
            }

        }

        for(int i=0;i<nrows;i++){
            int pid = -999, parent_pindex=-999, grandparent_pindex=-999;
            if(ip_pid>=0) pid = lund->getInt(ip_pid,i);
            float px = 0, py = 0, pz = 0, E = 0,type=-999;
            if(ip_px>=0) px = lund->getFloat(ip_px,i);
            if(ip_py>=0) py = lund->getFloat(ip_py,i);
            if(ip_pz>=0) pz = lund->getFloat(ip_pz,i);
            if(ip_E>=0)  E  = lund->getFloat(ip_E,i);
            if(ip_parent>=0) parent_pindex = lund->getInt(ip_parent,i);
            if(ip_parent>=0) grandparent_pindex = lund->getInt(ip_parent,parent_pindex-1);
            if(ip_type>=0) type = lund->getInt(ip_type,i);

            if(type!=1) continue; // only consider final state particles
            if(pid==22 && lund->getInt(ip_pid,parent_pindex-1)==111 && lund->getInt(ip_pid,grandparent_pindex-1)==213){ // photon from pi0 decay which came from a Rho
                TLorentzVector g; g.SetPxPyPzE(px,py,pz,E);
                photons.push_back(g);
            } else if(pid==211 && lund->getInt(ip_pid,parent_pindex-1)==213){ // pi+ from rho+ decay
                pion.SetPxPyPzE(px,py,pz,E);
                hasPion = true;
            } else if(pid==2112 && !hasNeutron){
                neutron.SetPxPyPzE(px,py,pz,E);
                neutron_parent_pindex = parent_pindex;
                neutron_parent_pid = lund->getInt(ip_pid,parent_pindex-1);
                hasNeutron = true;
            }
        }

        // Require at least one pi+, at least two photons, and one neutron
        if(hasPion && hasNeutron && photons.size() >= 2){
            // pick the two highest-energy photons
            // std::cout << "num photons: " << photons.size() << std::endl;
            std::sort(photons.begin(), photons.end(), [](const TLorentzVector &a, const TLorentzVector &b){ return a.E() > b.E(); });
            TLorentzVector diphoton = photons[0] + photons[1];
            TLorentzVector Ph = pion + diphoton;
            pion_M = pion.M();
            Mh = Ph.M();
            t = (q - Ph).M2(); 
            Mx = (q + target - Ph).M();
            W = (q + target).M();
            y = (target.Dot(q)) / (target.Dot(ebeam)); // inelasticity y = (P.q)/(P.k)
            xF = kin.xF(q,diphoton,target,W);
            z = kin.z(target,Ph,q);
            // t = (ebeam - escatter - Ph).M2(); 
            Mdiphoton = (diphoton).M();
            outTree->Fill();
            ++writtenEvents;
        }
    }

    fOut->cd();
    outTree->Write();
    fOut->Close();
    std::cout << "Written " << writtenEvents << " events to " << outputFile << std::endl;

    return 0;
}


    