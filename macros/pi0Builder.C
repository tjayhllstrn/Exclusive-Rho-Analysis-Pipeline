#include "../src/Constants.h"
#include "../src/Kinematics.C"
//  ./macros/pi0Builder.C  out/

int pi0Builder(const char *input_file="out/test/nSidis_005032.root"){
    
    // Read the TFile
    TFile *f = new TFile(input_file,"UPDATE");
    // Read the TTree
    TTree *EventTree = (TTree*)f->Get("EventTree");
    
    double x, Q2, W, Pol,y,nu;
    double truex, trueQ2, trueW;
    double tPol;
    int hel,run,A,_evnum,hwp,tSign,target;
    int Nmax=100;
    int isMC=0;
    double px[Nmax], py[Nmax], pz[Nmax], E[Nmax], vz[Nmax], chi2[Nmax], theta[Nmax], eta[Nmax], phi[Nmax];
    double truepx[Nmax] , truepy[Nmax] , truepz[Nmax], trueE[Nmax], truetheta[Nmax], trueeta[Nmax], truephi[Nmax];
    int trueparentid[Nmax],trueparentpid[Nmax],trueparentparentid[Nmax],trueparentparentpid[Nmax];
    int sector[Nmax];
    int pid[Nmax], truepid[Nmax];
    double p_gamma[Nmax];
    //link the TBranches to the variables
    EventTree->SetBranchAddress("A",&A);
    EventTree->SetBranchAddress("evnum",&_evnum);
    EventTree->SetBranchAddress("run",&run);
    EventTree->SetBranchAddress("Pol",&Pol);
    EventTree->SetBranchAddress("tPol",&tPol);
    EventTree->SetBranchAddress("target",&target);
    EventTree->SetBranchAddress("tSign",&tSign);
    EventTree->SetBranchAddress("hwp",&hwp);
    EventTree->SetBranchAddress("hel",&hel);
    EventTree->SetBranchAddress("x",&x);
    EventTree->SetBranchAddress("Q2",&Q2);
    EventTree->SetBranchAddress("truex",&truex);
    EventTree->SetBranchAddress("trueQ2",&trueQ2);
    EventTree->SetBranchAddress("trueW",&trueW);
    EventTree->SetBranchAddress("W",&W); 
    EventTree->SetBranchAddress("y",&y); 
    EventTree->SetBranchAddress("nu",&nu); 
    EventTree->SetBranchAddress("Nmax",&Nmax);
    EventTree->SetBranchAddress("px",px);
    EventTree->SetBranchAddress("py",py);
    EventTree->SetBranchAddress("pz",pz);
    EventTree->SetBranchAddress("E",E);
    EventTree->SetBranchAddress("sector",sector);
    EventTree->SetBranchAddress("vz",vz);
    EventTree->SetBranchAddress("chi2",chi2);
    EventTree->SetBranchAddress("pid",pid);
    EventTree->SetBranchAddress("theta",theta);
    EventTree->SetBranchAddress("eta",eta);
    EventTree->SetBranchAddress("phi",phi);
    EventTree->SetBranchAddress("trueE",&trueE);
    EventTree->SetBranchAddress("truepx",&truepx);
    EventTree->SetBranchAddress("truepy",&truepy);
    EventTree->SetBranchAddress("truepz",&truepz);
    EventTree->SetBranchAddress("truetheta",&truetheta);
    EventTree->SetBranchAddress("trueeta",&trueeta);
    EventTree->SetBranchAddress("truephi",&truephi);
    EventTree->SetBranchAddress("truepid",truepid);
    EventTree->SetBranchAddress("trueparentid", &trueparentid);
    EventTree->SetBranchAddress("trueparentpid", &trueparentpid);
    EventTree->SetBranchAddress("trueparentparentid", &trueparentparentid);
    EventTree->SetBranchAddress("trueparentparentpid", &trueparentparentpid);
    EventTree->SetBranchAddress("p_gamma",&p_gamma);
    TString treename = "";
    if (f->Get("pi0")) f->Delete("pi0*;*");
    treename = "pi0";
    TTree *outtree = new TTree(treename.Data(),"Tree");
    double z, pT, phih, Mx, Mgg, xF, Mh,eps,gamma;
    double z_true, pT_true, phih_true, Mx_true, xF_true, Mh_true;
    Mh=0;
    Mh_true=0;
    // Branching kinematic variables for the electron
    outtree->Branch("hel",&hel,"hel/I");
    outtree->Branch("run",&run,"run/I");
    outtree->Branch("x",&x,"x/D");
    outtree->Branch("eps",&eps,"eps/D");
    outtree->Branch("gamma",&gamma,"gamma/D");
    outtree->Branch("x_true",&truex,"x_true/D");
    outtree->Branch("Q2",&Q2,"Q2/D");
    outtree->Branch("Q2_true",&trueQ2,"Q2_true/D");
    
    outtree->Branch("z", &z, "z/D");
    outtree->Branch("pT", &pT, "pT/D");
    outtree->Branch("xF", &xF, "xF/D");
    outtree->Branch("phi", &phih, "phi/D");
    outtree->Branch("Mx", &Mx, "Mx/D");
    outtree->Branch("Mh", &Mh, "Mh/D");
    outtree->Branch("Mgg",&Mgg,"Mgg/D");
    
    outtree->Branch("z_true", &z_true, "z_true/D");
    outtree->Branch("pT_true", &pT_true, "pT_true/D");
    outtree->Branch("xF_true", &xF_true, "xF_true/D");
    outtree->Branch("phi_true", &phih_true, "phi_true/D");
    outtree->Branch("Mx_true", &Mx_true, "Mx_true/D");
    outtree->Branch("Mh_true", &Mh_true, "Mh_true/D");

    // Initial particles
    double eBeam = 10.6041;
    double mE = 0.000511;
    TLorentzVector init_electron(0,0,sqrt(eBeam*eBeam-mE*mE),eBeam);
    TLorentzVector init_target(0,0,0,0.938272);
    
    TLorentzVector pi0, g1,g2;
    TLorentzVector truepi0, trueg1,trueg2;
    // for loop over all events
    int N = EventTree->GetEntries();
    Kinematics kin;

    for (int ev = 0; ev < N; ++ev) {
        int step = (N >= 100 ? N/100 : 1);
        if (ev % step == 0) {
          int progress = (100 * ev) / N;
          std::cout << progress << "% complete.\r" << std::flush;
        }
        
        EventTree->GetEntry(ev);

        if(abs(run)==11){
            isMC=1;
        }
        
        //Loop over all particles in the event to find electron
        TLorentzVector electron;
        TLorentzVector trueelectron;
        TLorentzVector q; // virtual photon
        TLorentzVector trueq;
        int idx_e=-1;
        double max_e=-1;
        for (int i=0; i<Nmax; i++){
            if(isMC==1){
                if(trueE[i]>max_e&&truepid[i]==11){
                    idx_e=i;
                    max_e=trueE[i];
                }
            }
            else
            {
                if(E[i]>max_e&&pid[i]==11){
                    idx_e=i;
                    max_e=E[i];
                }
            }
        }
        
        if(idx_e==-1) continue;
        
        electron.SetPxPyPzE(px[idx_e],py[idx_e],pz[idx_e],E[idx_e]);
        trueelectron.SetPxPyPzE(truepx[idx_e],truepy[idx_e],truepz[idx_e],trueE[idx_e]);
        q = init_electron-electron;
        trueq = init_electron-trueelectron;

        for(int j = 0; j<Nmax;j++){
            if(pid[j]!=22)continue;
            g1.SetPxPyPzE(px[j],py[j],pz[j],E[j]);
            trueg1.SetPxPyPzE(truepx[j],truepy[j],truepz[j],trueE[j]);
            for(int k = j+1; k<Nmax;k++){
                if(pid[k]!=22)continue;
                g2.SetPxPyPzE(px[k],py[k],pz[k],E[k]);
                trueg2.SetPxPyPzE(truepx[k],truepy[k],truepz[k],trueE[k]);



                pi0 = g1+g2;
                truepi0 = trueg1+trueg2;
                z = kin.z(init_target,pi0,q);
                z_true = kin.z(init_target,truepi0,trueq);

                pT = kin.Pt(q,pi0,init_target);
                pT_true = kin.Pt(trueq,truepi0,init_target);

                phih = kin.phi_h(q,init_electron,pi0);
                phih_true = kin.phi_h(trueq,init_electron,truepi0);

                Mx = (init_electron+init_target-electron-pi0).M();
                Mx_true = (init_electron+init_target-trueelectron-truepi0).M();

                xF = kin.xF(q,pi0,init_target,W);
                xF_true = kin.xF(trueq,truepi0,init_target,trueW);
                
                Mgg = (g1+g2).M();
                eps=(1-y-pow(y*gamma,2)/4)/(1-y+pow(y,2)/2+pow(y*gamma,2)/4);
                gamma = 2*0.938272*x/sqrt(Q2);
                // Test true pids
                if (isMC){
                    if (!((trueparentid[j]==trueparentid[k])&&(trueparentpid[j]==111))){continue;} // pi0
                    if (Mx_true<1.5){continue;}
                    if (trueg1.E()<0.2){continue;}
                    if (trueg2.E()<0.2){continue;}
                    if (trueg1.Angle(trueelectron.Vect())< 8 * TMath::DegToRad()){continue;}
                    if (trueg2.Angle(trueelectron.Vect())< 8 * TMath::DegToRad()){continue;}
                    if (xF_true<0){continue;}
                }

                if(!(electron.E()>0&&g1.E()>0&&g2.E()>0&&p_gamma[j]>0.78&&p_gamma[k]>0.78&&xF>0)){ //makes cut based on GBT prediction. Both photons must have an assigned weight of over 0.78 to be considered a photon
                    continue;
                }
                outtree->Fill();
            }
        }
    }

    
    cout << "Writing Total TTree with " << outtree->GetEntries() << " entries" << endl;
    outtree->Write();
    f->Close();
    cout << "Done" << endl;
    return 0;
}


