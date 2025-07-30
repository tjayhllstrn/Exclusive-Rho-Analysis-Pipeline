#include "../src/Constants.h"
#include "../src/Kinematics.C"


int pippimBuilder(const char *input_file="out/test/nSidis_005032.root"){
    
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
    TString treename = "";
    if (f->Get("pippim")) f->Delete("pippim*;*");
    treename = "pippim";
    TTree *outtree = new TTree(treename.Data(),"Tree");
    double z, pT, phih, Mx, xF, xF1, xF2, Mh,eps,gamma,th,cth;
    double z_true, pT_true, phih_true, Mx_true, xF_true, xF1_true, xF2_true, Mh_true,th_true,cth_true;
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
    outtree->Branch("xF1", &xF1, "xF1/D");
    outtree->Branch("xF2", &xF2, "xF2/D");
    outtree->Branch("phi", &phih, "phi/D");
    outtree->Branch("Mx", &Mx, "Mx/D");
    outtree->Branch("Mh", &Mh, "Mh/D");
    
    outtree->Branch("z_true", &z_true, "z_true/D");
    outtree->Branch("pT_true", &pT_true, "pT_true/D");
    outtree->Branch("xF_true", &xF_true, "xF_true/D");
    outtree->Branch("xF1_true", &xF1_true, "xF1_true/D");
    outtree->Branch("xF2_true", &xF2_true, "xF2_true/D");
    outtree->Branch("phi_true", &phih_true, "phi_true/D");
    outtree->Branch("Mx_true", &Mx_true, "Mx_true/D");
    outtree->Branch("Mh_true", &Mh_true, "Mh_true/D");

    outtree->Branch("th", &th, "th/D");
    outtree->Branch("th_true", &th_true, "th_true/D");

    outtree->Branch("cth", &cth, "cth/D");
    outtree->Branch("cth_true", &cth_true, "cth_true/D");

    // Initial particles
    double eBeam = 10.6041;
    double mE = 0.000511;
    TLorentzVector init_electron(0,0,sqrt(eBeam*eBeam-mE*mE),eBeam);
    TLorentzVector init_target(0,0,0,0.938272);
    
    TLorentzVector pip, pim, dihadron;
    TLorentzVector truepip, truepim, truedihadron;
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
            if(pid[j]!=211)continue;
            pip.SetPxPyPzE(px[j],py[j],pz[j],E[j]);
            truepip.SetPxPyPzE(truepx[j],truepy[j],truepz[j],trueE[j]);
            for(int k = 0; k<Nmax;k++){
                if(pid[k]!=-211)continue;
                pim.SetPxPyPzE(px[k],py[k],pz[k],E[k]);
                truepim.SetPxPyPzE(truepx[k],truepy[k],truepz[k],trueE[k]);
                dihadron = pip+pim;
                truedihadron = truepip+truepim;

                Mh = dihadron.M();
                Mh_true = truedihadron.M();
                
                z = kin.z(init_target,dihadron,q);
                z_true = kin.z(init_target,truedihadron,trueq);

                pT = kin.Pt(q,dihadron,init_target);
                pT_true = kin.Pt(trueq,truedihadron,init_target);

                phih = kin.phi_h(q,init_electron,dihadron);
                phih_true = kin.phi_h(trueq,init_electron,truedihadron);

                Mx = (init_electron+init_target-electron-dihadron).M();
                Mx_true = (init_electron+init_target-trueelectron-truedihadron).M();

                xF = kin.xF(q,dihadron,init_target,W);
                xF_true = kin.xF(trueq,truedihadron,init_target,trueW);

                xF1 = kin.xF(q,pip,init_target,W);
                xF1_true = kin.xF(trueq,truedihadron,init_target,trueW);
                
                xF2 = kin.xF(q,pim,init_target,W);
                xF2_true = kin.xF(trueq,truedihadron,init_target,trueW);
                
                eps=(1-y-pow(y*gamma,2)/4)/(1-y+pow(y,2)/2+pow(y*gamma,2)/4);
                gamma = 2*0.938272*x/sqrt(Q2);

                th = kin.com_th(pip,pim);
                th_true = kin.com_th(truepip,truepim);

                cth = cos(th);
                cth_true = cos(th);

                if(!(electron.E()>0&&pip.P()>1.25&&pim.P()>1.25&&xF1>0&&xF2>0)){
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


