#include "CLAS12Ana.h"
#include "Constants.h"
#include "Structs.h"
CLAS12Ana::CLAS12Ana(){}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
CLAS12Ana::CLAS12Ana(const std::unique_ptr<clas12::clas12reader>& _c12){
     hipoBankInterface = HipoBankInterface(_c12);
}
// ***********************************************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::set_run_config(const std::unique_ptr<clas12::clas12reader>& _c12){
    _idx_RUNconfig = _c12->addBank("RUN::config"); //finds the index of the RUN::config bank in the hipo file
    _irun = _c12->getBankOrder(_idx_RUNconfig,"run"); //finds the index of the run number in that config bank
    _ievnum = _c12->getBankOrder(_idx_RUNconfig,"event"); //finds index of event number in the bank
    _itorus = _c12->getBankOrder(_idx_RUNconfig,"torus"); //finds index of torus in/outbending in the bank
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//fills event_info with relevant properties of that event
void CLAS12Ana::get_event_info(const std::unique_ptr<clas12::clas12reader>& _c12, EVENT_INFO &event_info){
    
    auto event = _c12->event(); // gets the 
    
    int prev_run = event_info.run; //set a dummy variable for previous run number before changing
    event_info.run = _c12->getBank(_idx_RUNconfig)->getInt(_irun,0); //run number  
    event_info.torus = _c12->getBank(_idx_RUNconfig)->getFloat(_itorus,0); //torus info 
    
    //MC runs correspond to the rungroup, for example RGA_inbending is -11, RGA_outbending is 11. RGB inbending is -22 RGB outbending is 22
    if(event_info.run==11){ 
        event_info.run *= event_info.torus; //makes the run name in event_info negative or positive based on the torus polarity
    }
    
    event_info.evnum = _c12->getBank(_idx_RUNconfig)->getInt(_ievnum,0);//gets the event number 
    event_info.hel   = event->getHelicity(); //gets the helicity
     
    //flips the helicity if within specific ranges defined in constants.h 
    if(runHelicityFlip(event_info.run))
        event_info.hel*=-1;

    //determines the event polarization based on certain conditions in constants.h 
    event_info.Pol = runPolarization(event_info.run);

    //if this is an RGC run there are special things you must do
    if (event_info.run >= 16082 && event_info.run <= 17738 && event_info.run!=prev_run){ // RGC, fill this info whenever a new run appears
        event_info.tPol = get_RGC_Tpol(event_info.run);
        event_info.hwp = get_RGC_HWP(event_info.run);
        event_info.tSign = get_RGC_TpolSign(event_info.run);
        event_info.target = get_RGC_target(event_info.run);
    }
    //calculate s
     _electron_beam_energy = runBeamEnergy(event_info.run);
     init_electron.SetPxPyPzE(0,0,sqrt(_electron_beam_energy*_electron_beam_energy-Me*Me),_electron_beam_energy);
     target.SetPxPyPzE(0,0,0,Mp);
     s = init_electron.M2()+target.M2()+2*target.M()*_electron_beam_energy; //but not saved anywhere????????????????????
    
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::set_beams(TLorentzVector eIn, TLorentzVector pIn){
    init_electron = eIn;
    target = pIn;
    _electron_beam_energy = eIn.E();
    s = init_electron.M2()+target.M2()+2*target.M()*_electron_beam_energy;
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<part> CLAS12Ana::load_reco_particles(const std::unique_ptr<clas12::clas12reader>& _c12){
    
    std::vector<part> vec_particles;
    
    // Loop over reconstructed particles
    // -------------------------------------------------------
    auto particles=_c12->getDetParticles();
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Create new part struct
      part partstruct;
      // Extract each particle from event one-at-a-time
      // -------------------------------------------------------
      auto particle = particles.at(idx);
      partstruct.pid = particle->getPid(); 

      partstruct.chi2 = particle->getChi2Pid();
      partstruct.theta = particle->getTheta();
      partstruct.eta = _kin.eta(partstruct.theta); //pseudorapidity calculation
      partstruct.phi = particle->getPhi();
      partstruct.p = particle->getP();
      partstruct.px = _kin.Px(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.py = _kin.Py(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.pz = _kin.Pz(partstruct.p,partstruct.theta,partstruct.phi);
      partstruct.pt = _kin.Pt(partstruct.px,partstruct.py); //transverse momentum
      if(partstruct.pid!=22) //if the particle is not a photon, then get the mass
        partstruct.m = particle->getPdgMass();
      else
        partstruct.m = 0;
      partstruct.beta = particle->getBeta();
      partstruct.pindex = particle->getIndex();
      partstruct.vx = particle->par()->getVx();
      partstruct.vy = particle->par()->getVy();
      partstruct.vz = particle->par()->getVz();
      partstruct.status = particle->getStatus();
      partstruct.E = _kin.E(partstruct.m,partstruct.p);
    
      // Ensure hadrons are not in CD - why??????????????????????????????
      if (partstruct.pid == 2212 || partstruct.pid == -2212 ||
        partstruct.pid == 2112 ||
        partstruct.pid == -321 || partstruct.pid == -211 ||
        partstruct.pid == 211 || partstruct.pid == 321) {
          if(partstruct.status>=4000 && partstruct.status<5000)
              continue;
      }
    
      hipoBankInterface.loadBankData(_c12,partstruct);  
      vec_particles.push_back(partstruct);

    }
    
    return vec_particles;
}


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<part> CLAS12Ana::load_mc_particles(const std::unique_ptr<clas12::clas12reader>& _c12){
    
    
    std::vector<part> vec_mcparticles;
    
      
    // Loop over all Monte Carlo particles
    // -------------------------------------
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows(); idx++){
      part partstruct;
      if(mcparticles->getType(idx)!=1) // Reject non-final state
        {continue;} 
      partstruct.pid = mcparticles->getPid(idx);
      partstruct.truepid = mcparticles->getPid(idx);
      partstruct.truepx = mcparticles->getPx(idx);
      partstruct.truepy = mcparticles->getPy(idx);
      partstruct.truepz = mcparticles->getPz(idx);
      partstruct.truem = mcparticles->getMass(idx);
      partstruct.truept = _kin.Pt(partstruct.truepx,partstruct.truepy);
      partstruct.truep  = _kin.P(partstruct.truepx,partstruct.truepy,partstruct.truepz);
      partstruct.trueE  = _kin.E(partstruct.truem,partstruct.truep);

      partstruct.truetheta = _kin.th(partstruct.truept,partstruct.truepz);
      partstruct.trueeta = _kin.eta(partstruct.truetheta);
      partstruct.truephi   = _kin.phi(partstruct.truepx,partstruct.truepy);

      partstruct.truevx = mcparticles->getVx(idx);
      partstruct.truevy = mcparticles->getVy(idx);
      partstruct.truevz = mcparticles->getVz(idx);

      partstruct.trueparentid = mcparticles->getParent(idx)-1;
      partstruct.trueparentpid = mcparticles->getPid(partstruct.trueparentid);
      partstruct.trueparentparentid = mcparticles->getParent(partstruct.trueparentid)-1;
      if(partstruct.trueparentparentid==-1){
          partstruct.trueparentparentpid = -999;
      }else{
          partstruct.trueparentparentpid = mcparticles->getPid(partstruct.trueparentparentid);
      }
      // for loop over the idxs until we find if this particle came from CFR
      int parent_idx = mcparticles->getParent(idx)-1;
      int parent_pid = 0;
      while(parent_idx>=0){
          parent_pid = mcparticles->getPid(parent_idx);
          if(parent_pid==0){break;}
          if(6-abs(parent_pid)>=0){
              partstruct.is_CFR=1;
              break;
          }
          parent_idx = mcparticles->getParent(parent_idx)-1;
      }
      if(partstruct.is_CFR!=1) partstruct.is_CFR=0;
      if(partstruct.truepid==11 && partstruct.trueparentid==0){ // scattered electro
        partstruct.is_scattered_electron=1;
      }
      
      // Add particle to list
      vec_mcparticles.push_back(partstruct);
    }
    
    return vec_mcparticles;
    
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CLAS12Ana::reco_event_contains_final_state(std::vector<part> vec_particles, FS fs){
    int num_e=0;
    int num_h1=0;
    int num_h2=0;
    for(part particle : vec_particles){
      if(particle.pid==11 && particle.is_scattered_electron==1) num_e++;
      else if(particle.pid==fs.pid_h1) num_h1++;
      else if(particle.pid==fs.pid_h2) num_h2++;
    }
    if(num_e<1 || num_h1 < fs.num_h1 || num_h2 < fs.num_h2) {
      return false;
    }
    
    return true;

}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CLAS12Ana::reco_event_contains_scattered_electron(std::vector<part> vec_particles){
    for(part particle : vec_particles){
      if(particle.pid==11 && particle.is_scattered_electron==1) return true;
    }
    return false;
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CLAS12Ana::find_reco_scattered_electron(std::vector<part>& vec_particles){
    // Code for determine the scattered electron from REC::Particle
    // --> Find pid==11 particle with largest energy
    // -->   If no electron is found, skip
    // --> Check if the status of the maximum energy electron is in FD
    // -->   Skip if not (i.e. always skip events if the max energy electron was not in FD)
    // --> Set that particle as the scattered electron
    int idx_e=-1;
    double max_energy = -1; 
    for (int i = 0; i < vec_particles.size(); i++) {
      part partstruct = vec_particles[i];
      // check if the particle is an electron
      if (partstruct.pid == 11) {
        // compare energy with the current maximum and update if necessary
        if (partstruct.E > max_energy) {
          max_energy = partstruct.E;
          idx_e=i;
         }
       }
    }
       
    
    if(idx_e==-1) return -1;
    
    
    // Toss events where the REC::Particle scattered electron was not in the FD
    
    if((vec_particles[idx_e].status <= -3000 || vec_particles[idx_e].status > -2000))
            return -1;
    
    return idx_e;
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CLAS12Ana::find_mc_scattered_electron(std::vector<part>& vec_mcparticles){
    for (int i = 0 ; i < vec_mcparticles.size(); i++){
        part partstruct = vec_mcparticles[i];
        if(partstruct.is_scattered_electron)
            return i;
        
    }
    return -1;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::fill_reco_event_variables(EVENT &event, std::vector<part> parts){
    
    part scattered_electron;
    for(auto p: parts){
        if(p.is_scattered_electron){
            scattered_electron = p;
            break;
        }
    }
    

    // If none was found, print an error message
    if(scattered_electron.is_scattered_electron != 1){
        cout << "ERROR: No rec::scattered electron found..." << endl;
    }
    
    event.Q2=_kin.Q2(_electron_beam_energy,scattered_electron.E,
               _kin.cth(scattered_electron.px,scattered_electron.py,scattered_electron.pz));
    event.y=_kin.y(_electron_beam_energy,scattered_electron.E);     
    event.nu=_kin.nu(_electron_beam_energy,scattered_electron.E);
    event.W=_kin.W(event.Q2,Mp,event.nu);
    event.x=_kin.x(event.Q2,s,event.y);
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::fill_mc_event_variables(EVENT &event, std::vector<part> parts){
    
    part scattered_electron;
    for(auto p: parts){
        if(p.is_scattered_electron){
            scattered_electron = p;
            break;
        }
    }
    
    // If none was found, print an error message
    if(scattered_electron.is_scattered_electron != 1){
        cout << "ERROR: No mc::scattered electron found..." << endl;
    }
    
    event.trueQ2=_kin.Q2(_electron_beam_energy,scattered_electron.trueE,
               _kin.cth(scattered_electron.truepx,scattered_electron.truepy,scattered_electron.truepz));
    event.truey=_kin.y(_electron_beam_energy,scattered_electron.trueE);     
    event.truenu=_kin.nu(_electron_beam_energy,scattered_electron.trueE);
    event.trueW=_kin.W(event.trueQ2,Mp,event.truenu);
    event.truex=_kin.x(event.trueQ2,s,event.truey);

}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::clear_dihadron_variables(EVENT &event){
    event.i = -1;
    event.ii = -1;
    event.j = -1;
    event.jj = -1;
    
    event.MCmatch = 0;
    event.MC_2h_match = 0;
    
    event.gamma = -999;
    event.eps = -999;
    event.depolA = -999;
    event.depolB = -999;
    event.depolC = -999;
    event.depolV = -999;
    event.depolW = -999;
    event.Mh = -999;
    event.phi_h = -999;
    event.pT1 = -999;
    event.pT2 = -999;
    event.pTtot = -999;
    event.phi_R0 = -999;
    event.phi_R1 = -999;
    event.th = -999;
    event.xF1 = -999;
    event.xF2 = -999;
    event.xF = -999;
    event.z1 = -999;
    event.z2 = -999;
    event.z = -999;
    event.Mx = -999;

    event.trueMh = -999;
    event.truephi_h = -999;
    event.truepT1 = -999;
    event.truepT2 = -999;
    event.truepTtot = -999;
    event.truephi_R0 = -999;
    event.truephi_R1 = -999;
    event.trueth = -999;
    event.truexF1 = -999;
    event.truexF2 = -999;
    event.truexF = -999;
    event.truez1 = -999;
    event.truez2 = -999;
    event.truez = -999;
    event.trueMx = -999;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::fill_mc_reco_dihadron_variables(EVENT &event, TLorentzVector q, TLorentzVector trueq, TLorentzVector electron , TLorentzVector trueelectron, std::vector<part> vec_particles, std::vector<int> dihadron_idx, int pid_h1, int pid_h2){
    
    clear_dihadron_variables(event);
    
    int i=0;
    int ii=0;
    int j=0;
    int jj=0;
    if(pid_h1==111){
        i=dihadron_idx.at(0);
        ii=dihadron_idx.at(1);
    }else{
        i=dihadron_idx.at(0);
    }
    if(pid_h2==111&&pid_h1!=111){
        j=dihadron_idx.at(1);
        jj=dihadron_idx.at(2);
    }else if(pid_h2==111&&pid_h1==111){
        j=dihadron_idx.at(2);
        jj=dihadron_idx.at(3);
    }else if(pid_h1==111){
        j=dihadron_idx.at(2);
    }else{
        j=dihadron_idx.at(1);
    }

    event.i=i;
    event.ii=ii;
    event.j=j;
    event.jj=jj;
    
    TLorentzVector h1;
    TLorentzVector trueh1;
    TLorentzVector h2;
    TLorentzVector trueh2;
    TLorentzVector dihadron;
    TLorentzVector truedihadron;

    if (pid_h1 == 111) {
        h1.SetPxPyPzE(vec_particles[i].px + vec_particles[ii].px,
                       vec_particles[i].py + vec_particles[ii].py,
                       vec_particles[i].pz + vec_particles[ii].pz,
                       vec_particles[i].E + vec_particles[ii].E);
        trueh1.SetPxPyPzE(vec_particles[i].truepx + vec_particles[ii].truepx,
                           vec_particles[i].truepy + vec_particles[ii].truepy,
                           vec_particles[i].truepz + vec_particles[ii].truepz,
                           vec_particles[i].trueE + vec_particles[ii].trueE);
    } else {
        h1.SetPxPyPzE(vec_particles[i].px, vec_particles[i].py, vec_particles[i].pz, vec_particles[i].E);
        trueh1.SetPxPyPzE(vec_particles[i].truepx, vec_particles[i].truepy,
                           vec_particles[i].truepz, vec_particles[i].trueE);
    }

    if (pid_h2 == 111) {
        h2.SetPxPyPzE(vec_particles[j].px + vec_particles[jj].px,
                       vec_particles[j].py + vec_particles[jj].py,
                       vec_particles[j].pz + vec_particles[jj].pz,
                       vec_particles[j].E + vec_particles[jj].E);
        trueh2.SetPxPyPzE(vec_particles[j].truepx + vec_particles[jj].truepx,
                           vec_particles[j].truepy + vec_particles[jj].truepy,
                           vec_particles[j].truepz + vec_particles[jj].truepz,
                           vec_particles[j].trueE + vec_particles[jj].trueE);
    } else {
        h2.SetPxPyPzE(vec_particles[j].px, vec_particles[j].py, vec_particles[j].pz, vec_particles[j].E);
        trueh2.SetPxPyPzE(vec_particles[j].truepx, vec_particles[j].truepy,
                           vec_particles[j].truepz, vec_particles[j].trueE);
    }
    
    if(pid_h1==pid_h2){
        double z1 = _kin.z(target,h1,q);
        double z2 = _kin.z(target,h2,q);
        if(z1<z2){
            TLorentzVector temp = h1;
            h2=h1;
            h1=temp;
            TLorentzVector truetemp = h1;
            trueh2=trueh1;
            trueh1=truetemp;
        }
    }

    // Build the dihadron
    dihadron = h1+h2;
    truedihadron = trueh1+trueh2;
    // fill results

    event.gamma = 2*0.938272*event.x/sqrt(event.Q2);
    event.eps   = (1-event.y-pow(event.y*event.gamma,2)/4)/(1-event.y+pow(event.y,2)/2+pow(event.y*event.gamma,2)/4);
    event.depolA = pow(event.y,2)/(2*(1-event.eps));
    event.depolB = event.depolA*event.eps;
    event.depolC = event.depolA*sqrt(1-event.eps*event.eps);
    event.depolV = event.depolA*sqrt(2*event.eps*(1+event.eps));
    event.depolW = event.depolA*sqrt(2*event.eps*(1-event.eps));
    event.Mh = dihadron.M();
    event.phi_h = _kin.phi_h(q,init_electron,h1,h2);
    event.pT1 = _kin.Pt(q,h1,target);
    event.pT2 = _kin.Pt(q,h2,target);
    event.pTtot = _kin.Pt(q,dihadron,target);
    event.phi_R0 = _kin.phi_R(q,init_electron,h1,h2,0);
    event.phi_R1 = _kin.phi_R(q,init_electron,h1,h2,1);
    event.th     = _kin.com_th(h1,h2);
    event.xF1 = _kin.xF(q,h1,target,event.W);
    event.xF2 = _kin.xF(q,h2,target,event.W);
    event.xF     = _kin.xF(q,dihadron,target,event.W);
    event.z1 = _kin.z(target,h1,q);
    event.z2 = _kin.z(target,h2,q);
    event.z = event.z1+event.z2;
    event.Mx = (init_electron+target-electron-dihadron).M();

    event.trueMh = truedihadron.M();
    event.truephi_h = _kin.phi_h(trueq,init_electron,trueh1,trueh2);
    event.truepT1 = _kin.Pt(trueq,trueh1,target);
    event.truepT2 = _kin.Pt(trueq,trueh2,target);
    event.truepTtot = _kin.Pt(trueq,truedihadron,target);
    event.truephi_R0 = _kin.phi_R(trueq,init_electron,trueh1,trueh2,0);
    event.truephi_R1 = _kin.phi_R(trueq,init_electron,trueh1,trueh2,1);
    event.trueth     = _kin.com_th(trueh1,trueh2);
    event.truexF1 = _kin.xF(trueq,trueh1,target,event.trueW);
    event.truexF2 = _kin.xF(trueq,trueh2,target,event.trueW);
    event.truexF     = _kin.xF(trueq,truedihadron,target,event.trueW);
    event.truez1 = _kin.z(target,trueh1,trueq);
    event.truez2 = _kin.z(target,trueh2,trueq);
    event.truez = event.truez1+event.truez2;
    event.trueMx = (init_electron+target-trueelectron-truedihadron).M();
    
    if(trueelectron.E()>0&&trueh1.E()>0&&trueh2.E()>0) event.MCmatch=1;
    else event.MCmatch = 0;
    
    if(trueh1.E()>0&&trueh2.E()>0) event.MC_2h_match=1;
    else event.MC_2h_match = 0;
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::fill_reco_dihadron_variables(EVENT &event, TLorentzVector q, TLorentzVector electron ,  std::vector<part> vec_particles, std::vector<int> dihadron_idx, int pid_h1, int pid_h2){
    
    clear_dihadron_variables(event);
    
    int i=0;
    int ii=0;
    int j=0;
    int jj=0;
    if(pid_h1==111){
        i=dihadron_idx.at(0);
        ii=dihadron_idx.at(1);
    }else{
        i=dihadron_idx.at(0);
    }
    if(pid_h2==111&&pid_h1!=111){
        j=dihadron_idx.at(1);
        jj=dihadron_idx.at(2);
    }else if(pid_h2==111&&pid_h1==111){
        j=dihadron_idx.at(2);
        jj=dihadron_idx.at(3);
    }else if(pid_h1==111){
        j=dihadron_idx.at(2);
    }else{
        j=dihadron_idx.at(1);
    }

    event.i=i;
    event.ii=ii;
    event.j=j;
    event.jj=jj;
    
    TLorentzVector h1;
    TLorentzVector h2;
    TLorentzVector dihadron;

    if (pid_h1 == 111) {
        h1.SetPxPyPzE(vec_particles[i].px + vec_particles[ii].px,
                       vec_particles[i].py + vec_particles[ii].py,
                       vec_particles[i].pz + vec_particles[ii].pz,
                       vec_particles[i].E + vec_particles[ii].E);
    } else {
        h1.SetPxPyPzE(vec_particles[i].px, vec_particles[i].py, vec_particles[i].pz, vec_particles[i].E);
    }

    if (pid_h2 == 111) {
        h2.SetPxPyPzE(vec_particles[j].px + vec_particles[jj].px,
                       vec_particles[j].py + vec_particles[jj].py,
                       vec_particles[j].pz + vec_particles[jj].pz,
                       vec_particles[j].E + vec_particles[jj].E);
    } else {
        h2.SetPxPyPzE(vec_particles[j].px, vec_particles[j].py, vec_particles[j].pz, vec_particles[j].E);
    }
    
    if(pid_h1==pid_h2){
        double z1 = _kin.z(target,h1,q);
        double z2 = _kin.z(target,h2,q);
        if(z1<z2){
            TLorentzVector temp = h1;
            h2=h1;
            h1=temp;
        }
    }

    // Build the dihadron
    dihadron = h1+h2;
    // fill results

    event.Mh = dihadron.M();
    event.phi_h = _kin.phi_h(q,init_electron,h1,h2);
    event.pT1 = _kin.Pt(q,h1,target);
    event.pT2 = _kin.Pt(q,h2,target);
    event.pTtot = _kin.Pt(q,dihadron,target);
    event.phi_R0 = _kin.phi_R(q,init_electron,h1,h2,0);
    event.phi_R1 = _kin.phi_R(q,init_electron,h1,h2,1);
    event.th     = _kin.com_th(h1,h2);
    event.xF1 = _kin.xF(q,h1,target,event.W);
    event.xF2 = _kin.xF(q,h2,target,event.W);
    event.xF     = _kin.xF(q,dihadron,target,event.W);
    event.z1 = _kin.z(target,h1,q);
    event.z2 = _kin.z(target,h2,q);
    event.z = event.z1+event.z2;
    event.Mx = (init_electron+target-electron-dihadron).M();
    
    
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::match_mc_to_reco(std::vector<part>& vec_particles,
                                      std::vector<part>& vec_mcparticles){
    
    for (int i=0; i < vec_particles.size(); i++){
      for (int j=0; j < vec_mcparticles.size(); j++){
        float dth = abs(vec_particles[i].theta - vec_mcparticles[j].truetheta)*180/PI;
        float dphi = abs(vec_particles[i].phi - vec_mcparticles[j].truephi)*180/PI;
        float dE = abs(vec_particles[i].E - vec_mcparticles[j].trueE);
        
        if (dth<2 && (dphi<4 || abs(dphi-2*PI)<4) && dE<1){
	  // Perform Pairing
	  vec_particles[i].truepx = vec_mcparticles[j].truepx;
	  vec_particles[i].truepy = vec_mcparticles[j].truepy;
	  vec_particles[i].truepz = vec_mcparticles[j].truepz;
	  vec_particles[i].truep = vec_mcparticles[j].truep;
	  vec_particles[i].truept = vec_mcparticles[j].truept;
	  vec_particles[i].trueE = vec_mcparticles[j].trueE;
	  vec_particles[i].truem = vec_mcparticles[j].truem;
	  vec_particles[i].truetheta = vec_mcparticles[j].truetheta;
	  vec_particles[i].trueeta = vec_mcparticles[j].trueeta;
	  vec_particles[i].truephi = vec_mcparticles[j].truephi;
	  vec_particles[i].truevx = vec_mcparticles[j].truevx;
	  vec_particles[i].truevy = vec_mcparticles[j].truevy;
	  vec_particles[i].truevz = vec_mcparticles[j].truevz;
      vec_particles[i].is_CFR = vec_mcparticles[j].is_CFR;
	  vec_particles[i].truepid = vec_mcparticles[j].truepid;
	  vec_particles[i].trueparentid = vec_mcparticles[j].trueparentid;
	  vec_particles[i].trueparentpid = vec_mcparticles[j].trueparentpid;
      vec_particles[i].trueparentparentid = vec_mcparticles[j].trueparentparentid;
	  vec_particles[i].trueparentparentpid = vec_mcparticles[j].trueparentparentpid;
            
      vec_mcparticles[j].px = vec_particles[i].px;
      vec_mcparticles[j].py = vec_particles[i].py;
      vec_mcparticles[j].pz = vec_particles[i].pz;
      vec_mcparticles[j].p = vec_particles[i].p;
      vec_mcparticles[j].pt = vec_particles[i].pt;
      vec_mcparticles[j].E = vec_particles[i].E;
      vec_mcparticles[j].m = vec_particles[i].m;
      vec_mcparticles[j].theta = vec_particles[i].theta;
      vec_mcparticles[j].eta = vec_particles[i].eta;
      vec_mcparticles[j].phi = vec_particles[i].phi;
      vec_mcparticles[j].vx = vec_particles[i].vx;
      vec_mcparticles[j].vy = vec_particles[i].vy;
      vec_mcparticles[j].vz = vec_particles[i].vz;
      vec_mcparticles[j].is_CFR = vec_particles[i].is_CFR;
            
    vec_mcparticles[j].pcal_sector = vec_particles[i].pcal_sector;
    vec_mcparticles[j].pcal_e = vec_particles[i].pcal_e;
    vec_mcparticles[j].pcal_x = vec_particles[i].pcal_x;
    vec_mcparticles[j].pcal_y = vec_particles[i].pcal_y;
    vec_mcparticles[j].pcal_z = vec_particles[i].pcal_z;
    vec_mcparticles[j].pcal_lu = vec_particles[i].pcal_lu;
    vec_mcparticles[j].pcal_lv = vec_particles[i].pcal_lv;
    vec_mcparticles[j].pcal_lw = vec_particles[i].pcal_lw;
    vec_mcparticles[j].pcal_m2u = vec_particles[i].pcal_m2u;
    vec_mcparticles[j].pcal_m2v = vec_particles[i].pcal_m2v;
    vec_mcparticles[j].pcal_m2w = vec_particles[i].pcal_m2w;

    vec_mcparticles[j].ecin_sector = vec_particles[i].ecin_sector;
    vec_mcparticles[j].ecin_e = vec_particles[i].ecin_e;
    vec_mcparticles[j].ecin_x = vec_particles[i].ecin_x;
    vec_mcparticles[j].ecin_y = vec_particles[i].ecin_y;
    vec_mcparticles[j].ecin_z = vec_particles[i].ecin_z;
    vec_mcparticles[j].ecin_lu = vec_particles[i].ecin_lu;
    vec_mcparticles[j].ecin_lv = vec_particles[i].ecin_lv;
    vec_mcparticles[j].ecin_lw = vec_particles[i].ecin_lw;
    vec_mcparticles[j].ecin_m2u = vec_particles[i].ecin_m2u;
    vec_mcparticles[j].ecin_m2v = vec_particles[i].ecin_m2v;
    vec_mcparticles[j].ecin_m2w = vec_particles[i].ecin_m2w;

    vec_mcparticles[j].ecout_sector = vec_particles[i].ecout_sector;
    vec_mcparticles[j].ecout_e = vec_particles[i].ecout_e;
    vec_mcparticles[j].ecout_x = vec_particles[i].ecout_x;
    vec_mcparticles[j].ecout_y = vec_particles[i].ecout_y;
    vec_mcparticles[j].ecout_z = vec_particles[i].ecout_z;
    vec_mcparticles[j].ecout_lu = vec_particles[i].ecout_lu;
    vec_mcparticles[j].ecout_lv = vec_particles[i].ecout_lv;
    vec_mcparticles[j].ecout_lw = vec_particles[i].ecout_lw;
    vec_mcparticles[j].ecout_m2u = vec_particles[i].ecout_m2u;
    vec_mcparticles[j].ecout_m2v = vec_particles[i].ecout_m2v;
    vec_mcparticles[j].ecout_m2w = vec_particles[i].ecout_m2w;

    vec_mcparticles[j].sector = vec_particles[i].sector;
    vec_mcparticles[j].traj_x1 = vec_particles[i].traj_x1;
    vec_mcparticles[j].traj_y1 = vec_particles[i].traj_y1;
    vec_mcparticles[j].traj_z1 = vec_particles[i].traj_z1;
    vec_mcparticles[j].traj_x2 = vec_particles[i].traj_x2;
    vec_mcparticles[j].traj_y2 = vec_particles[i].traj_y2;
    vec_mcparticles[j].traj_z2 = vec_particles[i].traj_z2;
    vec_mcparticles[j].traj_x3 = vec_particles[i].traj_x3;
    vec_mcparticles[j].traj_y3 = vec_particles[i].traj_y3;
    vec_mcparticles[j].traj_z3 = vec_particles[i].traj_z3;

    vec_mcparticles[j].nphe_ltcc = vec_particles[i].nphe_ltcc;
    vec_mcparticles[j].nphe_htcc = vec_particles[i].nphe_htcc;

	  break;
	}
      }
    }
    
    
}


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CLAS12Ana::generate_combinations(std::vector<int>& input, int num, int start_idx, std::vector<int>& curr_combination, std::vector<std::vector<int>>& result) {
    if (num == 0) {
        result.push_back(curr_combination);
        return;
    }
    for (int i = start_idx; i <= input.size() - num; i++) {
        curr_combination.push_back(input[i]);
        generate_combinations(input, num - 1, i + 1, curr_combination, result);
        curr_combination.pop_back();
    }
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> CLAS12Ana::unique_combinations(std::vector<int> input, int num) {
    std::vector<std::vector<int>> result;
    std::vector<int> curr_combination;
    std::sort(input.begin(), input.end());

    generate_combinations(input, num, 0, curr_combination, result);
    return result;
}
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> CLAS12Ana::remove_duplicates(std::vector<std::vector<int>> input) {
    std::vector<std::vector<int>> output;

    // Store the original indices and sort each inner vector based on the values
    std::vector<std::vector<size_t>> indices(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        indices[i].resize(input[i].size());
        std::iota(indices[i].begin(), indices[i].end(), 0);

        std::sort(indices[i].begin(), indices[i].end(),
                  [&](size_t a, size_t b) { return input[i][a] < input[i][b]; });
        std::sort(input[i].begin(), input[i].end());
    }

    // Sort and remove duplicates from the outer vector
    std::sort(input.begin(), input.end());
    input.erase(std::unique(input.begin(), input.end()), input.end());

    // Restore the original order of the inner vectors
    for (size_t i = 0; i < input.size(); ++i) {
        std::vector<int> temp(input[i].size());
        for (size_t j = 0; j < input[i].size(); ++j) {
            temp[indices[i][j]] = input[i][j];
        }
        input[i] = temp;
    }

    return input;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> CLAS12Ana::dihadron_idxs(int pid_h1, int pid_h2, int pid[], int Nmax){
    
    std::vector<int> h1_idxs;
    std::vector<int> h2_idxs;
    std::vector<std::vector<int>> twoh_idxs;
    for(int i = 0; i<Nmax; i++){
        if(pid[i]==pid_h1 || (pid[i]==22 && pid_h1==111)) h1_idxs.push_back(i);
        if(pid[i]==pid_h2 || (pid[i]==22 && pid_h2==111)) h2_idxs.push_back(i);
    }
    

    //Now form all possible dihadron index pairs
    if(pid_h1==pid_h2 && pid_h1!=111){twoh_idxs=unique_combinations(h1_idxs,2); }
    else if(pid_h1==pid_h2 && pid_h1==111){twoh_idxs=unique_combinations(h1_idxs,4);}
    else if(pid_h1!=pid_h2 && pid_h1==111 && pid_h2 != 111){
        for(int i = 0 ; i < h2_idxs.size(); i++){
            for(int j = 0 ; j < h1_idxs.size(); j++){
                for(int k = j+1 ; k < h1_idxs.size(); k++){
                    std::vector<int> dihadron_idx = {h1_idxs.at(j),h1_idxs.at(k),h2_idxs.at(i)}; // 2 photons at start
                    twoh_idxs.push_back(dihadron_idx);
                }
            }
        }
    }
    else if(pid_h1!=pid_h2 && pid_h1!=111 && pid_h2 == 111){
        for(int i = 0 ; i < h1_idxs.size(); i++){
            for(int j = 0 ; j < h2_idxs.size(); j++){
                for(int k = j+1 ; k < h2_idxs.size(); k++){
                    std::vector<int> dihadron_idx = {h1_idxs.at(i),h2_idxs.at(j),h2_idxs.at(k)}; // 2 photons at end
                    twoh_idxs.push_back(dihadron_idx);
                }
            }
        }
    }
    else{
        for(int i = 0 ; i < h1_idxs.size(); i++){
            for(int j = 0 ; j < h2_idxs.size(); j++){
                std::vector<int> dihadron_idx = {h1_idxs.at(i), h2_idxs.at(j)};
                twoh_idxs.push_back(dihadron_idx);
            }
        }
    }

    // Remove any instance of duplicate dihadrons
    twoh_idxs = remove_duplicates(twoh_idxs);
    
    return twoh_idxs;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> CLAS12Ana::dihadron_idxs(int pid_h1, int pid_h2, std::vector<int> pid){

    int* pidArray = pid.data();
    int size = pid.size();
    
    return dihadron_idxs(pid_h1,pid_h2,pidArray,size);
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<std::vector<int>> CLAS12Ana::dihadron_idxs(int pid_h1, int pid_h2, std::vector<part> vec_particles){

    std::vector<int> pid;
    for(auto particle: vec_particles){
        pid.push_back(particle.pid);
    }
    
    return dihadron_idxs(pid_h1,pid_h2,pid);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
