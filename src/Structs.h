#ifndef Structs_h
#define Structs_h
struct FS{
  int pid_h1=0;
  int pid_h2=0;
  int num_h1=0;
  int num_h2=0;
};

// RG-C Lookup Table
struct Row {
    int Run;
    std::string Target;
    double Nevents;
    double HWP;
    double TpolSign;
    double Tpol;
};

FS get_FS(int pid_h1, int pid_h2){
  FS fs;
  if(pid_h1==111&&pid_h2==111){ //Pi0Pi0
    fs.pid_h1=22;
    fs.num_h1=4;
    fs.pid_h2=0;
    fs.num_h2=0;
  }
  else if(pid_h1==111&&pid_h2!=111){ // One Pi0
    fs.pid_h1=22;
    fs.num_h1=2;
    fs.pid_h2=pid_h2;
    fs.num_h2=1;
  }
  else if(pid_h1!=111&&pid_h2==111){ // One Pi0
    fs.pid_h1=pid_h1;
    fs.num_h1=1;
    fs.pid_h2=22;
    fs.num_h2=2;
  }
  else if(pid_h1==pid_h2 && pid_h1!=0){ // PiPi
    fs.pid_h1=pid_h1;
    fs.num_h1=2;
    fs.pid_h2=0;
    fs.num_h2=0;
  }
  else{
    fs.pid_h1=pid_h1;
    fs.num_h1=1;
    fs.pid_h2=pid_h2;
    fs.num_h2=1;
  }
  
  if(pid_h1==0){
    fs.pid_h1=0;
    fs.num_h1=0;
  }
  if(pid_h2==0){
    fs.pid_h2=0;
    fs.num_h2=0;
  }

  return fs;
} 

struct EVENT_INFO {
    int A = 0;
    int evnum = 0;
    int fgId=0;
    int run=0;
    float torus=0;
    double Pol=0;
    int hel=0;
    int uID=0;
    int hwp=0;
    int tSign=0;
    // Targets
    // 0 --> NH3
    // 1 --> ND3
    // 2 --> C
    // 3 --> CH2
    // 4 --> CD2
    // 5 --> Empty
    int target=-999;
    double tPol = 0;
};

struct EVENT {
    
    double x=-999;
    double y=-999;
    double Q2=-999;
    double nu=-999;
    double W=-999;
    double gamma=-999;
    double eps=-999;
    double depolA=-999;
    double depolB=-999;
    double depolC=-999;
    double depolV=-999;
    double depolW=-999;
    double M1=-999;
    double M2=-999;
    double M12=-999;
    double Mh=-999;
    double phi_h=-999;
    double phi_R0=-999;
    double phi_R1=-999;
    double th=-999;
    double z1=-999;
    double z2=-999;
    double z=-999;
    double xF1=-999;
    double xF2=-999;
    double xF=-999;
    double Mx=-999;
    double phi_h1=-999;
    double phi_h2=-999;
    double delta_phi_h=-999;
    double pT1=-999;
    double pT2=-999;
    double pTtot=-999;
    double P1=-999;
    double P2=-999;
    double Ptot=-999;
    
    double truex=-999;
    double truey=-999;
    double trueQ2=-999;
    double truenu=-999;
    double trueW=-999;
    double truegamma=-999;
    double trueeps=-999;
    double truedepolA=-999;
    double truedepolB=-999;
    double truedepolC=-999;
    double truedepolV=-999;
    double truedepolW=-999;
    double trueM1 = -999;
    double trueM2 = -999;
    double trueM12 = -999;
    double trueMh = -999;
    double truephi_h = -999;
    double truephi_R0 = -999;
    double truephi_R1 = -999;
    double trueth = -999;
    double truez1 = -999;
    double truez2 = -999;
    double truez = -999;
    double truexF1 = -999;
    double truexF2 = -999;
    double truexF = -999;
    double trueMx = -999;
    double truephi_h1 = -999;
    double truephi_h2 = -999;
    double truedelta_phi_h = -999;
    double truepT1 = -999;
    double truepT2 = -999;
    double truepTtot = -999;
    double trueP1 = -999;
    double trueP2 = -999;
    double truePtot = -999;
    
    int truepid_1=-999;
    int truepid_2=-999;
    int truepid_11=-999;
    int truepid_12=-999;
    int truepid_21=-999;
    int truepid_22=-999;
    int trueparentpid_1=-999;
    int trueparentpid_2=-999;
    int trueparentid_1=-999;
    int trueparentid_2=-999;
    int trueparentparentpid_1=-999;
    int trueparentparentpid_2=-999;
    int trueparentparentid_1=-999;
    int trueparentparentid_2=-999;
    
    int MCmatch=0;
    int MC_2h_match=0;
    int isGoodEventWithoutML=0;
    int is_CFR_1=0;
    int is_CFR_2=0;
    
    double p_11=-999;
    double p_12=-999;
    double p_21=-999;
    double p_22=-999;
    
    int i = -999; // dihadron indices
    int ii = -999;
    int j = -999;
    int jj = -999;
};

struct part{
  // Reconstructed Info
  int pindex=-999;
  int status=0;
  double px=-999;
  double py=-999;
  double pz=-999;
  double pt=-999;
  double p=-999;
  double E=-999;
  double m=-999;
  double theta=-999;
  double eta=-999;
  double phi=-999;
  int pid=-999;
  double vx=-999;
  double vy=-999;
  double vz=-999;
  double chi2=-999;
  double beta=-999;
  int is_scattered_electron=0;
    
  // MC Lund Info
  double truepx=-999;
  double truepy=-999;
  double truepz=-999;
  double truep=-999;
  double truept=-999;
  double trueE=-999;
  double truem=-999;
  double truetheta=-999;
  double trueeta=-999;
  double truephi=-999;
  double truevx=-999;
  double truevy=-999;
  double truevz=-999;
  int is_CFR=-999;
  int truepid=-999;
  int trueparentid=-999;
  int trueparentpid=-999;
  int trueparentparentid=-999;
  int trueparentparentpid=-999;
  // Calorimeter Info
  int    pcal_sector=-999;
  double pcal_e=-999;
  double pcal_x=-999;
  double pcal_y=-999;
  double pcal_z=-999;
  double pcal_lu=-999;
  double pcal_lv=-999;
  double pcal_lw=-999;
  double pcal_m2u=-999;
  double pcal_m2v=-999;
  double pcal_m2w=-999;
    
  int    ecin_sector=-999;
  double ecin_e=-999;
  double ecin_x=-999;
  double ecin_y=-999;
  double ecin_z=-999;
  double ecin_lu=-999;
  double ecin_lv=-999;
  double ecin_lw=-999;
  double ecin_m2u=-999;
  double ecin_m2v=-999;
  double ecin_m2w=-999;

  int    ecout_sector=-999;
  double ecout_e=-999;
  double ecout_x=-999;
  double ecout_y=-999;
  double ecout_z=-999;
  double ecout_lu=-999;
  double ecout_lv=-999;
  double ecout_lw=-999;
  double ecout_m2u=-999;
  double ecout_m2v=-999;
  double ecout_m2w=-999;
    
  // Drift Chamber Info
  int sector = -999;
  double traj_x1=-999;
  double traj_y1=-999;
  double traj_z1=-999;
  double traj_x2=-999;
  double traj_y2=-999;
  double traj_z2=-999;
  double traj_x3=-999;
  double traj_y3=-999;
  double traj_z3=-999;
  double edge_1=-999; //added pass2 edge variable
  double edge_2=-999; //added pass2 edge variable
  double edge_3=-999; //added pass2 edge variable
  
  // Cherenkov Info
  double nphe_ltcc=-999;
  double nphe_htcc=-999;
};
#endif
