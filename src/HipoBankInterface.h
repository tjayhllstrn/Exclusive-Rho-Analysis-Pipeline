#ifndef HipoBankInterface_h
#define HipoBankInterface_h
#include "Constants.h"
#include "Structs.h"

class HipoBankInterface{
 public:
  HipoBankInterface();
  virtual ~HipoBankInterface() = default;
  HipoBankInterface(const std::unique_ptr<clas12::clas12reader>&);


  // Using the pindex from the SIDISParticle, attach extra bank info
  bool loadBankData(const std::unique_ptr<clas12::clas12reader>&, part&);
 protected:
  bool importDataToParticle(part &);
  int determineSectorDC(float, float, float);
  void clear();
  
  int _idx_RECCal;
  int _ipindex_RECCal;
  int _ix_RECCal;
  int _iy_RECCal; 
  int _iz_RECCal;
  int _ilu_RECCal;
  int _ilv_RECCal;
  int _ilw_RECCal;
  int _idu_RECCal;
  int _idv_RECCal;
  int _idw_RECCal;
  int _im2u_RECCal;
  int _im2v_RECCal;
  int _im2w_RECCal;
  int _im3u_RECCal;
  int _im3v_RECCal;
  int _im3w_RECCal;
  int _ilayer_RECCal;
  int _isector_RECCal;
  int _itime_RECCal;
  int _ipath_RECCal;
  int _ienergy_RECCal;
    
  int _idx_RECCherenkov;
  int _ipindex_RECCherenkov;
  int _iz_RECCherenkov;
  int _inphe_RECCherenkov;
  int _idetector_RECCherenkov;
    
  float _Ele_PCAL_e, _Ele_ECIN_e, _Ele_ECOUT_e;
  int _sector_Cal[3]={0,0,0};
  float _time_Cal[3]={0,0,0};
  float _path_Cal[3]={0,0,0};
  float _x_Cal[3]={0,0,0};
  float _y_Cal[3]={0,0,0};
  float _z_Cal[3]={0,0,0};
  float _lu_Cal[3]={0,0,0};
  float _lv_Cal[3]={0,0,0};
  float _lw_Cal[3]={0,0,0};
  float _du_Cal[3]={0,0,0};
  float _dv_Cal[3]={0,0,0};
  float _dw_Cal[3]={0,0,0};
  float _m2u_Cal[3]={0,0,0};
  float _m2v_Cal[3]={0,0,0};
  float _m2w_Cal[3]={0,0,0};
  float _m3u_Cal[3]={0,0,0};
  float _m3v_Cal[3]={0,0,0};
  float _m3w_Cal[3]={0,0,0};
  
  int _idx_RECTraj;
  int _ipindex_RECTraj;
  int _ilayer_RECTraj;
  int _idet_RECTraj;
  int _ipath_RECTraj;
  int _ix_RECTraj;
  int _iy_RECTraj;
  int _iz_RECTraj;

  int _sector_DC; 

  int _det_DC[3]={0,0,0};  
  float _path_DC[3]={0,0,0};
  float _x_DC[3]={0,0,0};
  float _y_DC[3]={0,0,0};
  float _z_DC[3]={0,0,0};
  float _edge[3]={0,0,0}; //added pass2 edge variable
  

  float _nphe_ltcc=0;
  float _nphe_htcc=0;
};  

#endif
