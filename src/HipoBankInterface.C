#include "HipoBankInterface.h"

HipoBankInterface::HipoBankInterface(){}

HipoBankInterface::HipoBankInterface(const std::unique_ptr<clas12::clas12reader>& _c12){
  
  // Add REC::Cal info
  _idx_RECCal = _c12->addBank("REC::Calorimeter");
  _ipindex_RECCal = _c12->getBankOrder(_idx_RECCal,"pindex");
  _ix_RECCal  = _c12->getBankOrder(_idx_RECCal,"x");
  _iy_RECCal  = _c12->getBankOrder(_idx_RECCal,"y");
  _iz_RECCal  = _c12->getBankOrder(_idx_RECCal,"z");
  _ilu_RECCal = _c12->getBankOrder(_idx_RECCal,"lu");
  _ilv_RECCal = _c12->getBankOrder(_idx_RECCal,"lv");
  _ilw_RECCal = _c12->getBankOrder(_idx_RECCal,"lw");
  _im2u_RECCal = _c12->getBankOrder(_idx_RECCal,"m2u");
  _im2v_RECCal = _c12->getBankOrder(_idx_RECCal,"m2v");
  _im2w_RECCal = _c12->getBankOrder(_idx_RECCal,"m2w");
  _ilayer_RECCal = _c12->getBankOrder(_idx_RECCal,"layer");
  _isector_RECCal = _c12->getBankOrder(_idx_RECCal,"sector");
  _itime_RECCal = _c12->getBankOrder(_idx_RECCal,"time");
  _ipath_RECCal = _c12->getBankOrder(_idx_RECCal,"path");
  _ienergy_RECCal = _c12->getBankOrder(_idx_RECCal,"energy");

  // Add REC::Traj info
  _idx_RECTraj = _c12->addBank("REC::Traj");
  _ipindex_RECTraj = _c12->getBankOrder(_idx_RECTraj,"pindex");
  _ilayer_RECTraj = _c12->getBankOrder(_idx_RECTraj,"layer");
  _idet_RECTraj = _c12->getBankOrder(_idx_RECTraj,"detector");
  _ipath_RECTraj = _c12->getBankOrder(_idx_RECTraj,"path");
  _ix_RECTraj = _c12->getBankOrder(_idx_RECTraj,"x");
  _iy_RECTraj = _c12->getBankOrder(_idx_RECTraj,"y");
  _iz_RECTraj = _c12->getBankOrder(_idx_RECTraj,"z");

  // Add REC::Cherenkov info
  _idx_RECCherenkov = _c12->addBank("REC::Cherenkov");
  _ipindex_RECCherenkov = _c12->getBankOrder(_idx_RECCherenkov,"pindex");
  _inphe_RECCherenkov = _c12->getBankOrder(_idx_RECCherenkov,"nphe");
  _iz_RECCherenkov = _c12->getBankOrder(_idx_RECCherenkov,"z");
  _idetector_RECCherenkov = _c12->getBankOrder(_idx_RECCherenkov,"detector");
}





bool HipoBankInterface::loadBankData(const std::unique_ptr<clas12::clas12reader>& _c12 , part &particle){
  clear();

  // Grab necessary particle info
  // -------------------------------------------------------------
  int pindex = particle.pindex;

  // -------------------------------------------------------------
  // Parse the REC::Calorimeter
  // -------------------------------------------------------------

  for(auto i = 0 ; i < _c12->getBank(_idx_RECCal)->getRows() ; i++){
    // Continue loop if the pindex in the calo bank does not match
    if(_c12->getBank(_idx_RECCal)->getInt(_ipindex_RECCal,i)!=pindex)
      continue;
    
    int sectorCal = _c12->getBank(_idx_RECCal)->getInt(_isector_RECCal,i);
    float timeCal = _c12->getBank(_idx_RECCal)->getFloat(_itime_RECCal,i);
    float pathCal = _c12->getBank(_idx_RECCal)->getFloat(_ipath_RECCal,i);
    float lu = _c12->getBank(_idx_RECCal)->getFloat(_ilu_RECCal,i);
    float lv = _c12->getBank(_idx_RECCal)->getFloat(_ilv_RECCal,i);
    float lw = _c12->getBank(_idx_RECCal)->getFloat(_ilw_RECCal,i);
    float m2u = _c12->getBank(_idx_RECCal)->getFloat(_im2u_RECCal,i);
    float m2v = _c12->getBank(_idx_RECCal)->getFloat(_im2v_RECCal,i);
    float m2w = _c12->getBank(_idx_RECCal)->getFloat(_im2w_RECCal,i);
    float x = _c12->getBank(_idx_RECCal)->getFloat(_ix_RECCal,i);
    float y = _c12->getBank(_idx_RECCal)->getFloat(_iy_RECCal,i);
    float z = _c12->getBank(_idx_RECCal)->getFloat(_iz_RECCal,i);
    int layerCal = _c12->getBank(_idx_RECCal)->getInt(_ilayer_RECCal,i);
    int calidx = -1; //Array index for lu, lv, lw

    switch(layerCal){
    case 1: //PCal
      calidx = 0;
      _Ele_PCAL_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;
    case 4: //ECIN
      calidx = 1;
      _Ele_ECIN_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;
    case 7: //ECOUT
      calidx = 2;
      _Ele_ECOUT_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;   
    }
    
    // If there was a pindex attached to one of the calo layers...
    if(calidx!=-1){
      _x_Cal[calidx]=x;
      _y_Cal[calidx]=y;
      _z_Cal[calidx]=z;
    
      _lu_Cal[calidx]=lu;
      _lv_Cal[calidx]=lv;
      _lw_Cal[calidx]=lw;

      _m2u_Cal[calidx]=m2u;
      _m2v_Cal[calidx]=m2v;
      _m2w_Cal[calidx]=m2w;
        
      _sector_Cal[calidx]=sectorCal;
      _time_Cal[calidx]=timeCal;
      _path_Cal[calidx]=pathCal;
    
    }
  }

  // -------------------------------------------------------------
  // Parse the REC::Trajectory
  // -------------------------------------------------------------

  for(auto i = 0 ; i < _c12->getBank(_idx_RECTraj)->getRows() ; i++){
    // Continue loop if the pindex in the traj bank does not match
    if(_c12->getBank(_idx_RECTraj)->getInt(_ipindex_RECTraj,i)!=pindex)
      continue;
    
    if(_c12->getBank(_idx_RECTraj)->getInt(_ilayer_RECTraj,i)==6){
      _det_DC[0] = _c12->getBank(_idx_RECTraj)->getInt(_idet_RECTraj,i);
      _path_DC[0] = _c12->getBank(_idx_RECTraj)->getFloat(_ipath_RECTraj,i);
      _x_DC[0] = _c12->getBank(_idx_RECTraj)->getFloat(_ix_RECTraj,i);
      _y_DC[0] = _c12->getBank(_idx_RECTraj)->getFloat(_iy_RECTraj,i);
      _z_DC[0] = _c12->getBank(_idx_RECTraj)->getFloat(_iz_RECTraj,i);
      _edge[0] = _c12->getBank(_idx_RECTraj)->getFloat("edge", i);
    }else if(_c12->getBank(_idx_RECTraj)->getInt(_ilayer_RECTraj,i)==18){
      _det_DC[1] = _c12->getBank(_idx_RECTraj)->getInt(_idet_RECTraj,i);
      _path_DC[1] = _c12->getBank(_idx_RECTraj)->getFloat(_ipath_RECTraj,i);
      _x_DC[1] = _c12->getBank(_idx_RECTraj)->getFloat(_ix_RECTraj,i);
      _y_DC[1] = _c12->getBank(_idx_RECTraj)->getFloat(_iy_RECTraj,i);
      _z_DC[1] = _c12->getBank(_idx_RECTraj)->getFloat(_iz_RECTraj,i);
      _edge[1] = _c12->getBank(_idx_RECTraj)->getFloat("edge", i);
    }else if(_c12->getBank(_idx_RECTraj)->getInt(_ilayer_RECTraj,i)==36){
      _det_DC[2] = _c12->getBank(_idx_RECTraj)->getInt(_idet_RECTraj,i);
      _path_DC[2] = _c12->getBank(_idx_RECTraj)->getFloat(_ipath_RECTraj,i);
      _x_DC[2] = _c12->getBank(_idx_RECTraj)->getFloat(_ix_RECTraj,i);
      _y_DC[2] = _c12->getBank(_idx_RECTraj)->getFloat(_iy_RECTraj,i);
      _z_DC[2] = _c12->getBank(_idx_RECTraj)->getFloat(_iz_RECTraj,i);
      _edge[2] = _c12->getBank(_idx_RECTraj)->getFloat("edge", i);
    }
  }
  
  // Get the azimuthal sector # from the middle drift chamber
  _sector_DC = determineSectorDC(_x_DC[1], _y_DC[1], _z_DC[1]);
  
  // -------------------------------------------------------------
  // Parse the REC::Cherenkov
  // -------------------------------------------------------------
  for(auto i = 0 ; i < _c12->getBank(_idx_RECCherenkov)->getRows() ; i++){
      if(_c12->getBank(_idx_RECCherenkov)->getInt(_ipindex_RECCherenkov,i)!=pindex)
          continue;
      
      // See https://clasweb.jlab.org/wiki/index.php/CLAS12_DSTs#Special_Banks for detector ints
      if(_c12->getBank(_idx_RECCherenkov)->getByte(_idetector_RECCherenkov,i)==15){// htcc == 15
          _nphe_htcc = _c12->getBank(_idx_RECCherenkov)->getFloat(_inphe_RECCherenkov,i);
      } else { // ltcc == 16
          _nphe_ltcc = _c12->getBank(_idx_RECCherenkov)->getFloat(_inphe_RECCherenkov,i);
      }
  }
    
  importDataToParticle(particle);
  
  return true;
}

bool HipoBankInterface::importDataToParticle(part &particle)
{
  
  // -------------------------------------------------------------
  // Import the REC::Calorimeter data
  // -------------------------------------------------------------
  particle.pcal_sector = _sector_Cal[0];
  particle.ecin_sector = _sector_Cal[1];
  particle.ecout_sector = _sector_Cal[2];

  particle.pcal_e = _Ele_PCAL_e;
  particle.ecin_e = _Ele_ECIN_e;
  particle.ecout_e = _Ele_ECOUT_e;
  
  particle.pcal_m2u = _m2u_Cal[0];
  particle.pcal_m2v = _m2v_Cal[0];
  particle.pcal_m2w = _m2w_Cal[0];
    
  particle.pcal_x = _x_Cal[0];
  particle.ecin_x = _x_Cal[1];
  particle.ecout_x = _x_Cal[2];

  particle.pcal_y = _y_Cal[0];
  particle.ecin_y = _y_Cal[1];
  particle.ecout_y = _y_Cal[2];

  particle.pcal_z = _z_Cal[0];
  particle.ecin_z = _z_Cal[1];
  particle.ecout_z = _z_Cal[2];
  particle.pcal_lu = _lu_Cal[0];
  particle.ecin_lu = _lu_Cal[1];
  particle.ecout_lu = _lu_Cal[2];

  particle.pcal_lv = _lv_Cal[0];
  particle.ecin_lv = _lv_Cal[1];
  particle.ecout_lv = _lv_Cal[2];

  particle.pcal_lw = _lw_Cal[0];
  particle.ecin_lw = _lw_Cal[1];
  particle.ecout_lw = _lw_Cal[2];

  particle.sector = _sector_DC;

  particle.traj_x1 = _x_DC[0];
  particle.traj_y1 = _y_DC[0];
  particle.traj_z1 = _z_DC[0];
  particle.edge_1 = _edge[0];

  particle.traj_x2 = _x_DC[1];
  particle.traj_y2 = _y_DC[1];
  particle.traj_z2 = _z_DC[1];
  particle.edge_2 = _edge[1];

  particle.traj_x3 = _x_DC[2];
  particle.traj_y3 = _y_DC[2];
  particle.traj_z3 = _z_DC[2];
  particle.edge_3 = _edge[2];
    
  particle.nphe_ltcc = _nphe_ltcc;
  particle.nphe_htcc = _nphe_htcc;
  return true;
}

int HipoBankInterface::determineSectorDC(float x, float y, float z){
  float phi = 180 / PI * atan2(y / sqrt(pow(x,2) + pow(y,2) + pow(z,2)),
			       x /sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
  if(phi<30 && phi>=-30){return 1;}
  else if(phi<90 && phi>=30){return 2;}
  else if(phi<150 && phi>=90){return 3;}
  else if(phi>=150 || phi<-150){return 4;}
  else if(phi<-90 && phi>=-150){return 5;}
  else if(phi<-30 && phi>-90){return 6;}

  return 0;

}


void HipoBankInterface::clear(){
  _Ele_PCAL_e = 0.0;
  _Ele_ECIN_e = 0.0;
  _Ele_ECOUT_e = 0.0;
  _sector_DC = -1;
  _nphe_ltcc = 0.0;
  _nphe_htcc = 0.0;
  for(int i = 0 ; i < 3 ; i++){
    _sector_Cal[i]=0;
    _time_Cal[i]=0;
    _path_Cal[i]=0;
    _x_Cal[i]=0;
    _y_Cal[i]=0;
    _z_Cal[i]=0;
    _m2u_Cal[i]=0;
    _m2v_Cal[i]=0;
    _m2w_Cal[i]=0;
    _lu_Cal[i]=0;
    _lv_Cal[i]=0;
    _lw_Cal[i]=0;
    _det_DC[i]=0;
    _path_DC[i]=0;
    _x_DC[i]=0;
    _y_DC[i]=0;
    _z_DC[i]=0;
  }
}
