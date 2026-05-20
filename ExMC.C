#include <TLorentzVector.h>
#include "clas12reader.h"
// #include "HipoChain.h"
const char* RESET  = "\e[0m";
const char* GREEN  = "\e[1;32m";
const char* YELLOW = "\e[1;33m";
const char* RED    = "\e[1;31m";
void ExMC(){

   clas12root::HipoChain chain;
   chain.Add("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/clasdis_rga_fa18_inb_45nA_10604MeV-0672.hipo");
   auto c12=chain.GetC12Reader(); //this line is only in for backward compatability with older versions of clas12root

   TLorentzVector p4;
   int i = 0;
   while (chain.Next()){
      if(i>200) break;
      std::cout << GREEN << "Event " << i << RESET << std::endl;
      auto c12=chain.GetC12Reader();

      auto mceve=c12->mcevent(); 
      cout<<" beam energy  "<<mceve->getEbeam()<<" type "<<mceve->getBtype()<<endl;

      auto mcpbank=c12->mcparts(); 
      const Int_t  Ngen=mcpbank->getRows();
      
      for(Int_t j=0;j<Ngen;j++){
         mcpbank->setEntry(j);
         
         auto px=mcpbank->getPx();
         auto py=mcpbank->getPy();
         auto pz=mcpbank->getPz();
         auto pm=mcpbank->getMass();
         p4.SetXYZM(px,py,pz,pm);

         auto pid = mcpbank->getPid();
         auto pindex = mcpbank->getIndex();

         if(mcpbank->getType()!=1) // Reject non-final state
         {continue;} 

         cout<<" particle "<<pindex<<" pid: "<<pid<<" p4 = "<<p4.X()<<","<<p4.Y()<<","<<p4.Z()<<","<<p4.T()<<" and mass "<<p4.M()<<endl;
      }

      auto mcmatch = mcpbank->getMatch(); // pointer to MC::RecMatch if present
      if (!mcmatch) {
      std::cout << "No MC::RecMatch bank available\n";
      continue;
      }

      const int nmatch = mcmatch->getRows();
      for (int k = 0; k < nmatch; ++k) {
         mcmatch->setEntry(k);
         int rec_idx = mcmatch->getPindex();
         int mc_idx  = mcmatch->getMCindex();
         float q     = mcmatch->getQuality();

         std::cout << "match row " << k
                     << " rec pindex=" << rec_idx
                     << " mcindex=" << mc_idx
                     << " quality=" << q << "\n";
      }
      
      i++;
         
      }

   }
