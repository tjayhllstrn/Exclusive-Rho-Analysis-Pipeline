#ifndef Chi2_Fitter_h
#define Chi2_Fitter_h
#include <vector>
#include <string>

class Chi2_Fitter {

    public:
    //constructor  
    Chi2_Fitter(const char* treename, const char* out_dir,
                   const char* obs_s,
                   const char* obs2_s,
                   const char* in_file,
                   std::vector<double> bn_edgs,
                   std::vector<double> obs2bn,
                   std::string fit_type);
        : TREENAME(treename),
        OUT_DIR(out_dir),
        OBS(obs_s),
        OBS2(obs2_s),
        IN_FILE(in_file),
        BN_EDGS(bn_edgs),
        OBS2BN(obs2bn),
        FIT_TYPE(fit_type)
        {
            // suppress RooFit messages
            RooMsgService* rms = &RooMsgService::instance();
            rms->setSilentMode(true);
            rms->setGlobalKillBelow(RooFit::ERROR);  

            //load the file and save the reference to the tree
            std::cout<<"Loading data from file: " << in_file << std::endl;
            TFile* file = TFile::Open(in_file);
            if (!file || file->IsZombie()) {
            throw std::runtime_error("Error: could not open ROOT file: " + std::string(in_file));
            }
            RAW_TREE = (TTree*)file->Get(treename);

            //populate the BN_CENTERS vector for later plotting
            BN_CENTERS.reserve(bn_edgs.size() - 1);
            for (size_t i = 0; i < bn_edgs.size() - 1; ++i) {
                double center = 0.5 * (bn_edgs[i] + bn_edgs[i + 1]);
                BN_CENTERS.push_back(center);
            }
        }

    //member variables

    //data structs for managing results

    //data structs for managing plots

    //methods

};
#endif