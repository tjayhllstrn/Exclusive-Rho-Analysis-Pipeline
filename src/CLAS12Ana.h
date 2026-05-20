#ifndef CLAS12Ana_H
#define CLAS12Ana_H

#include "Constants.h"
#include "Structs.h"
#include "HipoBankInterface.h"
#include "Kinematics.h"

class CLAS12Ana {
public:
    // Constructor and destructor
    CLAS12Ana();
       CLAS12Ana(const std::unique_ptr<clas12::clas12reader>&);
    
//     // Member functions
       void set_run_config(const std::unique_ptr<clas12::clas12reader>&);
       void set_beams(TLorentzVector, TLorentzVector);
       std::vector<part> load_reco_particles(const std::unique_ptr<clas12::clas12reader>&);
       std::vector<part> load_mc_particles(const std::unique_ptr<clas12::clas12reader>&);
       int find_reco_scattered_electron(std::vector<part>&);
       int find_mc_scattered_electron(std::vector<part>&);
       bool reco_event_contains_final_state(std::vector<part>, FS);
       bool reco_event_contains_scattered_electron(std::vector<part>);
    
       void get_event_info(const std::unique_ptr<clas12::clas12reader>&, EVENT_INFO &event_info);
       void fill_reco_event_variables(EVENT &, std::vector<part>);
       void fill_mc_event_variables(EVENT &, std::vector<part>);
    
       void match_mc_to_reco(std::vector<part>&, std::vector<part>&);
    void match_mc_to_reco_mcmatch(std::vector<part>&, std::vector<part>&, const std::unique_ptr<clas12::clas12reader>&);
    
       std::vector<std::vector<int>> dihadron_idxs(int,int,int[],int);
       std::vector<std::vector<int>> dihadron_idxs(int,int,std::vector<int>);
       std::vector<std::vector<int>> dihadron_idxs(int,int,std::vector<part>);
       void clear_dihadron_variables(EVENT &);
       void fill_mc_reco_dihadron_variables(EVENT &, TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector, std::vector<part>, std::vector<int>, int,int);
       void fill_reco_dihadron_variables(EVENT &, TLorentzVector, TLorentzVector, std::vector<part>, std::vector<int>, int,int);
    protected:
        // Dihadron indexing code
        void generate_combinations(std::vector<int>& input, int num, int start_idx, std::vector<int>& curr_combination, std::vector<std::vector<int>>& result);
        std::vector<std::vector<int>> unique_combinations(std::vector<int> input, int num);
        std::vector<std::vector<int>> remove_duplicates(std::vector<std::vector<int>> input);

    private:
        // Private member variables
        TLorentzVector init_electron, target;
        HipoBankInterface hipoBankInterface;
        Kinematics _kin;
        double _electron_beam_energy;
        double s; // com energy
        int _idx_RUNconfig;
        int _irun;
        int _ievnum;
        int _itorus;
};

#endif  // CLAS12ANA_H
