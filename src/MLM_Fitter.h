#ifndef MLM_Fitter_h
#define MLM_Fitter_h
#include <vector>
#include <string>

class MLM_Fitter {
    
    public:
    //constructor  
    MLM_Fitter(const char* treename, const char* out_dir,
                   const char* obs_s,
                   const char* obs2_s,
                   const char* in_file,
                   std::vector<double> bn_edgs,
                   std::vector<double> obs2bn,
                   std::string fit_type);

    //Destructor
    ~MLM_Fitter();
                   
    //member variables
    std::string TREENAME;
    std::string OUT_DIR;
    std::string OBS;
    std::string OBS2;
    std::string FIT_TYPE;
    std::string IN_FILE;
    std::vector<double> BN_EDGS;
    std::vector<double> BN_CENTERS;
    std::vector<double> OBS2BN;
    TTree* RAW_TREE;
    double sigbkg_min;
    double sigbkg_max;
    double bkg_min;
    double bkg_max;
    int bn_idx; //keeps track of what bin is being fit at any given time
    //RooDataSet* ptr_binned_data;
    double u;
    double u_err; 

    //data structures for managing results
    std::vector<std::vector<std::pair<double, double>>> A_sigbkg; //this object holds all obs2 bins (indexed by results[i]) and all obs bins (indexed by results[][j]). Each entry is a pair for (A,A_err)
    std::vector<std::vector<std::pair<double, double>>> A_bkg;
    std::vector<std::vector<std::pair<double, double>>> A_sig;
    std::vector<std::pair<double,double>> Purities; //this object holds the purity results for each obs bin in the current obs2 bin being processed

    //data structures for managing plots
    std::vector<TH1F*> purity_data_hists;      // Store data histograms
    std::vector<TGraph*> purity_total_graphs;  // Store total fit graphs
    std::vector<TGraph*> purity_sig_graphs;    // Store signal component graphs
    std::vector<TGraph*> purity_bkg_graphs;    // Store background component graphs
    std::vector<TLegend*> purity_legends;      // Store legends
    std::vector<TLatex*> purity_texts;         // Store text annotations
    std::vector<TPaveText*> purity_param_boxes; // Store parameter boxes


    //Purity Fit Functions
    void PurityCalc(TTree* tree);
    void PurityCalc_Mh(TTree* tree);
    void PurityCalc_Mx(TTree* tree);
    void PurityFromFit(RooAbsPdf& sig, RooAbsPdf& background,
                        RooRealVar& N_sig, RooRealVar& N_bkg,RooRealVar& x,
                        RooFitResult* fit_results);
    void PlotPurityGraph(RooDataSet& binned_data, RooRealVar& x,
                        RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext);
    double CalculateChi2(TH1F* data_hist, TGraph* total_graph);

    //MLM Fit Functions
    void RunMhFitMLM(int obs2bn_idx);
    void RunMxFitMLM(int obs2bn_idx);
    std::vector<std::pair<double,double>> FitMLM(TTree* tree,TCut* bounds);
    void Calc_A_sig_from_A_sigbkg(int obs2bn_idx);

    //Plotting Functions
    void PlotToCanvas_PostageStamp();
    void PlotToCanvas_overlayed(std::vector<TGraph*> graphs, const char* canvas_name, const char* canvas_title);
    void MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
        const char* y_title, const char* x_title, 
        const char* title, Color_t color,
        std::pair<double,double>* bounds=nullptr, std::vector<TGraph*>* out_graph=nullptr);

    
};
#endif
