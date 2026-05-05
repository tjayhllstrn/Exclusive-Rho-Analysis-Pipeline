#ifndef PurityFitter_h
#define PurityFitter_h

class PurityFitter{
    public:
    PurityFitter(TTree* input_tree,std::string out_dir);
    ~PurityFitter();

    //Member variables
    TTree* input_data;
    double minPhoE;
    double u;
    double u_err;
    std::vector<std::pair<double,double>> Purities; //this object holds the purity
    std::string OUT_DIR;
    std::string pre_cut_name;

    // data structures for managing plots
    std::vector<TH1F*> purity_data_hists;      // Store data histograms
    std::vector<TGraph*> purity_total_graphs;  // Store total fit graphs
    std::vector<TGraph*> purity_sig_graphs;    // Store signal component graphs
    std::vector<TGraph*> purity_bkg_graphs;    // Store background component graphs
    std::vector<TLegend*> purity_legends;      // Store legends
    std::vector<TLatex*> purity_texts;         // Store text annotations
    std::vector<TPaveText*> purity_param_boxes; // Store parameter boxes

    //Member Functions
    void FitPurity();
    void PurityFromFit(RooAbsPdf& sig, RooAbsPdf& background,
                       RooRealVar& N_sig, RooRealVar& N_bkg,RooRealVar& x,
                       RooFitResult* fit_results);
    void PlotPurityGraph(RooDataSet& binned_data, RooRealVar& x,
                         RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext);
    double CalculateChi2(TH1F* data_hist, TGraph* total_graph);
    void PlotToCanvas_PostageStamp();
    void MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
        const char* y_title, const char* x_title, 
        const char* title, Color_t color,  
        std::pair<double,double>* bounds=nullptr,std::vector<TGraph*>* out_graph=nullptr);

};

#endif