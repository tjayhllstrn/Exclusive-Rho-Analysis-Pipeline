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
                    std::vector<double> phibn_edges,
                    std::vector<double> bn_edgs,
                    std::vector<double> obs2bn,
                    std::string fit_type);

    //destructor
    ~Chi2_Fitter();

    //member variables
    std::string TREENAME;
    std::string OUT_DIR;
    std::string OBS;
    std::string OBS2;
    std::string FIT_TYPE;
    std::string IN_FILE;
    std::vector<double> PHIBN_EDGES;
    std::vector<double> PHIBN_CENTERS;
    std::vector<double> BN_EDGS;
    std::vector<double> BN_CENTERS; 
    std::vector<double> OBS2BN;
    TTree* RAW_TREE;
    int obs_bin_idx;
    int phi_bin_idx;
    int obs2_bin_idx;

    //data structs for managing results
    std::vector<std::vector<std::pair<double, double>>> N_sig_pos; //vector indexed by [obs bin, phibin] holding pairs of (N_sig, N_sig_err)
    std::vector<std::vector<std::pair<double, double>>> N_sig_neg;
    std::vector<std::vector<std::pair<double, double>>> alpha;
    std::vector<std::pair<double, double>> A_sig;

    //data structs for managing plots
    std::vector<std::vector<std::pair<TH1F*, TH1F*>>> N_sig_fitting_datathist; //pair is (pos,neg) graphs
    std::vector<std::vector<std::pair<TGraph*, TGraph*>>> N_sig_fitting_totalgraph; 
    std::vector<std::vector<std::pair<TGraph*, TGraph*>>> N_sig_fitting_siggraph;
    std::vector<std::vector<std::pair<TGraph*, TGraph*>>> N_sig_fitting_bkggraph;
    std::vector<std::vector<std::pair<TLegend*, TLegend*>>> N_sig_fitting_legends;
    std::vector<std::vector<std::pair<TLatex*, TLatex*>>> N_sig_fitting_texts;
    std::vector<std::vector<std::pair<TPaveText*, TPaveText*>>> N_sig_fitting_paramboxes;
    TCanvas* SinCanvas;

    
    //fitting methods
    void RunMhChi2Fit(int obs2bn_idx);
    void RunMxChi2Fit(int obs2bn_idx);
    std::vector<std::vector<std::pair<double, double>>> FitChi2(TTree* filteredTree,TCut* neg_hel, int helicity);
    std::pair<double,double> Mh_sig_fit(TTree* filteredTree, TCut bin_cut, int helicity);
    std::pair<double,double> Mx_sig_fit(TTree* filteredTree, TCut bin_cut, int helicity);
    void CalcAlpha(TTree* filteredTree, int obs_bin_idx, int phi_bin_idx);
    void FitToSin(std::vector<double>& x_vals, std::vector<std::pair<double,double>>& y, int obs_bin_idx);

    //plotting methods
    void BinningSchemePlot(TTree* filteredTree);
    void PlotToCanvas_N_sig_BarHist();
    void PlotSigFitGraph(RooDataSet& binned_data, RooRealVar& x,
                                     RooAbsPdf& sig, RooAbsPdf& background,RooRealVar& N_sig,RooRealVar& N_bkg,RooAddPdf& model_ext, int helicity);
    void PlotToCanvas_PostageStamp(std::vector<TH1F*>& data_hists,
                                            std::vector<TGraph*>& total_graphs,
                                            std::vector<TGraph*>& sig_graphs,
                                            std::vector<TGraph*>& bkg_graphs,
                                            std::vector<TLegend*>& legends,
                                            std::vector<TLatex*>& texts,
                                            std::vector<TPaveText*>& param_boxes,
                                            std::string title);
    void MakeGraphLinePlot(std::vector<std::pair<double,double>>& y, std::vector<double>& x, 
                            const char* y_title, const char* x_title, 
                            const char* title, Color_t color,  
                            std::pair<double,double>* bounds=nullptr,std::vector<TGraph*>* out_graph=nullptr);

};
#endif