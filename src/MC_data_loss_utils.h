// Utilities for MC/data loss studies: filtering and plotting helpers.

#ifndef MC_DATA_LOSS_UTILS_H
#define MC_DATA_LOSS_UTILS_H

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TCut.h>
#include <TH1F.h>
#include <TTree.h>

class filter {
public:
	explicit filter(const char* phi_var = "phi", const char* mh_var = "Mh");

	void SetPhiVar(const char* name);
	void SetMhVar(const char* name);

	void SetMhRange(double min_val, double max_val);
	void ClearMhRange();

	void SetExtraCut(const std::string& cut);

	TCut BuildBaseCut() const;
	TCut BuildPhiCut(double phi_min, double phi_max) const;
	TCut BuildBinCut(double phi_min, double phi_max) const;

	const std::string& PhiVar() const;
	const std::string& MhVar() const;
	bool HasMhRange() const;
	double MhMin() const;
	double MhMax() const;

private:
	std::string phi_var_;
	std::string mh_var_;
	std::string extra_cut_;
	bool use_mh_range_;
	double mh_min_;
	double mh_max_;
};

struct Chi2Result {
	double chi2;
	int ndf;
	double chi2_ndf;
	double mc_scale;
};

class plotter {
public:
	explicit plotter(const std::string& mh_var = "Mh", int mh_bins = 100, double mh_min = 0.5, double mh_max = 1.25);

	TH1F* MakeMhHist(const std::string& name, const std::string& title) const;
	void FillHist(TTree* tree, TH1F* hist, const TCut& cut) const;
	double NormalizeToData(TH1F* mc_hist, const TH1F* data_hist) const;
	Chi2Result ComputeChi2(const TH1F* data_hist, const TH1F* mc_hist) const;

	TCanvas* MakePostageStamp(const std::vector<TH1F*>& data_hists,
							  const std::vector<TH1F*>& mc_hists,
							  const std::vector<double>& phi_edges,
							  const std::vector<Chi2Result>& chi2_vals,
							  const std::string& title,
							  const std::string& canvas_name) const;

	const std::string& MhVar() const;
	int MhBins() const;
	double MhMin() const;
	double MhMax() const;

private:
	void StyleDataHist(TH1F* hist) const;
	void StyleMcHist(TH1F* hist) const;

	std::string mh_var_;
	int mh_bins_;
	double mh_min_;
	double mh_max_;
};

std::vector<double> parse_csv_to_doubles(const std::string& s);
std::string GetBaseName(const std::string& path);
std::string StripExtension(const std::string& filename);
std::string BuildComboName(const std::string& mc_file, const std::string& data_file);
void EnsureDir(const std::string& path);

#endif
