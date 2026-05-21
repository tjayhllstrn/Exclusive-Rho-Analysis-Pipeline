#include "../src/MC_data_loss_utils.h"
#include "../src/MC_data_loss_utils.cpp"

#include <TEnv.h>
#include <TFile.h>
#include <TTree.h>

#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
//This filter_and_strip file is meant to filter and strip a TTree to implement all the neccesary cuts and get rid
// of unnecesary branches. You can modify the filter and which branches to keep by changing the config file.

static std::string Trim(const std::string& value) {
	const size_t begin = value.find_first_not_of(" \t\n\r");
	if (begin == std::string::npos) {
		return "";
	}
	const size_t end = value.find_last_not_of(" \t\n\r");
	return value.substr(begin, end - begin + 1);
}

static std::vector<std::string> ParseCsvStrings(const std::string& value) {
	std::vector<std::string> out;
	std::istringstream iss(value);
	std::string token;
	while (std::getline(iss, token, ',')) {
		token = Trim(token);
		if (!token.empty()) {
			out.push_back(token);
		}
	}
	return out;
}

void filter_and_strip_TTree(const char* input_root, const char* config_file) {
	if (!input_root || !config_file) {
		std::cerr << "Input ROOT file and config file are required." << std::endl;
		return;
	}

	TEnv env;
	if (env.ReadFile(config_file, kEnvLocal) != 0) {
		std::cerr << "Failed to read config file: " << config_file << std::endl;
		return;
	}

	std::string label = Trim(env.GetValue("label", "default"));
	std::string tree_name = Trim(env.GetValue("treename", "pippi0"));
	std::string cut_string = Trim(env.GetValue("cut", ""));
	std::string branches_csv = env.GetValue("branches", "");
	std::string phi_var = Trim(env.GetValue("phi_var", "phi"));
	std::string mh_var = Trim(env.GetValue("mh_var", "Mh"));

	std::vector<std::string> branches = ParseCsvStrings(branches_csv);
	if (branches.empty()) {
		std::cerr << "No branches specified in config (key: branches)." << std::endl;
		return;
	}

	TFile* in_file = TFile::Open(input_root, "READ");
	if (!in_file || in_file->IsZombie()) {
		std::cerr << "Failed to open input ROOT file: " << input_root << std::endl;
		return;
	}

	TTree* in_tree = dynamic_cast<TTree*>(in_file->Get(tree_name.c_str()));
	if (!in_tree) {
		std::cerr << "Could not find tree '" << tree_name << "' in input file." << std::endl;
		in_file->Close();
		return;
	}

	filter tree_filter(phi_var.c_str(), mh_var.c_str());
	if (!cut_string.empty()) {
		tree_filter.SetExtraCut(cut_string);
	}
	if (env.Defined("mh_min") && env.Defined("mh_max")) {
		double mh_min = env.GetValue("mh_min", 0.0);
		double mh_max = env.GetValue("mh_max", 0.0);
		tree_filter.SetMhRange(mh_min, mh_max);
	}
	TCut base_cut = tree_filter.BuildBaseCut();

	in_tree->SetBranchStatus("*", 0);
	int enabled_branches = 0;
	for (const std::string& branch : branches) {
		if (in_tree->GetBranch(branch.c_str())) {
			in_tree->SetBranchStatus(branch.c_str(), 1);
			++enabled_branches;
		} else {
			std::cerr << "Warning: branch not found in input tree: " << branch << std::endl;
		}
	}
	if (enabled_branches == 0) {
		std::cerr << "None of the requested branches exist in the input tree." << std::endl;
		in_file->Close();
		return;
	}

	std::filesystem::path input_path(input_root);
	std::string base_name = StripExtension(GetBaseName(input_root));
	std::string out_file_name = base_name + "_filtered_" + label + ".root";
	std::filesystem::path out_path = input_path.parent_path() / out_file_name;

	TFile* out_file = TFile::Open(out_path.string().c_str(), "RECREATE");
	if (!out_file || out_file->IsZombie()) {
		std::cerr << "Failed to create output file: " << out_path.string() << std::endl;
		in_file->Close();
		return;
	}

	out_file->cd();
	TTree* out_tree = in_tree->CopyTree(base_cut);
	if (!out_tree) {
		std::cerr << "Failed to copy filtered tree." << std::endl;
		out_file->Close();
		in_file->Close();
		return;
	}
	out_tree->Write();
	out_file->Close();
	in_file->Close();

	std::cout << "Filtered tree written to: " << out_path.string() << std::endl;
	std::cout << "Tree: " << tree_name << ", branches kept: " << enabled_branches << std::endl;
	std::cout << "Cut: " << base_cut.GetTitle() << std::endl;
}
