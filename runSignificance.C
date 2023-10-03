#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ranges.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "boost/filesystem.hpp"
#include "./plugin/sigmaEff.h"

using namespace std;
namespace bfs = boost::filesystem;


vector<string> cats = {
    "diMu9.0MuPho_EBHR9",
    "diMu9.0MuPho_EBLR9",
    "diMu9.0MuPho_EE",
    "diMu9.0MuPho_VBF",
    "diMu9.0MuPho_BST",

    "diMu25MuPho_EBHR9",
    "diMu25MuPho_EBLR9",
    "diMu25MuPho_EE",
    "diMu25MuPho_VBF",
    "diMu25MuPho_BST",
    "diMu50MuPho",

    "diMu9.0IsoMu_EBHR9",
    "diMu9.0IsoMu_EBLR9",
    "diMu9.0IsoMu_EE",
    "diMu9.0IsoMu_VBF",
    "diMu9.0IsoMu_BST",

    "diMu50IsoMu_EBHR9",
    "diMu50IsoMu_EBLR9",
    "diMu50IsoMu_EE",
    "diMu50IsoMu_VBF",
    "diMu50IsoMu_BST"
};


map<string, double> calc_sigmaEff(ROOT::RDF::RNode df_sig_cat, ROOT::RDF::RNode df_bkg_cat, string mass_column, string weight_column, vector<int> mass_range){
    map<string, double> effsig;
    auto arr = df_sig_cat.Take<float>(mass_column);
    vector<float> mass = *arr;
    vector<float> sigma1 = sigmaEff(mass, 0.683);
    vector<float> sigma2 = sigmaEff(mass, 0.955);

    string sig_range = fmt::format("{} > {} && {} < {}", mass_column, sigma2[0], mass_column, sigma2[1]);
    auto sig_stats = df_sig_cat.Stats(mass_column, weight_column);
    auto sig_stats_range = df_sig_cat.Filter(sig_range).Stats(mass_column, weight_column);
    const double sig_yield = sig_stats->GetW();
    const double sig_yield_range = sig_stats_range->GetW();
    effsig["Resolution"] = sigma1[2];
    effsig["Signal"] = sig_yield;
    effsig["Signal_2sigma"] = sig_yield_range;

    // non-resonant background is estimated from data sideband region scaled to the same mass window as signal(2 sigma)
    string sideband_range = fmt::format("({} > {} && {} < {}) || ({} < {} && {} > {})", mass_column, mass_range[0], mass_column, sigma2[0], mass_column, mass_range[1], mass_column, sigma2[1]);
    auto sideband_stats = df_bkg_cat.Stats(mass_column);
    auto sideband_stats_range = df_bkg_cat.Filter(sideband_range).Stats(mass_column);
    const double sideband_yield = sideband_stats->GetN();
    const double sideband_yield_range = sideband_stats_range->GetN() * (sigma2[1] - sigma2[0])/(mass_range[1] - mass_range[0]);
    effsig["Sideband"] = sideband_yield;
    effsig["Sideband_2sigma"] = sideband_yield_range; // %
    effsig["Ratio"] = sig_yield_range*100/sideband_yield_range;
    effsig["AMS"] = sqrt(2 * ( (sig_yield_range + sideband_yield_range) * log( 1 + (sig_yield_range/sideband_yield_range) ) - sig_yield_range));

    return effsig;
}


void runSignificance(int isIsoMu){
    vector<string> bkg_files;
    vector<string> sig_files;

    if (isIsoMu != 1){
        bkg_files = {
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016B_UL2016preVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016C_UL2016preVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016D_UL2016preVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016E_UL2016preVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016F1_UL2016preVFP.root",

            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016F2_UL2016postVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016G_UL2016postVFP.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016H_UL2016postVFP.root",

            "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018A_UL2018.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018B_UL2018.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018C_UL2018.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018D_UL2018.root",

            "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017D_UL2017.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017E_UL2017.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017F_UL2017.root"
        };
    }
    else{
        bkg_files = {
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_SingleMu_Run2017B_UL2017.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_SingleMu_Run2017C_UL2017.root"
        };
    }

    if (isIsoMu != 1){
        sig_files = {
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ggF_mmg_125_UL2018_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_VBF_mmg_125_UL2018_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_WH_mmg_125_UL2018_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ZH_mmg_125_UL2018_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ttH_mmg_125_UL2018_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_bbH_mmg_125_UL2018_MuPhoTrig.root",

            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ggF_mmg_125_UL2017_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_VBF_mmg_125_UL2017_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_WH_mmg_125_UL2017_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ZH_mmg_125_UL2017_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ttH_mmg_125_UL2017_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_bbH_mmg_125_UL2017_MuPhoTrig.root",

            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_ggF_mmg_125_UL2016preVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_VBF_mmg_125_UL2016preVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_WH_mmg_125_UL2016preVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_ZH_mmg_125_UL2016preVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_ttH_mmg_125_UL2016preVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_HDalitz_bbH_mmg_125_UL2016preVFP_MuPhoTrig.root",

            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_ggF_mmg_125_UL2016postVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_VBF_mmg_125_UL2016postVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_WH_mmg_125_UL2016postVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_ZH_mmg_125_UL2016postVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_ttH_mmg_125_UL2016postVFP_MuPhoTrig.root",
            "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_HDalitz_bbH_mmg_125_UL2016postVFP_MuPhoTrig.root",
        };
    }
    else{
        sig_files = {
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ggF_mmg_125_UL2017_SingleMuTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_VBF_mmg_125_UL2017_SingleMuTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_WH_mmg_125_UL2017_SingleMuTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ZH_mmg_125_UL2017_SingleMuTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ttH_mmg_125_UL2017_SingleMuTrig.root",
            "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_bbH_mmg_125_UL2017_SingleMuTrig.root"
        };
    }

    ROOT::RDF::RNode df_sig = ROOT::RDataFrame("miniTree", sig_files);
    ROOT::RDF::RNode df_bkg = ROOT::RDataFrame("miniTree", bkg_files);

    vector<int> mass_range = {110, 170};
    string mass_column = "CMS_higgs_mass";
    auto weight_column = "weight";
    // calculate the combined significance and format to a table
    TString info = TString::Format("%-20s|%10s|%10s|%10s|%10s|%10s|%10s|%12s|\n", "category", "AMS", "f95.5[%]","sigma[GeV]", "Nsig", "Nsig95.5", "Nbkg", "Nbkg95.5");
    info += TString::Format("====================================================================================================\n");

    // int cat = (isIsoMu != 1) ? 20 : 11;

    int cat_init = (isIsoMu != 1) ? 0: 11;
    int cat_end = (isIsoMu != 1) ? 11: cats.size();
    for (int i = cat_init; i < cat_end; i++){
        string cat_name = fmt::format("{}", cats[i]);
        string cat_filter = fmt::format("category == {}", i+1);
        auto df_sig_cat = df_sig.Filter(cat_filter);
        auto df_bkg_cat = df_bkg.Filter(cat_filter);

        auto effsig = calc_sigmaEff(df_sig_cat, df_bkg_cat, mass_column, weight_column, mass_range);
        info += TString::Format("%-20s|%10.2f|%10.2f|%10.2f|%10.2f|%10.2f|%10d|%12.2f|\n", cat_name.c_str(), effsig["AMS"], effsig["Ratio"], effsig["Resolution"], effsig["Signal"], effsig["Signal_2sigma"], (int)effsig["Sideband"], effsig["Sideband_2sigma"]);
    }
    printf("%s\n", info.Data());

    // save the results to the file
    // string outName = "./ams/AMS_UL2017_SingleMu.txt";
    string outName = (isIsoMu != 1) ? "./ams/AMS_fullRun2_mupho.txt" : "./ams/AMS_fullRun2_isomu.txt";
    printf("Save table in: %s\n", outName.c_str());
    bfs::path outpath = outName;
    if (!bfs::exists(outpath.parent_path()))
        system(Form("mkdir -p %s", outpath.parent_path().c_str()));
    std::ofstream outfile(outName);
    if (outfile.is_open()){
        outfile << info.Data();
        outfile.close();
    }
    else
        throw std::runtime_error(Form("Failed to open the file: %s", outName.c_str()));
}