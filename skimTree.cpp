#include <iostream>
#include <string>
#include <vector>
#include <limits> // numeric_limits
#include <algorithm>
#include <filesystem>
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include "TFormula.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCollection.h"
#include "TChain.h"

#define CGREEN "\033[0;32m"
#define CEND "\033[0m"
#define cprintf(X) printf("%s%s%s", CGREEN, X, CEND)

using namespace std;
namespace fs = std::filesystem;
// script to skim the signal mc ntuples 
// g++ skimTree.cpp -o ./skimTree -I$CONDA_PREFIX/include -std=c++17 -Wall -O3 $(root-config --glibs --cflags --libs)
// ./skimTree [era], era = UL2016preVFP, UL2016postVFP, UL2017, UL2018


//===============================================//
//        Set up the input and ouput files       //
//                      UL2018                   //
//===============================================//
vector<const char*> inpath_UL2018 = {
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_VBF_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_WH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ZH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ttH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_bbH_m125/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_VBF_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_WH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ZH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ttH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_bbH_m120/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_VBF_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_WH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ZH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_ttH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL18_Dalitz_mmg_bbH_m130/*.root"
};
vector<const char*> outpath_UL2018 = {
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_VBF_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m125.root",

    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_VBF_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m120.root",

    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_VBF_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m130.root"
};


//===============================================//
//        Set up the input and ouput files       //
//                      UL2017                   //
//===============================================//
vector<const char*> inpath_UL2017 = {
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_VBF_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_WH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ZH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ttH_m125/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_bbH_m125/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_VBF_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_WH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ZH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ttH_m120/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_bbH_m120/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_VBF_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_WH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ZH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_ttH_m130/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_mmg_bbH_m130/*.root"
};
vector<const char*> outpath_UL2017 = {
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_VBF_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_WH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ZH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ttH_m125.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_bbH_m125.root",

    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_VBF_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_WH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ZH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ttH_m120.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_bbH_m120.root",

    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_VBF_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_WH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ZH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_ttH_m130.root",
    "/data4/chenghan/muon/skimTree/job_UL17_Dalitz_mmg_bbH_m130.root"
};


//===============================================//
//        Set up the input and ouput files       //
//                   UL2016preVFP                //
//===============================================//
vector<const char*> inpath_UL2016preVFP = {
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m125_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m125_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m125_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m125_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m125_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m125_preVFP/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m120_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m120_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m120_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m120_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m120_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m120_preVFP/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m130_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m130_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m130_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m130_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m130_preVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m130_preVFP/*.root"
};
vector<const char*> outpath_UL2016preVFP = {
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m125_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m125_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m125_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m125_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m125_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m125_preVFP.root",

    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m120_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m120_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m120_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m120_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m120_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m120_preVFP.root",

    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m130_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m130_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m130_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m130_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m130_preVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m130_preVFP.root"
};


//===============================================//
//        Set up the input and ouput files       //
//                   UL2016postVFP                //
//===============================================//
vector<const char*> inpath_UL2016postVFP = {
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m125_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m125_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m125_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m125_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m125_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m125_postVFP/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m120_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m120_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m120_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m120_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m120_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m120_postVFP/*.root",

    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_m130_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_VBF_m130_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_WH_m130_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ZH_m130_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_ttH_m130_postVFP/*.root",
    "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_mmg_bbH_m130_postVFP/*.root"
};
vector<const char*> outpath_UL2016postVFP = {
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m125_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m125_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m125_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m125_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m125_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m125_postVFP.root",

    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m120_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m120_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m120_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m120_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m120_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m120_postVFP.root",

    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_m130_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_VBF_m130_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_WH_m130_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ZH_m130_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_ttH_m130_postVFP.root",
    "/data4/chenghan/muon/skimTree/job_UL16_Dalitz_mmg_bbH_m130_postVFP.root"
};

float CalcMLL(
    // MLL cut at 60 GeV
    const ROOT::RVec<int>& lhePID,
    const ROOT::RVec<float>& lhePx,
    const ROOT::RVec<float>& lhePy,
    const ROOT::RVec<float>& lhePz
){
    // MLL cut at 60 GeV
    auto LHE_Muons = abs(lhePID) == 13;
    float mass = std::numeric_limits<float>::max();
    if (ROOT::VecOps::Nonzero(LHE_Muons).size() == 2){
        auto LHEMuPx = lhePx[LHE_Muons];
        auto LHEMuPy = lhePy[LHE_Muons];
        auto LHEMuPz = lhePz[LHE_Muons];
        ROOT::Math::PxPyPzMVector mu1(LHEMuPx[0], LHEMuPy[0], LHEMuPz[0], 0.1057);
        ROOT::Math::PxPyPzMVector mu2(LHEMuPx[1], LHEMuPy[1], LHEMuPz[1], 0.1057);
        mass = (mu1 + mu2).M();
    }
    return mass;
}

bool IsPhoIntConv(
    // remove internal conversion photon
    const int nMC,
    const ROOT::RVec<int>& mcPID,
    const ROOT::RVec<int>& mcMomPID,
    const ROOT::RVec<unsigned short> mcStatusFlag
){
    int realpho = 0;
    for (int i = 0; i < nMC; i++){
        bool hardProc = (mcStatusFlag[i] >> 0 & 1) == 1;
        bool isPrompt = (mcStatusFlag[i] >> 1 & 1) == 1;
        bool isRealPho = mcPID[i] == 22 && mcMomPID[i] == 25 && hardProc && isPrompt;
        if (isRealPho)
            realpho += 1;
    }
    bool isPhoIntConv = (realpho == 0) ? true : false; // no real photon
    return isPhoIntConv;
}


int main(int argc, char** argv){
    TStopwatch time;
	time.Start();

    if(argc < 2){
        printf("Usage: ./skimTree [era]\n");
        return 1;
    }
    std::string era = argv[1];

    vector<const char*> inpath;
    if (era == "UL2016preVFP")  inpath = inpath_UL2016preVFP;
    if (era == "UL2016postVFP") inpath = inpath_UL2016postVFP;
    if (era == "UL2017")        inpath = inpath_UL2017;
    if (era == "UL2018")        inpath = inpath_UL2018;

    vector<const char*> outpath;
    if (era == "UL2016preVFP")  outpath = outpath_UL2016preVFP;
    if (era == "UL2016postVFP") outpath = outpath_UL2016postVFP;
    if (era == "UL2017")        outpath = outpath_UL2017;
    if (era == "UL2018")        outpath = outpath_UL2018;

    if (inpath.size() == outpath.size()){
        ROOT::EnableImplicitMT(20); // enable multithreading
        cprintf("Start to skim the ggNtuples!\n");
        cprintf("Skim the following ntuples sequentially\n");
        for (size_t i = 0; i < inpath.size(); i++){
            cprintf(Form(" %ld. %s\n", i+1, inpath[i]));
        }
        printf("\n");
        
        for (size_t i = 0; i < inpath.size(); i++){
            printf("********************************   %ld    ********************************\n", i+1);
            printf("Find input file(s): %s\n", inpath[i]);
            fs::path p = outpath[i];
            if (!fs::exists(p.parent_path()))
                system(Form("mkdir -p %s", p.parent_path().c_str()));
            
            printf("Save skim file    : %s\n", outpath[i]);

            TChain chain1("ggNtuplizer/EventTree");
            TChain chain2("CorrTree");

            TString InPathStr(inpath[i]);
            TString InDir = gSystem->GetDirName(InPathStr.Data());
            TString SSInDir = gSystem->GetDirName(InPathStr.Data());
            bool sspath_found = false;
            // find the path containing shower shape corrections
            // eg: ntuple : "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_eeg_VBF_m130"
            //     sscorr : "/data5/ggNtuples/V10_06_30_02_sscorr/job_UL17_Dalitz_eeg_VBF_m130"
            TString tok;
            int position = 0;
            while (InDir.Tokenize(tok, position, "/")){
                if (tok.Contains("V10")){
                    SSInDir.ReplaceAll(tok, tok+"_sscorr");
                    break;
                }
            }
            if (fs::exists(SSInDir.Data())){
                printf(" [+] Shower shape corrected path found:\n");
                printf("     - %s\n", SSInDir.Data());
                sspath_found = true;
            }
            else 
                throw std::runtime_error("The path which contains shower shape corrected files cannot be found");


            if (InPathStr.EndsWith("*.root")){ // all the files
                TSystemDirectory SD("SD", InDir.Data());
                TList* FileList = SD.GetListOfFiles();
                if (FileList){
                    TIter Next(FileList);
                    while (auto File = (TSystemFile*) Next()){
                        TString fName = File->GetName();
                        if (File->IsDirectory() || !fName.EndsWith(".root"))
                            continue;
                        chain1.Add(Form("%s/%s", InDir.Data(), fName.Data()));
                        if (sspath_found){
                            fName.ReplaceAll(".root", "_corr.root");
                            const char* fullSSPath = Form("%s/%s", SSInDir.Data(), fName.Data());
                            if (fs::exists(fullSSPath))
                                chain2.Add(fullSSPath);
                            else
                                throw std::runtime_error(Form("File doesn't exist: %s", fullSSPath));
                        }
                    }
                    delete FileList;
                }
                if (sspath_found)
                    chain1.AddFriend(&chain2); 
            }
            else{ // one file
                chain1.Add(InPathStr.Data());
                if (sspath_found){
                    TString SSPathStr(InPathStr);
                    SSPathStr.ReplaceAll(InDir.Data(), SSInDir.Data());
                    SSPathStr.ReplaceAll(".root", "_corr.root");
                    if (fs::exists(SSPathStr.Data()))
                        chain2.Add(SSPathStr.Data());
                    else
                        throw std::runtime_error(Form("File doesn't exist: %s", SSPathStr.Data()));
                    chain1.AddFriend(&chain2); 
                }
            }
            
            auto df = ROOT::RDataFrame(chain1);
            auto df1 = df.Define("diLepMCMass", CalcMLL,      {"lhePID", "lhePx", "lhePy", "lhePz"})
                         .Define("HasIntPho",   IsPhoIntConv, {"nMC", "mcPID", "mcMomPID", "mcStatusFlag"})
                         .Filter("diLepMCMass < 60 && HasIntPho == 0", "mc preslections");

            auto report = df1.Report(); 
            df1.Snapshot("skimTree", outpath[i]);
            report->Print();
            printf("\n");
        }
    }
    else {
        printf("Error: the size of inpath and outpath are different for %s", era.c_str());
        return 1;
    }
    
    time.Stop();
    time.Print();

    return 0;
}