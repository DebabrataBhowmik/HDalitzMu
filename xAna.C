#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH2F.h"
#include "TMVA/RReader.hxx"
// #include "boost/filesystem.hpp"
#include "./plugin/JsonSFsReader.h"
#include "./plugin/roccor.Run2.v5/RoccoR.cc"
#include "./plugin/Reader.h"
#include "./plugin/Writter.h"
#include "./plugin/CutCalculater.h"
#include "./plugin/Generator.h"
#include "./plugin/jetSelections.h"
#include "./plugin/puweicalc.h"

using namespace std;
// using namespace TMVA::Experimental;
// namespace bfs = boost::filesystem;

void print_progress(const Long64_t& progress, const Long64_t& total){
    cout << "[";
    float percent = (progress + 1) * 100./ total;
    float unit_length = 100. / 50;
    for(int ibar = 0; ibar < 50; ++ibar){
        if(ibar * unit_length < percent)
            cout << "=";
        else
            cout << " ";
    }
    cout << "]" << setw(3) << ceil(percent) << "% \r" << flush;
}


// Function to sort the index via the pT
//* Description:
//*     1) Calibrated energys, gsf pT are not sortted from large to small in ggNtuple.
//*     2) phoCalibEt = {10, 20, 30}, nPho =  3
//*     3) accInd = {1, 2}
//*     2) ArgSort => accInd = {2, 1}, max pT ind = 2 (accInd[0])
void ArgSort(vector<float> vec, vector<int>& accInd){
    for (size_t i = 0; i < accInd.size() - 1; i++){
        for (size_t j = 0; j < accInd.size() - i - 1; j++){
            if (vec[accInd[j + 1]] > vec[accInd[j]]){
                int tempidx = accInd[j];
                accInd[j] = accInd[j + 1];
                accInd[j + 1] = tempidx;
            }
        }
    }
}


TH1F* GetSFsHisto1D(TString fName, TString histName){
    TFile* fSFs = TFile::Open(fName.Data(), "READ");
    if (!fSFs || fSFs->IsZombie()){
        printf("TFile::Open failed %s\n", fName.Data());
        gSystem->Exit(1);
    }

    TH1F* hist = (TH1F*) fSFs->Get(histName.Data());
    if (!hist){
        printf("TFile::Get() failed: no histogram called  %s\n", histName.Data());
        gSystem->Exit(1);
    }

    hist->SetDirectory(0);
    fSFs->Close();

    return hist;
}


TH2F* GetSFsHisto2D(TString fName, TString histName){
    TFile* fSFs = TFile::Open(fName.Data(), "READ");
    if (!fSFs || fSFs->IsZombie()){
        printf("TFile::Open failed %s\n", fName.Data());
        gSystem->Exit(1);
    }

    TH2F* hist = (TH2F*) fSFs->Get(histName.Data());
    if (!hist){
        printf("TFile::Get() failed: no histogram called  %s\n", histName.Data());
        gSystem->Exit(1);
    }

    hist->SetDirectory(0);
    fSFs->Close();

    return hist;
}


int FindHistBin(TH2F* hist, float xval, float yval){
    int xbin = TMath::Max(1, TMath::Min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(xval)));
    int ybin = TMath::Max(1, TMath::Min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(yval)));
    int gbin = hist->GetBin(xbin, ybin); // global bin
    return gbin;
}



void xAna(TString inpath, TString outpath, double xs, double lumi, TString era, TString HLT, TString variation, bool isMC){
    TStopwatch time;
	time.Start();

    //========================================================//
    //Set up the branches which are needed to read from TTree//
    //========================================================//
    TString treeName = "ggNtuplizer/EventTree";
    if (isMC && !inpath.Contains("MuMuG") && !inpath.Contains("MuMuJet"))
        treeName = "skimTree";
    auto data = new TreeReader(inpath, treeName.Data());

    // event branches
    float rho = 0., rhoAll = 0., genWeight = 1.;
    double L1ECALPrefire = 1., L1ECALPrefireUp = 1., L1ECALPrefireDown = 1.;
    double MuonPrefire = 1., MuonPrefireUp = 1., MuonPrefireDown = 1.;
    Long64_t event = 0.;
    int run = 0., lumis = 0.;
    ULong64_t HLTEleMuX = 0.;
    bool isPVGood = false;
    auto puTrue                   = new vector<float>();
    data->SetBranchAddress("rho",           rho);
    data->SetBranchAddress("rhoAll",        rhoAll);
    data->SetBranchAddress("event",         event);
    data->SetBranchAddress("run",           run);
    data->SetBranchAddress("lumis",         lumis);
    data->SetBranchAddress("HLTEleMuX",     HLTEleMuX);
    data->SetBranchAddress("isPVGood",      isPVGood);
    if (data->HasMC()){
        data->SetBranchAddress("puTrue",            puTrue);
        data->SetBranchAddress("L1ECALPrefire",     L1ECALPrefire);
        data->SetBranchAddress("L1ECALPrefireUp",   L1ECALPrefireUp);
        data->SetBranchAddress("L1ECALPrefireDown", L1ECALPrefireDown);
        // data->SetBranchAddress("MuonPrefire",       MuonPrefire);
        // data->SetBranchAddress("MuonPrefireUp",     MuonPrefireUp);
        // data->SetBranchAddress("MuonPrefireDown",   MuonPrefireDown);
        data->SetBranchAddress("genWeight",         genWeight);
    }

    // muon branches
    int nMu = 0;
    auto muPt                    = new vector<float>();
    auto muEta                   = new vector<float>();
    auto muPhi                   = new vector<float>();
    auto muEn                    = new vector<float>();
    auto muD0                    = new vector<float>();
    auto muDz                    = new vector<float>();
    auto muBestTrkPtError        = new vector<float>();
    auto muBestTrkPt             = new vector<float>();
    auto muSIP                   = new vector<float>();
    auto muPFChIso03             = new vector<float>();
    auto muPFPhoIso03            = new vector<float>();
    auto muPFNeuIso03            = new vector<float>();
    auto muPFPUIso03             = new vector<float>();
    auto muCharge                = new vector<int>();
    auto muType                  = new vector<int>();
    auto muTrkLayers             = new vector<int>();
    auto muBestTrkType           = new vector<int>();
    auto muPixelHits             = new vector<int>();
    auto muStations              = new vector<int>();
    auto muMatches               = new vector<int>();
    data->SetBranchAddress("nMu",               nMu);
    data->SetBranchAddress("muPt",              muPt);
    data->SetBranchAddress("muEta",             muEta);
    data->SetBranchAddress("muPhi",             muPhi);
    data->SetBranchAddress("muEn",              muEn);
    data->SetBranchAddress("muD0",              muD0);
    data->SetBranchAddress("muDz",              muDz);
    data->SetBranchAddress("muBestTrkPtError",  muBestTrkPtError);
    data->SetBranchAddress("muBestTrkPt",       muBestTrkPt);
    data->SetBranchAddress("muSIP",             muSIP);
    data->SetBranchAddress("muPFChIso03",       muPFChIso03);
    data->SetBranchAddress("muPFPhoIso03",      muPFPhoIso03);
    data->SetBranchAddress("muPFNeuIso03",      muPFNeuIso03);
    data->SetBranchAddress("muPFPUIso03",       muPFPUIso03);
    data->SetBranchAddress("muCharge",          muCharge);
    data->SetBranchAddress("muType",            muType);
    data->SetBranchAddress("muTrkLayers",       muTrkLayers);
    data->SetBranchAddress("muBestTrkType",     muBestTrkType);
    data->SetBranchAddress("muPixelHits",       muPixelHits);
    data->SetBranchAddress("muStations",        muStations);
    data->SetBranchAddress("muMatches",         muMatches);

    // photon branches
    int nPho = 0;
    auto phoE                    = new vector<float>();
    auto phoEt                   = new vector<float>();
    auto phoCalibEt              = new vector<float>();
    auto phoEta                  = new vector<float>();
    auto phoPhi                  = new vector<float>();
    auto phoSCEta                = new vector<float>();
    auto phoSCPhi                = new vector<float>();
    auto phoIDMVA                = new vector<float>();
    auto phoEleVeto              = new vector<int>();
    // Hgg Photon ID 
    auto phoSCRawE               = new vector<float>();
    auto phoSigmaIEtaIEtaFull5x5 = new vector<float>();
    auto phoSCEtaWidth           = new vector<float>();
    auto phoSCPhiWidth           = new vector<float>();
    auto phoSigmaIEtaIPhiFull5x5 = new vector<float>();
    auto phoPFPhoIso             = new vector<float>();
    auto phoPFChIso              = new vector<float>();
    auto phoPFChWorstIso         = new vector<float>();
    auto phoE2x2Full5x5          = new vector<float>();
    auto phoE5x5Full5x5          = new vector<float>();
    auto phoESEffSigmaRR         = new vector<float>();
    auto phoESEnP1               = new vector<float>();
    auto phoESEnP2               = new vector<float>();
    auto phoTrkIsoHollowConeDR03 = new vector<float>();
    auto phoHoverE               = new vector<float>();
    // scale and smear uncertainties  
    auto phoScale_stat_up        = new vector<float>();
    auto phoScale_syst_up        = new vector<float>();
    auto phoScale_gain_up        = new vector<float>();
    auto phoScale_stat_dn        = new vector<float>();
    auto phoScale_syst_dn        = new vector<float>();
    auto phoScale_gain_dn        = new vector<float>();
    auto phoResol_phi_up         = new vector<float>();
    auto phoResol_rho_up         = new vector<float>();
    auto phoResol_rho_dn         = new vector<float>();
    // photon ss corrected branches
    auto phoCorrR9Full5x5        = new vector<float>();
    auto phoCorrHggIDMVA         = new vector<float>();
    data->SetBranchAddress("nPho",            nPho);
    data->SetBranchAddress("phoE",            phoE);
    data->SetBranchAddress("phoEt",           phoEt);
    data->SetBranchAddress("phoCalibEt",      phoCalibEt);
    data->SetBranchAddress("phoEta",          phoEta);
    data->SetBranchAddress("phoPhi",          phoPhi);
    data->SetBranchAddress("phoSCEta",        phoSCEta);
    data->SetBranchAddress("phoSCPhi",        phoSCPhi);
    data->SetBranchAddress("phoCalibEt",      phoCalibEt);
    data->SetBranchAddress("phoIDMVA",        phoIDMVA);
    // data->SetBranchAddress("phoR9Full5x5",    phoR9Full5x5);
    data->SetBranchAddress("phoEleVeto",      phoEleVeto);

    data->SetBranchAddress("phoSCRawE",                     phoSCRawE);
    data->SetBranchAddress("phoSigmaIEtaIEtaFull5x5",       phoSigmaIEtaIEtaFull5x5);
    data->SetBranchAddress("phoSCEtaWidth",                 phoSCEtaWidth);
    data->SetBranchAddress("phoSCPhiWidth",                 phoSCPhiWidth);
    data->SetBranchAddress("phoSigmaIEtaIPhiFull5x5",       phoSigmaIEtaIPhiFull5x5);
    data->SetBranchAddress("phoPFPhoIso",                   phoPFPhoIso);
    data->SetBranchAddress("phoPFChIso",                    phoPFChIso);
    data->SetBranchAddress("phoPFChWorstIso",               phoPFChWorstIso);
    data->SetBranchAddress("phoE2x2Full5x5",                phoE2x2Full5x5);
    data->SetBranchAddress("phoE5x5Full5x5",                phoE5x5Full5x5);
    data->SetBranchAddress("phoESEffSigmaRR",               phoESEffSigmaRR);
    data->SetBranchAddress("phoESEnP1",                     phoESEnP1);
    data->SetBranchAddress("phoESEnP2",                     phoESEnP2);
    data->SetBranchAddress("phoTrkIsoHollowConeDR03",       phoTrkIsoHollowConeDR03);
    data->SetBranchAddress("phoHoverE",                     phoHoverE);

    data->SetBranchAddress("phoScale_stat_up",phoScale_stat_up);
    data->SetBranchAddress("phoScale_syst_up",phoScale_syst_up);
    data->SetBranchAddress("phoScale_gain_up",phoScale_gain_up);
    data->SetBranchAddress("phoScale_stat_dn",phoScale_stat_dn);
    data->SetBranchAddress("phoScale_syst_dn",phoScale_syst_dn);
    data->SetBranchAddress("phoScale_gain_dn",phoScale_gain_dn);
    data->SetBranchAddress("phoResol_phi_up", phoResol_phi_up);
    data->SetBranchAddress("phoResol_rho_up", phoResol_rho_up);
    data->SetBranchAddress("phoResol_rho_dn", phoResol_rho_dn);
    if (data->HasMC() && !inpath.Contains("MuMuG") && !inpath.Contains("MuMuJet")){
        data->SetBranchAddress((variation == "PhoNoR9Corr") ? "phoR9Full5x5" : "phoCorrR9Full5x5", phoCorrR9Full5x5);
        data->SetBranchAddress("phoCorrHggIDMVA",  phoCorrHggIDMVA);
    }
    else{
        data->SetBranchAddress("phoR9Full5x5", phoCorrR9Full5x5);
        // data->SetBranchAddress("phoIDMVA", phoCorrHggIDMVA);
    }

    // jet branches
    int nJet;
    auto jetPt                  = new vector<float>();
    auto jetEta                 = new vector<float>();
    auto jetPhi                 = new vector<float>();
    auto jetEn                  = new vector<float>();
    auto jetNHF                 = new vector<float>(); // Neutral Hadron Energy Fraction
    auto jetNEF                 = new vector<float>(); // Neutral EM Energy Fraction
    auto jetID                  = new vector<float>();
    auto jetCHF                 = new vector<float>(); // Charged Hadron Energy Fraction
    auto jetCEF                 = new vector<float>(); // Charged EM Energy Fraction
    auto jetMUF                 = new vector<float>(); // Muon energy fraction
    auto jetNCH                 = new vector<int>();   // Charged Multplicity
    auto jetNNP                 = new vector<int>();   // Neutral Multplicity
    auto jetP4Smear             = new vector<float>();
    auto jetP4SmearUp           = new vector<float>();
    auto jetP4SmearDo           = new vector<float>();
    auto jetJECUnc              = new vector<float>();
    data->SetBranchAddress("nJet",          nJet);
    data->SetBranchAddress("jetPt",         jetPt);
    data->SetBranchAddress("jetEta",        jetEta);
    data->SetBranchAddress("jetPhi",        jetPhi);
    data->SetBranchAddress("jetEn",         jetEn);
    data->SetBranchAddress("jetNHF",        jetNHF);
    data->SetBranchAddress("jetNEF",        jetNEF);
    data->SetBranchAddress("jetNNP",        jetNNP);
    data->SetBranchAddress("jetID",         jetID);
    data->SetBranchAddress("jetCHF",        jetCHF);
    data->SetBranchAddress("jetCEF",        jetCEF);
    data->SetBranchAddress("jetMUF",        jetMUF);
    data->SetBranchAddress("jetNCH",        jetNCH);
    if (data->HasMC()){
        data->SetBranchAddress("jetP4Smear",        jetP4Smear);
        data->SetBranchAddress("jetP4SmearUp",      jetP4SmearUp);
        data->SetBranchAddress("jetP4SmearDo",      jetP4SmearDo);
        data->SetBranchAddress("jetJECUnc",         jetJECUnc);
    }

    // generator branches (only for MC)
    int nMC = 0;
    int nLHE = 0;
    auto mcPID                  = new vector<int>();
    auto mcMomPID               = new vector<int>();
    auto mcGMomPID              = new vector<int>();
    auto mcPt                   = new vector<float>();
    auto mcEta                  = new vector<float>();
    auto mcPhi                  = new vector<float>();
    auto mcMass                 = new vector<float>();
    auto mcStatusFlag           = new vector<unsigned short>();
    auto lheE                   = new vector<float>();
    auto lhePx                  = new vector<float>();
    auto lhePy                  = new vector<float>();
    auto lhePz                  = new vector<float>();
    auto lhePID                 = new vector<int>();
    if (data->HasMC()){
        data->SetBranchAddress("nMC",              nMC);
        data->SetBranchAddress("mcPID",            mcPID);
        data->SetBranchAddress("mcMomPID",         mcMomPID);
        data->SetBranchAddress("mcGMomPID",        mcGMomPID);
        data->SetBranchAddress("mcPt",             mcPt);
        data->SetBranchAddress("mcEta",            mcEta);
        data->SetBranchAddress("mcPhi",            mcPhi);
        data->SetBranchAddress("mcMass",           mcMass);
        data->SetBranchAddress("mcStatusFlag",     mcStatusFlag);
        data->SetBranchAddress("nLHE",             nLHE);
        data->SetBranchAddress("lheE",             lheE);
        data->SetBranchAddress("lhePx",            lhePx);
        data->SetBranchAddress("lhePy",            lhePy);
        data->SetBranchAddress("lhePz",            lhePz);
        data->SetBranchAddress("lhePID",           lhePID);
    }

    //========================================================//
    //         Set up the Rochester Muon PT Correction        //
    //========================================================//
    TString RochCor_file;
    if (era == "UL2016preVFP")  RochCor_file = "./plugin/roccor.Run2.v5/RoccoR2016aUL.txt";
    if (era == "UL2016postVFP") RochCor_file = "./plugin/roccor.Run2.v5/RoccoR2016bUL.txt";
    if (era == "UL2017")        RochCor_file = "./plugin/roccor.Run2.v5/RoccoR2017UL.txt";
    if (era == "UL2018")        RochCor_file = "./plugin/roccor.Run2.v5/RoccoR2018UL.txt";
    if (RochCor_file.IsNull()){
        printf("No suitable rochester file for era: %s\n", era.Data());
        gSystem->Exit(1);
    }
    auto rc = new RoccoR(RochCor_file.Data());
    gRandom->SetSeed(123456);

    //========================================================//
    //           Set up the Pileup Weights Calculater         //
    //========================================================//
    vector<PUWeightCalculator*> puCalc = {new PUWeightCalculator(), new PUWeightCalculator(), new PUWeightCalculator()};
    if (data->HasMC()){
        if (era == "UL2016preVFP"){
            puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_69200nb.root");
            puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_72383nb.root");
            puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_66016nb.root");
        }
        else if (era == "UL2016postVFP"){
            puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_69200nb.root");
            puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_72383nb.root");
            puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_66016nb.root");
        }
        else if (era == "UL2017"){
            puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_69200nb.root");
            puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_72383nb.root");
            puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_66016nb.root");
        }
        else if (era == "UL2018"){
            puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_69200nb.root");
            puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_72383nb.root");
            puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_66016nb.root");
        }
        else{
            printf("No suitable pileup file for era: %s\n", era.Data());
            gSystem->Exit(1);
        }
    }

    //========================================================//
    //      Set up the Scale factor files and histograms      //
    //========================================================//
    TString HZZID_file;
    if (era.Contains("UL2016")) HZZID_file = "/data4/chenghan/external/ScaleFactor/HZZID/final_HZZ_SF_2016UL_mupogsysts_newLoose.root";
    if (era == "UL2017")        HZZID_file = "/data4/chenghan/external/ScaleFactor/HZZID/final_HZZ_SF_2017UL_mupogsysts_newLoose.root";
    if (era == "UL2018")        HZZID_file = "/data4/chenghan/external/ScaleFactor/HZZID/final_HZZ_SF_2018UL_mupogsysts_newLoose.root";
    if (HZZID_file.IsNull()){
        printf("No suitable HZZID SFs file for era: %s\n", era.Data());
        gSystem->Exit(1);
    }
    TH2F* hHZZID    = GetSFsHisto2D(HZZID_file, "FINAL");
    TH2F* hHZZIDErr = GetSFsHisto2D(HZZID_file, "ERROR");

    // TString Fall17PhoID_file;
    // if (era == "UL2016preVFP")  Fall17PhoID_file = "/data4/chenghan/external/ScaleFactor/Fall17PhoID/egammaEffi.txt_EGM2D_Pho_wp90_UL16.root";
    // if (era == "UL2016postVFP") Fall17PhoID_file = "/data4/chenghan/external/ScaleFactor/Fall17PhoID/egammaEffi.txt_EGM2D_Pho_MVA90_UL16_postVFP.root";
    // if (era == "UL2017")        Fall17PhoID_file = "/data4/chenghan/external/ScaleFactor/Fall17PhoID/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root";
    // if (era == "UL2018")        Fall17PhoID_file = "/data4/chenghan/external/ScaleFactor/Fall17PhoID/egammaEffi.txt_EGM2D_Pho_wp90.root_UL18.root";
    // if (Fall17PhoID_file.IsNull()){
    //     printf("No suitable phoID SFs file for era: %s\n", era.Data());
    //     gSystem->Exit(1);
    // }
    // TH2F* hFall17PhoID = GetSFsHisto2D(Fall17PhoID_file, "EGamma_SF2D");

    // TString CSEV_file;
    // if (era == "UL2016preVFP")  CSEV_file = "/data4/chenghan/external/ScaleFactor/CSEV/CSEV_SummaryPlot_UL16_preVFP.root";
    // if (era == "UL2016postVFP") CSEV_file = "/data4/chenghan/external/ScaleFactor/CSEV/CSEV_SummaryPlot_UL16_postVFP.root";
    // if (era == "UL2017")        CSEV_file = "/data4/chenghan/external/ScaleFactor/CSEV/CSEV_SummaryPlot_UL17.root";
    // if (era == "UL2018")        CSEV_file = "/data4/chenghan/external/ScaleFactor/CSEV/CSEV_SummaryPlot_UL18.root";
    // if (CSEV_file.IsNull()){
    //     printf("No suitable CSEV SFs file for era: %s\n", era.Data());
    //     gSystem->Exit(1);
    // }
    // TH1F* hCSEV = GetSFsHisto1D(CSEV_file, "MVAID/SF_CSEV_MVAID");

    TString RECOMu_file;
    if (era == "UL2016preVFP")  RECOMu_file = "/data4/chenghan/external/ScaleFactor/RECOMu/UL2016preVFP/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt.json";
    if (era == "UL2016postVFP") RECOMu_file = "/data4/chenghan/external/ScaleFactor/RECOMu/UL2016postVFP/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt.json";
    if (era == "UL2017")        RECOMu_file = "/data4/chenghan/external/ScaleFactor/RECOMu/UL2017/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt.json";
    if (era == "UL2018")        RECOMu_file = "/data4/chenghan/external/ScaleFactor/RECOMu/UL2018/NUM_TrackerMuons_DEN_genTracks_Z_abseta_pt.json";
    if (RECOMu_file.IsNull()){
        printf("No suitable RECOMu SFs file for era: %s\n", era.Data());
        gSystem->Exit(1);
    }
    auto jr = new JsonSFsReader(RECOMu_file.Data(), "NUM_TrackerMuons_DEN_genTracks");

    TString MuPhoTrig_file;
    if (era == "UL2016preVFP")  MuPhoTrig_file = "/data4/chenghan/external/ScaleFactor/MuPhoTrig/Mu17Pho30SF_2016preVFP.root";
    if (era == "UL2016postVFP") MuPhoTrig_file = "/data4/chenghan/external/ScaleFactor/MuPhoTrig/Mu17Pho30SF_2016postVFP.root";
    if (era == "UL2017")        MuPhoTrig_file = "/data4/chenghan/external/ScaleFactor/MuPhoTrig/Mu17Pho30SF_2017.root";
    if (era == "UL2018")        MuPhoTrig_file = "/data4/chenghan/external/ScaleFactor/MuPhoTrig/Mu17Pho30SF_2018.root";
    if (MuPhoTrig_file.IsNull()){
        printf("No suitable MuPhoTrig SFs file for era: %s\n", era.Data());
        gSystem->Exit(1);
    }
    TH2F* hMuPhoTrigEB;
    TH2F* hMuPhoTrigEE;
    if (data->HasMC() && era.Contains("UL2016")){
        hMuPhoTrigEB = GetSFsHisto2D(MuPhoTrig_file, "Mu17Pho30_DataEff_EB");
        hMuPhoTrigEE = GetSFsHisto2D(MuPhoTrig_file, "Mu17Pho30_DataEff_EE");
    }
    else{
        hMuPhoTrigEB = GetSFsHisto2D(MuPhoTrig_file, "Mu17Pho30_SFs_EB");
        hMuPhoTrigEE = GetSFsHisto2D(MuPhoTrig_file, "Mu17Pho30_SFs_EE");
    }

    TString IsoMuTrig_file = "/data4/chenghan/external/ScaleFactor/IsoMuTrig/NUM_IsoMu27_DEN_HZZid_abseta_pt.root";
    TH2F* hIsoMuTrig = GetSFsHisto2D(IsoMuTrig_file, "NUM_IsoMu27_DEN_HZZid_abseta_pt");

    //========================================================//
    //          Set up the output file and minitree           //
    //========================================================//
    auto dirName = gSystem->DirName(outpath.Data());
    gSystem->Exec(Form("mkdir -p %s", dirName));

    TreeWritter* outTree = nullptr;
    if (variation == "Nominal"){
        outTree = new TreeWritter(outpath, data->HasMC());
        outTree->SetTree("miniTree");
    }
    else{
        outTree = new TreeWritter(outpath, data->HasMC(), "UPDATE");
        outTree->SetTree(Form("miniTree_%s", variation.Data()));
    }

    //========================================================//
    //    Set up the TMVA reader for Hgg Photon ID (data)     //
    //========================================================//
    // vector<TMVA::Experimental::RReader*> reader = {NULL, NULL};
    string phoID_path_EB = "";
    if (era == "UL2016preVFP")  phoID_path_EB = "./HggPhoID/PhoID_barrel_UL16_SA_TMVA_BDTG.weights.xml";
    if (era == "UL2016postVFP") phoID_path_EB = "./HggPhoID/PhoID_barrel_UL16_SA_TMVA_BDTG.weights.xml";
    if (era == "UL2017")        phoID_path_EB = "./HggPhoID/PhoID_barrel_UL2017_GJetMC_SATrain_nTree2k_LR_0p1_13052020_BDTG.weights.xml";
    if (era == "UL2018")        phoID_path_EB = "./HggPhoID/PhoID_barrel_UL18_GJetMC_SATrain_BDTG_nTree2k_BDTG.weights.xml";
    string phoID_path_EE = "";
    if (era == "UL2016preVFP")  phoID_path_EE = "./HggPhoID/PhoID_endcap_UL16_SA_TMVA_BDTG.weights.xml";
    if (era == "UL2016postVFP") phoID_path_EE = "./HggPhoID/PhoID_endcap_UL16_SA_TMVA_BDTG.weights.xml";
    if (era == "UL2017")        phoID_path_EE = "./HggPhoID/PhoID_endcap_UL2017_GJetMC_SATrain_nTree2k_LR_0p1_13052020_BDTG.weights.xml";
    if (era == "UL2018")        phoID_path_EE = "./HggPhoID/PhoID_endcap_UL18_GJetMC_SATrain_BDTG_nTree2k_BDTG.weights.xml";

    TMVA::Experimental::RReader readerEB(phoID_path_EB);
    TMVA::Experimental::RReader readerEE(phoID_path_EE);

    map<string, vector<vector<float>>> HggPreselSFs;
    HggPreselSFs["UL2016preVFP"]  = {{0.9984, 0.9978}, {1.0064 , 1.0049}};
    HggPreselSFs["UL2016postVFP"] = {{0.9984, 0.9978}, {1.0064 , 1.0049}};
    HggPreselSFs["UL2017"]        = {{0.9961, 0.9981}, {1.0054 , 1.0061}};
    HggPreselSFs["UL2018"]        = {{1.0017, 0.9973}, {1.0030, 1.0031}};

    map<string, vector<vector<float>>> HggPreselSFsUnc;
    HggPreselSFsUnc["UL2016preVFP"]  = {{0.0258, 0.0039}, {0.0192, 0.0024}};
    HggPreselSFsUnc["UL2016postVFP"] = {{0.0258, 0.0039}, {0.0192, 0.0024}};
    HggPreselSFsUnc["UL2017"]        = {{0.0307, 0.0057}, {0.0129, 0.0018}};
    HggPreselSFsUnc["UL2018"]        = {{0.0237, 0.0067}, {0.0210, 0.0024}};

    map<string, vector<vector<float>>> HggPhoIDSFs;
    HggPhoIDSFs["UL2016preVFP"]  = {{0.9992, 0.9989}, {1.0000 , 1.0005}};
    HggPhoIDSFs["UL2016postVFP"] = {{0.9992, 0.9989}, {1.0000 , 1.0005}};
    HggPhoIDSFs["UL2017"]        = {{1.0021, 1.0001}, {1.0061 , 1.0016}};
    HggPhoIDSFs["UL2018"]        = {{1.0022, 1.0005}, {1.0058 , 1.0013}};

    map<string, vector<vector<float>>> HggPhoIDSFsUnc;
    HggPhoIDSFsUnc["UL2016preVFP"]  = {{0.0025, 0.0034}, {0.0033, 0.0029}};
    HggPhoIDSFsUnc["UL2016postVFP"] = {{0.0025, 0.0034}, {0.0033, 0.0029}};
    HggPhoIDSFsUnc["UL2017"]        = {{0.0014, 0.0038}, {0.0014, 0.0023}};
    HggPhoIDSFsUnc["UL2018"]        = {{0.0012, 0.0033}, {0.0020, 0.0028}};

    map<string, vector<vector<float>>> HggCSEVSFs;
    HggCSEVSFs["UL2016preVFP"]  = {{1.0004, 0.9976}, {0.9882 , 0.9971}};
    HggCSEVSFs["UL2016postVFP"] = {{1.0004, 0.9976}, {0.9882 , 0.9971}};
    HggCSEVSFs["UL2017"]        = {{0.9838, 0.9913}, {0.9777 , 0.9777}};
    HggCSEVSFs["UL2018"]        = {{0.9830, 0.9939}, {0.9603 , 0.9754}};

    map<string, vector<vector<float>>> HggCSEVSFsUnc;
    HggCSEVSFsUnc["UL2016preVFP"]  = {{0.0021, 0.0005}, {0.0072, 0.0016}};
    HggCSEVSFsUnc["UL2016postVFP"] = {{0.0021, 0.0005}, {0.0072, 0.0016}};
    HggCSEVSFsUnc["UL2017"]        = {{0.0024, 0.0009}, {0.0180, 0.0026}};
    HggCSEVSFsUnc["UL2018"]        = {{0.0021, 0.0005}, {0.0077, 0.0017}};


    //========================================================//
    //  calculate the mc weight (corresponding to luminosity) //
    //========================================================//
    Long64_t totalEvents = 0;
    double mcwei = 1.;
    if (data->HasMC()){
        for(Long64_t ev = 0; ev < data->GetEntries(); ev++){
            data->GetEntry(ev);
            if (genWeight > 0)
                totalEvents++;
            else
                totalEvents--;
        }
        mcwei = xs * lumi / totalEvents;
        printf("INFO: mcwei = %f with xs = %f and lumi = %f\n", mcwei, xs, lumi);
    }

    //========================================================//
    //                     START: Analysis                    //
    //========================================================//
    Long64_t tev = data->GetEntries();
    // Long64_t tev = 100000;
    Long64_t showev = tev / 100;
    auto cutflow = new CutCalculater(tev);
    for (Long64_t ev = 0; ev < tev; ev++){
        data->GetEntry(ev);

        if (ev % showev == 0 && showev > 0)
            print_progress(ev, tev);

        if (!isPVGood)
            continue;
        cutflow->AddCut("Good vtx");

        // reconstruction selections
        bool passTrig = false;
        if (HLT == "SingleMuTrig")
            passTrig = ((HLTEleMuX >> 19) & 1) == 1;
        else
            passTrig = ((HLTEleMuX >> 8) & 1) == 1;

        if (!data->HasMC() && era.Contains("UL2016")){ // 2016 data
            if (!passTrig)
                continue;
        }
        else if (!era.Contains("UL2016")){ // 2017 2018
            if (!passTrig)
                continue;
        }
        cutflow->AddCut("HLT cut");

        if (nMu < 2 || nPho < 1)
            continue;
        cutflow->AddCut("Non zero particle");

        // https://github.com/CMS-LUMI-POG/ZCounting/blob/master/ZUtils/plugins/MuonPATUserDataRochesterCorrectionAdder.cc
        // https://gitlab.cern.ch/akhukhun/roccor
        vector<float> muCalibPt(nMu);
        for (int i = 0; i < nMu; i++){
            double muSF = 1.;
            double fRand = gRandom->Rndm();
            if (data->HasMC()){
                muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand);
                if (variation == "MuCalibStat")   muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand, 1, 99);
                if (variation == "MuCalibZpt")    muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand, 2, 0);
                if (variation == "MuCalibEwk")    muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand, 3, 0);
                if (variation == "MuCalibdeltaM") muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand, 4, 0);
                if (variation == "MuCalibEwk2")   muSF = rc->kSmearMC(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), muTrkLayers->at(i), fRand, 5, 0);
            }
            else
                muSF = rc->kScaleDT(muCharge->at(i), muPt->at(i), muEta->at(i), muPhi->at(i), 0, 0);
            
            // muSF = 1.; // for syncronization
            muCalibPt[i] = muSF * muPt->at(i);
        }

        vector<int> pass_muonID;
        pass_muonID.clear();
        for (int i = 0; i < nMu; i++){
            if (muCalibPt.at(i) < 4. || fabs(muEta->at(i)) > 2.4)
                continue;

            // HZZ muon ID
            // muType >> 1: Global muon
            // muType >> 2: Tracker muon
            // muType >> 5: Particle-flow muon
            bool HZZTightID = false;
            bool IsGlobalMuon = ((muType->at(i) >> 1) & 1) == 1;
            bool IsTrackerMuon = ((muType->at(i) >> 2) & 1) == 1;
            bool IsHZZLoose = (IsGlobalMuon || (IsTrackerMuon && muStations->at(i) >= 1)) &&
                                muBestTrkType->at(i) != 2 &&
                                fabs(muD0->at(i)) < 0.5 &&
                                fabs(muDz->at(i)) < 1.0 &&
                                muSIP->at(i) < 4.;
            // bool IsHZZLoose = (IsGlobalMuon || (IsTrackerMuon && muStations->at(i) >= 1)) &&
            //                     muBestTrkType->at(i) != 2
                                // fabs(muD0->at(i)) < 0.5 &&
                                // fabs(muDz->at(i)) < 1.0 &&
                                // muSIP->at(i) < 4.;
            bool IsPFMuon = ((muType->at(i) >> 5) & 1) == 1;
            if (muPt->at(i) < 200.)
                HZZTightID = (IsPFMuon && IsHZZLoose);
            else{
                bool IsTrkHighPtID = IsTrackerMuon &&
                                        muStations->at(i) > 1 &&
                                        muPixelHits->at(i) > 0 &&
                                        muTrkLayers->at(i) > 5 &&
                                        fabs(muD0->at(i)) < 0.2 &&
                                        fabs(muDz->at(i)) < 0.5 &&
                                        muBestTrkPtError->at(i) / muBestTrkPt->at(i) < 0.3;
                HZZTightID = (IsHZZLoose && (IsTrkHighPtID || IsPFMuon));
            }
            if (!HZZTightID)
                continue;

            // HZZ Iso
            float isolation = (muPFChIso03->at(i) + TMath::Max(0., muPFNeuIso03->at(i) + muPFPhoIso03->at(i) - 0.5 * muPFPUIso03->at(i))) / muCalibPt.at(i);
            if (isolation > 0.35)
                continue;

            pass_muonID.push_back(i);
        }

        if (pass_muonID.size() < 2)
            continue;
        ArgSort(muCalibPt, pass_muonID);
        cutflow->AddCut("Good muons");

        vector<int> pairMmm, idx_passmu1, idx_passmu2;
        pairMmm.clear();
        idx_passmu1.clear();
        idx_passmu2.clear();
        float mu1ptcut = 20., mu2ptcut = 4.; // default muon pt cut, except for 2017 Run B & C
        if (HLT == "SingleMuTrig"){
            mu1ptcut = 30.;
            mu2ptcut = 4.;
        }
        for (size_t i = 0; i < pass_muonID.size(); i++){
            if (muCalibPt.at(pass_muonID[i]) < mu1ptcut)
                continue;

            TLorentzVector tmpmu1;
            tmpmu1.SetPtEtaPhiM(muCalibPt.at(pass_muonID[i]), muEta->at(pass_muonID[i]), muPhi->at(pass_muonID[i]), 0.1057);
            for (size_t j = i + 1; j < pass_muonID.size(); j++){
                if (muPt->at(pass_muonID[j]) < mu2ptcut)
                    continue;

                // Require two opposite sign muons
                if (muCharge->at(pass_muonID[i]) == muCharge->at(pass_muonID[j]))
                    continue;

                TLorentzVector tmpmu2;
                tmpmu2.SetPtEtaPhiM(muCalibPt.at(pass_muonID[j]), muEta->at(pass_muonID[j]), muPhi->at(pass_muonID[j]), 0.1057);

                float tmp_Mmm = (tmpmu1 + tmpmu2).M();

                if (tmp_Mmm > 9.3 && tmp_Mmm < 9.7)
                    continue;
                if (tmp_Mmm > 2.9 && tmp_Mmm < 3.3)
                    continue;
                if (tmp_Mmm > 50.)
                    continue;

                pairMmm.push_back(tmp_Mmm);
                idx_passmu1.push_back(pass_muonID[i]);
                idx_passmu2.push_back(pass_muonID[j]);
            }
        }
        if (pairMmm.size() < 1)
            continue;
        cutflow->AddCut("Good muon pair");

        // Sort the Mmm
        for (size_t i = 0; i < pairMmm.size() - 1; i++){
            for (size_t j = 0; j < pairMmm.size() - i - 1; j++){
                if (pairMmm[j] > pairMmm[j + 1]){
                    float temp = pairMmm[j];
                    int temp1 = idx_passmu1[j];
                    int temp2 = idx_passmu2[j];

                    pairMmm[j] = pairMmm[j + 1];
                    idx_passmu1[j] = idx_passmu1[j + 1];
                    idx_passmu2[j] = idx_passmu2[j + 1];

                    pairMmm[j + 1] = temp;
                    idx_passmu1[j + 1] = temp1;
                    idx_passmu2[j + 1] = temp2;
                }
            }
        }

        vector<float> phoCalibEt_NEW(nPho);
        for (int i = 0; i < nPho; i++){
            phoCalibEt_NEW[i] = phoCalibEt->at(i);
            if (data->HasMC()){
                if (variation == "PhoScaleStatUp") phoCalibEt_NEW[i] = phoScale_stat_up->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoScaleSystUp") phoCalibEt_NEW[i] = phoScale_syst_up->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoScaleGainUp") phoCalibEt_NEW[i] = phoScale_gain_up->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoScaleStatDo") phoCalibEt_NEW[i] = phoScale_stat_dn->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoScaleSystDo") phoCalibEt_NEW[i] = phoScale_syst_dn->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoScaleGainDo") phoCalibEt_NEW[i] = phoScale_gain_dn->at(i) / TMath::CosH(phoEta->at(i));
                if (variation == "PhoSigmaPhiUp")  phoCalibEt_NEW[i] = phoResol_phi_up->at(i)  / TMath::CosH(phoEta->at(i));
                if (variation == "PhoSigmaRhoUp")  phoCalibEt_NEW[i] = phoResol_rho_up->at(i)  / TMath::CosH(phoEta->at(i));
                if (variation == "PhoSigmaRhoDo")  phoCalibEt_NEW[i] = phoResol_rho_dn->at(i)  / TMath::CosH(phoEta->at(i));
            }
        }

        // photon selection
        vector<int> acc_pho;
        acc_pho.clear();
        float phoEtcut = 33.; // default photon Et cut, except for 2017 Run B & C
        if (HLT == "SingleMuTrig")
            phoEtcut = 15.;

        for (int i = 0; i < nPho; i++){
            if (phoCalibEt_NEW.at(i) < phoEtcut)
                continue;
            if (phoEleVeto->at(i) == 0) // Conversion safe electron veto
                continue;

            bool isEB = fabs(phoSCEta->at(i)) < 1.4442;
            bool isEE = fabs(phoSCEta->at(i)) > 1.566 && fabs(phoSCEta->at(i)) < 2.5;
            if (!(isEB || isEE))
                continue;

            // Hgg preselection
            float phoEffArea = (fabs(phoSCEta->at(i)) > 0 && fabs(phoSCEta->at(i)) < 1.5) ? 0.16544 : 0.13212;
            float rhocorr = rhoAll * phoEffArea;
            float phoPFPhoIso_corr = TMath::Max(phoPFPhoIso->at(i) - rhocorr, (float) 0.);
            bool isHR9_EB = isEB && phoCorrR9Full5x5->at(i) >  0.85;
            bool isLR9_EB = isEB && phoCorrR9Full5x5->at(i) <= 0.85 && phoCorrR9Full5x5->at(i) > 0.5 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03->at(i) < 6 && phoSigmaIEtaIEtaFull5x5->at(i) < 0.015;
            bool isHR9_EE = isEE && phoCorrR9Full5x5->at(i) >  0.9;
            bool isLR9_EE = isEE && phoCorrR9Full5x5->at(i) <= 0.9  && phoCorrR9Full5x5->at(i) > 0.8 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03->at(i) < 6 && phoSigmaIEtaIEtaFull5x5->at(i) < 0.035;
            bool isAOD = phoHoverE->at(i) < 0.08 && (phoCorrR9Full5x5->at(i) > 0.8 || phoPFChIso->at(i) < 20. || (phoPFChIso->at(i)/phoCalibEt->at(i) < 0.3));
            bool isHgg = isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE);
            if (!isHgg)
                continue;

            // Hgg photon ID 
            float MVA = -999;
            if (!data->HasMC() || inpath.Contains("MuMuG") || inpath.Contains("MuMuJet")){
                // int iBE = (isEB) ? 0 : 1;
                vector<float> HggIDInputs = {
                    phoSCRawE->at(i),
                    phoCorrR9Full5x5->at(i),
                    phoSigmaIEtaIEtaFull5x5->at(i),
                    phoSCEtaWidth->at(i),
                    phoSCPhiWidth->at(i),
                    phoSigmaIEtaIPhiFull5x5->at(i),
                    phoE2x2Full5x5->at(i)/phoE5x5Full5x5->at(i),
                    phoPFPhoIso->at(i),
                    phoPFChIso->at(i),
                    phoPFChWorstIso->at(i),
                    phoSCEta->at(i),
                    rho
                };
                if (isEE){
                    HggIDInputs.push_back(phoESEffSigmaRR->at(i));
                    HggIDInputs.push_back((phoESEnP1->at(i) + phoESEnP2->at(i)) / phoSCRawE->at(i));
                } 
                MVA = (isEB) ? readerEB.Compute(HggIDInputs)[0] : readerEE.Compute(HggIDInputs)[0];
            }
            else
                MVA = phoCorrHggIDMVA->at(i);

            if (MVA < -0.9)
                continue;
            
            acc_pho.push_back(i);
        }

        if (acc_pho.size() < 1)
            continue;
        ArgSort(phoCalibEt_NEW, acc_pho);
        cutflow->AddCut("Good photon");

        // kinematic selection
        TLorentzVector mu1, mu2, dimu, pho, mmg;
        int idx_mu1 = -1, idx_mu2 = -1, idx_pho = -1;

        float ptratiocut = 33. / 110.;
        // if (HLT == "SingleMuTrig")
        //     ptratiocut = 30. / 110.;

        bool passkinematic = false;
        for (size_t i = 0; i < pairMmm.size(); i++){
            mu1.SetPtEtaPhiM(muCalibPt.at(idx_passmu1[i]), muEta->at(idx_passmu1[i]), muPhi->at(idx_passmu1[i]), 0.1057);
            mu2.SetPtEtaPhiM(muCalibPt.at(idx_passmu2[i]), muEta->at(idx_passmu2[i]), muPhi->at(idx_passmu2[i]), 0.1057);
            dimu = mu1 + mu2;

            for (size_t j = 0; j < acc_pho.size(); j++){
                pho.SetPtEtaPhiM(phoCalibEt_NEW.at(acc_pho[j]), phoEta->at(acc_pho[j]), phoPhi->at(acc_pho[j]), 0.);
                mmg = dimu + pho;

                if (mu1.DeltaR(pho) < 1.)
                    continue;
                if (mu2.DeltaR(pho) < 1.)
                    continue;
                if (mmg.M() < 110. || mmg.M() > 170.)
                    continue;
                if (phoCalibEt_NEW.at(acc_pho[j]) / mmg.M() < ptratiocut || dimu.Pt() / mmg.M() < ptratiocut)
                    continue;

                passkinematic = true;
                if (passkinematic){
                    idx_mu1 = idx_passmu1[i];
                    idx_mu2 = idx_passmu2[i];
                    idx_pho = acc_pho[j];
                    break;
                }
            }
            if (passkinematic)
                break;
        }
        if (!passkinematic)
            continue;
        cutflow->AddCut("kinematics selections");

        // VBF TAGGING
        vector<int> pass_JetID;
        pass_JetID.clear();
        vector<float> jetCalibPt(nJet);
        vector<float> jetCalibEn(nJet);
        for (int i = 0; i < nJet; i++){
            float jetPt_tmp = jetPt->at(i);
            float jetEn_tmp = jetEn->at(i);
            if (data->HasMC()){
                jetPt_tmp = jetP4Smear->at(i) * jetPt->at(i);
                jetEn_tmp = jetP4Smear->at(i) * jetEn->at(i);
                if (variation == "JERUp"){
                    jetPt_tmp = jetP4SmearUp->at(i) * jetPt->at(i);
                    jetEn_tmp = jetP4SmearUp->at(i) * jetEn->at(i);
                }
                if (variation == "JERDo"){
                    jetPt_tmp = jetP4SmearDo->at(i) * jetPt->at(i);
                    jetEn_tmp = jetP4SmearDo->at(i) * jetEn->at(i);
                }
                if (variation == "JECUp"){
                    jetPt_tmp = jetP4Smear->at(i) * jetPt->at(i) + jetJECUnc->at(i);
                    jetEn_tmp = jetP4Smear->at(i) * jetEn->at(i) + jetJECUnc->at(i);
                }
                if (variation == "JECDo"){
                    jetPt_tmp = jetP4Smear->at(i) * jetPt->at(i) - jetJECUnc->at(i);
                    jetEn_tmp = jetP4Smear->at(i) * jetEn->at(i) - jetJECUnc->at(i);
                }
            }
            jetCalibPt[i] = jetPt_tmp;
            jetCalibEn[i] = jetEn_tmp;

            if (jetPt_tmp < 30.)
                continue;
            if (fabs(jetEta->at(i)) > 4.7)
                continue;
            bool isJetID = JetIDUL(jetEta->at(i), jetNHF->at(i), jetNEF->at(i), jetNNP->at(i), jetNCH->at(i), jetCHF->at(i), jetCEF->at(i), jetMUF->at(i), era);
            if (!isJetID)
                continue;

            pass_JetID.push_back(i);
        }

        int vbf_tag = 0;
        TLorentzVector jet1, jet2, jj;
        vector<int> accJet;
        accJet.clear();
        float zepen_ = -999;
        if (pass_JetID.size() > 1){
            for (size_t i = 0; i < pass_JetID.size(); i++){
                jet1.SetPtEtaPhiE(jetCalibPt.at(pass_JetID[i]), jetEta->at(pass_JetID[i]), jetPhi->at(pass_JetID[i]), jetEn->at(pass_JetID[i]));

                if (jet1.DeltaR(mu1) < 0.4 || jet1.DeltaR(mu2) < 0.4 || jet1.DeltaR(pho) < 0.4)
                    continue;

                for (size_t j = i + 1; j < pass_JetID.size(); j++){
                    jet2.SetPtEtaPhiE(jetCalibPt.at(pass_JetID[j]), jetEta->at(pass_JetID[j]), jetPhi->at(pass_JetID[j]), jetEn->at(pass_JetID[j]));

                    if (jet2.DeltaR(mu1) < 0.4 || jet2.DeltaR(mu2) < 0.4 || jet2.DeltaR(pho) < 0.4)
                        continue;

                    jj = jet1 + jet2;
                    float zepen = fabs(mmg.Eta() - ((jet1.Eta() + jet2.Eta()) * 0.5));
                    if (fabs(zepen) > 2.5)
                        continue;
                    if (fabs(jet1.Eta() - jet2.Eta()) < 3.5)
                        continue;
                    if (fabs(jj.DeltaPhi(mmg)) < 2.4)
                        continue;

                    // if (jj.M() <= 500. && jj.M() > 360){
                    //     vbf_tag = 2;
                    //     accJet.push_back(pass_JetID[i]);
                    //     accJet.push_back(pass_JetID[j]);
                    //     zepen_ = zepen;
                    //     break;
                    // }
                    if (jj.M() > 500){
                        vbf_tag = 1;
                        accJet.push_back(pass_JetID[i]);
                        accJet.push_back(pass_JetID[j]);
                        zepen_ = zepen;
                        break;
                    }
                }
                if (accJet.size() > 1)
                    break;
            }
        }

        // int hovbfTag_ = (vbf_tag == 1) ? 1 : 0;
        // int lovbfTag_ = (vbf_tag == 2) ? 1 : 0;
        int vbfTag_ = (vbf_tag == 1) ? 1 : 0;
        int bstTag_ = (mmg.Pt() > 60.) ? 1 : 0;

        int isEB_ = (fabs(phoSCEta->at(idx_pho)) < 1.4442) ? 1 : 0;
        int isEE_ = ((fabs(phoSCEta->at(idx_pho)) > 1.566) && (fabs(phoSCEta->at(idx_pho)) < 2.5)) ? 1 : 0;
        int isHR9_ = (phoCorrR9Full5x5->at(idx_pho) > 0.96) ? 1 : 0;
        int isLR9_ = (phoCorrR9Full5x5->at(idx_pho) < 0.96) ? 1 : 0;

        int cat_ = 0;
        if(HLT != "SingleMuTrig"){
            // if(dimu.M() < 0.7){
            //     if(isEB_ == 1 && isHR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 1;
            //     if(isEB_ == 1 && isLR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 2;
            //     if(isEE_ == 1                && vbfTag_ == 0 && bstTag_ == 0) cat_ = 3;
            //     if(                             vbfTag_ == 1                ) cat_ = 4;
            //     if(                             vbfTag_ == 0 && bstTag_ == 1) cat_ = 5;
            // }
            if(dimu.M() < 9.){
                if(isEB_ == 1 && isHR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 1;
                if(isEB_ == 1 && isLR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 2;
                if(isEE_ == 1                && vbfTag_ == 0 && bstTag_ == 0) cat_ = 3;
                if(                             vbfTag_ == 1                ) cat_ = 4;
                if(                             vbfTag_ == 0 && bstTag_ == 1) cat_ = 5;
            }
            else if(dimu.M() < 25.){
                if(isEB_ == 1 && isHR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 6;
                if(isEB_ == 1 && isLR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 7;
                if(isEE_ == 1                && vbfTag_ == 0 && bstTag_ == 0) cat_ = 8;
                if(                             vbfTag_ == 1                ) cat_ = 9;
                if(                             vbfTag_ == 0 && bstTag_ == 1) cat_ = 10;
            }
            else if(dimu.M() < 50.) cat_ = 11;
            else cat_ = -1;
        }
        else {
            if(dimu.M() < 9.){
                if(isEB_ == 1 && isHR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 12;
                if(isEB_ == 1 && isLR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 13;
                if(isEE_ == 1                && vbfTag_ == 0 && bstTag_ == 0) cat_ = 14;
                if(                             vbfTag_ == 1                ) cat_ = 15;
                if(                             vbfTag_ == 0 && bstTag_ == 1) cat_ = 16;
            }
            else if(dimu.M() < 50.){
                if(isEB_ == 1 && isHR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 17;
                if(isEB_ == 1 && isLR9_ == 1 && vbfTag_ == 0 && bstTag_ == 0) cat_ = 18;
                if(isEE_ == 1                && vbfTag_ == 0 && bstTag_ == 0) cat_ = 19;
                if(                             vbfTag_ == 1                ) cat_ = 20;
                if(                             vbfTag_ == 0 && bstTag_ == 1) cat_ = 21;
            }
            else cat_ = -1;
        }

        // Fill the value to outTree
        outTree->category               = cat_;
        outTree->CMS_higgs_mass_        = mmg.M();

        // dataset variables
        outTree->event_                 = event;
        outTree->run_                   = run;
        outTree->lumis_                 = lumis;

        // physics objects
        outTree->mu1_                   = mu1;
        outTree->mu2_                   = mu2;
        outTree->pho_                   = pho;
        outTree->jet1_                  = jet1;
        outTree->jet2_                  = jet2;
        outTree->jj_                    = jj;
        outTree->H_                     = mmg;
        outTree->dimu_                  = dimu;

        // photon variables
        outTree->phoSCEta_              = phoSCEta->at(idx_pho);
        outTree->phoSCPhi_              = phoSCPhi->at(idx_pho);
        outTree->phoCorrR9Full5x5_      = phoCorrR9Full5x5->at(idx_pho);

        // kinematics:
        float dRmu1pho = mu1.DeltaR(pho), dRmu2pho = mu2.DeltaR(pho);
        outTree->dRmumu_                = mu1.DeltaR(mu2);
        outTree->dRmax_                 = (dRmu1pho > dRmu2pho) ? dRmu1pho : dRmu2pho;
        outTree->dRmin_                 = (dRmu1pho > dRmu2pho) ? dRmu2pho : dRmu1pho;
        outTree->muPtRatio_             = mu1.Pt() / mu2.Pt();
        outTree->phoEtRatio_            = dimu.Pt() / pho.Pt();

        // jet variables
        if (accJet.size() > 1){
            outTree->dPhi_              = jj.DeltaPhi(mmg);
            outTree->dEta_              = fabs(jet2.Eta() - jet1.Eta());
            outTree->jjphi_             = jj.Phi();
            outTree->zep_               = zepen_;
            outTree->mjj_               = jj.M();
            outTree->dRjgmax_           = (jet1.DeltaR(pho) > jet2.DeltaR(pho)) ? jet1.DeltaR(pho) : jet2.DeltaR(pho);
            outTree->dRjgmin_           = (jet1.DeltaR(pho) < jet2.DeltaR(pho)) ? jet1.DeltaR(pho) : jet2.DeltaR(pho);

            double dRjm[4] = {jet1.DeltaR(mu1), jet1.DeltaR(mu2), jet2.DeltaR(mu1), jet2.DeltaR(mu2)};
            outTree->dRjmumax_          = (float) TMath::MaxElement(4, dRjm);
            outTree->dRjmumin_          = (float) TMath::MinElement(4, dRjm);
        }

        // weight evaluation
        if (data->HasMC()){
            outTree->mcwei_             = mcwei;
            outTree->puwei_             = (float) puCalc[0]->GetWeight(run, puTrue->at(1));
            outTree->puwei_up_          = (float) puCalc[1]->GetWeight(run, puTrue->at(1));
            outTree->puwei_down_        = (float) puCalc[2]->GetWeight(run, puTrue->at(1));
            outTree->genwei_            = (genWeight > 0) ? 1. : -1.;
            outTree->L1PF_              = L1ECALPrefire;
            outTree->L1PFUp_            = L1ECALPrefireUp;
            outTree->L1PFDo_            = L1ECALPrefireDown;
            outTree->MuPF_              = MuonPrefire;
            outTree->MuPFUp_            = MuonPrefireUp;
            outTree->MuPFDo_            = MuonPrefireDown;
            outTree->RECOMu1SF_         = jr->GetScaleFactor(fabs(mu1.Eta()), mu1.Pt());
            outTree->RECOMu1SFUp_       = outTree->RECOMu1SF_ + jr->GetScaleFactorErr(fabs(mu1.Eta()), mu1.Pt());
            outTree->RECOMu1SFDo_       = outTree->RECOMu1SF_ - jr->GetScaleFactorErr(fabs(mu1.Eta()), mu1.Pt());
            outTree->RECOMu2SF_         = jr->GetScaleFactor(fabs(mu2.Eta()), mu2.Pt());
            outTree->RECOMu2SFUp_       = outTree->RECOMu2SF_ + jr->GetScaleFactorErr(fabs(mu2.Eta()), mu2.Pt());
            outTree->RECOMu2SFDo_       = outTree->RECOMu2SF_ - jr->GetScaleFactorErr(fabs(mu2.Eta()), mu2.Pt());

            if (isEB_ == 1 && HLT != "SingleMuTrig"){
                int MuPhoTrigEBSF_bin = FindHistBin(hMuPhoTrigEB, mu1.Pt(), pho.Pt());
                outTree->HLTSF_         = hMuPhoTrigEB->GetBinContent(MuPhoTrigEBSF_bin);
                outTree->HLTSFUp_       = hMuPhoTrigEB->GetBinContent(MuPhoTrigEBSF_bin) + hMuPhoTrigEB->GetBinError(MuPhoTrigEBSF_bin);
                outTree->HLTSFDo_       = hMuPhoTrigEB->GetBinContent(MuPhoTrigEBSF_bin) - hMuPhoTrigEB->GetBinError(MuPhoTrigEBSF_bin);
            }
            if (isEB_ == 0 && HLT != "SingleMuTrig"){
                int MuPhoTrigEESF_bin = FindHistBin(hMuPhoTrigEE, mu1.Pt(), pho.Pt());
                outTree->HLTSF_         = hMuPhoTrigEE->GetBinContent(MuPhoTrigEESF_bin);
                outTree->HLTSFUp_       = hMuPhoTrigEE->GetBinContent(MuPhoTrigEESF_bin) + hMuPhoTrigEE->GetBinError(MuPhoTrigEESF_bin);
                outTree->HLTSFDo_       = hMuPhoTrigEE->GetBinContent(MuPhoTrigEESF_bin) - hMuPhoTrigEE->GetBinError(MuPhoTrigEESF_bin);
            }
            if (HLT == "SingleMuTrig"){
                int IsoMuTrigSF_bin = FindHistBin(hIsoMuTrig, fabs(mu1.Eta()), mu1.Pt());
                outTree->HLTSF_         = hIsoMuTrig->GetBinContent(IsoMuTrigSF_bin);
                outTree->HLTSFUp_       = hIsoMuTrig->GetBinContent(IsoMuTrigSF_bin) + hIsoMuTrig->GetBinError(IsoMuTrigSF_bin);
                outTree->HLTSFDo_       = hIsoMuTrig->GetBinContent(IsoMuTrigSF_bin) - hIsoMuTrig->GetBinError(IsoMuTrigSF_bin);
            }

            int HZZMu1IDSF_bin = FindHistBin(hHZZID, fabs(mu1.Eta()), mu1.Pt());
            int HZZMu1IDSFErr_bin = FindHistBin(hHZZIDErr, fabs(mu1.Eta()), mu1.Pt());
            outTree->HZZMu1IDSF_        = hHZZID->GetBinContent(HZZMu1IDSF_bin);
            outTree->HZZMu1IDSFUp_      = hHZZID->GetBinContent(HZZMu1IDSF_bin) + hHZZIDErr->GetBinContent(HZZMu1IDSFErr_bin);
            outTree->HZZMu1IDSFDo_      = hHZZID->GetBinContent(HZZMu1IDSF_bin) - hHZZIDErr->GetBinContent(HZZMu1IDSFErr_bin);

            int HZZMu2IDSF_bin = hHZZID->FindBin(fabs(mu2.Eta()), mu2.Pt());
            int HZZMu2IDSFErr_bin = hHZZIDErr->FindBin(fabs(mu2.Eta()), mu2.Pt());
            outTree->HZZMu2IDSF_        = hHZZID->GetBinContent(HZZMu2IDSF_bin);
            outTree->HZZMu2IDSFUp_      = hHZZID->GetBinContent(HZZMu2IDSF_bin) + hHZZIDErr->GetBinContent(HZZMu2IDSFErr_bin);
            outTree->HZZMu2IDSFDo_      = hHZZID->GetBinContent(HZZMu2IDSF_bin) - hHZZIDErr->GetBinContent(HZZMu2IDSFErr_bin);

            // int Fall17PhoIDSF_bin = hFall17PhoID->FindBin(phoSCEta->at(idx_pho), phoEt->at(idx_pho));
            // outTree->Fall17PhoIDSF_     = hFall17PhoID->GetBinContent(Fall17PhoIDSF_bin);
            // outTree->Fall17PhoIDSFUp_   = hFall17PhoID->GetBinContent(Fall17PhoIDSF_bin) + hFall17PhoID->GetBinError(Fall17PhoIDSF_bin);
            // outTree->Fall17PhoIDSFDo_   = hFall17PhoID->GetBinContent(Fall17PhoIDSF_bin) - hFall17PhoID->GetBinError(Fall17PhoIDSF_bin);

            // float SF = 1., SF_err = 1.;
            int eta_bin = -1;
            int r9_bin = -1;
            if (fabs(phoSCEta->at(idx_pho)) > 0 && fabs(phoSCEta->at(idx_pho)) < 1.5)
                eta_bin = 0;
            if (fabs(phoSCEta->at(idx_pho)) > 1.5 && fabs(phoSCEta->at(idx_pho)) < 999)
                eta_bin = 1;
            if (eta_bin == 0 && phoCorrR9Full5x5->at(idx_pho) < 0.85)
                r9_bin = 0;
            if (eta_bin == 0 && phoCorrR9Full5x5->at(idx_pho) > 0.85)
                r9_bin = 1;
            if (eta_bin == 1 && phoCorrR9Full5x5->at(idx_pho) < 0.9)
                r9_bin = 0;
            if (eta_bin == 1 && phoCorrR9Full5x5->at(idx_pho) > 0.9)
                r9_bin = 1;

            string era_str = era.Data();

            outTree->HggPhoIDSF_   = HggPhoIDSFs[era_str][eta_bin][r9_bin];
            outTree->HggPhoIDSFUp_ = HggPhoIDSFs[era_str][eta_bin][r9_bin] + HggPhoIDSFsUnc[era_str][eta_bin][r9_bin];
            outTree->HggPhoIDSFDo_ = HggPhoIDSFs[era_str][eta_bin][r9_bin] + HggPhoIDSFsUnc[era_str][eta_bin][r9_bin];

            outTree->HggPreselSF_   = HggPreselSFs[era_str][eta_bin][r9_bin];
            outTree->HggPreselSFUp_ = HggPreselSFs[era_str][eta_bin][r9_bin] + HggPreselSFsUnc[era_str][eta_bin][r9_bin];
            outTree->HggPreselSFDo_ = HggPreselSFs[era_str][eta_bin][r9_bin] + HggPreselSFsUnc[era_str][eta_bin][r9_bin];
            
            outTree->HggCSEVSF_   = HggCSEVSFs[era_str][eta_bin][r9_bin];
            outTree->HggCSEVSFUp_ = HggCSEVSFs[era_str][eta_bin][r9_bin] + HggCSEVSFsUnc[era_str][eta_bin][r9_bin];
            outTree->HggCSEVSFDo_ = HggCSEVSFs[era_str][eta_bin][r9_bin] + HggCSEVSFsUnc[era_str][eta_bin][r9_bin];


            outTree->weight_                = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_MuIDUp_         = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SFUp_   * outTree->RECOMu2SFUp_    * outTree->HZZMu1IDSFUp_   * outTree->HZZMu2IDSFUp_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_MuIDDo_         = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SFDo_   * outTree->RECOMu2SFDo_    * outTree->HZZMu1IDSFDo_   * outTree->HZZMu2IDSFDo_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_PhoIDUp_         = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSFUp_  * outTree->HggPreselSFUp_   * outTree->HggPreselSFUp_;

            outTree->weight_PhoIDDo_         = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSFDo_  * outTree->HggPreselSFDo_   * outTree->HggPreselSFDo_;

            outTree->weight_L1PFUp_         = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PFUp_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_L1PFDo_          = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PFDo_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_MuPFUp_          = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPFUp_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_MuPFDo_          = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPFDo_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;
            
            outTree->weight_HLTUp_           = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSFUp_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_HLTDo_           = outTree->mcwei_ * outTree->puwei_      * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSFDo_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_puweiUp_        = outTree->mcwei_ * outTree->puwei_up_    * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;

            outTree->weight_puweiDo_        = outTree->mcwei_ * outTree->puwei_down_    * outTree->genwei_ * outTree->L1PF_   * outTree->MuPF_   * outTree->HLTSF_ * outTree->RECOMu1SF_   * outTree->RECOMu2SF_    * outTree->HZZMu1IDSF_   * outTree->HZZMu2IDSF_   * outTree->HggPhoIDSF_  * outTree->HggPreselSF_   * outTree->HggPreselSF_;
            
        }

        outTree->FillTree();
    }
    cout << endl;
    cutflow->Print();

    delete data;
    delete outTree;
    delete cutflow;
    delete rc;
    delete hHZZID;
    delete hHZZIDErr;
    delete jr;
    delete hMuPhoTrigEB;
    delete hMuPhoTrigEE;
    delete hIsoMuTrig;
    for (size_t i = 0; i < puCalc.size(); i++){
        delete puCalc[i];
    }

    time.Stop();
    time.Print();
    cout << endl;
}