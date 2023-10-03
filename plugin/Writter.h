#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TString.h"


class TreeWritter{
public:
    TreeWritter(TString fName_, bool isMC_, TString option = "RECREATE"){
        isMC = isMC_;
        fout = new TFile(fName_, option.Data());
        printf("Save_File(): %s\n", fName_.Data());
    }
    ~TreeWritter(){
        // fout->cd();
        // https://root.cern.ch/doc/master/classTTree.html#a4680b0dfd17292acc29ba9a8f33788a3
        fout = tout->GetCurrentFile();
        fout->Write("", TObject::kOverwrite);
        fout->Close();

        delete fout;
    };

    void SetTree(TString treeName_){
        tout = new TTree(treeName_, "miniTree for H->mmg");

        tout->Branch("category",            &category);
        tout->Branch("CMS_higgs_mass",      &CMS_higgs_mass_);

        tout->Branch("event",               &event_);
        tout->Branch("run",                 &run_);
        tout->Branch("lumis",               &lumis_);

        tout->Branch("mu1",                 &mu1_);
        tout->Branch("mu2",                 &mu2_);
        tout->Branch("pho",                 &pho_);
        tout->Branch("jet1",                &jet1_);
        tout->Branch("jet2",                &jet2_);
        tout->Branch("dimu",                &dimu_);
        tout->Branch("jj",                  &jj_);
        tout->Branch("H",                   &H_);

        tout->Branch("phoSCEta",            &phoSCEta_);
        tout->Branch("phoSCPhi",            &phoSCPhi_);
        tout->Branch("phoCorrR9Full5x5",    &phoCorrR9Full5x5_);

        tout->Branch("dRmumu",              &dRmumu_);
        tout->Branch("dRmax",               &dRmax_);
        tout->Branch("dRmin",               &dRmin_);
        tout->Branch("muPtRatio",           &muPtRatio_);
        tout->Branch("phoEtRatio",          &phoEtRatio_);

        tout->Branch("dPhi",                &dPhi_);
        tout->Branch("dEta",                &dEta_);
        tout->Branch("mjj",                 &mjj_);
        tout->Branch("jjphi",               &jjphi_);
        tout->Branch("zep",                 &zep_);
        tout->Branch("dRjgmax",             &dRjgmax_);
        tout->Branch("dRjgmin",             &dRjgmin_);
        tout->Branch("dRjmumax",            &dRjmumax_);
        tout->Branch("dRjmumin",            &dRjmumin_);

        tout->Branch("weight",              &weight_);

        if (isMC){
            tout->Branch("mcwei",           &mcwei_);
            tout->Branch("puwei",           &puwei_);
            tout->Branch("puwei_up",        &puwei_up_);
            tout->Branch("puwei_down",      &puwei_down_);
            tout->Branch("genwei",          &genwei_);
            tout->Branch("L1PF",            &L1PF_);
            tout->Branch("L1PFUp",          &L1PFUp_);
            tout->Branch("L1PFDo",          &L1PFDo_);
            tout->Branch("MuPF",            &MuPF_);
            tout->Branch("MuPFUp",          &MuPFUp_);
            tout->Branch("MuPFDo",          &MuPFDo_);
            tout->Branch("HZZMu1IDSF",      &HZZMu1IDSF_);
            tout->Branch("HZZMu1IDSFUp",    &HZZMu1IDSFUp_);
            tout->Branch("HZZMu1IDSFDo",    &HZZMu1IDSFDo_);
            tout->Branch("HZZMu2IDSF",      &HZZMu2IDSF_);
            tout->Branch("HZZMu2IDSFUp",    &HZZMu2IDSFUp_);
            tout->Branch("HZZMu2IDSFDo",    &HZZMu2IDSFDo_);
            tout->Branch("RECOMu1SF",       &RECOMu1SF_);
            tout->Branch("RECOMu1SFUp",     &RECOMu1SFUp_);
            tout->Branch("RECOMu1SFDo",     &RECOMu1SFDo_);
            tout->Branch("RECOMu2SF",       &RECOMu2SF_);
            tout->Branch("RECOMu2SFUp",     &RECOMu2SFUp_);
            tout->Branch("RECOMu2SFDo",     &RECOMu2SFDo_);
            tout->Branch("HggPhoIDSF",      &HggPhoIDSF_);
            tout->Branch("HggPhoIDSFUp",    &HggPhoIDSFUp_);
            tout->Branch("HggPhoIDSFDo",    &HggPhoIDSFDo_);
            tout->Branch("HggPreselSF",      &HggPreselSF_);
            tout->Branch("HggPreselSFUp",    &HggPreselSFUp_);
            tout->Branch("HggPreselSFDo",    &HggPreselSFDo_);
            tout->Branch("HggCSEVSF",      &HggCSEVSF_);
            tout->Branch("HggCSEVSFUp",    &HggCSEVSFUp_);
            tout->Branch("HggCSEVSFDo",    &HggCSEVSFDo_);
            tout->Branch("HLTSF",           &HLTSF_);
            tout->Branch("HLTSFUp",         &HLTSFUp_);
            tout->Branch("HLTSFDo",         &HLTSFDo_);
            tout->Branch("weight_MuIDUp",      &weight_MuIDUp_);
            tout->Branch("weight_MuIDDo",      &weight_MuIDDo_);
            tout->Branch("weight_L1PFUp",      &weight_L1PFUp_);
            tout->Branch("weight_L1PFDo",      &weight_L1PFDo_);
            tout->Branch("weight_PhoIDUp",      &weight_PhoIDUp_);
            tout->Branch("weight_PhoIDDo",      &weight_PhoIDDo_);
            tout->Branch("weight_MuPFUp",      &weight_MuPFUp_);
            tout->Branch("weight_MuPFDo",      &weight_MuPFDo_);
            tout->Branch("weight_HLTUp",       &weight_HLTUp_);
            tout->Branch("weight_HLTDo",       &weight_HLTDo_);
            tout->Branch("weight_puweiUp",      &weight_puweiUp_);
            tout->Branch("weight_puweiDo",      &weight_puweiDo_);

        }
    }

    void FillTree(){
        tout->Fill();
    }

    int category = 0;

    // Higgs mass
    float CMS_higgs_mass_ = -999;

    // dataset variables
    int run_ = -999, lumis_ = -999;
    Long64_t event_ = -999;
    float puTrue_ = -999;

    // physics objects
    TLorentzVector mu1_, mu2_, pho_, jet1_, jet2_, jj_;
    TLorentzVector H_, dimu_;

    // photon variables
    float phoSCEta_ = -999, phoSCPhi_ = -999;
    float phoCorrR9Full5x5_ = -999;

    // kinematics
    float dRmumu_ = -999, dRmax_ = -999, dRmin_ = -999;
    float muPtRatio_ = -999, phoEtRatio_ = -999;

    // jet variables
    float dPhi_ = -999, dEta_ = -999;
    float mjj_ = -999, jjphi_ = -999;
    float zep_ = -999;
    float dRjgmax_ = -999, dRjgmin_ = -999, dRjmumax_ = -999, dRjmumin_ = -999;

    // weights
    float mcwei_ = -999.;
    float puwei_ = -999, puwei_up_ = -999, puwei_down_ = -999;
    float genwei_ = -999;
    float L1PF_ = -999, L1PFUp_ = -999, L1PFDo_ = -999;
    float MuPF_ = -999, MuPFUp_ = -999, MuPFDo_ = -999;
    float RECOMu1SF_ = -999., RECOMu1SFUp_ = -999., RECOMu1SFDo_ = -999.;
    float RECOMu2SF_ = -999., RECOMu2SFUp_ = -999., RECOMu2SFDo_ = -999.;
    float HZZMu1IDSF_ = -999, HZZMu1IDSFUp_ = -999, HZZMu1IDSFDo_ = -999;
    float HZZMu2IDSF_ = -999, HZZMu2IDSFUp_ = -999, HZZMu2IDSFDo_ = -999;
    float HggPhoIDSF_ = -999, HggPhoIDSFUp_ = -999, HggPhoIDSFDo_ = -999;
    float HggPreselSF_ = -999, HggPreselSFUp_ = -999, HggPreselSFDo_ = -999;
    float HggCSEVSF_ = -999, HggCSEVSFUp_ = -999, HggCSEVSFDo_ = -999;
    float HLTSF_ = -999, HLTSFUp_ = -999, HLTSFDo_ = -999;
    

    float weight_ = 1.;
    float weight_MuIDUp_ = 1., weight_MuIDDo_ = 1.;
    float weight_L1PFUp_ = 1., weight_L1PFDo_ = 1.;
    float weight_PhoIDUp_ = 1., weight_PhoIDDo_ = 1.;
    float weight_MuPFUp_ = 1., weight_MuPFDo_ = 1.;
    float weight_HLTUp_ = 1., weight_HLTDo_ = 1.;
    float weight_puweiUp_ = 1., weight_puweiDo_ = 1.;


private:
    TFile* fout = nullptr;
    TTree* tout = nullptr;
    bool isMC = false;
};