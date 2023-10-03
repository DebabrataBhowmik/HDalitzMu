#include <iostream>
#include <string>
#include <vector>
#include "boost/filesystem.hpp"
#include "yaml-cpp/yaml.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "./plugin/tdrstyle.h"
#include "./plugin/CMS_lumi_modification.h"


using namespace std;
namespace bfs = boost::filesystem;

void draw_var(ROOT::RDF::RNode df_data, ROOT::RDF::RNode df_MuMuG, ROOT::RDF::RNode df_MuMuJet, string varName, string sel, vector<double> bins, string legName){
    // preprocess mc dataframe
    auto nf_MuMuG = df_MuMuG;
    if (!nf_MuMuG.HasColumn(varName))
        nf_MuMuG = nf_MuMuG.Define(varName, sel);
    
    auto nf_MuMuJet = df_MuMuJet;
    if (!nf_MuMuJet.HasColumn(varName))
        nf_MuMuJet = nf_MuMuJet.Define(varName, sel);

    // preprocess data dataframe
    auto nf_data = df_data;
    if (!nf_data.HasColumn(varName))
        nf_data = nf_data.Define(varName, sel);
    // nf_data = nf_data.Define(varName, sel).Filter("(higgsMass > 110 && higgsMass <= 120) || (higgsMass >= 130 && higgsMass < 170)");

    auto h_MuMuG = nf_MuMuG.Histo1D({Form("h_%s_MuMuG", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_MuMuJet = nf_MuMuJet.Histo1D({Form("h_%s_MuMuJet", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_data = nf_data.Histo1D({Form("h_%s_data", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName);
    // cout << h_data->Integral()/h_mc->Integral() << endl;
    // h_mc->Scale(h_data->Integral()/h_mc->Integral());
    

    setTDRStyle();
    auto canv = new TCanvas("c", "c", 800, 800);
    canv->cd();

    TPad* pad1 = new TPad("pad1", " ", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.05);
    pad1->SetTopMargin(0.08);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.03);
    pad1->Draw();    
    pad1->SetLogy();         
    pad1->cd();

    // h_data->GetXaxis()->SetTitle(legName.c_str());
    // h_data->GetYaxis()->SetTitle(Form("Events / %.2f", (bins[2]-bins[1])/bins[0]));
    // h_data->GetXaxis()->SetTitleSize(0.04);
    // h_data->GetXaxis()->SetLabelSize(0.04);
    // h_data->GetXaxis()->SetLabelOffset(0.015);
    // h_data->GetXaxis()->SetTitleOffset(1.25);
    // h_data->GetYaxis()->SetTitleSize(0.04);
    // h_data->GetYaxis()->SetLabelSize(0.04);
    // h_data->GetYaxis()->SetTitleOffset(1.8);
    auto hs = new THStack("hs", " ");
    h_MuMuG->SetFillColor(TColor::GetColor(248, 206, 104));
    h_MuMuG->SetLineColor(TColor::GetColor(248, 206, 104));
    hs->Add(h_MuMuG.GetPtr());
    h_MuMuJet->SetFillColor(TColor::GetColor("#00A88F"));
    h_MuMuJet->SetLineColor(TColor::GetColor("#00A88F"));
    hs->Add(h_MuMuJet.GetPtr());

    h_data->GetXaxis()->SetLabelOffset(0.05);
    h_data->GetXaxis()->SetTickSize(0.03);
    h_data->GetYaxis()->SetTitle(Form("Events / %.2f", (bins[2]-bins[1])/bins[0]));
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetRangeUser(1, h_data->GetBinContent(h_data->GetMaximumBin()) * 20);
    h_data->GetYaxis()->SetTickSize(0.03);
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetLabelSize(0.055);
    h_data->GetYaxis()->SetTitleOffset(0.8);

    h_data->SetMarkerColor(TColor::GetColor("#202020"));
    h_data->SetMarkerSize(1.1);
    h_data->SetMarkerStyle(20);
    h_data->SetLineColor(TColor::GetColor("#202020"));
    h_data->SetLineWidth(2);
    h_data->Draw("EP");
    hs->Draw("hist same");
    h_data->Draw("EP same");

    TString lumi_text("59.82 fb^{-1}");
    int year = 2018;
    TString extra_text("Work-in-progress");
    TString proc_text("H #rightarrow #gamma* #gamma #rightarrow #mu#mu#gamma");
    CMS_lumi(pad1 , 5, 10, lumi_text, year, true, extra_text, proc_text, "");
    pad1->RedrawAxis();
    pad1->Update();

    auto leg = new TLegend(0.69, 0.7, 0.9, 0.88);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->AddEntry(h_data.GetPtr(),  "data", "LE1P");
    leg->AddEntry(h_MuMuG.GetPtr(), "#mu#mu + #gamma", "f");
    leg->AddEntry(h_MuMuJet.GetPtr(), "#mu#mu + Jet", "f");
    leg->Draw("same");
    canv->cd();

    //---------- 2nd Pad's Setting(Draw the Scale Factors) ----------//
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad2->SetGridy();
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();
    
    auto h_ratio = (TH1F*) h_data.GetPtr()->Clone();
    auto hs_hist = (TH1F*) h_MuMuG.GetPtr()->Clone();
    hs_hist->Add(h_MuMuJet.GetPtr(), 1);
    h_ratio->Divide(hs_hist);

    h_ratio->GetXaxis()->SetTitle(legName.c_str());
    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetYaxis()->SetRangeUser(0, 4);

    h_ratio->SetMarkerColor(TColor::GetColor("#202020"));
    h_ratio->SetMarkerSize(1.1);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetLineColor(TColor::GetColor("#202020"));
    h_ratio->SetLineWidth(2);
    h_ratio->GetXaxis()->SetTitleSize(0.16);
    h_ratio->GetXaxis()->SetTitleOffset(1);
    h_ratio->GetXaxis()->SetLabelSize(0.12);
    h_ratio->GetXaxis()->SetLabelOffset(0.03);
    h_ratio->GetYaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetTitleOffset(0.4);
    h_ratio->GetYaxis()->SetLabelSize(0.12);
    h_ratio->GetYaxis()->SetNdivisions(502);
    h_ratio->Draw("EP");

    // TString lumi_text("86.9 fb^{-1}");
    // int year = 2017;
    // TString extra_text("Work-in-progress");
    // TString proc_text("H #rightarrow #gamma* #gamma #rightarrow #mu#mu#gamma");
    // CMS_lumi((TCanvas*) canv, 5, 10, lumi_text, year, true, extra_text, proc_text, "");

    // string plot_path = "./plots/UL2018_mupho_bkg";
    // if (!bfs::exists(plot_path)){
    //     system(Form("mkdir -p %s", plot_path.c_str()));
    // }
    // canv->Print(Form("%s/%s.pdf", plot_path.c_str(), varName.c_str()));
    canv->Print("./test2.pdf");
    canv->Print("./test2.png");
    // delete canv;
    // delete leg;
}

void drawBkgComp(){
     vector<string> data_files = {
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018A_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018B_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018C_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018D_UL2018.root",

        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017D_UL2017.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017E_UL2017.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017F_UL2017.root"
    };

    vector<string> mc_files = {
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ggF_mmg_125_UL2018_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_VBF_mmg_125_UL2018_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_WH_mmg_125_UL2018_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ZH_mmg_125_UL2018_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ttH_mmg_125_UL2018_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_bbH_mmg_125_UL2018_MuPhoTrig.root",

        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ggF_mmg_125_UL2017_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_VBF_mmg_125_UL2017_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_WH_mmg_125_UL2017_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ZH_mmg_125_UL2017_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_ttH_mmg_125_UL2017_MuPhoTrig.root",
        // "/data4/chenghan/muon/miniTree/UL2017/miniTree_HDalitz_bbH_mmg_125_UL2017_MuPhoTrig.root"
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuG_mmg_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuJet_mmg_UL2018_MuPhoTrig.root"
    };

    ROOT::RDF::RNode df_data    = ROOT::RDataFrame("miniTree", data_files);
    ROOT::RDF::RNode df_MuMuG   = ROOT::RDataFrame("miniTree", "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuG_mmg_UL2018_MuPhoTrig.root");
    ROOT::RDF::RNode df_MuMuJet = ROOT::RDataFrame("miniTree", "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuJet_mmg_UL2018_MuPhoTrig.root");
    
    ROOT::RDF::RNode df_mc = ROOT::RDataFrame("miniTree", mc_files);

    map<string, string> var_sel = {
        // {"mu1Pt",               "mu1.Pt()"},
        // {"mu2Pt",               "mu2.Pt()"},
        // {"dimuPt",              "(mu1+mu2).Pt()"},
        {"dimuM",               "(mu1+mu2).M()"},
        // {"higgsMass",               "(mu1+mu2+pho).M()"},
        // {"mu1Eta",              "mu1.Eta()"},
        // {"mu2Eta",              "mu2.Eta()"},
        // {"phoEt",               "pho.Pt()"},
        // {"dRmumu",              "dRmumu"},
        // {"dRmax",               "dRmax"},
        // {"dRmin",               "dRmin"},
        // {"muPtRatio",           "muPtRatio"},
        // {"phoEtRatio",          "phoEtRatio"},
        // {"phoCorrR9Full5x5",    "phoCorrR9Full5x5"},
    };

    map<string, vector<double>> var_bins = {
        // {"mu1Pt",               {35, 0, 140}},
        // {"mu2Pt",               {35, 0, 140}},
        // {"dimuPt",              {35, 0, 140}},
        // {"dimuM",               {50, 0, 5}},
        {"dimuM",               {50, 0, 50}},
        // {"higgsMass",               {60, 110, 170}},
        // {"mu1Eta",              {60, -3, 3}},
        // {"mu2Eta",              {60, -3, 3}},
        // {"phoEt",               {35, 0, 140}},
        // {"dRmumu",              {80, 0, 4}},
        // {"dRmax",               {80, 0, 4}},
        // {"dRmin",               {80, 0, 4}},
        // {"muPtRatio",           {40, 0, 20}},
        // {"phoEtRatio",          {50, 0, 5}},
        // {"phoCorrR9Full5x5",    {50, 0, 1}},
    };

    map<string, string> var_legend = {
        // {"mu1Pt",               "p^{#mu^{1}}_{T} (GeV)"},
        // {"mu2Pt",               "p^{#mu^{2}}_{T} (GeV)"},
        // {"dimuPt",              "p^{#mu#mu}_{T} (GeV)"},
        {"dimuM",              "m_{#mu#mu} (GeV)"},
        // {"higgsMass",              "m_{#mu#mu#gamma} (GeV)"},
        // {"mu1Eta",              "#eta^{#mu^{1}}"},
        // {"mu2Eta",              "#eta^{#mu^{2}}"},
        // {"phoEt",               "E^{#gamma}_{T} (GeV)"},
        // {"dRmumu",              "#DeltaR(#mu^{1},#mu^{2})"},
        // {"dRmax",               "#DeltaR_{max}(#mu,#gamma)"},
        // {"dRmin",               "#DeltaR_{min}(#mu,#gamma)"},
        // {"muPtRatio",           "p^{#mu^{1}}_{T} / p^{#mu^{2}}_{T}"},
        // {"phoEtRatio",          "p^{#mu#mu}_{T} / p^{#gamma}_{T}"},
        // {"phoCorrR9Full5x5",    "Full5x5 R_{9}"},
    };

    for (auto it = var_sel.begin(); it != var_sel.end(); it++){
        const string varName        = it->first;
        const string sel            = it->second;
        const vector<double> bins   = var_bins[varName];
        const string legName        = var_legend[varName];
        draw_var(df_data, df_MuMuG, df_MuMuJet, varName, sel, bins, legName);
    }
}