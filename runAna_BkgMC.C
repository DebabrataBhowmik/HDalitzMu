{
    TStopwatch time_all;
	time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C+O");

    vector<TString> mode = {
        "MuMuG",
        "MuMuJet"
    };

    vector<double> xs = {
        (7.425  * 1000),
        (25.78  * 1000)
    };

    vector<TString> inpath = {
        "/data5/ggNtuples/V10_06_30_00/job_UL18_MuMuG_mll_0to60/*.root",
        "/data5/ggNtuples/V10_06_30_00/job_UL18_MuMuJet_mll_0to60/*.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m125.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m125.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m125.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m125.root",

        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m120.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_VBF_m120.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m120.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m120.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m120.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m120.root",

        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m130.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_VBF_m130.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_WH_m130.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ZH_m130.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_ttH_m130.root",
        // "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_bbH_m130.root"
    };

    vector<TString> outpath = {
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuG_mmg_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuMuJet_mmg_UL2018_MuPhoTrig.root"
    };

    vector<TString> variation = {
        "Nominal",
        // "JERUp",
        // "JERDo",
        // "JECUp",
        // "JECDo",
        // "MuCalibStat",
        // "MuCalibZpt",
        // "MuCalibEwk",
        // "MuCalibdeltaM",
        // "MuCalibEwk2",
        // "PhoNoR9Corr",
        // "PhoScaleStatUp",
        // "PhoScaleSystUp",
        // "PhoScaleGainUp",
        // "PhoScaleStatDo",
        // "PhoScaleSystDo",
        // "PhoScaleGainDo",
        // "PhoSigmaPhiUp",
        // "PhoSigmaRhoUp",
        // "PhoSigmaRhoDo",
    };

    for (size_t i = 0; i < inpath.size(); i++){
        printf("\033[0;32m Higgs production mode: %s \033[0m\n", mode[i].Data());
        for (size_t j = 0; j < variation.size(); j++){
            printf ("=========================== variation: %s ==============================\n", variation[j].Data());
            xAna(inpath[i], outpath[i], xs[i], 59.82, "UL2018", "MuPhoTrig", variation[j], true);
        }
    }

    time_all.Stop();
    time_all.Print();
    
    cout << "all done!" << endl;
}