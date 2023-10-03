{
    TStopwatch time_all;
	time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C+O");

    vector<TString> mode = {
        "ggF",
        "VBF",
        "WH",
        "ZH",
        "ttH",
        "bbH",

        "ggF",
        "VBF",
        "WH",
        "ZH",
        "ttH",
        "bbH",

        "ggF",
        "VBF",
        "WH",
        "ZH",
        "ttH",
        "bbH"
    };

    vector<double> xs = {
        (48.58  * 1000) * 3.92E-5,
        (3.782  * 1000) * 3.92E-5,
        (1.373  * 1000) * 3.92E-5,
        (0.8839 * 1000) * 3.92E-5,
        (0.5071 * 1000) * 3.92E-5,
        (0.4880 * 1000) * 3.92E-5,

        (52.22  * 1000) * 3.80E-5,
        (3.935  * 1000) * 3.80E-5,
        (1.565  * 1000) * 3.80E-5,
        (0.9939 * 1000) * 3.80E-5,
        (0.5697 * 1000) * 3.80E-5,
        (0.5534 * 1000) * 3.80E-5,

        (45.31  * 1000) * 3.90E-5,
        (3.637  * 1000) * 3.90E-5,
        (1.209  * 1000) * 3.90E-5,
        (0.7899 * 1000) * 3.90E-5,
        (0.4539 * 1000) * 3.90E-5,
        (0.4304 * 1000) * 3.90E-5
    };

    vector<TString> inpath = {
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

    vector<TString> outpath = {
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ggF_mmg_125_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_VBF_mmg_125_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_WH_mmg_125_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ZH_mmg_125_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ttH_mmg_125_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_bbH_mmg_125_UL2018_MuPhoTrig.root",

        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ggF_mmg_120_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_VBF_mmg_120_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_WH_mmg_120_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ZH_mmg_120_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ttH_mmg_120_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_bbH_mmg_120_UL2018_MuPhoTrig.root",

        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ggF_mmg_130_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_VBF_mmg_130_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_WH_mmg_130_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ZH_mmg_130_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_ttH_mmg_130_UL2018_MuPhoTrig.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_HDalitz_bbH_mmg_130_UL2018_MuPhoTrig.root"
    };

    vector<TString> variation = {
        "Nominal",
        "JERUp",
        "JERDo",
        "JECUp",
        "JECDo",
        "MuCalibStat",
        "MuCalibZpt",
        "MuCalibEwk",
        "MuCalibdeltaM",
        "MuCalibEwk2",
        "PhoNoR9Corr",
        "PhoScaleStatUp",
        "PhoScaleSystUp",
        "PhoScaleGainUp",
        "PhoScaleStatDo",
        "PhoScaleSystDo",
        "PhoScaleGainDo",
        "PhoSigmaPhiUp",
        "PhoSigmaRhoUp",
        "PhoSigmaRhoDo",
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