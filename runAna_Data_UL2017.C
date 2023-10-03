{
    TStopwatch time_all;
	time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C+O");

    vector<TString> inpath = {
        "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017B_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017C_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2017D_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2017E_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2017F_UL/*.root"
    };

    vector<TString> outpath = {
        "/data4/chenghan/muon/miniTree/UL2017/miniTree_SingleMu_Run2017B_UL2017.root",
        "/data4/chenghan/muon/miniTree/UL2017/miniTree_SingleMu_Run2017C_UL2017.root",
        "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017D_UL2017.root",
        "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017E_UL2017.root",
        "/data4/chenghan/muon/miniTree/UL2017/miniTree_MuEG_Run2017F_UL2017.root"
    };

    for (size_t i = 0; i < inpath.size(); i++){
        if (i < 2)
            xAna(inpath[i], outpath[i], 1, 14.42, "UL2017", "SingleMuTrig", "Nominal", false);
        else
            xAna(inpath[i], outpath[i], 1, 27.07, "UL2017", "MuPhoTrig", "Nominal", false);
    }

    time_all.Stop();
    time_all.Print();
    cout << "all done!" << endl;
}