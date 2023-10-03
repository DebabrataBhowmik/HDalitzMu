{
    TStopwatch time_all;
    time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C+O");

    vector<TString> inpath = {
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2018A_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2018B_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2018C_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2018D_UL/*.root"
    };

    vector<TString> outpath = {
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018A_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018B_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018C_UL2018.root",
        "/data4/chenghan/muon/miniTree/UL2018/miniTree_MuEG_Run2018D_UL2018.root"
    };

    for (size_t i = 0; i < inpath.size(); i++){
        xAna(inpath[i], outpath[i], 1, 59.82, "UL2018", "MuPhoTrig", "Nominal", false);
    }

    time_all.Stop();
    time_all.Print();
    cout << "all done!" << endl;
}