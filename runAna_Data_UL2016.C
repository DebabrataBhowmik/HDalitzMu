{
    TStopwatch time_all;
	time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C+O");

    vector<TString> inpath = {
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016B_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016C_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016D_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016E_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016F1_UL/*.root",

        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016F2_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016G_UL/*.root",
        "/data3/ggNtuples/V10_06_30_00/job_MuEG_Run2016H_UL/*.root"
    };

    vector<TString> outpath = {
        "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016B_UL2016preVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016C_UL2016preVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016D_UL2016preVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016E_UL2016preVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016preVFP/miniTree_MuEG_Run2016F1_UL2016preVFP.root",

        "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016F2_UL2016postVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016G_UL2016postVFP.root",
        "/data4/chenghan/muon/miniTree/UL2016postVFP/miniTree_MuEG_Run2016H_UL2016postVFP.root"
    };

    for (size_t i = 0; i < inpath.size(); i++){
        if (i < 5)
            xAna(inpath[i], outpath[i], 1, 19.52, "UL2016preVFP", "MuPhoTrig", "Nominal", false);
        else
            xAna(inpath[i], outpath[i], 1, 16.81, "UL2016postVFP", "MuPhoTrig", "Nominal", false);
    }

    time_all.Stop();
    time_all.Print();
    cout << "all done!" << endl;
}