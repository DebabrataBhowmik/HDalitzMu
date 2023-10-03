{
    TStopwatch time_all;
	time_all.Start();

    gSystem->SetBuildDir("build", kTRUE);
    gROOT->ProcessLine(".L xAna.C++");

    TString inpath = "/data4/chenghan/muon/skimTree/job_UL18_Dalitz_mmg_m125.root";
    TString outpath = "/data4/chenghan/muon/miniTree/test_mmg_ch.root";
    xAna(inpath, outpath, (48.58 * 1000) * 3.92E-5, 59.82, "UL2018", "MuPhoTrig", "Nominal", true);

    time_all.Stop();
    time_all.Print();
}