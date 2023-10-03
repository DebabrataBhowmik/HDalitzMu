#include <iostream>
#include <fstream>
#include <vector>
#include <string>
// #include <algorithm>
#include <cmath>
#include "TString.h"
#include "nlohmann/json.hpp"

class JsonSFsReader{
public:
    JsonSFsReader(std::string fileName, std::string keyName){
        std::ifstream f(fileName);
        data = nlohmann::json::parse(f)[keyName];
        f.close();
        abseta_bins = data["abseta_pt"]["binning"][0]["binning"].get<std::vector<float>>();
        pt_bins     = data["abseta_pt"]["binning"][1]["binning"].get<std::vector<float>>();
    }

    float GetScaleFactor(const float abseta, const float pt){
        int abseta_bin = FindBins(abseta_bins, abseta);
        int pt_bin     = FindBins(pt_bins, pt);
        if (abseta_bin!= -1 && pt_bin != -1){
            const char* etabinName = Form("abseta:[%.1f,%.1f]", abseta_bins[abseta_bin], abseta_bins[abseta_bin+1]);
            const char* ptbinName = Form("pt:[%d,%d]", (int)pt_bins[pt_bin], (int)pt_bins[pt_bin+1]);
            return data["abseta_pt"][etabinName][ptbinName]["value"];
        }
        return 1.;
    }

    float GetScaleFactorErr(const float abseta, const float pt){
        int abseta_bin = FindBins(abseta_bins, abseta);
        int pt_bin     = FindBins(pt_bins, pt);
        if (abseta_bin!= -1 && pt_bin != -1){
            const char* etabinName = Form("abseta:[%.1f,%.1f]", abseta_bins[abseta_bin], abseta_bins[abseta_bin+1]);
            const char* ptbinName = Form("pt:[%d,%d]", (int)pt_bins[pt_bin], (int)pt_bins[pt_bin+1]);
            float syst = data["abseta_pt"][etabinName][ptbinName]["syst"];
            float stat = data["abseta_pt"][etabinName][ptbinName]["stat"];
            return hypot(syst, stat);
        }
        return 0.;
    }

private:
    int FindBins(const std::vector<float> bin, const float var){
        // std::sort(bin.begin(), bin.end());
        auto lower = std::lower_bound(bin.begin(), bin.end(), var);
        int pos = std::distance(bin.begin(), lower) - 1;
        if (pos < 0 || pos > (int) bin.size()-2) // underflow or overflow
            return (int) -1;
        return pos;
    }
    nlohmann::json data;
    std::vector<float> abseta_bins;
    std::vector<float> pt_bins;
};
