#include <iostream>
#include <map>
#include <vector>
#include "TString.h"
#include "boost/algorithm/cxx11/any_of.hpp"

class CutCalculater{
public:
    CutCalculater(Long64_t totalev_){
        totalev = totalev_;
    }

    void AddCut(TString cutName_){
        if (!boost::algorithm::any_of_equal(name, cutName_)){
            name.push_back(cutName_);
            cutflow[cutName_] = 0;
        }
        cutflow[cutName_] += 1;
    }

    void Print(){
        vector<Long64_t> all = {totalev};
        vector<Long64_t> pass;

        for (size_t i = 0; i < name.size(); i++){
            if (i < name.size()-1)
                all.push_back(cutflow[name[i]]);
            pass.push_back(cutflow[name[i]]);
        }

        TString cutReport("Cut flow report: \n");
        for (size_t i = 0; i < all.size(); i++){
            cutReport += TString::Format("  - %-21s: pass = %-10lldall = %-10lldeff = %3.2f%%\n", name[i].Data(), pass[i], all[i], (float)pass[i]*100./all[i]);
        }
        printf("%s\n", cutReport.Data());
    }

private:
    Long64_t totalev = 0;
    std::map<TString, Long64_t> cutflow;
    std::vector<TString> name;
};