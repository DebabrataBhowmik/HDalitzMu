#include <TChain.h>
#include <TTree.h>
#include <vector>
#include "TSystem.h"
#include "TString.h"

class TreeReader{
public:
    TreeReader(TString _fName, TString _tName){
        ch = new TChain(_tName);
        int nfiles = ch->Add(_fName);
        if (nfiles < 1){
            printf("TChain::Add return 0 file were found: %s\n", _fName.Data());
            gSystem->Exit(1);
        }
        printf("Find_files(): %d files are found in %s\n", nfiles, _fName.Data());
        printf("Get_TTree(): %s (Total entries = %lli)\n", _tName.Data(), ch->GetEntries());

        isMC = ch->GetBranch("nMC") ? true : false;
        ch->SetBranchStatus("*", 0);
    }

    ~TreeReader(){delete ch;}

    template <typename T>
    void SetBranchAddress(TString _bName, T &obj){
        ch->SetBranchStatus(_bName, 1);
        int status = ch->SetBranchAddress(_bName, &obj);
        if (status == -1 || status == -2 || status == -3 || status == -4 || status == -5){
            printf("TTree::CheckBranchAddressType %s return error type %d\n", _bName.Data(), status);
            gSystem->Exit(1);
        }
    }

    Long64_t GetEntries(){
        return ch->GetEntries();
    }

    void GetEntry(Long64_t ev){
        ch->GetEntry(ev);
    }

    // specific function for ggNtuple (to see if there is "nMC" branch)
    bool HasMC(){
        return isMC;
    }

private:
    TChain *ch = nullptr;
    bool isMC = false;
};