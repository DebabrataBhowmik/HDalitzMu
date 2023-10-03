#include <vector>
#include <limits>
// #include "Math/Vector4D.h"
#include "TLorentzVector.h"

using namespace std;

bool IsPhoIntConv(
    const int nMC,
    const vector<int> mcPID,
    const vector<int> mcMomPID,
    const vector<unsigned short> mcStatusFlag
){
    int realpho = 0;
    for (int i = 0; i < nMC; i++){
        bool hardProc = (mcStatusFlag[i] >> 0 & 1) == 1;
        bool isPrompt = (mcStatusFlag[i] >> 1 & 1) == 1;
        bool isRealPho = mcPID[i] == 22 && mcMomPID[i] == 25 && hardProc && isPrompt;

        if (isRealPho)
            realpho += 1;
    }
    bool isPhoIntConv = (realpho == 0) ? true : false; // no real photon

    return isPhoIntConv;
}


float CalcLHEMmumu(
    const int nLHE,
    const vector<int> lhePID,
    const vector<float> lhePx,
    const vector<float> lhePy,
    const vector<float> lhePz,
    const vector<float> lheE
){
    std::vector<int> muIdx; // find lhe level muons
    for (int i = 0; i < nLHE; i++){
        if (abs(lhePID[i]) != 13)
            continue;
        muIdx.push_back(i);
    }

    float mass = std::numeric_limits<float>::max(); // assign a large value that will be filtered out
    if (muIdx.size() > 1){
        
        TLorentzVector mu1, mu2;
        mu1.SetPxPyPzE(lhePx[muIdx[0]], lhePy[muIdx[0]], lhePz[muIdx[0]], lheE[muIdx[0]]);
        mu2.SetPxPyPzE(lhePx[muIdx[1]], lhePy[muIdx[1]], lhePz[muIdx[1]], lheE[muIdx[1]]);
        // ROOT::Math::PxPyPzMVector mu1(lhePx[muIdx[0]], lhePy[muIdx[0]], lhePz[muIdx[0]], 0.1057);
        // ROOT::Math::PxPyPzMVector mu2(lhePx[muIdx[1]], lhePy[muIdx[1]], lhePz[muIdx[1]], 0.1057);
        mass = (mu1 + mu2).M();
    }
    return mass;
}