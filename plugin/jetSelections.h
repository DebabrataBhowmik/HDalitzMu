#include <vector>
#include "TString.h"

bool JetIDUL(
    float jetEta,
    float jetNHF,
    float jetNEF,
    int jetNNP,
    int jetNCH,
    float jetCHF,
    float jetCEF,
    float jetMUF,
    TString era
){
    bool cutID = false;
    if (era.Contains("2016")){
        const bool eta1 = fabs(jetEta) <= 2.4;
        const bool eta2 = fabs(jetEta) > 2.4 && fabs(jetEta) <= 2.7;
        const bool eta3 = fabs(jetEta) > 2.7 && fabs(jetEta) <= 3.;
        const bool eta4 = fabs(jetEta) > 3. && fabs(jetEta) <= 5.;

        int numConst = jetNCH + jetNNP; // Number of Constituents
        if (eta1)
            cutID = jetNHF < 0.9 && jetNEF < 0.9 && numConst > 1 && /*jetMUF < 0.8 &&*/ jetCHF > 0 && jetNCH > 0 /*&& jetCEF < 0.8*/;
        else if (eta2)
            cutID = jetNHF < 0.9 && jetNEF < 0.9;
        else if (eta3)
            cutID = jetNHF < 0.9 && jetNEF < 0.99 && jetNEF > 0. && jetNNP > 1;
        else if (eta4)
            cutID = jetNHF > 0.2 && jetNEF < 0.9 && jetNNP > 10;
    }
    else{
        const bool eta1 = fabs(jetEta) <= 2.6;
        const bool eta2 = fabs(jetEta) > 2.6 && fabs(jetEta) <= 2.7;
        const bool eta3 = fabs(jetEta) > 2.7 && fabs(jetEta) <= 3.;
        const bool eta4 = fabs(jetEta) > 3. && fabs(jetEta) <= 5.;

        int numConst = jetNCH + jetNNP; // Number of Constituents
        if (eta1)
            cutID = jetNHF < 0.9 && jetNEF < 0.9 && numConst > 1 && jetCHF > 0. && jetNCH > 0 /*&& jetMUF < 0.8 && jetCEF < 0.8*/;
        else if (eta2)
            cutID = jetNHF < 0.9 && jetNEF < 0.99 && jetNCH > 0 /*&& jetMUF < 0.8 && jetCEF < 0.8*/;
        else if (eta3)
            cutID = jetNEF <  0.99 && jetNEF > 0.01 && jetNNP > 1;
        else if (eta4)
            cutID = jetNHF > 0.2 && jetNEF < 0.9 && jetNCH > 10;
    }

    return cutID;
}