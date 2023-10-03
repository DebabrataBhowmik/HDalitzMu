
echo "Submit runAna_MC.C to process UL2018 signal MC" 
nohup root -l -b -q runAna_MC.C &> ./logger/UL2018_SignalMC.txt &   
echo "---> Logger file: ./logger/UL2018_SignalMC.txt"
echo "---> Wait for the compilation (60s)..." 
sleep 60
echo "" 

echo "Submit runAna_MC_UL2017.C to process UL2017 signal MC" 
nohup root -l -b -q runAna_MC_UL2017.C &> ./logger/UL2017_SignalMC.txt &   
echo "---> Logger file: ./logger/UL2017_SignalMC.txt"
sleep 20

echo "Submit runAna_MC_UL2017_IsoMu.C to process UL2017 signal MC" 
nohup root -l -b -q runAna_MC_UL2017_IsoMu.C &> ./logger/UL2017_SignalMC.txt &   
echo "---> Logger file: ./logger/UL2017_SignalMC_IsoMu.txt"
sleep 20

echo "Submit runAna_MC_UL2016preVFP.C to process UL2016preVFP signal MC" 
nohup root -l -b -q runAna_MC_UL2016preVFP.C &> ./logger/UL2016preVFP_SignalMC.txt &   
echo "---> Logger file: ./logger/UL2016preVFP_SignalMC.txt"
sleep 20

echo "Submit runAna_MC_UL2016postVFP.C to process UL2016postVFP signal MC" 
nohup root -l -b -q runAna_MC_UL2016postVFP.C &> ./logger/UL2016postVFP_SignalMC.txt &   
echo "---> Logger file: ./logger/UL2016postVFP_SignalMC.txt"
sleep 20

echo "Submit runAna_Data.C to process UL2018 data" 
nohup root -l -b -q runAna_Data.C &> ./logger/UL2018_Data.txt &   
echo "---> Logger file: ./logger/UL2018_Data.txt"
sleep 20

echo "Submit runAna_Data_UL2017.C to process UL2017 data" 
nohup root -l -b -q runAna_Data_UL2017.C &> ./logger/UL2017_Data.txt &   
echo "---> Logger file: ./logger/UL2017_Data.txt"
sleep 20

echo "Submit runAna_Data_UL2016.C to process UL2016 data" 
nohup root -l -b -q runAna_Data_UL2016.C &> ./logger/UL2016_Data.txt &   
echo "---> Logger file: ./logger/UL2016_Data.txt"
echo "---> All jobs sumitted" 