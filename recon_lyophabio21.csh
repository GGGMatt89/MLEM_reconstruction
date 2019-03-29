#!/bin/bash


echo
echo "starting reconstruction..." 

energy=(245 555 1099 1524)
Emeas="keV"
logDir="/home/collat/fontana/reconstruction_MLEM/reconstruction_output/logdir/"
resSubDir="3Siplanes/"
resDir="/home/collat/fontana/reconstruction_MLEM/reconstruction_output/Results/geometryStudy/"$resSubDir
saveDir="/Users/mattia.fontana/PhD/Work/Simulations/Geant4/ComptonCamera_Simulation/CC_Nuclear_Medicine/reconstruction_MLEM/results/geometryStudy/"

for i in 245 555 1099 1524
do
    mkdir $resDir$i$Emeas
    
    echo "launching reconstruction for $i..." 

    $PWD/reconstruction $PWD/NM_3Siplanes_$i.m > $logDir/log_$i.txt

    echo "Reconstruction over"

    echo "Transferring files to lyophabio15"

    scp -r $resDir$i$Emeas/* mattia.fontana@lyophabio15:$saveDir$resSubDir$i$Emeas/

    scp -r $logDir/log_$i.txt mattia.fontana@lyophabio15:$saveDir$resSubDir$i$Emeas/ 
done

echo " DONE - Bye! "
