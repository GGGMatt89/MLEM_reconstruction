#!/bin/bash


echo
echo "starting reconstruction..."

energy=(245 555 1099 1524)
Emeas="keV"
logDir="/Users/mattia.fontana/Desktop/reconstruction_MLEM/results/logdir/"
resSubDir="4Siplanes/"
resDir="/Users/mattia.fontana/Desktop/reconstruction_MLEM/results/geometryStudy/"$resSubDir
saveDir="/Users/mattia.fontana/Desktop/reconstruction_MLEM/results/geometryStudy/"

for i in 245 555 1099 1524
do
    mkdir $resDir$i$Emeas

    echo "launching reconstruction for $i..."

    $PWD/reconstruction $PWD/NM_4Siplanes_$i.m > $logDir/log_$i.txt

    echo "Reconstruction over"

done

echo " DONE - Bye! "
