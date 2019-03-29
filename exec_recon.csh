#!/bin/tcsh

echo
echo "starting qsub..." 

set logDir="/afs/in2p3.fr/home/m/mfontana/reconstruction_output/logdir/"

echo
echo

echo "launching qsub..." 

qsub -P P_hadronth  -l fsize=1G,ct=10:00:00,vmem=2G,sps=1 -e $logDir -o $logDir -j y $PWD/recon_submit.csh
