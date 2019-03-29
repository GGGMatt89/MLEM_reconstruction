#$ -j y
#$ -l sps=1
#$ -l ct=10:00:00
#$ -l s_rss=2G
#$ -e /sps/hep/hadronth/mfontana/SimulationOutput/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/logdir
#$ -o /sps/hep/hadronth/mfontana/SimulationOutput/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/logdir

echo "\nMLEM reconstruction started!\n"
echo "\nI want an output\n"
echo "\nPlease give me an output\n"

/afs/in2p3.fr/home/m/mfontana/Simulations/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/reconstruction_code/reconstruction /afs/in2p3.fr/home/m/mfontana/Simulations/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/reconstruction_code/NM_ReducedAbs_EnSel.m

echo "\nReconstruction over\n"