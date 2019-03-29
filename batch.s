#version 09-2016 fontana
#!/bin/csh 

### exec_path : folder of the executable file 
### logdir : folder where all the infos about the job will be saved : stopped job, problems, etc.
### macro_path : folder where is set the simulation macro
### exec_name : name of the simulation executable file
### macro_name : name of the macro
### outfile_prefix : name of the output root file 
### job_prefix : prefix for the job name (not so important)
### jobid : user defined ID for the job - here it's the date of the job execution in the format (yy-mm-dd_HH.MM)	

set exec_path=/afs/in2p3.fr/home/m/mfontana/Simulations/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/reconstruction_code
set logdir=/sps/hep/hadronth/mfontana/SimulationOutput/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/logdir
set job_path=/afs/in2p3.fr/home/m/mfontana/Simulations/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/job
set exec_name="reconstruction"

set argfile="NM_ReducedAbs_EnSel.m"
set job_prefix="MLEM_rec"
set today=`date +%y-%m-%d_%H.%M`

### Printf to check the information set are ok before the job starts

printf "Executable path is:\n"
printf "$exec_path\n"
printf "Executable name is $exec_name\n"
printf "The execution command will be $exec_path"/"$exec_name $exec_path"/"$argfile\n"	

set jobid="$job_prefix"_"$today"   
printf "Job ID: $jobid\n"
set jobname="$job_path"/"$jobid"
     
### Store the command line used to send a job in the variable "jobname". This variable is used as a script to send the reconstruction on the CC.
printf "$exec_path"/"$exec_name $exec_path"/"$argfile\n\n" > $jobname
	
### Change the job permissions to make it executable	
chmod u+x $jobname	
		
# Command line who sends the job. ct is the time estimate for the job and s_rss the job size
qsub -V -l sps=1 -l ct=10:00:00 -l s_rss=2G -e $logdir -o $logdir -j y -N $jobid $jobname

echo "\nReconstruction starts!\n"

