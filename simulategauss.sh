#!/bin/bash
#This code is the automated pipeline to recover dm error for different numbers of frequency channels
#We first acquire the inputs
#Full pulsar Jname.
PSR=$1
#Number of frequency channels we want to use.
chns=$2

#Activate conda environment
#conda activate test_pss_jacob

#We cd into a working directory
#cd workingdir
#Create a file to recover the variables later in the python scripts
#touch tempvars.txt
#echo "$psr" > tempvars.txt
#echo "$chns" >> tempvars.txt

par=$(ls -t /home/jdc0059/SimulatorProj/templatefiles/${PSR}/*.par | head -1)

#subprocess.run('source activate environment-pint_psr && python /home/jdc0059/SimulatorProj/simcode.py && source deactivate', shell=True)
#Run the simulator (Remember to change file paths/names and CHNS numbers) 
python /home/jdc0059/SimulatorProj/simgauss.py $PSR $chns

#To fix any issued with the .tim file we want to do the following:
pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt  -C rcvr:name -C be:name -f tempo2 IPTA -s "/home/jdc0059/SimulatorProj/templatefiles/J0000+0000/fits_Lband_ASP0.fits.sm" "/home/jdc0059/SimulatorProj/sim_folder/fits_Lband_ASP0.fits" > "/home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.tim"

cd sim_folder

sed 's/ ao / @ /g' fits_Lband_0.tim > fits_Lband_0_2.tim

sed 's/-rcvr:name lbw/-fe L-wide/g' fits_Lband_0_2.tim > fits_Lband_0_3.tim

#Does Puppi or Guppi need to change?
sed 's/-be:name PUPPI/-f L-wide_PUPPI/g' fits_Lband_0_3.tim > fits_Lband_0_4.tim

cd ..

#Here we delete unecessairy files and move the others
rm /home/jdc0059/SimulatorProj/sim_folder/*0.fits.tim
rm /home/jdc0059/SimulatorProj/sim_folder/*.fits_ALL.tim
rm /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_1.tim
rm /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_2.tim
rm /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_3.tim
rm /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_ALL.tim
rm /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.tim

mkdir /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/
mkdir /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_4.tim  /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /home/jdc0059/SimulatorProj/sim_folder/${PSR}.tim /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /home/jdc0059/SimulatorProj/sim_folder/Combined_tim_file.tim /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /home/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.fits /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /home/jdc0059/SimulatorProj/sim_folder/*.fits /home/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/

#Create the output file for pint program
touch PINT_output.txt

#Activate second conda environment that works with PINT
source /home/jdc0059/anaconda3/etc/profile.d/conda.sh
conda activate pint_psr

#Run pint on our data (Remember to change file paths/names)
python /home/jdc0059/SimulatorProj/PINT_gauss.py $PSR $chns $par


#Deactive PINT conda environment
conda deactivate







#rm tempvars.txt
