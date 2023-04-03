#!/bin/bash
#This code is the automated pipeline to recover dm error for different numbers of frequency channels
#We first acquire the inputs
#Full pulsar Jname.
PSR=$1
#Number of frequency channels we want to use.
chns=$2

template=$(ls -t /bowser/jdc0059/SimulatorProj/templatefiles/${PSR}/${PSR}.Rcvr1_2.*.sum.sm | head -1)
PSRnoj="${PSR:1}"
#puppior=$(ls -t /nanograv/timing/data/${PSRnoj}/puppi/2016/rawdata/puppi_*5581*_0001.fits | head -1)
#cp $puppior /bowser/jdc0059/SimulatorProj/templatefiles/$PSR/.
#puppi=$(ls -t /bowser/jdc0059/SimulatorProj/templatefiles/$PSR/puppi_* | head -1)
par=$(ls -t /bowser/jdc0059/SimulatorProj/templatefiles/${PSR}/*.par | head -1)

echo "template:"
echo $template

#cd templatefiles/${PSR}/.
#rm ${PSR}.par
#touch ${PSR}.par
#touch temp1.par
#sed '/^DMX/d' ${PSR}*.par >> temp1.par
#sed 's/RNAMP//' temp1.par >> temp1.par
#sed 's/RNIDX//' temp1.par >> temp1.par
#sed 's/T2EFAC//' temp1.par >> temp1.par
#sed 's/T2EQUAD//' temp1.par >> temp1.par 
#sed 's/ECORR//' temp1.par >> temp1.par
#sed 's/JUMP//' temp1.par >> temp1.par
#sed -i 's/ 1 / 0 /' temp1.par >> temp1.par
#sed '/DM/s/$/  1/' temp1.par >> ${PSR}.par
#rm temp1.par
#cd ..
#cd ..

#par=$(ls -t /bowser/jdc0059/SimulatorProj/templatefiles/${PSR}/${PSR}.par | head -1)
echo "parfile:"
echo $par
#Activate conda environment
#conda activate test_pss_jacob

#We cd into a working directory
#cd workingdir
#Create a file to recover the variables later in the python scripts
#touch tempvars.txt
#echo "$psr" > tempvars.txt
#echo "$chns" >> tempvars.txt

#subprocess.run('source activate environment-pint_psr && python /home/jdc0059/SimulatorProj/simcode.py && source deactivate', shell=True)
#Run the simulator (Remember to change file paths/names and CHNS numbers) 
python /bowser/jdc0059/SimulatorProj/simcode.py ${PSR} ${chns} ${template} ${par}

echo "Start of fix:"
#To fix any issued with the .tim file we want to do the following: (Check to make sure second entry name is right, this will change depending on the telescope either ASP0 (GBT) or 0 (AO) or something)
pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt  -C rcvr:name -C be:name -f tempo2 IPTA -s "/bowser/jdc0059/SimulatorProj/templatefiles/J0740+6620/J0740+6620.Rcvr1_2.GUPPI.12y.x.sum.sm" "/bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_ASP0.fits" > "/bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.tim"
echo "end of fix:"

echo "Going to modify the file"
cd sim_folder

sed 's/ ao / @ /g' fits_Lband_0.tim > fits_Lband_0_2.tim
#(the line below is used for AO)
#sed 's/-rcvr:name lbw/-fe L-wide/g' fits_Lband_0_2.tim > fits_Lband_0_3.tim # what it was originally
#The line below is used for GBT (Some moficiations might still be needed)
sed 's/-rcvr:name Rcvr1_2/-fe L-wide/g' fits_Lband_0_2.tim > fits_Lband_0_3.tim
#A third option is this:
#sed 's/-rcvr:name Rcvr_800/-fe L-wide/g' fits_Lband_0_2.tim > fits_Lband_0_3.tim

#Does Puppi or Guppi needs to change depending on if using AO or GBT.
sed 's/-be:name GUPPI/-f L-wide_GUPPI/g' fits_Lband_0_3.tim > fits_Lband_0_4.tim
echo "Done modifying the file"

cd ..

#Here we delete unecessairy files and move the oth
rm /bowser/jdc0059/SimulatorProj/sim_folder/*0.fits.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/*.fits_ALL.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_1.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_2.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_3.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_ALL.tim
rm /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.tim

mkdir /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/
mkdir /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0_4.tim  /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /bowser/jdc0059/SimulatorProj/sim_folder/${PSR}.tim /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /bowser/jdc0059/SimulatorProj/sim_folder/Combined_tim_file.tim /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /bowser/jdc0059/SimulatorProj/sim_folder/fits_Lband_0.fits /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/
mv /bowser/jdc0059/SimulatorProj/sim_folder/*.fits /bowser/jdc0059/SimulatorProj/sim_folder/results/${PSR}/${chns}_chns/

#Create the output file for pint program
touch PINT_output.txt

#Activate second conda environment that works with PINT
#source /bowser/jdc0059/anaconda3/etc/profile.d/conda.sh
#conda activate pint_psr

#Run pint on our data (Remember to change file paths/names)
/minish/jdc0059/.conda/envs/pint_psr/bin/python /bowser/jdc0059/SimulatorProj/PINT_processing.py ${PSR} ${chns} ${par}


#Deactive PINT conda environment
#conda deactivate







#rm tempvars.txt 
