#!/usr/bin/env python
# coding: utf-8

# # Simulating Data Notebook on bowser
# 
# This notebook will be used to generate simulated data using the pulsar signal simulator. This notebook will use the new api to generate the data with various effects, such as DM variations, FD parameters, pulse profile variation, and scattering tails. This notebook will be used only for generating the data, and is designed only to work with the structure set up on bowser. A separate notebook will be used to for analyzing the data.
# 
# We start by importing all the packages we need

# In[1]:


# These are all standard packages
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import ast
# These are more specialized pulsar packages, necessary for this notebook
import pypulse as pp # not required for psrsigsim
# and import the signal simulator
import psrsigsim as pss
import psrsigsim.telescope.receiver as rcvr
import psrsigsim.telescope.backend as bcknd
from psrsigsim.utils import make_quant




# Define a function for easier plotting later on/throughout the testing
def plotsignal(signals, nbins=2048):
    # signals can be a list of multiple signals to overplot
    for ii in range(len(signals)):
        # Define the x axis
        phases = np.linspace(0.0, len(signals[ii]), len(signals[ii]))/nbins
        # now plot it
        plt.plot(phases, signals[ii], label="signal %s" % (ii))
    plt.xlim([0.0, np.max(phases)])
    plt.xlabel("Pulse Phase")
    plt.ylabel("Arb. Flux")
    plt.show()
    plt.close()



# Functions that can be used in-notebook to run command-line commands, not strictly necessary for the PsrSigSim
def call(x):
    subprocess.call(x,shell=True)


def callvar(x):
    variable = subprocess.check_output(x, shell=True)
    return variable


# We define a function to create a variety of DM variations (at least simply)
# Same DM variation function as in the simulate data notebook. Used for comparison with recovered values.
def DMvars(dm, epochs, epoch_length = 30.0, variation = 'both', slope = 1.0, amplitude = 1.0, period = 365.0):
    """
    This version of the function will take in a slope (10^-3 pc cm^-3 yr^-1), 
    amplitude (10^-4 pc cm^-3), and period (days) as reported as in Jones et al. 2017
    to determine what the DM variations should look like.
    
    The inputs here are the DM of the pulsar and the number of epochs, or simulated files
    or values of DM we want to input. For now these are all done with linear spacings using
    numpy.linspace. Future versions could add different sampling. The other inputs are:
    epoch_length - time between epochs, e.g. in days, e.g. the default is one month,
                   meaning that there will be one month between epochs.
    variation - can be 'none', linear', 'periodic' (sinusoidal), or 'both'. Gives the trend 
                of DM variations to simulate. If none will return a list with epoch number
                of the input DM.
                NOTE: for sinusiodal variations, we will have one cycle for every 4 epochs
                      as a default for now.
    slope - slope of the linear trend of DM variations over one year (10^-3 pc cm^-3)
    amplitude - The amplitude of the sinusoidal trend (10^-4 pc cm^-3)
    period - period of the sinusoidal trend (days)
    """
    # Find the total time the DM is varying over in days:
    if isinstance(epochs, int):
        set_length = epochs*epoch_length # total days the data set spans
    # Otherwise we assume it's a list of MJDs to calculate the DM at
    else: 
        set_length = np.max(epochs)-np.min(epochs)
    
    if variation == 'linear':
        # find the overall scale of the variations depending on the slope
        set_slope = slope*(set_length/365.0) # overall linear change over the data set
        if isinstance(epochs, int):
            DM_values = np.linspace(dm-(set_slope*10**-3/2.0), dm+(set_slope*10**-3/2.0), epochs)
        else:
            print(epochs-np.min(epochs))
            DM_values = slope*((epochs-np.min(epochs))/set_length)*10**-3 + (dm-(set_slope*10**-3/2.0))
    elif variation == 'periodic':
        cycles = set_length/period # number of cycles
        if isinstance(epochs, int):
            DM_values = np.sin(np.linspace(0, cycles*2*np.pi, epochs))*amplitude*10**-4+dm
        else:
            DM_values = np.sin((cycles*2*np.pi)*((epochs-np.min(epochs))/set_length))*amplitude*10**-4+dm
    elif variation == 'both':
        # find the overall scale of the variations depending on the slope
        if isinstance(epochs, int):
            set_slope = slope*(set_length/365.0) # overall linear change over the data set
            DM_linear = np.linspace(dm-(set_slope*10**-3/2.0), dm+(set_slope*10**-3/2.0), epochs)
            # Now get the periodic component
            cycles = set_length/period # number of cycles
            DM_periodic = np.sin(np.linspace(0, cycles*2*np.pi, epochs))*amplitude*10**-4
        else:
            set_slope = slope*(set_length/365.0) # overall linear change over the data set
            DM_linear = slope*((epochs-np.min(epochs))/set_length)*10**-3 + (dm-(set_slope*10**-3/2.0))
            # Now get the periodic component
            cycles = set_length/period # number of cycles
            DM_periodic = np.sin((cycles*2*np.pi)*((epochs-np.min(epochs))/set_length))*amplitude*10**-4
        # Get the total DM variation
        DM_values = DM_linear+DM_periodic
    elif variation == 'none':
        if isinstance(epochs, int):
            DM_values = np.repeat(dm, epochs)
        else:
            DM_values = np.repeat(dm, len(epochs))
    else:
        print("ERROR: Incompatible DM variation specified. Returning input DM value.")
        DM_values = np.repeat(dm, epochs)
    return DM_values





"""
BJS: May need some updates here. Also recommmend using new branch with PSRFTIS saving issues fixed.
"""
"""
JCT: I removed all the "test" in fron of the file names here.
"""
# We define a flexible function to run the full simulations given some input parameters
def simulate_full(sim_folder, increments, f0, bw, Nf, f_samp, fold, subintlen, period, flux, psr_name,\
                  ObsTime, DMs, tempfits, parfile, prof_array, twoD = True, scattering = False, FD = False,
                  ref_MJD = 56000.0, MJD_start = 55999.9861, alpha = 0.0, ref_freq = 1400.0):
    # Make the folder to put the simulated data if it doesn't already exist
    if not os.path.exists(sim_folder):
        print("1")
        call("mkdir %s" % (sim_folder))
    # Make a list of the new saved file names
    fitsfiles = []
    # Now we loop through and make the data files; we start with the L-band data
    for ii in range(len(increments)):
        print("2")
        #start by defining the signal
        #sim_signal = pss.signal.FilterBankSignal(fcent = f0, bandwidth = bw, Nsubband=Nf,\
        #                                      sample_rate=f_samp, fold=fold, sublen=subintlen)
        if f0 > 500.0 and 'puppi' in tempfits:
            print("2.1")
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_Lband_'+str(ii)+'.fits'
        elif f0 < 500.0 and 'puppi' in tempfits:
            print("2.2")
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_430_'+str(ii)+'.fits'
        elif f0 > 500.0 and 'puppi' not in tempfits:
            print("2.3")
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_Lband_ASP'+str(ii)+'.fits'
        elif f0 < 500.0 and 'puppi' not in tempfits:
            print("2.4")
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_430_ASP'+str(ii)+'.fits'
           
        pfit = pss.io.PSRFITS(path=fitspath, template=tempfits, fits_mode='copy', \
                              obs_mode='PSR')
        print("3")
        #try:
        #    sim_signal = pfit.make_signal_from_psrfits()
        #except:
        sim_signal = pss.signal.FilterBankSignal(fcent = f0, bandwidth = bw, Nsubband=Nf,\
                                sample_rate=f_samp, fold=fold, sublen=subintlen)
        # Need to edit some parameters
        sim_signal._sublen = make_quant(subintlen, 's')
        
        if f0 > 500:
            print("4")
            sim_signal._dat_freq = np.flipud(sim_signal._dat_freq)
       
        if twoD:
            print("5.1")
            # Define the profile - this is a 2D wideband profile array
            prof = pss.pulsar.DataPortrait(prof_array, phases = None)
        else:
            print("5.2")
            # Define the profile - this is a single 1D profile
            prof = pss.pulsar.DataProfile(prof_array, phases = None, \
                                                Nchan=sim_signal.Nchan)
       
        # Define the pulsar
        pulsar = pss.pulsar.Pulsar(period=period, Smean=flux, profiles=prof, name=psr_name,\
                                   specidx = alpha, ref_freq = ref_freq)
        print("6")
    # Define the ISM object
        print("7")
        ism_ob = pss.ism.ISM()
        # Add scattering here if desired; input should be a list: [tau_d (s), ref_freq (MHz)], assumes convolve=True
        if scattering != False:
            ism_ob.scatter_broaden(sim_signal, scattering[0][ii], scattering[1], convolve=True, pulsar=pulsar)
        # make the pulses
        pulsar.make_pulses(sim_signal, tobs = ObsTime)
        # disperse the data, assumes DMs are in a list, even if constant
        ism_ob.disperse(sim_signal, DMs[ii])
        # Do FD shift directly if desired, assumes a list of parameters in units of seconds
        if FD != False:
            ism_ob.FD_shift(sim_signal, FD)
        # Define telescope and add radiometer noise
        tscope = pss.telescope.telescope.Arecibo()
        # Check frequency for noise and savefile name
        print("8")
        print(sim_signal.data)
        if f0 > 500.0 and 'puppi' in tempfits:
            tscope.observe(sim_signal, pulsar, system="Lband_PUPPI", noise=True)
            # Need to rescale the data due to bit overflow issues?
            sim_signal._data = sim_signal.data[:,:]/np.sqrt(2000.0)
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_Lband_'+str(ii)+'.fits'
            fitsfiles.append(fitspath)
        elif f0 < 500.0 and 'puppi' in tempfits:
            tscope.observe(sim_signal, pulsar, system="430_PUPPI", noise=True)
            sim_signal._data = sim_signal.data[:,:]/np.sqrt(2000.0)
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_430_'+str(ii)+'.fits'
            fitsfiles.append(fitspath)
        elif f0 > 500.0 and 'puppi' not in tempfits:
            tscope.observe(sim_signal, pulsar, system="Lband_ASP", noise=True)
            # Need to rescale the data due to bit overflow issues?
            sim_signal._data = sim_signal.data[:,:]/np.sqrt(2000.0)
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_Lband_ASP'+str(ii)+'.fits'
            fitsfiles.append(fitspath)
        elif f0 < 500.0 and 'puppi' not in tempfits:
            tscope.observe(sim_signal, pulsar, system="430_ASP", noise=True)
            sim_signal._data = sim_signal.data[:,:]/np.sqrt(2000.0)
            # Get our template PSRFITS file
            fitspath = sim_folder+'/'+'fits_430_ASP'+str(ii)+'.fits'
            fitsfiles.append(fitspath)
       
        pfit = pss.io.PSRFITS(path=fitspath, template=tempfits, fits_mode='copy', \
                              obs_mode='PSR')
        pfit.get_signal_params(signal = sim_signal)
        print("9")
        # Now save the data
        pfit.save(sim_signal, pulsar, parfile = parfile, \
                  MJD_start = MJD_start+increments[ii],segLength = 60.0,\
                  inc_len = increments[ii], ref_MJD = ref_MJD, usePint = True)
        print("10")
    # And return the list of saved fits files
    return fitsfiles




# Now we need to load in some functions from Tim Pennucci's pulse portrature python package so that we can read and process the pickled wideband pulse profiles.

# ## Setting up the Simulation
# 
# Now we will set up the different parameters we need in order to simulate our signals. These will include different input parameters (observing time, period, dm, telescope, etc.) as well as the pulse profiles, and any other convenience functions or parameters we need before running the simulation.

#Gaussian Test Pulsar
# We define all variables needed for the simulation at the begining; start with L-band values
dm = 35.6 # pc cm^-3
F0 = np.double(223.713647) # pulsar frequency, Hz
period = np.double(1.0)/F0 # pulsar period, seconds
bw_Lband = 800.0 # bandwidth, MHz
#WE WILL WANT TO EDIT THIS NUMBER
Nf_Lband = 1 # number of frequency channels
f0_Lband = 1380 # central observing frequency, MHz
telescope = "GBT" # Telescope name (for default telescopes)
psr_name = "J0000+0000" # pulsar name
flux_Lband = 0.00119 # mean flux (AO) for GB it is 0.00202 Jy (from Alam et al. 2020)
f_samp = F0*2048*10**-6 # sampling rate, chosen to be 2048 bins per profile
# This was for a single file, want to simulate full observation
ObsTime_Lband = 1200 #Before it was 1269.6189 do we have the right values? # total observing time, seconds
subintlen_Lband = ObsTime_Lband # length of subintegration, seconds




### Loading the Profiles
# 
# Now we load the NANOGrav 11-yr data profiles. These are 1-D referenced just to a single frequency

# start with 1400 MHz profile
#Might need these profiles?
#template_Lband = "templatefiles/J0340+4130/J0340+4130.Rcvr1_2.GUPPI.12y.x.sum.sm"
#temp_Lband = pp.archive.Archive(template_Lband, lowmem=True)
#temp_Lband_ar = temp_Lband.getData()

#Let's try creating the gaussian profile here:
gauss_prof = pss.pulsar.GaussProfile(peak = 0.25, width = 0.01, amp = 1.0)
gauss_prof.init_profiles(2048, Nchan = 1)
#gauss_prof.init_profiles(2048, Nchan = 4)
#Do I then need to ass the other lines which are added to real profile?



# Now we want to load the 2-D wideband NANOGrav profiles so that we can model the profile evolution.

"""We can only extrapolate within the bounds of the profile frequencies.
We will need to figure out how to do placeholder frequencies for other channels..."""
#real_data_highSN_profs_file = \
#"templatefiles/J1713+0747/puppi_J1713+0747_1400_profile.TF.FR.calibP"
#real_data_highSN_profs = np.load(real_data_highSN_profs_file)
#real_data_highSN_profs = real_data_highSN_profs/np.max(real_data_highSN_profs, axis=1).reshape(512,1)
#prof_models_interp_Lband = real_data_highSN_profs

# Define the pulse profiles used to get TOAs with PSRCHIVE
#May need this?
#temp_Lband = "templatefiles/J0340+4130/J0340+4130.Rcvr1_2.GUPPI.12y.x.sum.sm"
temp_Lband = pss.pulsar.GaussProfile(peak = 0.25, width = 0.01, amp = 1.0)
#Do we need to change these Chns?
temp_Lband.init_profiles(2048, Nchan = 1)



# We define a function to get TOAs from the simulated data
def getSimTOAs(fitsfiles, tempfile, scrunch = False, nchan = 512, nsubint = 1, npol = 1, nbins = 2048):
    """
    This function will take a single or list of fitsfiles of simulated data and run 
    PSRCHIVE calls to get TOAs from the simulated fits files, and scrunch the data to a 
    number of frequency channels and subintegrations is desired. It then also barycenters
    all of the TOAs (e.g. replaces the observatory code with '@') and also saves all the
    TOAs as one big file with "_ALL.tim" at the end. The inputs are as follows:
    fitsfiles [string] - a single or list of fits file names (and paths to) to get TOAs 
                        from
    tempfile [string] - Profile template fits file to make TOAs with (since we are only
                        using PSRCHIVE for now)
    scrunch [bool] - if False, don't manipulate the data, just get TOAs. If True, we will
              first scrunch the fits files to the number of frequency channels (nchan),
              subintegrations (nsubint), and polarizations (npol), given as input.
              The new fits files will be stored in the same place as the originals with
              the extension: .f'nchan't'nsubint'p'npol'
    """
    # If a single fitsfile string, put it into a list
    if isinstance(fitsfiles, str):
        fitsfiles = [fitsfiles]
    # Now check if we want to scrunch the files and if so do it
    if scrunch:
        # figure out the factor we need to scrunch; start with frequency
        if nchan == 1:
            freq_flag = " -F "
        else:
            freq_flag = " --setnchn %s " % (nchan)
        # then polarization
        if npol == 1:
            pol_flag = " -p "
        else:
            pol_flag = "" # don't thing we can scrunch to 2 from 4, not sure though
        # and subintegrations
        if nsubint == 1:
            sub_flag = " -T "
        else:
            sub_flag = " --setnsub %s " % (nsubint)
        # Set the bins flag
        bin_flag = " --setnbin %s " % (nbins)
        # Now put it all together
        scrunchfits = []
        ext = "f%st%sp%s" % (nchan, nsubint, npol)
        for ff in fitsfiles:
            scrunchcall = "pam -e " + ext + freq_flag + pol_flag + sub_flag + bin_flag + ff
            call(scrunchcall)
            if 'ASP' in ff:
                scrunchfits.append(ff.split(".")[0]+"."+ext)
            else:
                scrunchfits.append(ff.split(".")[0]+"."+ext)
        # Then reassign the fits files if we needed to scrunch the,
        fitsfiles = scrunchfits
    # Now once we've scrunched we get the TOAs # IPTA does not work with tempo2 on this machine for some reason
    TOAcall = "pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt  -C rcvr:name -C be:name -f tempo2 IPTA -s %s " % (tempfile)
    timfiles = []
    for ff in fitsfiles:
        print(ff)
        call(TOAcall+"%s > %s.tim" % (ff, ff))
        timfiles.append(ff+".tim")
    # now we need to tell the tim files that the TOAs are barycentred
    alltimlines = ["FORMAT 1 \n"]
    for t_f in timfiles:
        # Go through the tim file line by line
        lines = []
        with open(t_f, 'r') as tf:
            for line in tf:
                # Get other necessary lines
                if "FORMAT" in line:
                    lines.append(line)
                # replace the observatory codes and receiver name and add the lines
                elif 'ao' in line:
                    newline = line.replace('ao', '@')
                    newline = newline.replace("-rcvr:name lbw", "-fe L-wide")
                    newline = newline.replace("-rcvr:name 430", "-fe 430")
                    newline = newline.replace("-rcvr:name L-wide", "-fe L-wide")
                    if "-fe L-wide" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f L-wide_PUPPI")
                    elif "-fe 430" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f 430_PUPPI")
                    elif "-fe L-wide" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f L-wide_ASP")
                    elif "-fe 430" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f 430_ASP")
                    lines.append(newline)
                    alltimlines.append(newline)
                elif 'gbt' in line:
                    newline = line.replace('gbt', '@')
                    newline = newline.replace("-rcvr:name lbw", "-fe L-wide")
                    newline = newline.replace("-rcvr:name 430", "-fe 430")
                    newline = newline.replace("-rcvr:name L-wide", "-fe L-wide")
                    if "-fe L-wide" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f L-wide_PUPPI")
                    elif "-fe 430" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f 430_PUPPI")
                    elif "-fe L-wide" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f L-wide_ASP")
                    elif "-fe 430" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f 430_ASP")
                    lines.append(newline)
                    alltimlines.append(newline)
                elif line.split()[4] == '0':
                    split_line = line.split()
                    split_line[4] = '@'
                    newline = ""
                    for val in split_line:
                        if val == split_line[0]:
                            newline += val
                        else:
                            newline += " %s" % (val)
                    newline += " \n"
                    newline = newline.replace("-rcvr:name lbw", "-fe L-wide")
                    newline = newline.replace("-rcvr:name 430", "-fe 430")
                    newline = newline.replace("-rcvr:name L-wide", "-fe L-wide")
                    if "-fe L-wide" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f L-wide_PUPPI")
                    elif "-fe 430" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f 430_PUPPI")
                    elif "-fe L-wide" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f L-wide_ASP")
                    elif "-fe 430" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f 430_ASP")
                    lines.append(newline)
                    alltimlines.append(newline)
                else:
                    newline = line.replace("-rcvr:name lbw", "-fe L-wide")
                    newline = newline.replace("-rcvr:name 430", "-fe 430")
                    newline = newline.replace("-rcvr:name L-wide", "-fe L-wide")
                    if "-fe L-wide" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f L-wide_PUPPI")
                    elif "-fe 430" in newline and "PUPPI" in newline:
                        newline = newline.replace("-be:name PUPPI", "-f 430_PUPPI")
                    elif "-fe L-wide" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f L-wide_ASP")
                    elif "-fe 430" in newline and "ASP" in newline:
                        newline = newline.replace("-be:name ASP", "-f 430_ASP")
                    lines.append(newline)
                    alltimlines.append(newline)
            tf.close()
        # now we write out the new file
        with open(t_f, 'w') as nt:
            nt.writelines(lines)
            nt.close()
    # Now write one big file with all the TOAs in it
    fulltim = timfiles[0].split('.ti')[0]+'_ALL.tim'
    print("All simulated TOAs can be found in:", fulltim)
    with open(fulltim, 'w') as ft:
        ft.writelines(alltimlines)
        ft.close()
    print("11")
    print(fulltim)
    # and return the name of this last fits file
    return fulltim



# We define a function to combine multiple tim files generated by get_SimTOAs above
def combine_tim(timfiles, outfile = "J0000.tim"):
    """
    This function will take a list of tim files (e.g. the files generated by the 
    gen_SimTOAs function above) and combine them into one big tim file with the name 
    of the output file designated by the 'outfile' input. Inputs;
    
    timfiles [list of strings] : a list of tim files to combine
    outfile [string] : name of output file constiting of all combined TOAs
    """
    tim_lines = ["FORMAT 1 \n"]
    for timfile in timfiles:
        with open(timfile, 'r') as tf:
            for line in tf:
                if 'FORMAT' in line:
                    pass
                else:
                    tim_lines.append(line)
            tf.close()
    # Now save the file
    with open(outfile, 'w') as of:
        of.writelines(tim_lines)
        of.close()
    print("All simulated TOAs can be found in:", outfile)
    return outfile


# ### Assign Simulation Files and Values
# 
# Now we assign the template files and par files that will be used to generate the simualted data. We also will define the increment lengths we want for the simulated data sets.

"""
BJS: You'll need different appropriate template files for J1713, as well as a correctly edited par file.
"""
# Now lets make the L-band simulated data; We will loop through a few times
#Not sure if we need these lines?
tempfits_Lband = "/home/jdc0059/SimulatorProj/templatefiles/J0000+0000/guppi_58908_J0340+4130_0007_0001.fits"

parfile = "templatefiles/J0000+0000/J0000+0000.par"



"""
BJS: I don't think you'll want all epochs of J1713 to simulate, especially right at the beginning.
You should probably hard code a few things here, especially as these will result in the wrong observing
epochs anyways.
I recomment hardcoding the following for now:
    
epochs = 1
ref_MJD = 56000.0 
MJD_start = 56000.0

# Either
inc_len = 30.0
increments_Lband_all = np.linspace(MJD_start, MJD_start+inc_len*(epochs-1), 1)
# or for now just:
increments_Lband_all = [0.0]

All_Lband_DMs_const = np.repeat(dm, epochs)
 
BJS: Them comment out from here down to where I say stop:
"""
epochs = 1
ref_MJD = 56000.0
MJD_start = 56000.0

inc_len = 30.0
increments_Lband_all = [0.0]

All_Lband_DMs_const = np.repeat(dm, epochs)



# ## Simulating the Data
# 
# Now we have assigned basically all of the values we need to assign in order to simulate the data set. We will now run through and simulate the data set, saving the files to the appropriate places. 
# 
# If the folder we want to save the simulations in does not exist, we want to make the folder.
"""
Here we define all of the different simulation runs we want to do so we can easily just comment and uncomment
the different simulations. Note, we should be able to run multiple at once this way.
"""

"""
------------------------------------------------------
This first set are for no variations
------------------------------------------------------
"""
"""
"""
sim_folder = "sim_folder"

# Full (short) simulations for L-band observations
#fitsfiles_Lband = simulate_full(sim_folder, increments_Lband_all, f0_Lband, bw_Lband, Nf_Lband, f_samp,True, subintlen_Lband, period, flux_Lband, psr_name, ObsTime_Lband, All_Lband_DMs_const, tempfits_Lband, parfile, temp_Lband_ar, twoD = False, ref_MJD=56000.0, MJD_start=56000.0)
fitsfiles_Lband = simulate_full(sim_folder, increments_Lband_all, f0_Lband, bw_Lband, Nf_Lband, f_samp,True, subintlen_Lband, period, flux_Lband, psr_name, ObsTime_Lband, All_Lband_DMs_const, tempfits_Lband, parfile, gauss_prof.profiles[0], twoD = False, ref_MJD=56000.0, MJD_start=56000.0)
#fitsfiles_Lband = simulate_full(gauss_prof)

# Get the 1400 MHz TOAs
#timfile_Lband = getSimTOAs(fitsfiles_Lband, temp_Lband.profiles[0], scrunch = False, nchan = 64, nsubint = 1, npol = 1)
#full_timfile = combine_tim([timfile_Lband], outfile=sim_folder+"/Combined_tim_file.tim")
#print(full_timfile)

sim_folder = "DM_Vars"
