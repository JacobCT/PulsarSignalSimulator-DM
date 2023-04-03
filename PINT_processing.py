import pint.toa as toa
import pint.models as models
import pint.fitter as fitter
import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path

#t = toa.get_TOAs("/home/jdc0059/SimulatorProj/sim_folder/results/" + sys.argv[1] + "/" + sys.argv[2] + "_chns/fits_Lband_0_4.tim", ephem="DE436") # put tim file name in here

t = "/bowser/jdc0059/SimulatorProj/sim_folder/results/" + sys.argv[1] + "/" + sys.argv[2] + "_chns/fits_Lband_0_4.tim" #Put tim file name here

def NANOfreqs(t, badchans=[[380,423], [443,480], [980,1150], [1618,1630]],\
              timfile = True):
    if timfile == False:
        print("timfile is False")
        tim_freqs = t.get_freqs().value
        all_idxs = np.arange(0,len(tim_freqs))
        good_idxs = []
        for i in all_idxs:
            ADD = True
            for chanband in badchans:
                if tim_freqs[i] >= chanband[0] and tim_freqs[i] <= chanband[1]:
                    ADD = False
            if ADD == True:
                good_idxs.append(i)
        return good_idxs
    elif timfile == True:
        # assign the new tim file
        print("timefile is True")
        newtimfile = "/bowser/jdc0059/SimulatorProj/sim_folder/results/" + sys.argv[1] + "/" + sys.argv[2] + "_chns/fits_Lband_clean.tim"
        newlines = ["FORMAT 1 \n"]
        with open(t, 'r') as oldtf:
            next(oldtf)
            for line in oldtf:
                f = float(line.split()[1])
                rcvr = line.split()[-3]
                if rcvr == "430":
                    if f >= 423.0 and f <= 443.0:
                        newlines.append(line)
                    else:
                        newlines.append("# "+line)
                elif rcvr == "L-wide":
                    if f >= 1150.0 and f <= 1618.0:
                        newlines.append(line)
                    elif f >= 1630.0:
                        newlines.append(line)
                    else:
                        newlines.append("# "+line)
                else:
                    print("Unknown Reciever: %s" % (rcvr))
            oldtf.close()
        # Save the new file
        with open(newtimfile, 'w') as ntf:
            ntf.writelines(newlines)
            ntf.close()
            print("did output newtim?")
        return newtimfile #newtimfile

NANOidxs = NANOfreqs(t, timfile=True)
#t.select(NANOidxs)
#freqs = t.get_freqs() #Adding this to check.
#print(freqs)

#This like takes the file that was cleaned through the NANOfreqs function and re-inputs the data as a new variable (ntf) which now replaces the old variable (t).
ntf = toa.get_TOAs("/bowser/jdc0059/SimulatorProj/sim_folder/results/" + sys.argv[1] + "/" + sys.argv[2] + "_chns/fits_Lband_clean.tim")

m = models.get_model(sys.argv[3]) # put par file name in here
f = fitter.WLSFitter(ntf, m)  #(t, m)
f.model.DM.frozen = False # if all not fit (0's) in par file, do this (Brent Originaly set this to False)
print(f.get_fitparams())

f.fit_toas()

rs = f.resids.time_resids.to('us').value
mjds = f.toas.get_mjds().value
re_errs = f.toas.get_errors().value

#plt.errorbar(mjds, rs, yerr=re_errs, fmt = 'o')
#plt.xlabel("MJD")
#plt.ylabel("Residual [us]")
#plt.show()
#plt.close()

# Print just DM
print(f.model.DM)

# Print whole par file
print(f.model.as_parfile(), file=open("PINT_output.txt", "a"))
