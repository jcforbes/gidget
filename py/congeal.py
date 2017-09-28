import glob
import os
import pdb
import numpy as np

andirec = os.environ['GIDGETDIR']+'/analysis/' 

SIs = glob.glob( andirec + 'broadDistr24*_sampleInfo.txt' )

nchar = 4+4+7
ncol = 1000+24+1+200

arr = np.zeros((ncol, len(SIs) ))

for k,si in enumerate(SIs):
    thisSI = np.loadtxt(si)
    assert len(thisSI)==ncol-1000
    direc = si[:-nchar]
    name = direc[direc.rfind('/')+1:]
    input_fn = direc+'/'+name+'_inputRandomFactors.txt'
    inputs = np.zeros(1000)
    with open(input_fn, 'r') as f:
        f.readline()
        for i in range(1000):
            line = f.readline()
            if len(line.split())==5:
	        inputs[i] = float(line.split()[0])
            else:
                print "WARNING: not adding line with length ", len(line), line
    arr[:len(thisSI),k] = thisSI[:]
    arr[len(thisSI):,k] = inputs[:]

counter = np.ones(len(SIs))
successes = np.sum(counter[ arr[-1003,:]>0 ])
print "Success rate: ", successes, " of ",len(counter)
    
np.savetxt( 'broad24partial_to_lasso.txt', arr.T )


