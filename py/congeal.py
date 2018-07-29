import glob
import os
import pdb
import numpy as np

andirec = os.environ['GIDGETDIR']+'/analysis/' 

#bn = 'broadDistr133'
bn = 'broadPrior162a'
#bn = 'broadPrior160'
SIs = glob.glob( andirec + bn+'*_sampleInfo.txt' )

nchar = 4+4+7
ncol = 1000+24+1+200

arr = np.zeros((ncol, len(SIs) ))

for k,si in enumerate(SIs):
    thisSI = np.loadtxt(si)
    try:
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
    
        if k%3000==10:
            print "Progress! ",k,"/",len(SIs)
            np.savetxt( bn+'_to_lasso.txt', arr[:,:k+1].T )
    
    	    successes = np.sum( arr[-1003,:]>0 )
    	    print "Success rate: ", successes, " of ",k+1
    except:
        print "Failed at k=",k
        
np.savetxt( bn+'_to_lasso.txt', arr.T )


