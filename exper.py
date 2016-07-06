import copy
import os, argparse
import sys,glob
import shutil
import subprocess
import time
import math
import numpy as np
import pdb
import random
from bolshoireader import bolshoireader

# This is a script to allow you to run the GIDGET code with
# a wide variety of systematically varied parameters. The functions here are described
# in some detail, but if you skip to the end of the script, there are plenty of 
# examples of what you can run from the command line, e.g.
#  $ python exper.py ex1
# will execute the example used in the README file.


t0=time.time()
allModels={} # a global dictionary of all experiments created.
allProcs=[] # a global list of all processes we've started
allStartTimes=[] # a global list of the start times of all of these processes


globalBolshoiReader = bolshoireader('rf_registry3.txt',3.0e11,3.0e12, '/Users/jforbes/bolshoi/')
bolshoiSize =  len(globalBolshoiReader.keys)

def HowManyStillRunning(procs):
    ''' Given a list of processes created with subprocess.Popen, (each of which
          in this case will be a single run of the gidget code), count up how
          many of those runs are still going.'''
    nStillRunning=0
    for proc in procs:
        if(proc.poll()==None):
            nStillRunning+=1
    return nStillRunning 

def successCode(filename):
    ''' Given a file assumed to contain the standard error output of
    a run of the gidget code, return 0 for success (file is empty)
    1 for a time step below floor result, 2 for a problem with Qst,
    3 for a rejected accr. history, or 4 for a miscellaneous error code.  '''
    if (not os.path.exists(filename)):
      return -1 # if the file *stde.txt doesn't exist, this particular run has not been started.

    if os.stat(filename).st_size == 0:
      return 0 # the *_stde.txt file is empty, i.e. there were no errors (so far) for this run

    with open(filename,'r') as f:
        contents = f.read()
        if 'floor' in contents:
            return 1 # the time step has gotten very small. Typically
                     # this is caused by 1 cell in the gas column density
                     # distribution becoming very small.
        if 'Qst' in contents:  
            return 2 # a spike in the stellar column density distribution
                     # which in turn causes an interpolation error
        if 'merger' in contents:
            return 3 # accretion history rejected because of a merger
    return 4 # none of the most common types of failures
    

class experiment:
    ''' This class is a container to generate and store the necessary information to run
         a series of GIDGET models with parameters varied in a systematic way.'''
    def __init__(self,name):
        ''' All we need here is a name. This will be the directory containing and the prefix for
             all files created in all the GIDGET runs produced as part of this experiment.'''
        # fiducial model
        self.p=[name,200,1.5,.01,4.0, \
                1,1,.01,1,10, \
                220.0,20.0,7000.0,2.5,.5, \
                1.0,2.0,5000.0,int(1e9),1.0e-3, \
                1.0,-2.0/3.0,0,.5,2.0, \
                2.0,0,0.0,1.5,1, \
                2.0,1.0,1.0e12,5.0,3.0, \
                0,2.0,-1.0,2.5,0.0, \
                0.54,0.1,200,.30959,0.38, \
                -0.25,1.0,1.0,0.3,1.0, \
                0,0.0,.1,.03,2.0, \
                .002,.054,0.0,0.0, bolshoiSize/2]
        self.p_orig=self.p[:] # store a copy of p, possibly necessary later on.
        self.pl=[self.p[:]] # define a 1-element list containing a copy of p.
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat', \
                'analyticQ','cosmologyOn','xmin','NActive','NPassive', \
                'vphiR','R','gasTemp','Qlim','fg0', \
                'phi0','zstart','tmax','stepmax','TOL', \
                'muNorm','muColScaling','b','innerPowerLaw','softening', \
                'diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive', \
                'fixedQ','kappaMetals','Mh0','minSigSt','NChanges', \
                'dbg','accScaleLength','zquench','zrelax','xiREC', \
                'RfREC','deltaOmega','Noutputs','accNorm','accAlphaZ', \
                'accAlphaMh','accCeiling','fscatter','invMassRatio','fcool', \
                'whichAccretionProfile','alphaAccretionProfile','widthAccretionProfile','fH2Min','tDepH2SC', \
                'ZIGM','yREC','concentrationRandomFactor','muFgScaling', 'bolshoiWeight']
        assert len(self.p)==len(self.names)
        self.keys={}
        ctr=0
        for n in self.names:
            self.keys[n]=ctr
            ctr=ctr+1
        self.expName=name
	self.covariables=np.zeros(len(self.p),int)

        # store the location of various expected subdirectories in the gidget distribution.
        #self.base=os.getcwd() # Assume we are in the base directory - alter this to /path/to/gidget/directory if necessary
        self.base = os.environ['GIDGETDIR']
        self.src=self.base+'/src'
        self.analysis=self.base+'/analysis'
        self.bin=self.base+'/bin'
        self.out=self.base+'/output'
        self.xgrid=self.base+'/xgrid'

        allModels[name] = self
       
    def changeName(self,newName):
        ''' If this object is copied to form the starting point for a separate experiment,
            the self.expName, the names of each GIDGET run, and the entry in allModels will
            need to be updated. This method takes care of all of these tasks.'''
        oldName = self.expName

        self.p[0] = newName
        self.p_orig[0] = newName
        for p in self.pl:
            p[0] = newName + p[0][len(oldName):] # new prefix + old suffixes
        self.expName = newName
        allModels[newName] = self
 
    def vary(self,name,mmin,mmax,n,log,cov=0):
        '''Set up to vary a particular parameter, name, over the
        range [mmin,mmax] with n steps spaced linearly if log=0 or
        logarithmically if log=1. The cov flag means that all variables
        with the same value of cov, e.g. accScaleLength and whichAccretionHistory,
        will be varied together rather than independently. '''
        if(n>1 and log==0):
            self.p[self.keys[name]]=[mmin + (mmax-mmin)*type(self.p_orig[self.keys[name]])(i)/type(self.p_orig[self.keys[name]])(n-1) for i in range(n)]
        elif(n>1 and log==1 and mmin!=0 and mmin!=0.0):
            rat = (float(mmax)/float(mmin))
            typ = type(self.p_orig[self.keys[name]])
            self.p[self.keys[name]]=[typ(float(mmin) * (rat)**(float(i)/float(n-1))) for i in range(n)]
        elif(n==1):
            self.p[self.keys[name]]=mmin # a mechanism for changing the parameter from its default value.
    	self.covariables[self.keys[name]] = cov
    	return self.p[self.keys[name]]

    def multiply(self,name,multiple):
        keynum = self.keys[name]
        if(type(self.p[keynum]) == type([])): # we're dealing with a list
            self.p[keynum] = [a*multiple for a in self.p[keynum]]
        else:
            self.p[keynum] = [self.p[keynum] * multiple]
        return self.p[keynum]

    def irregularVary(self,name,values,cov=0):
        ''' Similar to vary, except, rather than generating the
        array of parameter values to use for parameter "name",
        this list is specified explicitly in "values"'''
        keyNum = self.keys[name] # the element of the list self.p which we're editing today.
        typ = type(self.p_orig[keyNum]) # what type of variable is this in a single run of gidget?
        if(type(values) == type([1])): # if values is a list..
            # replace the given index of p with the list given by values.
            self.p[keyNum]=[typ(value) for value in values]
        else:
           # even though we weren't given a list, put our new value in a single-parameter list.
            self.p[keyNum]=[typ(values)]
        # Next, we want to pass along the information we've been given on whether this parameter
        # should be covaried:
    	self.covariables[self.keys[name]] = cov
        # and finally, give back what we've done.
        return self.p[keyNum]
    
    def ConsistencyCheck(self):
        ''' Here we check whether all parameters in a covarying set have
	the same length. Otherwise it makes little sense to vary them together.
	'''
        consistent = True
        distinctCovSetIndices = set(self.covariables)
	for ind in distinctCovSetIndices:
		covIndicesForThisSet = np.where(self.covariables == ind)
		lengths=[]
		for var in covIndicesForThisSet:
			lengths.append(len(self.p[var]))
		if(len(set(lengths)) != 1):
			consistent=False
			print "Problem found in the ", ind, " set of covarying variables."
			print "This corresponds to the variables ", covIndicesForThisSet
			print "A set of variables you have asked to covary do not have the same lengths!"

    def generatePl(self):
        '''vary() will change a certain element of self.p into a list
        of values for that variable. This method assumes that some
        (perhaps many) variables have been varied in this manner, and
        expands p into a list of well-defined parameter lists, each
        one of which corresponds to a run of the program. This list can
        then be sent to the xgrid or run on your machine. See other
        methods for details on how to do that.'''
        self.pl=[self.p_orig[:]] # in case generatePl() has already been run, reset pl.
	# for each separate parameter (e.g. number of cells, dissipation rate, etc.)
        for j in range(len(self.p)): 
            param = self.p[j] 
	    # if instead of a single value of the parameter, we have a whole list, create a set of runs
	    # such that for every previous run we had, we now have len(param) more, each with a different
	    # value from the list param.
            if(type(param)==list):
		covaryOtherVarsWithJ=False
		varyThisVariable=True # i.e. the variable corresponding to j.
		cov = self.covariables[j] # the flag for covarying variables.
		if (cov!=0): # if this parameter (j) is part of a set of covarying parameters...
			covIndices = np.where(self.covariables == cov)[0]
			# if our current variable, j, is the first in a set of variables we're covarying,
			# let's go ahead and vary other the other variables in the same set.
                        
			if( covIndices[0] == j ): 
				covaryOtherVarsWithJ = True
			# whereas if this is not the case, then this variable has already been varied
			else:
				varyThisVariable = False
		
		if(varyThisVariable): 
	                # make len(param) copies of the current pl
	                pl2=[]
	                for i in range(len(param)):
	                    f=copy.deepcopy(self.pl)
	                    pl2.append(f) # pl2 ~ [pl,pl,pl]~[[p,p,p],[p,p,p],[p,p,p]]
	                    # in each copy, set the jth parameter in p to param[i]
	                    for a_p in pl2[i]: # each element is a list of parameters for
        	                a_p[j]=param[i]
				if(covaryOtherVarsWithJ): 
					for covIndex in covIndices:
						if(covIndex != j): # already taken care of with a_p[j]=...
							print covIndex, i, len(a_p), len(self.p), len(self.p[covIndex]) 
							a_p[covIndex] = self.p[covIndex][i]
                	        # in each copy, append to the name a...z corresponding to 
                            # which copy is currently being edited
	                        if(i<=25):
        	                    base='a'
                	            app=''
                                # If there are more than 26 variations, add a numeral for each
                                # time we have looped through the alphabet when this number is
                                # greater than zero.
                        	else:
	                            base='a'
        	                    app=str((i-(i%26))/26)
                                # avoid adding more letters to distinguish individual models
                                # if such a distinction is unnecessary.
                                if(len(param) > 1):
                                    # but if such a distinction is necessary, by all means:
                	            a_p[self.keys['name']]+=chr(ord(base)+(i%26))+app
	                    # end loop over p's in pl2[i]
        	        # end loop over range of this varied parameter
		
			# collapse pl to be a 1d array of p's.
	                self.pl=[element for sub in pl2 for element in sub]
	    # end of param-being-a-list contingency
            else : #just a regular non-varying parameter
                # If the user has used vary() but with N=1, the new argument will
                # not be a list! In this case, we set this parameter in each element 
		# of pl to be the value the user picked, instead of the default value.
                for i in range(len(self.pl)):
                    if(self.pl[i][j]!=param):
                        self.pl[i][j]=param
                    else:
                        pass # nothing to do- just use default value
            # end of non-varying parameter contingency 
        # end of for loop over each parameter.

            
    def write(self,name):
        '''Write a text file which can be handled by GridStuffer to
        run the experiment. This consists of a list of lists, each
        one being a set of command line arguments to give to the
        executable. Note that before running this, you need to run
        generatePl() to generate this list of lists, and before that,
        you should run vary() to set up which parameter(s) to change
        in each run. Once you have written this file, to use it on
        your local xgrid, you use GridStuffer to create a new metajob
        with the "Input File" as the text file, and the "Output Folder"
        as gidget/output'''
        self.generatePl()
        with open(name,'w') as f:
            for a_p in self.pl:
                f.write('xgrid -job submit -in '+self.bin+' -out '+self.analysis+'/'+self.expName+' ./gidget')
                for param in a_p:
                    f.write(' '+repr(param))
                if(self.pl[-1]!=a_p):
                    f.write('\n') # no return after the last line
        if(os.path.exists(self.analysis+self.expName)):
            print "Warning: this experiment already exists! It would be a good idea to rename your experiment or delete the pre-existing one manually."
        os.rename(name,self.xgrid+'/'+name) # move the text file to gidget/xgrid
        os.mkdir(self.xgrid+'/output/'+self.expName) # make gidget/xgrid/output/name to store stde and stdo files.
        print "Prepared file ",name," for submission to xGrid."

    def ExamineExperiment(self):
	self.generatePl()
        varied=[] # which elements of p are lists
	Nvaried=[] # how many variations for each such list.
	theLists=[] # a copy of each list.
        for j in range(len(self.p)): 
            param = self.p[j] 
            if(type(param)==list ):
		varied.append(j)
		Nvaried.append(len(param))
		theLists.append(param[:])
	successTable = np.ndarray(shape=tuple(Nvaried[:]),dtype=np.int64)
	successTableKeys = np.zeros(tuple(Nvaried[:]))
	for a_p in self.pl: # for each run..
	    successTableKey =[]
	    for i in range(len(varied)): # for each parameter that was varied
                j = varied[i]
                successTableKey.append(theLists[i].index(a_p[j]))
            successTable[tuple(successTableKey)] = successCode(self.analysis+'/'+self.expName+'/'+a_p[0]+"_stde.txt")
            successTableKeys[tuple(successTableKey)] = 1#tuple(successTableKey)
	print
	print "Success Table:"
        print successTable.T
	print 
        return successTable.T



    def makeList(self):
        """ This function is under construction. The idea is that when we want to run a large
          number of experiments in series, we don't want to wait for the first experiment to
          finish if, say, it is still using up 2 processors but there are 2 other processors free.
        """
        self.generatePl()
        ctr=0
        binary = self.bin+'/gidget'
        expDir = self.analysis+'/'+self.expName
        if(os.path.exists(expDir) and startAt == 0):
          print "Cancelling this run! ",self.expName
          return ([],[],[])
        elif(os.path.exists(expDir) and startAt != 0):
          pass
        else:
          os.mkdir(expDir)
        stdo=[]
        stde=[]
        cmds=[]
        expDirs=[]
        for a_p in self.pl[startAt:]:
            ctr += 1
            tmpap = a_p[:]
            stdo.append(expDir+'/'+a_p[self.keys['name']]+'_stdo.txt')
            stde.append(expDir+'/'+a_p[self.keys['name']]+'_stde_aux.txt')
            cmds.append([binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]])
            expDirs.append(expDir)

        return (cmds,stdo,stde,expDirs)


    def localRun(self,nproc,startAt,maxTime=36000,overwrite=False):
        ''' Run the specified experiment on this machine,
        using no more than nproc processors, and starting
        at the "startAt"th run in the experiment '''
        self.generatePl()
        procs=[]
        startTimes=[]
        ctr=0
        global allProcs
        binary=self.bin+'/gidget'
        expDir=self.analysis+'/'+self.expName #directory for output files
        #if(os.path.exists(expDir) and startAt == 0 and not overwrite):
        #print "dbg: ",startAt,overwrite,glob.glob(expDir+'/*comment.txt'),os.path.exists(expDir)
        if( len(glob.glob(expDir+'/*comment.txt'))!=0 and startAt == 0 and not overwrite ):
            print "************"
            print "This directory already contains output. CANCELLING this run!"
            print
            return
        elif(os.path.exists(expDir) and (startAt != 0 or overwrite or len(glob.glob(expDir+'/*comment.txt'))==0)):
            pass # let the user overwrite whatever runs they so desire.
        else:
            os.mkdir(expDir)

        for a_p in self.pl[startAt:]:
            ctr+=1
            tmpap=a_p[:]
            globalBolshoiReader.copynearest(expDir, a_p[self.keys['name']]+'_inputRandomFactors.txt', a_p[self.keys['bolshoiWeight']])
            with open(expDir+'/'+a_p[self.keys['name']]+'_stdo.txt','w') as stdo:
                with open(expDir+'/'+a_p[self.keys['name']]+'_stde_aux.txt','w') as stde:
                    print "Sending run #",ctr,"/",len(self.pl[startAt:])," , ",tmpap[0]," to a local core."
                    #print "Parameters: "
                    #print [binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]]
                    os.chdir(expDir)
                    procs.append(subprocess.Popen([binary]+tmpap[:1]+[repr(el) for el in tmpap[1:-1]],stdout=stdo,stderr=stde))
                    cmd = ''
                    for el in tmpap:
                        cmd+=str(el)+' '
                    print "Using cmd: ",cmd
                    allProcs.append(procs[-1])
                    startTimes.append(time.time())
                    allStartTimes.append(time.time())
                    nPrinted=True

            # we've started a process off and running, but there are
            # probably more processes waiting to go. We do not want
            # to exceed the number of processors nproc the user
            # is willing to let run at once, so let's check what
            # our processes are up to:
            while True:

                # Check whether each process has exceeded the max time:
                for k,proc in enumerate(procs):
                    if proc.poll()==None:
                        # The process is still running.
                        if(time.time() - startTimes[k] > maxTime):
                            print "WARNING: process has reached maximum allowed time of ",maxTime," seconds! Sending the kill signal."
                            print " If you're sure the process should be running longer, increase the maxTime argument of localRun."
                            print " Usually a run this long means something has gone horribly wrong."
                            with open(expDir+'/'+a_p[self.keys['name']]+'_stde_aux.txt','a') as stde:
                                stde.write('Reached max allowed time '+str(maxTime)+'\n')
                            proc.kill()


                nStillRunning=HowManyStillRunning(allProcs)
                nStillRunning2=HowManyStillRunning(procs)
                #print "dbg2",nStillRunning,nStillRunning2,expDir
    
                # if our processors are all booked, wait a minute and try again
                if(nStillRunning >= nproc):
                    time.sleep(10) # wait a little bit
                    if(nPrinted):
                        print "Waiting for a free processor..."
                        nPrinted=False # only print the above message if we haven't seen it since 
					# a new job has been submitted
                # whereas if we have free processors, loop back around to submit the next run. 
                # If there are no more runs, we're done!
                else:
                    break


def LocalRun(runBundle,nproc):
    cmds,stdo,stde,expDirs = runBundle
    for i in range(len(cmds)):
        with open(stdo[i],'w') as stdo_file:
            with open(stde[i],'w') as stde_file: 
                os.chdir(expDirs[i])
                procs.append(subprocess.Popen(cmds[i],stdout=stdo_file,stderr=stde_file))
                nPrinted=True
        while True:
            nStillRunning = HowManyStillRunning(procs)
            # if our processors are all booked, wait a minute and try again
            if(nStillRunning >= nproc):
                time.sleep(10) # wait a little bit
                if(nPrinted):
                    print "Waiting for a free processor..."
                    nPrinted=False # only print the above message if we haven't seen it since 
    	                            # a new job has been submitted
                   # whereas if we have free processors, loop back around to submit the next run. 
                   # If there are no more runs, we're done!
            else:
                break
    # now all of our processes have been sent off
    nPrev = 0
    while True:
        nStillRunning=HowManyStillRunning(procs)
        # has anything changed since the last time we checked?
        if(nStillRunning == nPrev and nStillRunning != 0):
            # do nothing except wait a little bit
            time.sleep(5)
        else:
            nPrev=nStillRunning
            if(nPrev == 0):
               break # we're done!
            print "Still waiting for ",nPrev, " processes to finish; I'll check every few seconds for changes."
    print "Local run complete!"


def GetScaleLengths(N,median=0.045,scatter=0.5,Mh0=1.0e12,sd=100,lower=0.0,upper=1.0e3,multiple=1.0):
    ''' Return a vector of scale lengths in kpc such that the haloes will have
        spin parameters of median, with a given scatter in dex. These 
        are also scaled by halo mass, since the Virial radius goes as M_h^(1/3).
        sd seeds the random number generator so that this function is deterministic.'''
    random.seed(sd)
    scLengths=[]
    while len(scLengths) < N:
	spinParameter = median  * (10.0 ** random.gauss(0.0,scatter)) # lognormal w/ scatter in dex
        if(type(Mh0) == type([1,2,3])): # if Mh0 is a list
            if(len(Mh0) == N):
                Mh=Mh0[len(scLengths)]
            else:
                print "You've given GetScaleLengths the wrong number of halo masses"
                Mh=1.0e12
        else:
            Mh=Mh0
        length = multiple * (311.34 * (Mh/1.0e12)**(1.0/3.0) * spinParameter/(2.0**.5))
        #length = 2.5*(Mh/1.0e12)**(1./3.) * (spinParameter / .045) 
        if(length > lower and length < upper):
	        scLengths.append(length)
    return scLengths



def PrintSuccessTables(successTables):
    ''' Take the results of all experiments named when this script was called
     and sum them up such that it's obvious which runs succeeded and which failed.
     A schematic example output for a situation where 5 different experiments were
     run, each of which varied the same 2 variables, 2 values for 1, 3 for the other:

       [ [ 00000,  00000 ], [ 00010, 00010], [00011, 00011] ]
	
     indicates that, as far as whether the runs crash or succeed,
     the 2-value variable doesn't make a difference, but for increasing
     values of the 3-value variable, fewer runs succeed. The numbers which
     appear here are explained in the successCode function above. '''

    sumOfSuccessTables = np.zeros(shape=successTables[0].shape,dtype=np.int64)
    for i in range(len(successTables)):
        table = successTables[i]
        sumOfSuccessTables += (table * (10**i))

    strSumOfSuccessTables = np.array(sumOfSuccessTables, dtype=str)
    it = np.nditer(sumOfSuccessTables, flags=['multi_index'])
    while not it.finished:
        strSumOfSuccessTables[it.multi_index] = str(it[0]).zfill(len(successTables)-1)
        it.iternext()

    print
    print "Here is the sum of all success tables, in homogenized form"
    print
    print strSumOfSuccessTables

def letter(i):
    return chr(ord("a")+i)

def NewSetOfExperiments(copyFrom, name, N=1):
    if(type(copyFrom)==type([1])):
#        if(N!=len(copyFrom)):
#            print "WARNING: you asked for ",N,"experiments copied from ",\
#                    name," containing ",len(copyFrom),"experiments. Ignoring ",N,"."
        theList = [copy.deepcopy(copyFrom[i]) for i in range(len(copyFrom))]
    else:
        theList=[copy.deepcopy(copyFrom) for i in range(N)]

    [theList[i].changeName(name+letter(i)) for i in range(len(theList))]
    return theList

if __name__ == "__main__":

    # The structure here is straightforward:
    #   1) First, parse the arguments
    #   2) Then define a bunch of experiments
    #   3) Run the experiments the user told us to run in the command line arguments

    # Read in the arguments. To see the results of this bit of code, just try to run
    # this script without any arguments, i.e. $ python exper.py
    parser = argparse.ArgumentParser(description='Analyze data and/or run experiments associated with a list of experiment names.')
    parser.add_argument('models', metavar='experiment', type=str, nargs='+',
                   help='a string contained in any experiment to be run, e.g. if rs04a rs04b rs05a are the models defined, then exper.py rs04 will run the first two, exper.py a will run the first and last, exper.py rs will run all three, exper.py rs04b rs05 will run the last two, etc.')
    parser.add_argument('--nproc',type=int,help="maximum number of processors to use (default: 16)",default=16)
    parser.add_argument('--start',metavar='startingModel',type=int,
                   help='The number of the model in the experiment (as ordered by GeneratePl) at which we will start sending experiments to the processors to run. (default: 0)',default=0)
    parser.add_argument('--xgrid',type=bool,help="run on an xGrid (requires the user to submit the generated file to an xGrid (default: False)",default=False)
    args = parser.parse_args()
    
    # Store the names of the experiments the user has told us to run.
    modelList=[]
    for modelName in args.models:
        modelList.append(modelName)

    # Begin defining experiments!
    re01=experiment("re01")
    re01.irregularVary("R",40)
    re01.irregularVary('accScaleLength', .05) # now a fraction of r200
    re01.irregularVary('NPassive',3)
    re01.irregularVary('dbg',2)
    re01.irregularVary('xmin',.01)
    re01.irregularVary('alphaMRI',.01)
    re01.irregularVary('fcool',0.6)
    re01.irregularVary('nx',200)
    re01.irregularVary('kappaMetals',1.0)
    re01.irregularVary('xiREC',0.0)
    re01.irregularVary('muNorm',1)

    re02=NewSetOfExperiments(re01,'re02')
    re02[0].irregularVary('muNorm',[.1,1,10])

    re03=NewSetOfExperiments(re01,'re03')
    re03[0].irregularVary('accScaleLength', [.01,.03,.1])

    re04=NewSetOfExperiments(re01,'re044')
    re04[0].irregularVary('NPassive',1)
    re04[0].irregularVary('zstart',[4.8],3)
    re04[0].irregularVary('zrelax',[5.0],3)
    re04[0].irregularVary('accScaleLength',[.01,.031,.1,.25],4)
    re04[0].irregularVary('R',[15,20,35,40],4)
    re04[0].irregularVary('alphaMRI',0)
    re04[0].irregularVary('muNorm',0.0)

    re05=NewSetOfExperiments(re01,'re05')
    re05[0].irregularVary('zstart',[4.9],3)
    re05[0].irregularVary('zrelax',[5.0],3)
    re05[0].irregularVary('accScaleLength',[.01,.031,.1],4)
    re05[0].irregularVary('R',[15,30,45],4)
    re05[0].irregularVary('muNorm',[.1,1,10])
    re05[0].irregularVary('eta',[.5,4.5])
    re05[0].irregularVary('fixedQ',[1.3,2.0,3.0])
    re05[0].irregularVary('fcool',[.1,.9])
    re05[0].irregularVary('fg0',[.3,.6,.9])
    re05[0].irregularVary('Mh0',[3.0e11,1.0e12,3.0e12])

    # Playing around with mass dependence
    re06=NewSetOfExperiments(re01,'re06666')
    re06[0].irregularVary('NPassive',1)
    re06[0].irregularVary('zstart',[4.99],3)
    re06[0].irregularVary('zrelax',[5.0],3)
    re06[0].irregularVary('accScaleLength',.05)
    re06[0].vary('Mh0',1.0e10,1.0e13,9,1,4)
    re06[0].vary('R',5,50,9,1,4)
    re06[0].irregularVary('fg0',.9999)
    re06[0].irregularVary('fcool',.05)
    re06[0].irregularVary('alphaMRI',0)
    re06[0].irregularVary('muNorm',0.0)
    re06[0].irregularVary('ZIGM',.00002)
    re06[0].irregularVary('whichAccretionProfile', 2)
    re06[0].irregularVary('widthAccretionProfile', 0.6)
    re06[0].irregularVary('dbg',2**4+2**1+2**0)

    # Adding some stochasticity before messing around with previous expt.
    re07=NewSetOfExperiments(re01,'re07')
    re07[0].irregularVary('NPassive',1)
    re07[0].irregularVary('zstart',[4.8],3)
    re07[0].irregularVary('zrelax',[5.0],3)
    re07[0].irregularVary('accScaleLength',.05)
    re07[0].vary('Mh0',1.0e10,1.0e13,200,1,4)
    re07[0].vary('R',5,50,200,1,4)
    re07[0].vary('whichAccretionHistory',-300,-101,0,4)
    re07[0].irregularVary('alphaMRI',0)
    re07[0].irregularVary('muNorm',0.0)

    # Use params from toy model.
    re08=NewSetOfExperiments(re01,'re08')
    re08[0].irregularVary('NPassive',1)
    re08[0].irregularVary('zstart',[4.99],3)
    re08[0].irregularVary('zrelax',[5.0],3)
    re08[0].irregularVary('accScaleLength',.05)
    re08[0].vary('Mh0',1.0e10,1.0e13,9,1,4)
    re08[0].vary('R',5,50,9,1,4)
    re08[0].irregularVary('fg0',.9999)
    re08[0].irregularVary('fcool',.05)
    re08[0].irregularVary('alphaMRI',0)
    re08[0].vary('muNorm',22,.21,9,1,4)
    re08[0].irregularVary('ZIGM',.00002)
    re08[0].irregularVary('dbg',2**4+2**1+2**0)

    re09=NewSetOfExperiments(re08[0],'re09')
    re09[0].irregularVary('dbg',2**4+2**1+2**0+2**12)

    # Modified the code to use mu = .5 (Mh/10^12)^(-2/3)
    re10=NewSetOfExperiments(re08[0],'re10')
    re10[0].irregularVary('accCeiling',0.6)
    re10[0].vary('R',15,100,9,1,4)
    re10[0].irregularVary('accScaleLength',0.1)
    re10[0].vary('fcool',0.05/2,0.05,9,1,4)
    #re10[0].irregularVary('accNorm',.1)
    #re10[0].irregularVary('dbg',2**4+2**1+2**0+2**

    # sample to study MS.
    re11=NewSetOfExperiments(re01,'re11')
    re11[0].vary('whichAccretionHistory',-300,-101,200,0,4)
    re11[0].vary('R', 15, 100, 200, 1, 4)
    re11[0].vary('fcool', 0.05/2, 0.05, 200, 1, 4)
    re11[0].vary('Mh0',1.0e10, 1.0e13, 200, 1, 4)
    re11[0].vary('muNorm', 22, .21, 200, 1, 4)
    re11[0].irregularVary('fscatter', .45)
    re11[0].irregularVary('NChanges', list(np.random.binomial(40,.25,size=200)) ,4)
    re11[0].irregularVary('accCeiling',0.6)
    re11[0].irregularVary('accScaleLength',0.1)
    re11[0].irregularVary('NPassive',1)
    re11[0].irregularVary('zstart',[4.99],3)
    re11[0].irregularVary('zrelax',[5.0],3)
    re11[0].irregularVary('fg0',.9999)
    re11[0].irregularVary('alphaMRI',0)
    re11[0].irregularVary('ZIGM',.00002)
    re11[0].irregularVary('dbg',2**4+2**1+2**0)

    # more reasonable ZIGM - are galaxies more in equilibrium?
    re12=NewSetOfExperiments(re10[0], 're12')
    re12[0].irregularVary('ZIGM',.002)
    re12[0].irregularVary('accCeiling',0.4)
    re12[0].irregularVary('accScaleLength',0.1)

    # MW progenitors
    #re13=NewSetOfExperiments(re01[0], 're13')

    # 
    re14=NewSetOfExperiments(re01, 're14')
    re14[0].irregularVary('xmin', [.00031, .001, .0031, .01, .031])
    re14[0].irregularVary('zstart',[4.95],3)
    re14[0].irregularVary('zrelax',[5.0],3)
    re14[0].irregularVary('dbg',2**4+2**1+2**0)


    re15 = NewSetOfExperiments(re01,'re15')
    re15[0].irregularVary('zstart',[4.95],3)
    re15[0].irregularVary('zrelax',[5.0],3)
    re15[0].irregularVary('dbg',2**4+2**1+2**0)
    re15[0].irregularVary('Mh0',[1.0e11,1.0e12,1.0e13,1.0e14])
    re15[0].irregularVary('accScaleLength',[.01,.03,.05,.1,.2])
    re15[0].irregularVary('concentrationRandomFactor',[-.3,-.1,0.0,.1,.3])


    #from mcmc_individual import emceeParameterSpaceToGidgetExperiment
    #favoriteParams = [  2.83269124e-01,   9.30881877e-03,   1.82351732e-01,   1.50469698e+00,
    #           3.24674422e-01,   1.85339104e+00,   7.43781644e-03,   1.94868075e-01,
    #              4.05944370e+12,   1.38101233e-01,  -2.36399004e-01,  -5.50757004e-01,
    #                 5.94584202e-01,  -6.06652441e-01,   8.88202578e+00,   3.92808177e-01,
    #                   -2.84901649e+00]

    #re16 = emceeParameterSpaceToGidgetExperiment(favoriteParams,'re16')[0]
    #re16.irregularVary("Noutputs",200)
    #allModels['re16']=re16

    # Barro plot...
    re17=NewSetOfExperiments(re01,'re17')
    re17[0].irregularVary('zstart',4.95)
    re17[0].irregularVary('zrelax',5.0)
    re17[0].irregularVary('dbg',2**4+2**1+2**0)
    re17[0].vary('whichAccretionHistory',-300,-101,200,0,4)
    re17[0].irregularVary('R', 10)
    re17[0].irregularVary('fcool', .8)
    re17[0].vary('Mh0',1.0e12, 1.0e13, 200, 1, 4)
    re17[0].irregularVary('muNorm', 0.5)
    re17[0].irregularVary('fscatter', .45)
    re17[0].irregularVary('NChanges', list(np.random.binomial(40,.25,size=200)) ,4)
    re17[0].irregularVary('accCeiling',0.6)
    re17[0].irregularVary('accScaleLength',0.1)
    re17[0].irregularVary('NPassive',1)
    re17[0].irregularVary('eta',0.5)

    re18=NewSetOfExperiments(re17,'re18')
    re18[0].irregularVary('R',30)
    re18[0].irregularVary('Noutputs',300)
    re18[0].vary('Mh0',3.0e12, 1.0e13, 200, 1, 4)
    re18[0].irregularVary('yREC',.02)
    re18[0].irregularVary('kappaMetals',.5)

    re19a=NewSetOfExperiments(re18,'re19a')
    re19a[0].irregularVary('accScaleLength', .01)
    re19a[0].irregularVary('R',10)

    re19b=NewSetOfExperiments(re18,'re19b')
    re19b[0].irregularVary('accScaleLength', .04)
    re19b[0].irregularVary('R',20)

    re19c=NewSetOfExperiments(re18,'re19c')
    re19c[0].irregularVary('accScaleLength', .09)
    re19c[0].irregularVary('R',30)

    re20=NewSetOfExperiments(re01,'re20')
    re20[0].irregularVary('Noutputs',600)
    re20[0].irregularVary('zstart',3.95)
    re20[0].irregularVary('zrelax',4.0)
    re20[0].irregularVary('dbg',2**4+2**1+2**0 )
    asls = np.array([.012,.027,.042,.057,.072])
    re20[0].irregularVary('accScaleLength',list(asls), 5)
    re20[0].irregularVary('R', list(asls * 20/0.012), 5)
    re20[0].irregularVary('fcool', .8)
    re20[0].irregularVary('Mh0', [1.0e13])
    re20[0].irregularVary('muNorm', 1.5)
    re20[0].irregularVary('muFgScaling', 0.4)
    re20[0].irregularVary('muColScaling', 0)
    re20[0].irregularVary('fscatter', .45)
    re20[0].irregularVary('accCeiling',0.6)
    re20[0].irregularVary('NPassive',1)
    re20[0].irregularVary('eta',0.5)

    re21=NewSetOfExperiments(re20,'re21')
    re21[0].irregularVary('dbg', 2**4+2**1+2**0 + 2**12 )

    re22=NewSetOfExperiments(re20,'re22')
    re22[0].irregularVary('accScaleLength', list(asls*4), 5)
    re22[0].irregularVary('fixedQ',1.5)
    re22[0].irregularVary('whichAccretionProfile', 2)
    re22[0].irregularVary('widthAccretionProfile', 0.5)
    re22[0].irregularVary('TOL',1.0e-4)

    re23=NewSetOfExperiments(re20,'re23')
    re23[0].irregularVary('Qlim',.05)
    re23[0].irregularVary('fixedQ',.05)

    re30=NewSetOfExperiments(re20,'re30')
    re30[0].irregularVary('Noutputs',600)
    re30[0].irregularVary('zstart',3.95)
    re30[0].irregularVary('zrelax',4.0)
    re30[0].irregularVary('dbg',2**4+2**1+2**0 + 2**14 )
    asls = np.array([.012,.027,.042,.057,.072])
    re30[0].irregularVary('accScaleLength',.042 )
    mhl = np.power(10.0, np.linspace(10,13,200))
    re30[0].irregularVary('R', list(.042 * 20/0.012 * np.power(mhl/1.0e12,1.0/3.0) ), 6)
    re30[0].irregularVary('fcool', .8)
    re30[0].irregularVary('Mh0', list(mhl), 6)
    re30[0].irregularVary('muNorm', list(1.5*np.power(mhl/1.0e12, -2.0/3.0)), 6)
    re30[0].irregularVary('muFgScaling', 0.4)
    re30[0].irregularVary('muColScaling', 0)
    re30[0].irregularVary('fscatter', .45)
    re30[0].irregularVary('accCeiling',0.6)
    re30[0].irregularVary('NPassive',1)
    re30[0].irregularVary('NChanges',30)
    re30[0].irregularVary('eta',0.5)
    re30[0].vary('whichAccretionHistory',-300,-101,200,0,6)


    # Same as re30 (scatter in accr. rate) PLUS scatter in scale length
    re31=NewSetOfExperiments(re20,'re31')
    re31[0].irregularVary('Noutputs',600)
    re31[0].irregularVary('zstart',3.95)
    re31[0].irregularVary('zrelax',4.0)
    re31[0].irregularVary('dbg',2**4+2**1+2**0 + 2**14 )
    asls = np.random.normal( np.random.normal(size=200) )*(0.042-0.027) + 0.042
    re31[0].irregularVary('accScaleLength', list(asls), 6)
    mhl = np.power(10.0, np.linspace(10,13,200))
    re31[0].irregularVary('R', list(asls * 25/0.012 * np.power(mhl/1.0e12,1.0/3.0) ), 6)
    re31[0].irregularVary('fcool', .8)
    re31[0].irregularVary('Mh0', list(mhl), 6)
    re31[0].irregularVary('muNorm', list(1.5*np.power(mhl/1.0e12, -2.0/3.0)), 6)
    re31[0].irregularVary('muFgScaling', 0.4)
    re31[0].irregularVary('muColScaling', 0)
    re31[0].irregularVary('fscatter', .45)
    re31[0].irregularVary('accCeiling',0.6)
    re31[0].irregularVary('NPassive',1)
    re31[0].irregularVary('NChanges',30)
    re31[0].irregularVary('eta',0.5)
    re31[0].vary('whichAccretionHistory',-300,-101,200,0,6)
    successTables=[]


    ## Try to actually get 1st order galaxy scaling relations right.
    ## Use low-ish time resolution and mass spacing, no scatter.
    ## Use relations from appendix of Hayward & Hopkins
    def Moster(Mh, mparams):
        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
        zti = 4.0
        logM1z = M10 + M11*zti/(zti+1.0)
        Nz = N10 + N11*zti/(zti+1.0)
        betaz = beta10 + beta11*zti/(zti+1.0)
        gammaz = gamma10 + gamma11*zti/(zti+1.0)
        M1 = np.power(10.0, logM1z)
        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
        return eff
    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329])
    Mhz0 = np.power(10.0, np.linspace(9.5, 12.3, 12))
    Mhz4 = Mhz0/(10.0 * np.power(Mhz0/1.0e12,0.14)) # very roughly...
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re40 = NewSetOfExperiments(re01,'re40')
    re40[0].irregularVary('fg0', list(fgz4), 5)
    re40[0].irregularVary('Noutputs',100)
    re40[0].irregularVary('nx', 300)
    re40[0].irregularVary('zstart',3.95)
    re40[0].irregularVary('zrelax',4.0)
    re40[0].irregularVary('dbg',2**4+2**1+2**0 )
    re40[0].irregularVary('accScaleLength',0.042)
    re40[0].irregularVary('R', list(reff4*51.0 ), 5)
    re40[0].irregularVary('xmin', .00003)
    re40[0].irregularVary('Mh0', list(Mhz0), 5) 
    re40[0].irregularVary('muNorm', list(.3*np.power(Mhz0/1.0e12,-0.6)), 5)
    fcools = 1.0/(1.0-fgz4) * mst/(0.17*Mhz4)
    re40[0].irregularVary('fcool', list(fcools), 5)
    re40[0].irregularVary('muFgScaling', 0.5)
    re40[0].irregularVary('muColScaling', 0.6)
    re40[0].irregularVary('fscatter', .45)
    re40[0].irregularVary('accCeiling',0.8)
    re40[0].irregularVary('NPassive',1)
    re40[0].irregularVary('eta',0.5)
    re40[0].irregularVary('yREC',0.04)
    re40[0].irregularVary('fixedQ', 1.25)
    re40[0].irregularVary('Qlim', 1.75)
    re40[0].irregularVary('concentrationRandomFactor', 0.7) ## units of dex! 
    re40[0].irregularVary('kappaMetals', list(.1*np.power(Mhz0/3.0e12,2.0/3.0)), 5)
    re40[0].irregularVary('ZIGM', list(ZHayward), 5)



    Mhz0 = np.power(10.0, np.linspace(9.5, 12.3, 240))
    Mhz4 = Mhz0/(10.0 * np.power(Mhz0/1.0e12,0.14)) # very roughly...
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re41 = NewSetOfExperiments(re01,'re41')
    re41[0].irregularVary('fg0', list(fgz4), 5)
    re41[0].irregularVary('Noutputs',200)
    re41[0].irregularVary('zstart',3.95)
    re41[0].irregularVary('zrelax',4.0)
    re41[0].irregularVary('dbg',2**4+2**1+2**0 )
    re41[0].irregularVary('accScaleLength',0.052)
    re41[0].irregularVary('R', list(np.power(reff4/reff4[-1],-0.1)*25), 5)
    re41[0].irregularVary('Mh0', list(Mhz0), 5) 
    re41[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.6)), 5)
    re41[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re41[0].irregularVary('muFgScaling', -0.3)
    re41[0].irregularVary('muColScaling', 0.2)
    re41[0].irregularVary('fscatter', .45)
    re41[0].irregularVary('accCeiling',0.5)
    re41[0].irregularVary('NPassive',1)
    re41[0].irregularVary('eta',0.5)
    re41[0].irregularVary('yREC',0.03)
    re41[0].irregularVary('fixedQ', 1.25)
    re41[0].irregularVary('Qlim', 1.75)
    re41[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re41[0].irregularVary('ZIGM', list(ZHayward), 5)
    asls = np.random.normal( np.random.normal(size=240) )*(0.042-0.027) + 0.052
    re41[0].irregularVary('accScaleLength', list(asls), 5)
    re41[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re41[0].vary('whichAccretionHistory',-340,-101,240,0,5)



    ## Try to set up some galaxies that look sort of like the enzo sims!
    ## Best guesses for many parameters based on gidgetmcmc fits, e.g. eta=0.5
    ## For the mass loading factor and xi, choose values based on the results
    ## of the simulations themselves. Note that these are for the 5 kpc run.
    dw01 = NewSetOfExperiments(re01, 'dw01')
    dw01[0].irregularVary('fg0', .95)
    dw01[0].irregularVary('Noutputs', 200)
    dw01[0].irregularVary('zstart', 3.95)
    dw01[0].irregularVary('zrelax', 4.0)
    dw01[0].irregularVary('dbg', 2**4+2**1+2**0)
    dw01[0].irregularVary('accScaleLength',0.052)
    dw01[0].irregularVary('R', 20)
    dw01[0].irregularVary('Mh0', 1.0e10) 
    dw01[0].irregularVary('muNorm', 6)
    dw01[0].irregularVary('fcool', 5.0e-4/.18)
    dw01[0].irregularVary('muFgScaling', 0)
    dw01[0].irregularVary('muColScaling', 0)
    dw01[0].irregularVary('accCeiling',0.5)
    dw01[0].irregularVary('NPassive',1)
    dw01[0].irregularVary('eta',0.5)
    dw01[0].irregularVary('yREC',0.03)
    dw01[0].irregularVary('xiREC',[.7,.8,.9])
    dw01[0].irregularVary('fixedQ', 2)
    dw01[0].irregularVary('Qlim', 3)
    dw01[0].irregularVary('kappaMetals', .1)
    dw01[0].irregularVary('ZIGM', 0.0006)
    dw01[0].irregularVary('whichAccretionHistory',0)
    dw01[0].irregularVary('alphaMRI',1.0e-4)

    # Try to tweak the above set of disks to look a bit more like the simulations.
    dw02=NewSetOfExperiments(dw01, 'dw02')
    dw02[0].irregularVary('xiREC',[0,.7,.8,.9])
    dw02[0].irregularVary('accScaleLength', [.026, .052, .104])
    dw02[0].irregularVary('accCeiling',.95)
    dw02[0].irregularVary('fcool', 2.0e-4/0.18)

    # A few more tweaks, plus... In dw02 above, I find that the SFR is essentially the same for all 12 simulations.
    # Examining them in a bit more detail suggests that this is because hte inner regions of the disks are gravitationally unstable.
    # Assuming this is true, perhaps the enzo sims haven't had enough time to adjust their gas profiles accordingly. To make the comparison
    # fair and see if we can get different rates of star formation, how about essentially turning off GI? Just to check what will happen
    dw03 = NewSetOfExperiments(dw02, 'dw03')
    dw03[0].irregularVary('accCeiling', 0.1)
    dw03[0].irregularVary('fixedQ', 0.1 )
    dw03[0].irregularVary('Qlim', 0.2 )
    # Some results: lower sfr. With this I realized that fH2min is still active in the code(!).

    # Correct the fH2min thing.
    dw04 = NewSetOfExperiments(dw03, 'dw04')
    dw04[0].irregularVary('accCeiling', 0.5)
    dw04[0].irregularVary('fH2Min', 0)
    dw04[0].irregularVary('ZIGM', 0.001)
    dw04[0].irregularVary('accScaleLength', [.013, .026, .052, .104, .208] )
    dw04[0].irregularVary('xiREC', [0, .3, .6] )
    dw04[0].irregularVary('muNorm', [1, 6])
    # Result: now we get quite large differences in SFR as we expected based on the 3D sims. How to proceed? Consider turning GI back on, and try to set up two different models that end up at about the same place as the sims at z=0

    # These models are:
    #  Weak ej feedback +       perfect mixing      + preventative feedback
    #  Intermediate feedback +  imperfect mixing    + no preventative 
    #  Strong feedback +        perfect mixing      + no preventative 
    dw05=NewSetOfExperiments(dw04, 'dw05')
    dw05[0].irregularVary( 'muNorm', [0.5, 6, 60], 3)
    dw05[0].irregularVary( 'accScaleLength', .2 )
    dw05[0].irregularVary( 'xiREC', [0, .6, 0] , 3)
    dw05[0].irregularVary( 'accCeiling', [0.2, 1, 1] , 3)



    # re52 and re53 are the same, but 53 is fed with multidark
    Mhz0 = np.power(10.0, np.linspace(11.5, 13.3, 300))
    Mhz4 = Mhz0/(10.0 * np.power(Mhz0/1.0e12,0.14)) # very roughly...
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re51 = NewSetOfExperiments(re01,'re53')
    re51[0].irregularVary('fg0', list(fgz4), 5)
    re51[0].irregularVary('Noutputs',200)
    re51[0].irregularVary('zstart',3.98)
    re51[0].irregularVary('zrelax',4.0)
    re51[0].irregularVary('dbg',2**4+2**1+2**0 )
    re51[0].irregularVary('accScaleLength',0.052)
    re51[0].irregularVary('R', list(np.power(reff4/reff4[-1],-0.1)*25), 5)
    re51[0].irregularVary('Mh0', list(Mhz0), 5) 
    re51[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.6)), 5)
    re51[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re51[0].irregularVary('muFgScaling', -0.3)
    re51[0].irregularVary('muColScaling', 0.2)
    re51[0].irregularVary('fscatter', 1.0)
    re51[0].irregularVary('accCeiling',0.5)
    re51[0].irregularVary('NPassive',1)
    re51[0].irregularVary('eta',0.5)
    re51[0].irregularVary('yREC',0.03)
    re51[0].irregularVary('fixedQ', 1.25)
    re51[0].irregularVary('Qlim', 1.75)
    re51[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re51[0].irregularVary('ZIGM', list(ZHayward), 5)
    asls = np.random.normal( np.random.normal(size=len(Mhz0)) )*(0.042-0.027) + 0.052
    re51[0].irregularVary('accScaleLength', list(asls), 5)
    #re51[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re51[0].irregularVary('NChanges', 301)
    #re51[0].vary('whichAccretionHistory',-340,-101,240,0,5)
    re51[0].irregularVary('whichAccretionHistory',-413)
    #re51[0].irregularVary('bolshoiWeight', bolshoiSize/4)
    bolweights = bolshoiSize * np.random.random(size=len(Mhz0))
    re51[0].irregularVary('bolshoiWeight', [int(bolweights[i]) for i in range(len(bolweights))], 5)

    # re54 is very similar to those above, but now fed with the ordinary Bouche accretion history ----------- this is also a good experiment: old-school variations in accretion histories to verify that our method of constructing them from multidark is reasonable.
    Mhz0 = np.power(10.0, np.linspace(11.5, 13.3, 300))
    Mhz4 = Mhz0/(10.0 * np.power(Mhz0/1.0e12,0.14)) # very roughly...
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re54 = NewSetOfExperiments(re01,'re54')
    re54[0].irregularVary('fg0', list(fgz4), 5)
    re54[0].irregularVary('Noutputs',200)
    re54[0].irregularVary('zstart',3.98)
    re54[0].irregularVary('zrelax',4.0)
    re54[0].irregularVary('dbg',2**4+2**1+2**0 )
    re54[0].irregularVary('accScaleLength',0.052)
    re54[0].irregularVary('R', list(np.power(reff4/reff4[-1],-0.1)*25), 5)
    re54[0].irregularVary('Mh0', list(Mhz0), 5) 
    re54[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.6)), 5)
    re54[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re54[0].irregularVary('muFgScaling', -0.3)
    re54[0].irregularVary('muColScaling', 0.2)
    re54[0].irregularVary('fscatter', 1.0)
    re54[0].irregularVary('accCeiling',0.5)
    re54[0].irregularVary('NPassive',1)
    re54[0].irregularVary('eta',0.5)
    re54[0].irregularVary('yREC',0.03)
    re54[0].irregularVary('fixedQ', 1.25)
    re54[0].irregularVary('Qlim', 1.75)
    re54[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re54[0].irregularVary('ZIGM', list(ZHayward), 5)
    asls = np.random.normal( np.random.normal(size=len(Mhz0)) )*(0.042-0.027) + 0.052
    re54[0].irregularVary('accScaleLength', list(asls), 5)
    #re54[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re54[0].irregularVary('NChanges', 301)
    #re54[0].vary('whichAccretionHistory',-340,-101,240,0,5)
    re54[0].irregularVary('whichAccretionHistory',0)
    #re54[0].irregularVary('bolshoiWeight', bolshoiSize/4)
    bolweights = bolshoiSize * np.random.random(size=len(Mhz0))
    re54[0].irregularVary('bolshoiWeight', [int(bolweights[i]) for i in range(len(bolweights))], 5)





    # re55 - fix up the run a bit. First, try to get rid of the high redshift cases where r_* is unresolved.
    Ngal = 15
    Mhz0 = np.power(10.0, np.linspace(11.5, 13.3, Ngal))
    Mhz4 = Mhz0/(.75*33.3 * np.power(Mhz0/1.5e13,0.359)) # very roughly... We know for Mh0 = 3e11, Mhz4 =4e10, and same for Mh0=2e13->Mhz4=6e11. Using those 4 numbers and assuming a powerlaw, we get this relation.
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re55 = NewSetOfExperiments(re01,'re55')
    re55[0].irregularVary('fg0', list(fgz4), 5)
    re55[0].irregularVary('Noutputs',200)
    re55[0].irregularVary('zstart',3.98)
    re55[0].irregularVary('zrelax',4.0)
    re55[0].irregularVary('dbg',2**4+2**1+2**0 )
    #re55[0].irregularVary('accScaleLength',0.104) # twice the canonical value I've been using
    re55[0].irregularVary('R', list(np.power(reff4/reff4[-1],0.2)*40), 5)
    re55[0].irregularVary('Mh0', list(Mhz0), 5) 
    re55[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.2)), 5)
    re55[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re55[0].irregularVary('muFgScaling', -0.1)
    re55[0].irregularVary('muColScaling', 0.1)
    re55[0].irregularVary('fscatter', 1.0)
    re55[0].irregularVary('accCeiling',0.5)
    re55[0].irregularVary('NPassive',1)
    re55[0].irregularVary('eta',0.5)
    re55[0].irregularVary('xmin',0.005)
    re55[0].irregularVary('yREC',0.03)
    re55[0].irregularVary('fixedQ', 1.25)
    re55[0].irregularVary('Qlim', 1.75)
    re55[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re55[0].irregularVary('ZIGM', list(ZHayward), 5)
    width = (0.042-0.027)/2.0
    asls = np.random.normal( np.random.normal(size=len(Mhz0)) )*width + 0.104
    re55[0].irregularVary('accScaleLength', list(asls), 5)
    #re51[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re55[0].irregularVary('NChanges', 301)
    #re51[0].vary('whichAccretionHistory',-340,-101,240,0,5)
    re55[0].irregularVary('whichAccretionHistory',-413)
    #re51[0].irregularVary('bolshoiWeight', bolshoiSize/4)
    bolweights = bolshoiSize * np.random.random(size=len(Mhz0))
    re55[0].irregularVary('bolshoiWeight', [int(bolweights[i]) for i in range(len(bolweights))], 5)




    # re56 - fix up the run a bit. First, try to get rid of the high redshift cases where r_* is unresolved.
    Ngal = 15
    Mhz0 = np.power(10.0, np.linspace(11.5, 13.3, Ngal))
    Mhz4 = Mhz0/(.75*33.3 * np.power(Mhz0/1.5e13,0.359)) # very roughly... We know for Mh0 = 3e11, Mhz4 =4e10, and same for Mh0=2e13->Mhz4=6e11. Using those 4 numbers and assuming a powerlaw, we get this relation.
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re56 = NewSetOfExperiments(re01,'re56')
    re56[0].irregularVary('fg0', list(fgz4), 5)
    re56[0].irregularVary('Noutputs',200)
    re56[0].irregularVary('zstart',3.98)
    re56[0].irregularVary('zrelax',4.0)
    re56[0].irregularVary('dbg',2**4+2**1+2**0+2**13 )
    #re56[0].irregularVary('accScaleLength',0.104) # twice the canonical value I've been using
    re56[0].irregularVary('R', list(np.power(reff4/reff4[-1],0.2)*40), 5)
    re56[0].irregularVary('Mh0', list(Mhz0), 5) 
    re56[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.2)), 5)
    re56[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re56[0].irregularVary('muFgScaling', -0.1)
    re56[0].irregularVary('muColScaling', 0.1)
    re56[0].irregularVary('fscatter', 1.0)
    re56[0].irregularVary('accCeiling',0.5)
    re56[0].irregularVary('NPassive',1)
    re56[0].irregularVary('eta',0.5)
    re56[0].irregularVary('xmin',0.005)
    re56[0].irregularVary('yREC',0.03)
    re56[0].irregularVary('fixedQ', 1.25)
    re56[0].irregularVary('Qlim', 1.75)
    re56[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re56[0].irregularVary('ZIGM', list(ZHayward), 5)
    width = (0.042-0.027)/2.0
    asls = np.random.normal( np.random.normal(size=len(Mhz0)) )*width + 0.104
    re56[0].irregularVary('accScaleLength', list(asls), 5)
    #re51[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re56[0].irregularVary('NChanges', 301)
    #re56[0].vary('whichAccretionHistory',-340,-101,240,0,5)
    re56[0].irregularVary('whichAccretionHistory',-413)
    #re56[0].irregularVary('bolshoiWeight', bolshoiSize/4)
    bolweights = bolshoiSize * np.random.random(size=len(Mhz0))
    re56[0].irregularVary('bolshoiWeight', [int(bolweights[i]) for i in range(len(bolweights))], 5)


    # re57 - continue to debug newly-discovered issue with rotation curves
    Ngal = 15
    Mhz0 = np.power(10.0, np.linspace(11.5, 13.3, Ngal))
    Mhz4 = Mhz0/(.75*33.3 * np.power(Mhz0/1.5e13,0.359)) # very roughly... We know for Mh0 = 3e11, Mhz4 =4e10, and same for Mh0=2e13->Mhz4=6e11. Using those 4 numbers and assuming a powerlaw, we get this relation.
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.25)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    ZHayward = -8.69 + 9.09*np.power(1.0+4.0,-0.017) - 0.0864*np.power(np.log10(mst) - 11.07*np.power(1.0+4.0,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
    re57 = NewSetOfExperiments(re01,'re57')
    re57[0].irregularVary('fg0', list(fgz4), 5)
    re57[0].irregularVary('Noutputs',200)
    re57[0].irregularVary('zstart',3.98)
    re57[0].irregularVary('zrelax',4.0)
    re57[0].irregularVary('dbg',2**4+2**1+2**0+2**13 )
    #re57[0].irregularVary('accScaleLength',0.104) # twice the canonical value I've been using
    re57[0].irregularVary('R', list(np.power(reff4/reff4[-1],0.2)*60), 5)
    re57[0].irregularVary('Mh0', list(Mhz0), 5) 
    re57[0].irregularVary('muNorm', list(.2*np.power(Mhz0/1.0e12,-0.2)), 5)
    re57[0].irregularVary('fcool', list(1.0/(1.0-fgz4) * mst/(0.17*Mhz4)), 5)
    re57[0].irregularVary('muFgScaling', -0.1)
    re57[0].irregularVary('muColScaling', 0.1)
    re57[0].irregularVary('fscatter', 1.0)
    re57[0].irregularVary('accCeiling',0.5)
    re57[0].irregularVary('NPassive',1)
    re57[0].irregularVary('eta',0.5)
    re57[0].irregularVary('xmin',0.005)
    re57[0].irregularVary('yREC',0.03)
    re57[0].irregularVary('fixedQ', 1.25)
    re57[0].irregularVary('Qlim', 1.75)
    re57[0].irregularVary('kappaMetals', list(np.power(Mhz0/3.0e12,1.0/3.0)), 5)
    re57[0].irregularVary('ZIGM', list(ZHayward), 5)
    width = (0.042-0.027)/2.0
    asls = np.random.normal( np.random.normal(size=len(Mhz0)) )*width + 0.052
    re57[0].irregularVary('accScaleLength', list(asls), 5)
    #re51[0].irregularVary('NChanges', list(np.random.binomial(19,.53,size=240)), 5)
    re57[0].irregularVary('NChanges', 301)
    #re56[0].vary('whichAccretionHistory',-340,-101,240,0,5)
    re57[0].irregularVary('whichAccretionHistory',-413)
    #re56[0].irregularVary('bolshoiWeight', bolshoiSize/4)
    bolweights = bolshoiSize * np.random.random(size=len(Mhz0))
    re57[0].irregularVary('bolshoiWeight', [int(bolweights[i]) for i in range(len(bolweights))], 5)


    # re58: try turning off GI - confirm that it's the cause of the oscillations.
    re58= NewSetOfExperiments(re57,'re58')
    re58[0].irregularVary( 'fixedQ', 0.1)

    # re59
    re59 = NewSetOfExperiments(re58,'re59')
    re59[0].irregularVary( 'fixedQ', .01)
    re59[0].irregularVary( 'Qlim', 0.015 )
    
    # re60 - turn off stars basically. Confirm that the instability still occurs in a gas-only disk. If it does, as we expect, then we can check whether it can be killed with a large alpha viscosity.
    re60 = NewSetOfExperiments(re57, 're60')
    re60[0].irregularVary( 'Qlim', 0 )
    re60[0].irregularVary( 'fixedQ', 1.5 )
    re60[0].irregularVary( 'fg0', .999999 )
    re60[0].irregularVary( 'epsff', 1.0e-7 )
    re60[0].irregularVary( 'muNorm', 0)

    # re61 - now check that we can kill the oscillation with high alpha.
    re61 = NewSetOfExperiments(re60, 're61')
    re61[0].irregularVary( 'alphaMRI', 0.2 )

    re62 = NewSetOfExperiments( re57, 're62' )
    re62[0].irregularVary( 'nx', 256 ) # FFT o'clock.


    re63 = NewSetOfExperiments( re62, 're63')
    re63[0].irregularVary('dbg',2**4+2**1+2**0+2**13 + 2**7)



    for inputString in modelList: # aModelName will therefore be a string, obtained from the command-line args
        # Get a list of all defined models (allModels.keys())
        # for each such key (aModel) check whether this inputString is contained in its name
        matches = [aModel for aModel in sorted(allModels.keys()) if inputString in aModel]
        if(len(matches) != 0): 
            for model in matches: #if(model in allModels): 
                print "dbg0",inputString,matches
                if(not args.xgrid): #local run
                    allModels[model].localRun(args.nproc,args.start)
                else: # write a file to run on the xgrid
                    allModels[model].write('runExperiment_'+model+'.txt')
        else:
            print "You asked me to run ",inputString," but did not define it in the script."


    # now all of our processes have been sent off
    nPrev = 0
    while True:
        nStillRunning=HowManyStillRunning(allProcs)
        # has anything changed since the last time we checked?
        if(nStillRunning == nPrev and nStillRunning != 0):
            # do nothing except wait a little bit
            time.sleep(2)
        else:
            nPrev=nStillRunning
            if(nPrev == 0):
                break # we're done!
            print "Still waiting for ",nPrev, " processes to finish; I'll check every few seconds for changes."
    print "All local runs complete!"


    print "Time elapsed (seconds): ", time.time()-t0
