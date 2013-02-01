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

# This is a script to allow you to run the GIDGET code with
# a wide variety of systematically varied parameters. The functions here are described
# in some detail, but if you skip to the end of the script, there are plenty of 
# examples of what you can run from the command line, e.g.
#  $ python exper.py ex1
# will execute the example used in the README file.


t0=time.time()
allModels={} # a global dictionary of all experiments created.
allProcs=[] # a global list of all processes we've started

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
        self.p=[name,200,1.5,.01,4.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,1.0,2.0,50.0,int(1e9),1.0e-3,1.0,0,.5,2.0,2.0,0,0.0,1.5,1,2.0,.001,1.0e12,5.0,3.0,0,2.0,-1.0,2.5,1.0,0.46,0.1,200,.30959,0.38,-0.25,1.0,1.0,0.3,1.0,0,0.0,.1]
        self.p_orig=self.p[:] # store a copy of p, possibly necessary later on.
        self.pl=[self.p[:]] # define a 1-element list containing a copy of p.
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat','analyticQ','cosmologyOn','xmin','NActive','NPassive','vphiR','R','gasTemp','Qlim','fg0','phi0','zstart','tmax','stepmax','TOL','mu','b','innerPowerLaw','softening','diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive','fixedQ','kappaMetals','Mh0','minSigSt','NChanges','dbg','accScaleLength','zquench','zrelax','zetaREC','RfREC','deltaOmega','Noutputs','accNorm','accAlphaZ','accAlphaMh','accCeiling','fscatter','invMassRatio','fcool','whichAccretionProfile','alphaAccretionProfile','widthAccretionProfile']
        self.keys={}
        ctr=0
        for n in self.names:
            self.keys[n]=ctr
            ctr=ctr+1
        self.expName=name
	self.covariables=np.zeros(len(self.p),int)

        # store the location of various expected subdirectories in the gidget distribution.
        self.base=os.getcwd() # Assume we are in the base directory - alter this to /path/to/gidget/directory if necessary
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


    def localRun(self,nproc,startAt):
        ''' Run the specified experiment on this machine,
        using no more than nproc processors, and starting
        at the "startAt"th run in the experiment '''
        self.generatePl()
        procs=[]
        ctr=0
        binary=self.bin+'/gidget'
        expDir=self.analysis+'/'+self.expName #directory for output files
        if(os.path.exists(expDir) and startAt == 0):
            print "This directory already contains output. CANCELLING this run!"
            print
            print "This is fine if you've already run the desired models and are"
            print "just interested in looking at the results of ExamineExperiment."
            print "Otherwise, just delete the appropriate directory and run this"
            print "script again, or specify a nonzero starting run number."
            return
        elif(os.path.exists(expDir) and startAt != 0):
            pass # let the user overwrite whatever runs they so desire.
        else:
            os.mkdir(expDir)

        for a_p in self.pl[startAt:]:
            ctr+=1
            tmpap=a_p[:]
            with open(expDir+'/'+a_p[self.keys['name']]+'_stdo.txt','w') as stdo:
                with open(expDir+'/'+a_p[self.keys['name']]+'_stde_aux.txt','w') as stde:
                    print "Sending run #",ctr,"/",len(self.pl[startAt:])," , ",tmpap[0]," to a local core."
                    #print "Parameters: "
                    #print [binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]]
                    os.chdir(expDir)
                    procs.append(subprocess.Popen([binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]],stdout=stdo,stderr=stde))
                    allProcs.append(procs[-1])
                    nPrinted=True

            # we've started a process off and running, but there are
            # probably more processes waiting to go. We do not want
            # to exceed the number of processors nproc the user
            # is willing to let run at once, so let's check what
            # our processes are up to:
            while True:
                nStillRunning=HowManyStillRunning(allProcs)
    
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
##        # now all of our processes have been sent off
##        nPrev = 0
##        while True:
##            nStillRunning=HowManyStillRunning(allProcs)
##            # has anything changed since the last time we checked?
##            if(nStillRunning == nPrev and nStillRunning != 0):
##                # do nothing except wait a little bit
##                time.sleep(5)
##            else:
##                nPrev=nStillRunning
##                if(nPrev == 0):
##                    break # we're done!
##                print "Still waiting for ",nPrev, " processes to finish; I'll check every few seconds for changes."
##        print "Local run complete!"


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


    # A vanilla setup: 11 models equally spaced in log Mh0.
    rs01=experiment("rs01")
    Mh0s = rs01.vary("Mh0",1.0e10,1.0e12,11,1,3)
    Mh0s11 = Mh0s
    rs01.irregularVary('R',40)
    rs01.irregularVary('accScaleLength',GetScaleLengths(11,Mh0=Mh0s,scatter=1.0e-10),3)
    rs01.irregularVary('diskScaleLength',GetScaleLengths(11,Mh0=Mh0s,scatter=1.0e-10),3)
    rs01.irregularVary('mu',list(1.0*(np.array(Mh0s)/1.0e12)**(-1./3.)),3)
    rs01.irregularVary('vphiR',list(220.0*(np.array(Mh0s)/1.0e12)**(1./3.)),3)
    rs01.irregularVary('NPassive',1)

    # Another vanilla setup: 101 models equally spaced in log Mh0 with varying accr. histories.
    # I'm hoping this will be a minimal set of models which can test the functionality of variability3.
    rs02=experiment("rs02")
    Mh0s = rs02.vary("Mh0",1.0e10,1.0e12,101,1,3)
    Mh0s101 = Mh0s
    rs02.irregularVary("R",40)
    rs02.irregularVary('accScaleLength',GetScaleLengths(101,Mh0=Mh0s,scatter=1.0e-10),3)
    rs02.irregularVary('diskScaleLength',GetScaleLengths(101,Mh0=Mh0s,scatter=1.0e-10),3)
    rs02.vary('whichAccretionHistory',1000,1100,101,0,3)
    rs02.irregularVary('mu',list(1.0*(np.array(Mh0s)/1.0e12)**(-1./3.)),3)
    rs02.irregularVary('vphiR',list(220.0*(np.array(Mh0s)/1.0e12)**(1./3.)),3)
    rs02.irregularVary('NPassive',1)


    # What I'd like to do here is compare Bouche+ to Neistein+ w/ fscatter = 0
    rs03=[copy.deepcopy(rs01) for i in range(4)]
    [rs03[i].changeName("rs03"+letter(i)) for i in range(len(rs03))]
    # no scatter - should be the median case
    rs03[1].irregularVary('fscatter',0.0)
    rs03[1].vary('whichAccretionHistory',1000,1010,11,0,3)
    # small scatter
    rs03[2].irregularVary('fscatter',.03)
    rs03[2].vary('whichAccretionHistory',1000,1010,11,0,3)
    # no scatter again, but with smaller Delta omega.
    rs03[3].irregularVary('fscatter',0.0)
    rs03[3].vary('whichAccretionHistory',1000,1010,11,0,3)
    rs03[3].irregularVary('deltaOmega',.01)

    # A simple comparison between low and high angular momentum systems w/ realistic accr scatter
    rs04=[copy.deepcopy(rs02) for i in range(2)]
    values=[0.5, 5.0]
    [rs04[i].changeName("rs04"+letter(i)) for i in range(len(rs04))]
    for i in range(len(values)):
        rs04[i].irregularVary('accScaleLength', \
	  list(np.array(GetScaleLengths(101,Mh0=Mh0s101,scatter=1.0e-10))*values[i]),3)

    # some runs in rs04b were crashing because the gas added per time step at very large radii
    # was very large compared to the initial column density there. rs05 implements a fix wherein
    # the simulation executes the first 2 time steps without altering dt.
    rs05=[copy.deepcopy(rs04[i]) for i in range(len(rs04))]
    [rs05[i].changeName("rs05"+letter(i)) for i in range(len(rs05))]

    # see what happens when we take our default model and artificially alter fscatter.
    fscatters = [0.2,1.0,1.5]
    rs06=[copy.deepcopy(rs02) for i in range(len(fscatters))]
    [rs06[i].changeName("rs06"+letter(i)) for i in range(len(rs06))]
    [rs06[i].irregularVary('fscatter',fscatters[i]) for i in range(len(rs06))]

    # With a smooth accretion history, see what effect we get from Q-dep forcing terms.
    rs07=[copy.deepcopy(rs01) for i in range(2)]
    [rs07[i].changeName("rs07"+letter(i)) for i in range(len(rs07))]
    rs07[1].irregularVary("dbg",2**4)

    # Look at a different scaling of mu with halo mass
    rs08 = copy.deepcopy(rs02)
    rs08.changeName("rs08")
    rs08.irregularVary('mu',1.0)

    # vary kappa.
    kappa = [3.0e-4, 3.0e-3]
    dbgs = [0, 2**15]
    rs09 = [copy.deepcopy(rs02) for i in range(len(kappa)*len(dbgs))]
    [rs09[i].changeName("rs09"+letter(i)) for i in range(len(rs09))]
    [rs09[i].irregularVary('kappaMetals',kappa[i/2]) for i in range(len(rs09))]
    [rs09[i].irregularVary('dbg',dbgs[1]) for i in [1,3]]

    # vary metal loading efficiency yields
    rs10 = copy.deepcopy(rs02)
    rs10.changeName("rs10")
    rs10.irregularVary('zetaREC',2.0)

    # vary kappa, assuming high angular momentum inflows
    kappas=[1.0e-5,1.0e-4,1.0e-3,3.0e-3]
    rs11 = [copy.deepcopy(rs04[1]) for i in range(len(kappas))]
    [rs11[i].changeName("rs11"+letter(i)) for i in range(len(rs11))]
    [rs11[i].irregularVary("kappaMetals",kappas[i]) for i in range(len(rs11))]

    # just decay based on the initial conditions!
    rs12 = [copy.deepcopy(rs04[i]) for i in range(len(rs04))]
    [rs12[i].changeName("rs12"+letter(i)) for i in range(len(rs12))]
    [rs12[i].irregularVary("zrelax",2.01) for i in range(len(rs12))]
    [rs12[i].irregularVary("zquench",1.8) for i in range(len(rs12))]

    # try another round of varying fscatter. This time do not restrict the mass ratio so much, and try smaller values.
    fscatters = [0.1,0.3,0.6,1.0]
    rs13 = [copy.deepcopy(rs02) for i in range(len(fscatters))]
    [rs13[i].changeName("rs13"+letter(i)) for i in range(len(rs13))]
    [rs13[i].irregularVary('fscatter',fscatters[i]) for i in range(len(rs13))]
    [rs13[i].irregularVary('invMassRatio',1.0) for i in range(len(rs13))]

    # Try lognormal variation
    rs14= [copy.deepcopy(rs02) for i in range(len(fscatters))]
    [rs14[i].changeName("rs14"+letter(i)) for i in range(len(rs14))]
    [rs14[i].irregularVary('fscatter',fscatters[i]) for i in range(len(rs14))]
    [rs14[i].irregularVary('NChanges',33) for i in range(len(rs14))]
    [rs14[i].irregularVary('accNorm',1.0e4) for i in range(len(rs14))]
    [rs14[i].vary('whichAccretionHistory',-1100,-1000,101,0,3) for i in range(len(rs14))]

    rs15=[copy.deepcopy(rs14[i]) for i in range(len(rs14))]
    [rs15[i].changeName("rs15"+letter(i)) for i in range(len(rs15))]
    [rs15[i].irregularVary('NChanges',5) for i in range(len(rs15))]

    rs16=[copy.deepcopy(rs15[i]) for i in range(len(rs15))]
    [rs16[i].changeName("rs16"+letter(i)) for i in range(len(rs16))]
    [rs16[i].irregularVary('NChanges',100) for i in range(len(rs16))]
    
    # The next few are analogous to the previous few, with much reduced accretion efficiencies.
    rs17=[copy.deepcopy(rs14[i]) for i in range(len(rs14))]
    [rs17[i].changeName("rs17"+letter(i)) for i in range(len(rs17))]
    [rs17[i].irregularVary("accCeiling",.05) for i in range(len(rs17))]

    rs18=[copy.deepcopy(rs15[i]) for i in range(len(rs15))]
    [rs18[i].changeName("rs18"+letter(i)) for i in range(len(rs18))]
    [rs18[i].irregularVary("accCeiling",.05) for i in range(len(rs18))]

    rs19=[copy.deepcopy(rs16[i]) for i in range(len(rs16))]
    [rs19[i].changeName("rs19"+letter(i)) for i in range(len(rs19))]
    [rs19[i].irregularVary("accCeiling",.05) for i in range(len(rs19))]

    rs21=copy.deepcopy(rs14[2])
    rs21.changeName("rs21")

    # try exp thing, didn't work
    rs22=copy.deepcopy(rs21)
    rs22.changeName("rs22")
    rs22.irregularVary('dbg',2**4)

    rs23=copy.deepcopy(rs21)
    rs23.changeName("rs23")
    rs23.irregularVary("tauHeat",.5)

    ### OK, we've fixed some bugs. Let's rerun some things.
    rs20 = [copy.deepcopy(rs05[i]) for i in range(len(rs05))]
    [rs20[i].changeName("rs20"+letter(i)) for i in range(len(rs20))]


    # do a study of different input angular momentum:
    scLengthMult=[1.0,4.0]
    rs30=[copy.deepcopy(rs02) for i in range(len(scLengthMult))]
    [rs30[i].changeName("rs30"+letter(i)) for i in range(len(rs30))]
    [rs30[i].multiply('accScaleLength',scLengthMult[i]) for i in range(len(rs30))]

    # another ang mom study, but explicitly set torque to zero everywhere
    rs31=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs31[i].changeName("rs31"+letter(i)) for i in range(len(rs31))]
    [rs31[i].irregularVary('dbg',2**12) for i in range(len(rs31))]

    # a high-mass-only case and a low-mass-only case
    Mh0s=[1.0e9,1.0e10,1.0e11,1.0e12,1.0e13]
    rs32=[experiment("rs32"+letter(i)) for i in range(len(Mh0s))]
    [rs32[i].irregularVary("R",40) for i in range(len(rs32))]
    [rs32[i].irregularVary("Mh0",Mh0s[i]) for i in range(len(rs32))]
    [rs32[i].irregularVary("accScaleLength",GetScaleLengths(1,Mh0=[Mh0s[i]],scatter=1.0e-10)) for i in range(len(rs32))]
    [rs32[i].irregularVary('diskScaleLength',GetScaleLengths(1,Mh0=[Mh0s[i]],scatter=1.0e-10)) for i in range(len(rs32))]
    [rs32[i].vary('whichAccretionHistory',1000,1100,101,0) for i in range(len(rs32))]
    [rs32[i].irregularVary('mu',(Mh0s[i]/1.0e12)**(-1./3.)) for i in range(len(rs32))]
    [rs32[i].irregularVary('vphiR',220.0*(Mh0s[i]/1.0e12)**(1./3.)) for i in range(len(rs32))]
    [rs32[i].irregularVary("NPassive",1) for i in range(len(rs32))]

    rs33=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs33[i].changeName("rs33"+letter(i)) for i in range(len(rs33))]
    [rs33[i].irregularVary('fscatter',.5) for i in range(len(rs33))]

    rs34=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs34[i].changeName("rs34"+letter(i)) for i in range(len(rs34))]
    [rs34[i].irregularVary('fscatter',.1) for i in range(len(rs34))]

    rs35=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs35[i].changeName("rs35"+letter(i)) for i in range(len(rs35))]
    [rs35[i].irregularVary('deltaOmega',.05) for i in range(len(rs35))]

    rs36=copy.deepcopy(rs32[1])
    rs36.changeName('rs36')
    rs36.irregularVary('b',2)

    rs37=copy.deepcopy(rs32[1])
    rs37.changeName('rs37')
    rs37.irregularVary("R",10)
    rs37.irregularVary('xmin',.005)

    rs38=copy.deepcopy(rs37)
    rs38.changeName('rs38')
    rs38.irregularVary('b',2)

    rs39=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs39[i].changeName("rs39"+letter(i)) for  i in range(len(rs39))]
    [rs39[i].irregularVary('epsff',.005) for i in range(len(rs39))]

    rs40=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs40[i].changeName("rs40"+letter(i)) for i in range(len(rs40))]
    [rs40[i].irregularVary('eta',.5) for i in range(len(rs40))]

    rs41=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs41[i].changeName("rs41"+letter(i)) for i in range(len(rs41))]
    [rs41[i].irregularVary('vphiR',list(150.0*(np.array(Mh0s101)/1.0e12)**(1./3.)),3) for i in range(len(rs41))]

    rs42=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs42[i].changeName("rs42"+letter(i)) for  i in range(len(rs42))]
    [rs42[i].irregularVary('mu',list(.1*(np.array(Mh0s101)/1.0e12)**(-1./3.)),3) for i in range(len(rs42))]

    rs43=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs43[i].changeName("rs43"+letter(i)) for i in range(len(rs43))]
    [rs43[i].irregularVary("kappaMetals",1.0e-4) for i in range(len(rs43))]

    rs44=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs44[i].changeName("rs44"+letter(i)) for i in range(len(rs44))]
    [rs44[i].irregularVary("kappaMetals",1.0e-2) for i in range(len(rs44))]   
    
    rs45=[copy.deepcopy(rs32[i]) for i in range(len(rs32))]
    [rs45[i].changeName("rs45"+letter(i)) for i in range(len(rs45))]
    [rs45[i].irregularVary("accScaleLength",GetScaleLengths(1,Mh0=[Mh0s[i]],scatter=1.0e-10)[0]*4.0) for i in range(len(rs45))]

    rs46=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs46[i].changeName("rs46"+letter(i)) for i in range(len(rs46))]
    [rs46[i].irregularVary("alphaMRI",0.1) for i in range(len(rs46))]
    
    # repeate rs32 and rs35 but simulate to smaller radii
    # High J
    rs47=[copy.deepcopy(rs45[i]) for i in range(len(rs45))]
    [rs47[i].changeName("rs47"+letter(i)) for i in range(len(rs47))]
    [rs47[i].irregularVary("R",50.0*(Mh0s[i]/1.0e12)**(1.0/3.0)) for i in range(len(rs47))]
    [rs47[i].irregularVary("xmin",.005) for i in range(len(rs47))]
    [rs47[i].irregularVary("dbg",2**1) for i in range(len(rs47))]
    
    # Low J
    rs48=[copy.deepcopy(rs32[i]) for i in range(len(rs32))]
    [rs48[i].changeName("rs48"+letter(i)) for i in range(len(rs48))]
    [rs48[i].irregularVary("R",50.0*(Mh0s[i]/1.0e12)**(1.0/3.0)) for i in range(len(rs48))]
    [rs48[i].irregularVary("xmin",.005) for i in range(len(rs48))]
    [rs48[i].irregularVary("dbg",2**1) for i in range(len(rs48))]
    
    fscatters=[0.1,0.5,1.0]
    # High J, vary fscatter
    rs49=[copy.deepcopy(rs47[3]) for i in range(len(fscatters))]
    [rs49[i].changeName("rs49"+letter(i)) for i in range(len(rs49))]
    [rs49[i].irregularVary('fscatter',fscatters[i]) for i in range(len(rs49))]

    # Low J, vary fscatter
    rs50=[copy.deepcopy(rs48[3]) for i in range(len(fscatters))]
    [rs50[i].changeName("rs50"+letter(i)) for i in range(len(rs50))]
    [rs50[i].irregularVary('fscatter',fscatters[i]) for i in range(len(rs50))]

    # High J, turnover radii
    rs51=[copy.deepcopy(rs47[i]) for i in range(len(rs47))]
    [rs51[i].changeName("rs51"+letter(i)) for i in range(len(rs51))]
    [rs51[i].irregularVary("b",2.0) for i in range(len(rs51))]

    # low J
    rs52=[copy.deepcopy(rs48[i]) for i in range(len(rs48))]
    [rs52[i].changeName("rs52"+letter(i)) for i in range(len(rs52))]
    [rs52[i].irregularVary("b",2.0) for i in range(len(rs52))]


    # High J, Qf
    rs53=[copy.deepcopy(rs47[i]) for i in range(len(rs47))]
    [rs53[i].changeName("rs53"+letter(i)) for i in range(len(rs53))]
    [rs53[i].irregularVary("fixedQ",1.0) for i in range(len(rs53))]

    # low J
    rs54=[copy.deepcopy(rs48[i]) for i in range(len(rs48))]
    [rs54[i].changeName("rs54"+letter(i)) for i in range(len(rs54))]
    [rs54[i].irregularVary("fixedQ",1.0) for i in range(len(rs54))]

    #High J Vary kappaMetals for MW-mass halos
    kappaMetals = [1.0e-4,3.0e-4,1.0e-3,3.0e-3]
    rs55=[copy.deepcopy(rs47[3]) for i in range(len(kappaMetals))]
    [rs55[i].changeName("rs55"+letter(i)) for i in range(len(rs55))]
    [rs55[i].irregularVary("kappaMetals",kappaMetals[i]) for i in range(len(rs55))]
    [rs55[i].irregularVary("dbg",2**1) for i in range(len(rs55))]

    # Low J
    rs56=[copy.deepcopy(rs48[3]) for i in range(len(kappaMetals))]
    [rs56[i].changeName("rs56"+letter(i)) for i in range(len(rs56))]
    [rs56[i].irregularVary("kappaMetals",kappaMetals[i]) for i in range(len(rs56))]
    [rs56[i].irregularVary("dbg",2**1) for i in range(len(rs56))]

    #High J Vary kappaMetals for MW-mass halos, with time-varying kappa
    rs57=[copy.deepcopy(rs47[3]) for i in range(len(kappaMetals))]
    [rs57[i].changeName("rs57"+letter(i)) for i in range(len(rs57))]
    [rs57[i].irregularVary("kappaMetals",kappaMetals[i]) for i in range(len(rs57))]
    [rs57[i].irregularVary("dbg",2**15+2**1) for i in range(len(rs57))]

    # Low J
    rs58=[copy.deepcopy(rs48[3]) for i in range(len(kappaMetals))]
    [rs58[i].changeName("rs58"+letter(i)) for i in range(len(rs58))]
    [rs58[i].irregularVary("kappaMetals",kappaMetals[i]) for i in range(len(rs58))]
    [rs58[i].irregularVary("dbg",2**15+2**1) for i in range(len(rs58))]

    # High J, Full mass range for tau=0
    rs59=[copy.deepcopy(rs47[i]) for i in range(len(rs47))]
    [rs59[i].changeName("rs59"+letter(i)) for i in range(len(rs59))]
    [rs59[i].irregularVary("dbg",2**12) for i in range(len(rs59))]

    #Low J
    rs60=[copy.deepcopy(rs48[i]) for i in range(len(rs48))]
    [rs60[i].changeName("rs60"+letter(i)) for i in range(len(rs60))]
    [rs60[i].irregularVary("dbg",2**12) for i in range(len(rs60))]

    # A little experiment to try large zstarts
    rs61=[copy.deepcopy(rs48[i]) for i in range(len(rs48))]
    [rs61[i].changeName("rs61"+letter(i)) for i in range(len(rs61))]
    [rs61[i].irregularVary("fg0",.999) for i in range(len(rs61))]
    [rs61[i].irregularVary("zrelax",9) for i in range(len(rs61))]
    [rs61[i].irregularVary("zstart",2) for i in range(len(rs61))]
    [rs61[i].irregularVary("invMassRatio",2) for i in range(len(rs61))]

    rs62=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs62[i].changeName("rs62"+letter(i)) for i in range(len(rs62))]
    [rs62[i].irregularVary("fg0",.999) for i in range(len(rs62))]
    [rs62[i].irregularVary("zrelax",9) for i in range(len(rs62))]
    [rs62[i].irregularVary("zstart",2) for i in range(len(rs62))]
    [rs62[i].irregularVary("invMassRatio",2) for i in range(len(rs62))]

    # try turning NPassive>1 back on
    rs63=[copy.deepcopy(rs30[i]) for i in range(len(rs30))]
    [rs63[i].changeName("rs63"+letter(i)) for i in range(len(rs63))]
    [rs63[i].irregularVary("NPassive",10) for i in range(len(rs63))]



    l045 = GetScaleLengths(1,Mh0=1.0e12,scatter=1.0e-10)[0]

    lengthFactors=[ 0.3, 0.7, 1.1 ]
    ru01=[experiment("ru01"+letter(i)) for i in range(3)]
    [ru01[i].irregularVary("R",40) for i in range(len(ru01))]
    [ru01[i].irregularVary('diskScaleLength',l045/3.0) for i in range(len(ru01))]
    [ru01[i].irregularVary('accScaleLength',l045*lengthFactors[i]) for i in range(len(ru01))]
    [ru01[i].irregularVary('mu',1.0) for i in range(len(ru01))]
    [ru01[i].irregularVary('vphiR',220.0) for i in range(len(ru01))]
    [ru01[i].irregularVary('NPassive',10) for i in range(len(ru01))]
    [ru01[i].irregularVary('invMassRatio',1.0) for i in range(len(ru01))]
    [ru01[i].irregularVary('dbg',2**1) for i in range(len(ru01))]
    [ru01[i].irregularVary('xmin',.01) for i in range(len(ru01))]

    ru02=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru02[i].changeName('ru02'+letter(i)) for i in range(len(ru02))]
    [ru02[i].vary('whichAccretionHistory',1000,1100,101,0) for i in range(len(ru02))]

    ru03=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru03[i].changeName("ru03"+letter(i)) for i in range(len(ru03))]
    [ru03[i].vary('whichAccretionHistory',-1100,-1000,101,0) for i in range(len(ru03))]
    [ru03[i].irregularVary('NChanges',10) for i in range(len(ru03))]
    [ru03[i].irregularVary('fscatter',.3) for i in range(len(ru03))]

    ru04=[copy.deepcopy(ru03[i]) for i in range(len(ru03))]
    [ru04[i].changeName("ru04"+letter(i)) for i in range(len(ru04))]
    [ru04[i].irregularVary('NChanges',5) for i in range(len(ru04))]

    ru05=[copy.deepcopy(ru03[i]) for i in range(len(ru03))]
    [ru05[i].changeName("ru05"+letter(i)) for i in range(len(ru05))]
    [ru05[i].irregularVary('NChanges',100) for i in range(len(ru05))]

    # Artificial alterations in the SF law:
    # c=5/3, c=15, fH2=1, fH2=.03, Zgrad, fH2->0, sig*=sig
    sfAlter = [2+2**6, 2+2**8, 2+2**9, 2+2**10, 2+2**13, 2+2**17, 2+2**16]
    ru06=[copy.deepcopy(ru01[1]) for i in range(len(sfAlter))]
    [ru06[i].changeName("ru06"+letter(i)) for i in range(len(sfAlter))]
    [ru06[i].irregularVary('dbg',sfAlter[i]) for i in range(len(ru06))] # clumping = 5/3

    # low-high extermes:
    ru07=[copy.deepcopy(ru01[1]) for i in range(2)]
    [ru07[i].changeName("ru07"+letter(i)) for i in range(len(ru07))]
    [ru07[i].irregularVary("fixedQ",1.0+2*i) for i in range(len(ru07))]

    ru08=[copy.deepcopy(ru07[i]) for i in range(len(ru07))]
    [ru08[i].changeName("ru08"+letter(i)) for i in range(len(ru08))]
    [ru08[i].irregularVary("dbg",2+2**4) for i in range(len(ru08))]

    ru09=[copy.deepcopy(ru01[1]) for i in range(2)]
    [ru09[i].changeName("ru09"+letter(i)) for i in range(len(ru09))]
    [ru09[i].irregularVary('vphiR',180.0+70*i) for i in range(len(ru09))]

    ru10=[copy.deepcopy(ru01[1]) for i in range(2)]
    [ru10[i].changeName("ru10"+letter(i)) for i in range(len(ru10))]
    [ru10[i].irregularVary("mu",0.1+1.9*i) for i in range(len(ru10))]

    ru11=[copy.deepcopy(ru01[1]) for i in range(2)]
    [ru11[i].changeName("ru11"+letter(i)) for i in range(len(ru11))]
    [ru11[i].irregularVary('b',1.0+4.0*i) for i in range(len(ru11))]

    ru12=[copy.deepcopy(ru01[1]) for i in range(2)]
    [ru12[i].changeName("ru12"+letter(i)) for i in range(len(ru12))]
    [ru12[i].irregularVary("eta",.3+4.7*i) for i in range(len(ru12))]

    ru13=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru13[i].changeName("ru13"+letter(i)) for i in range(len(ru13))]
    [ru13[i].irregularVary("alphaMRI",0.1) for i in range(len(ru13))]

    ru14=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru14[i].changeName("ru14"+letter(i)) for i in range(len(ru14))]
    [ru14[i].irregularVary("dbg",2+2**12) for i in range(len(ru14))]

    ru15=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru15[i].changeName("ru15"+letter(i)) for i in range(len(ru15))]
    [ru15[i].vary('kappaMetals',1.0e-5,1.0e-2,20,1) for i in range(len(ru15))]
    
    ru16=[copy.deepcopy(ru01[i]) for i in range(len(ru01))]
    [ru16[i].changeName("ru16"+letter(i)) for i in range(len(ru16))]
    [ru16[i].vary('fixedQ',1.0,3.0,20,0) for i in range(len(ru16))]



    


    # A series of experiments where we concentrate only on Mh0=10^12
    rt01 = experiment("rt01")
    rt01.irregularVary("R",40)
    # Scale length for lambda = 0.045:
    rt01.irregularVary('accScaleLength',l045/3.0)
    rt01.irregularVary('diskScaleLength',l045/3.0)
    rt01.irregularVary('mu',1.0)
    rt01.irregularVary('vphiR',220.0)
    rt01.irregularVary('NPassive',10)
    rt01.irregularVary('invMassRatio',1.0)
    rt01.irregularVary('dbg',2**1)

    # effect of realistic accretion history variations
    rt02=copy.deepcopy(rt01)
    rt02.changeName("rt02")
    rt02.vary('whichAccretionHistory',1000,1100,101,0)

    # different assumptions for accretion history variation
    rt03=copy.deepcopy(rt01)
    rt03.changeName("rt03")
    rt03.vary('whichAccretionHistory',-1100,-1000,101,0)
    rt03.irregularVary('NChanges',10) # by default, const. in redshift 
    rt03.irregularVary("fscatter",.3) # half a dex of scatter

    rt04=copy.deepcopy(rt01)
    rt04.changeName("rt04")
    rt04.irregularVary("whichAccretionProfile",0) # exponential. Vary scale length
    rt04.irregularVary('accScaleLength',[l045*(i+1)/10.0 for i in range(10)])

    rt05=copy.deepcopy(rt01)
    rt05.changeName("rt05")
    rt05.irregularVary('whichAccretionProfile',1) # const accretion profile
    rt05.irregularVary('accScaleLength',[l045*(i+1)*3.0/10.0 for i in range(10)]) # same ang mom distr as rt04

    rt06=copy.deepcopy(rt01)
    rt06.changeName("rt06")
    rt06.irregularVary('whichAccretionProfile',2) # normal distr.
    rt06.irregularVary('accScaleLength',[l045*(i+1)/10.0 for i in range(10)]) # hopefully roughly the same ang mom as rt04? Still need to do that integral...

    rt07=copy.deepcopy(rt02)
    rt07.changeName("rt07")
    rt07.irregularVary('deltaOmega',.4)

    # vary clumping factor: 5/3 or 15
    rt08=copy.deepcopy(rt01)
    rt08.changeName('rt08')
    rt08.irregularVary('dbg',[2+2**6,2+2**8])

    # artificially set fH2=1 or .03
    rt09=copy.deepcopy(rt01)
    rt09.changeName('rt09')
    rt09.irregularVary('dbg',[2+2**9,2+ 2**10])

    rt10=copy.deepcopy(rt01)
    rt10.changeName("rt10")
    rt10.vary('kappaMetals',1.0e-5,1.0e-2,20,1)

    rt11=copy.deepcopy(rt01)
    rt11.changeName("rt11")
    rt11.vary('fixedQ',1.0,3.0,20,0)

    rt12=copy.deepcopy(rt01)
    rt12.changeName("rt12")
    rt12.vary('Qlim',2.0,4.0,20,0)

    rt13=copy.deepcopy(rt11)
    rt13.changeName("rt13")
    rt13.irregularVary("dbg",2+2**4) # exp dQ/dt

    rt14=copy.deepcopy(rt12)
    rt14.changeName("rt14")
    rt14.irregularVary("dbg",2+2**4) # exp dQ/dt

    rt15=copy.deepcopy(rt01)
    rt15.changeName("rt15")
    rt15.vary('b',0.0,5.0,30,0)

    rt16=copy.deepcopy(rt01)
    rt16.changeName("rt16")
    rt16.vary('alphaMRI',1.0e-4,0.1,20,1)

    rt17=copy.deepcopy(rt01)
    rt17.changeName("rt17")
    rt17.vary("eta",.3,4.5,20,1)

    rt18=copy.deepcopy(rt01)
    rt18.changeName("rt18")
    rt18.irregularVary('dbg',2+2**12)




    # Guess for a reasonable model.
    l045 = GetScaleLengths(1,Mh0=1.0e12,scatter=1.0e-10)[0]
    rv01=experiment("rv01")
    rv01.irregularVary("R",40)
    rv01.irregularVary('diskScaleLength',l045*.7)
    rv01.irregularVary('accScaleLength',l045*.7)
    rv01.irregularVary('mu',.5)
    rv01.irregularVary('vphiR',220.0)
    rv01.irregularVary('NPassive',10)
    rv01.irregularVary('invMassRatio',1.0)
    rv01.irregularVary('dbg',2)
    rv01.irregularVary('xmin',.004)
    rv01.irregularVary('alphaMRI',.01)

    # Vary the exponential scale length of the IC and accretion together.
    rv02=[copy.deepcopy(rv01) for i in range(2)]
    [rv02[i].changeName("rv02"+letter(i)) for i in range(len(rv02))]
    rv02[0].vary('diskScaleLength',l045*.10,l045*.67,10,0,3)
    rv02[1].vary('diskScaleLength',l045*.73,l045*2.1,20,0,3)
    rv02[0].vary('accScaleLength',l045*.10,l045*.67,10,0,3)
    rv02[1].vary('accScaleLength',l045*.73,l045*2.1,20,0,3)

    # Vary the metal diffusion constant 
    rv03=[copy.deepcopy(rv01) for i in range(2)]
    [rv03[i].changeName('rv03'+letter(i)) for i in range(len(rv03))]
    rv03[0].vary('kappaMetals',1.0e-4,9.5e-4,10,1)
    rv03[1].vary('kappaMetals',1.1e-3,3.0e-3,5,1)

    # Vary the metal diffusion constant normalization, but this time let it vary in proportion to the sigma H
    rv04=[copy.deepcopy(rv03[i]) for i in range(len(rv03))]
    [rv04[i].changeName('rv04'+letter(i)) for i in range(len(rv04))]
    [rv04[i].irregularVary('dbg',2+2**15) for i in range(len(rv04))]

    # Vary the Q below which the disk will be unstable
    rv05=NewSetOfExperiments(rv01,"rv05",N=2)
    rv05[0].vary('fixedQ',1.3,1.9,7,0)
    rv05[1].vary('fixedQ',2.1,3.0,10,0)

    # Vary the Q below which the disk will be unstable, but relax the timescale on which turbulence will stabilize the disk
    rv06=NewSetOfExperiments(rv05,"rv06")
    [rv06[i].irregularVary('dbg',2+2**4) for i in range(len(rv06))]

    # Vary MRI torques
    rv07=NewSetOfExperiments(rv01,"rv07",N=2)
    rv07[0].vary('alphaMRI',0.0,.009, 5,0)
    rv07[1].vary('alphaMRI',.02, 1.0,15,1)

    # Vary mass loading factor
    rv08=NewSetOfExperiments(rv01,"rv08",N=2)
    rv08[0].vary('mu',.1,.4,4,0)
    rv08[1].vary('mu',.6,2.0,15,0)

    # Vary dissipation rate
    rv09=NewSetOfExperiments(rv01,"rv09",N=2)
    rv09[0].vary('eta',.5,1.5,10,1)
    rv09[1].vary('eta',1.5,4.5,10,1)

    # Vary turnover radius of rot. curve with an inner power law of 0.5 and 1
    rv10=NewSetOfExperiments(rv01,"rv10",N=2)
    rv10[0].vary('b',1,10,10,0)
    rv10[1].vary('b',1,10,10,0)
    rv10[1].irregularVary('innerPowerLaw',1.0)

    # Vary initial gas fraction
    rv11=NewSetOfExperiments(rv01,"rv11",N=2)
    rv11[0].vary('fg0',0.1,.45,8,0)
    rv11[1].vary('fg0',.55,.99,10,0)

    # Vary circular velocity
    rv12=NewSetOfExperiments(rv01,"rv12",N=2)
    rv12[0].vary('vphiR',160,215,12,0)
    rv12[1].vary('vphiR',225,250, 6,0)

    # Vary accretion scale, but this time use a narrow gaussian profile
    rv13=NewSetOfExperiments(rv02,"rv13")
    [rv13[i].irregularVary('whichAccretionProfile',2) for i in range(len(rv13))]

    # Vary star formation efficiency
    rv14=NewSetOfExperiments(rv01,"rv14",N=2)
    rv14[0].vary('epsff',.0031,.01,10,1)
    rv14[1].vary('epsff',.01,.031,10,1)

    # Vary metal loading factor
    rv15=NewSetOfExperiments(rv01,"rv15",N=2)
    rv15[0].vary('zetaREC',.5,1,10,1)
    rv15[1].vary('zetaREC',1,2,10,1)

    # Vary efficiency of accretion at high redshift
    rv16=NewSetOfExperiments(rv01,"rv16",N=2)
    rv16[0].vary('accAlphaZ',.1,.38,4,0)
    rv16[1].vary('accAlphaZ',.4,.8,5,0)

    # Vary accretion scale, this time with a wide gaussian profile
    rv17=NewSetOfExperiments(rv13,"rv17")
    [rv17[i].irregularVary('widthAccretionProfile',0.4) for i in range(len(rv17))]

    # Vary recycling fraction
    rv18=NewSetOfExperiments(rv01,"rv18",N=2)
    rv18[0].vary("RfREC",.22,.45,10,0)
    rv18[1].vary("RfREC",.47,.7,10,0)

    # Vary recycling fraction, with large RfREC leading to non-instantaneous recycling
    rv19=NewSetOfExperiments(rv18,"rv19")
    [rv19[i].irregularVary('dbg',2+2**6) for i in range(len(rv19))]

    # Vary delta omega
    rv20=NewSetOfExperiments(rv01,"rv20",N=3)
    [rv20[i].vary('whichAccretionHistory',1000,1400,401,0) for i in range(len(rv20))]
    rv20[0].irregularVary('deltaOmega',.3)
    rv20[1].irregularVary('deltaOmega',.5)
    rv20[2].irregularVary('deltaOmega',.8)

    # Lognormal acc history variation, with different coherence redshift intervals
    rv21=NewSetOfExperiments(rv01,"rv21",N=2)
    [rv21[i].vary('whichAccretionHistory',-1010,-1000,11,0) for i in range(len(rv21))]
    [rv21[i].irregularVary("fscatter",0.3) for i in range(len(rv21))]
    rv21[0].irregularVary("NChanges",5)
    rv21[1].irregularVary("NChanges",30)

    # Lognormal acc history variation, with different scatter amplitude
    rv22=NewSetOfExperiments(rv01,"rv22",N=2)
    [rv22[i].vary("whichAccretionHistory",-1010,-1000,11,0) for i in range(len(rv22))]
    [rv22[i].irregularVary("NChanges",10) for i in range(len(rv22))]
    rv22[0].irregularVary("fscatter",0.1)
    rv22[1].irregularVary("fscatter",0.5)

    # Gaussian profile with non-inst. recycling. Can we fill those holes?
    rv23=NewSetOfExperiments(rv17,"rv23")
    [rv23[i].irregularVary("RfREC",0.8) for i in range(len(rv23))]
    [rv23[i].irregularVary('dbg',2+2**6) for i in range(len(rv23))]

    # Gaussian again with delayed recycling, but this time have an evolving radius.
    rv24=NewSetOfExperiments(rv23,"rv24")
    [rv24[i].irregularVary('alphaAccretionProfile',0.5) for i in range(len(rv24))]

    # Build on the previous model, but significantly reduce the z=0 acc rate
    # tuned so that efficiency is the same as before at z=2, but evolution goes as (1+z)^2
    rv25=NewSetOfExperiments(rv24, "rv25")
    [rv25[i].irregularVary("accAlphaZ",2) for i in range(len(rv25))]
    [rv25[i].irregularVary("accNorm",.0522216) for i in range(len(rv25))]

    successTables=[]

    for inputString in modelList: # aModelName will therefore be a string, obtained from the command-line args
        # Get a list of all defined models (allModels.keys())
        # for each such key (aModel) check whether this inputString is contained in its name
        matches = [aModel for aModel in sorted(allModels.keys()) if inputString in aModel]
        if(len(matches) != 0): 
	    for model in matches: #if(model in allModels): 
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
