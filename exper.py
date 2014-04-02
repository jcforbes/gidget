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
        self.p=[name,200,1.5,.01,4.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,1.0,2.0,5000.0,int(1e9),1.0e-3,1.0,0,.5,2.0,2.0,0,0.0,1.5,1,2.0,1.0,1.0e12,5.0,3.0,0,2.0,-1.0,2.5,0.0,0.54,0.1,200,.30959,0.38,-0.25,1.0,1.0,0.3,1.0,0,0.0,.1,.03,2.0,.002,.054]
        self.p_orig=self.p[:] # store a copy of p, possibly necessary later on.
        self.pl=[self.p[:]] # define a 1-element list containing a copy of p.
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat','analyticQ','cosmologyOn','xmin','NActive','NPassive','vphiR','R','gasTemp','Qlim','fg0','phi0','zstart','tmax','stepmax','TOL','mu','b','innerPowerLaw','softening','diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive','fixedQ','kappaMetals','Mh0','minSigSt','NChanges','dbg','accScaleLength','zquench','zrelax','xiREC','RfREC','deltaOmega','Noutputs','accNorm','accAlphaZ','accAlphaMh','accCeiling','fscatter','invMassRatio','fcool','whichAccretionProfile','alphaAccretionProfile','widthAccretionProfile','fH2Min','tDepH2SC','ZIGM','yREC']
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
            print "************"
            print "This directory already contains output. CANCELLING this run!"
            print
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


    # Guess for a reasonable model.
    l045 = GetScaleLengths(1,Mh0=1.0e12,scatter=1.0e-10)[0]
    print "l045 = ",l045," kpc"
    rg01=experiment("rg01")
    rg01.irregularVary("R",40)
    rg01.irregularVary('accScaleLength',l045*.7)
    rg01.irregularVary('mu',.5)
    rg01.irregularVary('vphiR',220.0)
    rg01.irregularVary('NPassive',20)
    rg01.irregularVary('invMassRatio',1.0)
    rg01.irregularVary('dbg',2)
    rg01.irregularVary('xmin',.002)
    rg01.irregularVary('alphaMRI',.01)
    rg01.irregularVary('fcool',0.6)
    rg01.irregularVary('innerPowerLaw',0.5)
    rg01.irregularVary('b',3.0)
    rg01.irregularVary('nx',200)
    rg01.irregularVary('kappaMetals',1.0)
    rg01.irregularVary('xiREC',0.0)
    rg01.irregularVary('alphaAccretionProfile',1./3.)
    rg01.irregularVary('deltaOmega',.1)

    rg01t=NewSetOfExperiments(rg01,'rg01t',N=1)
    rg01t[0].irregularVary('dbg',2+2**4)

    rg01t2=NewSetOfExperiments(rg01t,'rg01t2',N=1)
    rg01t2[0].irregularVary('TOL',1.0e-3)

    rg01t3=NewSetOfExperiments(rg01t2,'rg01t3',N=1)
    rg01t3[0].vary('mu',0.1,10.0,3,1)

    rg01t4=NewSetOfExperiments(rg01t3,'rg01t4',N=1) # turn off disk contrib to rot curve, fix stellar scale height error
    rg01t4[0].irregularVary('fH2Min',1.0e-10)
    rg01t4[0].irregularVary('TOL',1.0e-4)
    rg01t4[0].irregularVary('NPassive',3)

    rg01t5=NewSetOfExperiments(rg01t4,'rg01t5',N=1) # now with exponential forcing

    rg01t6=NewSetOfExperiments(rg01t5,'rg01t6',N=1) # now using approx. prescription for calculating fH2 and SFR.

    rg01t7=NewSetOfExperiments(rg01t6,'rg01t7',N=1) # now using no updates to rot. curve. In theory this should be very similar
    # to the results we had before, in terms of speed, etc.

    rg01t8=NewSetOfExperiments(rg01t7,'rg01t8',N=1) # add in some smoothing out of the torques because the sims are very slow
    rg01t8[0].irregularVary('dbg',2+2**4+2**18)

    rg01t9=NewSetOfExperiments(rg01t8,'rg01t9',N=1) # allow the stars a longer time to reach Qlim
    rg01t9[0].irregularVary('tauHeat',30)

    rg01t10=NewSetOfExperiments(rg01t7,'rg01t10',N=1) # allow the stars a long time to reach Qlim, AND do not include smoothing
    rg01t10[0].irregularVary('tauHeat',50)

    rg01t11=NewSetOfExperiments(rg01t7,'rg01t11',N=1) # try the opposite -- immediate equlibration
    rg01t11[0].irregularVary('tauHeat',.0001)

    rg01t12=NewSetOfExperiments(rg01t7,'rg01t12',N=1)
    rg01t12[0].irregularVary('tauHeat',100000)

    rg01t13=NewSetOfExperiments(rg01t7,'rg01t13',N=1) # try taking out the exponential approach to Qfixed, return to exact forcing eq.
    rg01t13[0].irregularVary('dbg',2)

    rg01t14=NewSetOfExperiments(rg01t13,'rg01t14',N=1) # exact Q>=Qfixed. Only update rot. curve @ first time step, then leave it there the whole sim.
    rg01t14[0].irregularVary('dbg',2)


    rg01t15=NewSetOfExperiments(rg01t14,'rg01t15',N=1) # exact Q>=Qfixed. Try again update rotation curve frequently.
    rg01t15[0].irregularVary('dbg',2)

    rg01t16=NewSetOfExperiments(rg01t15,'rg01t16',N=1) # exact Q>=Qfixed. Bug fixed in rot. curve
    rg01t16[0].irregularVary('dbg',2)


    # Vary the scale length of the accretion.
    rg02=NewSetOfExperiments(rg01,"rg02",N=2)
    rg02[0].vary('accScaleLength',l045*.10,l045*.67,5,0,3)
    rg02[0].vary('R',l045*.10*7,l045*.67*7,5,0,3)
    rg02[1].vary('accScaleLength',l045*.73,l045*2.1,10,0,3)
    rg02[1].vary('R',l045*.73*7,l045*2.1*7,10,0,3)


    rh02=NewSetOfExperiments(rg01,'rh02')[0]
    rh02.vary('accScaleLength',1,10,12,1,3)
    rh02.vary('R',5,50,12,1,3)
    rh02.irregularVary('mu',1.0)
    rh02.irregularVary('alphaMRI',0.02)

    rh03=NewSetOfExperiments(rh02,'rh03')[0]
    rh03.vary('xmin',.02/(1.0*l045),.02/(9*l045),12,1,3)
    #rh03.vary('nx',200,200*(40/3),12,1,3)


    rh04=NewSetOfExperiments(rg01,'rh04')[0]
    rh04.vary('alphaMRI',.01,1,10,1)

    rh05=NewSetOfExperiments(rg01,'rh05')[0]
    rh05.vary('fixedQ',1,3,10,0)

    # Vary mass and accretion history
    rh06=NewSetOfExperiments(rg01,'rh06')[0]
    rh06.vary('vphiR',220.0*(.01)**(1.0/3.0), 220.0*(10.0)**(1.0/3.0), 300, 1, 17)
    rh06.vary('accScaleLength',7.0*(.01)**(1.0/3.0), 7.0*(10.0)**(1./3.), 300, 1, 17)
    rh06.vary('mu',1.0*(.01)**(-1.0/3.0),1.0*(10.0)**(-1./3.),300,1,17)
    rh06.vary('R',50.0*(.01)**(1.0/3.0), 60.0*(10.0)**(1./3.), 300, 1, 17)
    rh06Mh0 = rh06.vary('Mh0',1.0e10,1.0e13,300,1,17)
    rh06.vary('whichAccretionHistory',3000,3299,300,0,17)
    rh06.irregularVary('alphaMRI',0.02)
    rh06.irregularVary('dbg',2+2**3) 
    rh06.irregularVary('deltaOmega',0.1)
    rh06.irregularVary('b',0)

    # Vary mass, accretion history, and scale length
    rh07=NewSetOfExperiments(rh06,'rh07')[0]
    rh07.irregularVary('accScaleLength',GetScaleLengths(300,Mh0=rh06Mh0,scatter=0.3,multiple=7.0/l045), 17)
    rh07.irregularVary('R',GetScaleLengths(300,Mh0=rh06Mh0,scatter=0.3,multiple=60.0/l045), 17)
    rh07.irregularVary('b',0.0) # formerly b=3 kpc! That setup caused 68/300 models to fail (!)

    rh07x=NewSetOfExperiments(rh07,'rh07x')[0]
    rh07x.irregularVary('fg0',.999999999999)

    rh07y=NewSetOfExperiments(rh07x,'rh07y')[0]
    rh07y.irregularVary('ZIGM',.00002)

    rh07z=NewSetOfExperiments(rh07y,'rh07z')[0]
    rh07z.irregularVary('fH2Min',.001)

    rh07q=NewSetOfExperiments(rh07y,'rh07q')[0]
    rh07q.irregularVary('fH2Min',.0001)


    # Vary mass but nothing else
    rh08=NewSetOfExperiments(rh06,'rh08')[0]
    rh08.irregularVary('whichAccretionHistory',0)

    rh09=NewSetOfExperiments(rg01,'rh09',N=3)
    # no star formation! - rg36d
    rh09[1].irregularVary('epsff',0)
    rh09[1].irregularVary('tDepH2SC',1000000000.0)
    rh09[1].irregularVary('fH2Min',0.03)
    rh09[2].irregularVary('dbg',2+2**12) # no GI - rg36e

    # Vary mass and scale length, but not accretion history
    rh10=NewSetOfExperiments(rh07,'rh10')[0]
    rh10.irregularVary('whichAccretionHistory',0)

    # Vary mu
    rh11=NewSetOfExperiments(rg01,"rh11")[0]
    rh11.vary('mu',.5,10,12,1)

    rh12=NewSetOfExperiments(rg01,'rh12')[0]
    rh12.vary('whichAccretionHistory',1200,1499,300,0)
    rh12.irregularVary('alphaMRI',.05)

    rh13=NewSetOfExperiments(rg01,'rh13')[0]
    rh13.vary('zquench',.1,1.9,5,0)
    rh13.irregularVary('alphaMRI',.05)

    rh14=NewSetOfExperiments(rh12,'rh14')[0]
    rh14.irregularVary('b',0)

    # Make some attempts to bring sims into aggreement with R Genzel.
    rh16=NewSetOfExperiments(rg01,'rh16')[0]
    rh16.vary('fixedQ',1.3,2.5,3,0)
    rh16.vary('b',0,3,3,0)
    rh16.vary('vphiR',150,300,4,0)
    rh16.vary('eta',.1,10.0,3,1)
    rh16.vary('innerPowerLaw',-.5,.5,2,0,14)
    rh16.irregularVary('softening',[-2,2],14)
    rh16.vary('alphaAccretionProfile',1./3.,3,3,1)


    # Vary the metal diffusion constant 
    rg03=NewSetOfExperiments(rg01,"rg03",N=2)
    rg03[0].vary('kappaMetals',1.0e-1,1.0,6,1)
    rg03[1].vary('kappaMetals',1.0,1.0e1,6,1)


    # rg04 and rg04y study various prescriptions for the accretion history
    rg04=NewSetOfExperiments(rg01,"rg04",N=4)
    [rg04[i].irregularVary('deltaOmega',0.5) for i in range(len(rg04))]
    # rg04a Dekel13 prescription for WMAP5
    rg04[0].irregularVary('dbg',2+2**8) 
    # rg04b Average of NMD w/ standard params
    rg04[1].irregularVary('dbg',2+2**15) 
    # rg04c Average of NMD w/ no merger cut
    rg04[2].irregularVary('dbg',2+2**15)
    rg04[2].irregularVary('invMassRatio',100.0)
    # rg04d Average of NMD w/ strict merger cut
    rg04[3].irregularVary('dbg',2+2**15)
    rg04[3].irregularVary('invMassRatio',0.2)

    rg04y=NewSetOfExperiments(rg01,"rg04y",N=4)
    [rg04y[i].irregularVary('deltaOmega',0.1) for i in range(len(rg04y))]
    # rg04ya Dekel13 prescription for WMAP5
    rg04y[0].irregularVary('dbg',2+2**8+2**3) 
    # rg04yb Average of NMD w/ standard params
    rg04y[1].irregularVary('dbg',2+2**15+2**3) 
    # rg04yc Average of NMD w/ no merger cut
    rg04y[2].irregularVary('dbg',2+2**15+2**3)
    rg04y[2].irregularVary('invMassRatio',100.0)
    # rg04yd Average of NMD w/ strict merger cut
    rg04y[3].irregularVary('dbg',2+2**15+2**3)
    rg04y[3].irregularVary('invMassRatio',0.2)



    # Vary the Q below which the disk will be unstable
    rg05=NewSetOfExperiments(rg01,"rg05",N=2)
    rg05[0].vary('fixedQ',1.3,1.9,4,0)
    rg05[1].vary('fixedQ',2.1,3.0,6,0)

    # Vary the Q below which the disk will be unstable, but relax the timescale on which turbulence will stabilize the disk
    rg06=NewSetOfExperiments(rg05,"rg06")
    [rg06[i].irregularVary('dbg',2+2**4+2**15) for i in range(len(rg06))]

    # Vary MRI torques
    rg07=NewSetOfExperiments(rg01,"rg07",N=2)
    rg07[0].vary('alphaMRI',0.0,.009, 6,0)
    rg07[1].vary('alphaMRI',.011, .5,12,0)

    # Vary mass loading factor
    rg08=NewSetOfExperiments(rg01,"rg08",N=2)
    rg08[0].vary('mu',.1,.4,4,0)
    rg08[1].vary('mu',.6,1.5,10,0)

    # Vary dissipation rate
    rg09=NewSetOfExperiments(rg01,"rg09",N=2)
    rg09[0].vary('eta',.5,1.5,5,1)
    rg09[1].vary('eta',1.5,4.5,5,1)

    # Vary turnover radius of rot. curve with an inner power law of 0.5 and 1
    rg10=NewSetOfExperiments(rg01,"rg10",N=2)
    rg10[0].vary('b',0,2.5,3,0)
    rg10[1].vary('b',3.5,10,7,0)


    # Vary initial gas fraction
    rg11=NewSetOfExperiments(rg01,"rg11",N=2)
    rg11[0].vary('fg0',0.1,.45,8,0)
    rg11[1].vary('fg0',.55,.7,4,0)

    # Vary circular velocity
    rg12=NewSetOfExperiments(rg01,"rg12",N=2)
    rg12[0].vary('vphiR',160,215, 9,0)
    rg12[1].vary('vphiR',225,250, 5,0)

    # Vary accretion scale, but this time use a narrow gaussian profile
    # This one has some problems
    rg13=NewSetOfExperiments(rg02,"rg13")
    [rg13[i].irregularVary('whichAccretionProfile',2) for i in range(len(rg13))]

    # Vary star formation efficiency
    rg14=NewSetOfExperiments(rg01,"rg14",N=2)
    rg14[0].vary('epsff',.0031,.0095,6,1)
    rg14[1].vary('epsff',.0105,.031,6,1)

    # Vary metal ejecta mixing factor
    rg15=NewSetOfExperiments(rg01,"rg15",N=2)
    rg15[0].vary('xiREC',-0.3,-0.1,3,0)
    rg15[1].vary('xiREC',0.1,1.0,5,0)

    rh15=NewSetOfExperiments(rg01,"rh15")[0]
    rh15.vary('xiREC',0,.6,7,0)

    # Vary Qlim
    rg16=NewSetOfExperiments(rg01,"rg16",N=2)
    rg16[0].vary("Qlim",1.8,2.4,7,0)
    rg16[1].vary("Qlim",2.6,3.0,5,0)

    # Vary accretion scale, this time with a wide gaussian profile
    # similar problems as rg13.
    rg17=NewSetOfExperiments(rg13,"rg17")
    [rg17[i].irregularVary('widthAccretionProfile',0.75) for i in range(len(rg17))]

    # Vary recycling fraction
    rg18=NewSetOfExperiments(rg01,"rg18",N=2)
    rg18[0].vary("RfREC",.4,.5,4,0)
    rg18[1].vary("RfREC",.6,.7,4,0)

    # Vary recycling fraction, with large RfREC leading to non-instantaneous recycling
    # Note that non-inst. recycling is not well-tested
    rg19=NewSetOfExperiments(rg18,"rg19")
    [rg19[i].irregularVary('dbg',2+2**6) for i in range(len(rg19))]

    # Vary delta omega - this is not guaranteed to be physically meaningful.
    rg20=NewSetOfExperiments(rg01,"rg20",N=2)
    [rg20[i].vary('whichAccretionHistory',1001,1400,400,0,3) for i in range(len(rg20))]
    rg20[0].irregularVary('deltaOmega',.2)
    rg20[1].irregularVary('deltaOmega',.5)
#    rg20[2].irregularVary('deltaOmega',.8)

    # Lognormal acc history variation, with different coherence redshift intervals (set by NChanges -
    # every new draw from the lognormal distribution is equally distributed in redshift.
    rg21=NewSetOfExperiments(rg01,"rg21",N=2)
    [rg21[i].vary('whichAccretionHistory',-1200,-1000,201,0) for i in range(len(rg21))]
    [rg21[i].irregularVary("fscatter",0.2) for i in range(len(rg21))]
    rg21[0].irregularVary("NChanges",5)
    rg21[1].irregularVary("NChanges",30)

    # Lognormal acc history variation, with different scatter amplitude
    rg22=NewSetOfExperiments(rg01,"rg22",N=2)
    [rg22[i].vary("whichAccretionHistory",-1050,-1000,51,0) for i in range(len(rg22))]
    [rg22[i].irregularVary("NChanges",10) for i in range(len(rg22))]
    rg22[0].irregularVary("fscatter",0.1)
    rg22[1].irregularVary("fscatter",0.5)

    # Gaussian profile with non-inst. recycling. 
    # similar problems as w/ rg13 rg17
    rg23=NewSetOfExperiments(rg17,"rg23")
    [rg23[i].irregularVary("RfREC",0.8) for i in range(len(rg23))]
    [rg23[i].irregularVary('dbg',2+2**6) for i in range(len(rg23))]

    # Vary the sharpness of the turnover in the rotation curve
    rg24=NewSetOfExperiments(rg01,"rg24",N=2)
    rg24[0].vary('softening',1.0,1.9,5,0)
    rg24[1].vary('softening',2.1,3.0,5,0)

    # Vary the fraction of baryons which have `cooled' into a disk at z=zrelax
    rg25=NewSetOfExperiments(rg01,"rg25",N=2)
    rg25[0].vary('fcool',.20,.55,8,0)
    rg25[1].vary('fcool',.65,1.0,8,0)

    # Vary the inner power law.
    rg26=NewSetOfExperiments(rg01,"rg26",N=3)
    rg26[0].vary('innerPowerLaw',0,.45,5,0)
    rg26[1].vary('innerPowerLaw',.55,.95, 5,0)
    rg26[2].vary('innerPowerLaw',-.5,-.05, 5, 0)
    rg26[2].irregularVary('softening',-2)

    # Sanity check. Very strong efficiency evolution
    rg27=NewSetOfExperiments(rg01,"rg27",N=2)
    eps0Low = [.30959 * 1.0e-2 * 10.0**(2.0*i/10.0) for i in range(10)]
    eps0High= [.30959 + .05 *i for i in range(3)]
    eps2 = .30959 * 3.0**0.38
    rg27[0].irregularVary('accNorm',eps0Low,3)
    rg27[0].irregularVary('accAlphaZ',[math.log(eps2/eps0Low[i])/math.log(3.0) for i in range(len(eps0Low))],3)
    rg27[1].irregularVary('accNorm',eps0High,3)
    rg27[1].irregularVary('accAlphaZ',[math.log(eps2/eps0High[i])/math.log(3.0) for i in range(len(eps0High))],3)

    # Very strong efficiency evolution with wind recycling.
    rg28=NewSetOfExperiments(rg27,"rg28")
    [rg28[i].irregularVary('dbg',2+2**6) for i in range(len(rg28))]
    [rg28[i].irregularVary('RfREC',.8) for i in range(len(rg28))]

    # Vary initial ratio of stellar : gas velocity dispersion
    rg29=NewSetOfExperiments(rg01,"rg29",N=2)
    rg29[0].vary('phi0',.5,.95,6,1)
    rg29[1].vary('phi0',1.05,2.0,6,1)

    # Vary the floor of the H2 fraction.
    rg30=NewSetOfExperiments(rg01,"rg30",N=2)
    rg30[0].vary('fH2Min',.003,.029,6,1)
    rg30[1].vary('fH2Min',.031,.3,6,1)

    # Vary the gas temperature.
    rg31=NewSetOfExperiments(rg01,"rg31",N=2)
    rg31[0].vary('gasTemp',100,6900,8,1)
    rg31[1].vary('gasTemp',7100,30000,4,1)

    # Exp accretion with a stricter merger cut
    rg32=NewSetOfExperiments(rg20,"rg32")
    [rg32[i].irregularVary('invMassRatio',0.3) for i in range(len(rg32))]

    # Vary the single-cloud star formation timescale
    rg33=NewSetOfExperiments(rg01,"rg33",N=2)
    rg33[0].vary("tDepH2SC",1.0,1.9,5,1)
    rg33[1].vary('tDepH2SC',2.1,4.0,5,1)

    # A set of experimental parameters, including the `No GI' and `No SF' models from the most recent paper.
    rg36=NewSetOfExperiments(rg01,"rg36",N=9)
    rg36[0].irregularVary('dbg',2+2**4) # exp Delta Q
    rg36[1].irregularVary('dbg',2+2**17) # upstream
    rg36[2].irregularVary('dbg',2+2**18) # overshoot
    # no star formation! - rg36d
    rg36[3].irregularVary('epsff',0)
    rg36[3].irregularVary('tDepH2SC',1000000000.0)
#    rg36[3].irregularVary('ZIGM',10**-10)
    rg36[3].irregularVary('fH2Min',0.03)
    rg36[4].irregularVary('dbg',2+2**12) # no GI - rg36e
    rg36[5].irregularVary('dbg',2+1) # tau=0 when F<0
    rg36[5].irregularVary('kappaMetals',1.0e-7)
    rg36[5].irregularVary('xiREC',.9)
    rg36[6].irregularVary('innerPowerLaw',0) # g
    rg36[7].irregularVary('alphaAccretionProfile',1) #h
    rg36[8].irregularVary('dbg',2+2**12+2**13), #i


    # Same as rg36 but with a larger computational domain
    rg36x=NewSetOfExperiments(rg36,"rg36x")
    [rg36x[i].irregularVary("R",60) for i in range(len(rg36x))]

    
    # Vary ZIGM - this is the initial metallicity and the metallicity of incoming material
    rg38=NewSetOfExperiments(rg01,"rg38",N=2)
    rg38[0].vary("ZIGM",.0002,.002,5,1) # 1/100 - 1/10 solar
    rg38[1].vary("ZIGM",.002,.02,5,1) # 1/10 - 1 solar

    # Vary halo mass, but only as it affects the accretion history (vcirc, racc(z=0), etc. are unchanged)
    rg39=NewSetOfExperiments(rg01,"rg39",N=2)
    Mhlo = [1.0e10 * 10**(i/10.0) for i in range(20)]
    Mhhi = [1.0e12 * 10**((i+1.0)/10.0) for i in range(10)]
    rg39[0].irregularVary('Mh0',Mhlo,5)
    rg39[1].irregularVary('Mh0',Mhhi,5)


    rh40 = NewSetOfExperiments(rg01,"rh40")[0]
    Mh = [1.0e10 * 10**(i/10.0) for i in range(30)]
    rh40.irregularVary('Mh0',Mh,5)
    rh40.irregularVary('b',0)
    rh40.irregularVary("R",GetScaleLengths(30,Mh0=Mh,scatter=1.0e-10,multiple=4.1),5)
    rh40.irregularVary('vphiR',[220.0*(Mh[i]/1.0e12)**(1.0/3.0) for i in range(len(Mh))],5)
    rh40.irregularVary('accScaleLength',GetScaleLengths(30,Mh0=Mh,scatter=1.0e-10,multiple=0.7),5)
    rh40.irregularVary('mu',[0.5*(Mh[i]/1.0e12)**(-1.0/3.0) for i in range(len(Mh))],5)

    rh41=NewSetOfExperiments(rh40,"rh41")[0]
    rh41.irregularVary('alphaMRI',0)

    rh42=NewSetOfExperiments(rh41,"rh42")[0]
    rh42.irregularVary('xmin',.0005)

    rh43=NewSetOfExperiments(rh41,"rh43")[0]
    rh43.irregularVary('tauHeat',1.0e16)

    rh44=NewSetOfExperiments(rh41,"rh44")[0]
    rh44.irregularVary('Qlim',4.0)

    rh45=NewSetOfExperiments(rh40,"rh45")[0]
    rh45.irregularVary('alphaMRI',0)

    rh46=NewSetOfExperiments(rh40,"rh46")[0]
    rh46.irregularVary('b',1.0)
    rh46.irregularVary('innerPowerLaw',0.25)

    rh48=NewSetOfExperiments(rh40,"rh48")[0]
    rh48.irregularVary('mu',0.5)

    # Study dwarfs
    rh47=NewSetOfExperiments(rg01,"rh47")[0]
    Mh = 1.0e10
    rh47.irregularVary('Mh0',Mh)
    rh47.irregularVary('b',0)
    rh47.irregularVary("R",GetScaleLengths(1,Mh0=[Mh],scatter=1.0e-10,multiple=4.1))
    rh47.irregularVary('vphiR',220.0*(Mh/1.0e12)**(1.0/3.0))
    rh47.irregularVary('accScaleLength',GetScaleLengths(1,Mh0=[Mh],scatter=1.0e-10,multiple=0.7))
    rh47.irregularVary('mu',0.5*(Mh/1.0e12)**(-1.0/3.0))
    rh47.vary('whichAccretionHistory',1000,1019,20,0)



    # Now vary halo mass with some other parameters expected to scale along with it. 
    rg40=NewSetOfExperiments(rg39,"rg40x")
    #[rg40[i].irregularVary('TOL',1.0e-4) for i in range(len(rg40))]
    [rg40[i].irregularVary('b',0) for i in range(len(rg40))]
    rg40[0].irregularVary("R",GetScaleLengths(20,Mh0=Mhlo,scatter=1.0e-10,multiple=4.1),5)
    rg40[0].irregularVary("vphiR",[220.0*(Mhlo[i]/1.0e12)**(1.0/3.0) for i in range(20)],5)
    rg40[0].irregularVary("accScaleLength",GetScaleLengths(20,Mh0=Mhlo,scatter=1.0e-10,multiple=0.7),5)
    rg40[0].irregularVary("mu",[1.0*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(20)],5)
    #rg40[0].irregularVary('xmin',[.002*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(20)],5)

    rg40[1].irregularVary("R",GetScaleLengths(10,Mh0=Mhhi,scatter=1.0e-10,multiple=4.1),5)
    rg40[1].irregularVary("vphiR",[220.0*(Mhhi[i]/1.0e12)**(1.0/3.0) for i in range(10)],5)
    rg40[1].irregularVary("accScaleLength",GetScaleLengths(10,Mh0=Mhhi,scatter=1.0e-10,multiple=0.7),5)
    rg40[1].irregularVary("mu",[0.5*(Mhhi[i]/1.0e12)**(-1.0/3.0) for i in range(10)],5)

    rg40y=NewSetOfExperiments(rg40,"rg40y")
    [rg40y[i].irregularVary('dbg',2+2**12) for i in range(len(rg40y))]

    # The fiducial model with only one passive stellar population (saves time and disk space)
    rg41x=NewSetOfExperiments(rg01,"rg41x")
    rg41x[0].irregularVary("NPassive",1)

    # A mass study with stochastic accretion histories
    rg41=NewSetOfExperiments(rg41x[0],"rg41n",N=2)
    Mhlo = [1.0e10 * 10**(i/100.0) for i in range(200)]
    Mhhi = [1.0e12 * 10**((i+1.0)/100.0) for i in range(100)]

    rg41[0].irregularVary("R",GetScaleLengths(200,Mh0=Mhlo,scatter=1.0e-10,multiple=4.1),5)
    rg41[0].irregularVary("vphiR",[220.0*(Mhlo[i]/1.0e12)**(1.0/3.0) for i in range(200)],5)
    rg41[0].irregularVary("diskScaleLength",GetScaleLengths(200,Mh0=Mhlo,scatter=1.0e-10,multiple=0.35),5)
    rg41[0].irregularVary("accScaleLength",GetScaleLengths(200,Mh0=Mhlo,scatter=1.0e-10,multiple=0.7),5)
    rg41[0].irregularVary("b",0)
    rg41[0].irregularVary("mu",[0.5*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(200)],5)
    rg41[0].irregularVary('whichAccretionHistory',[i+1000 for i in range(200)],5)
    rg41[0].irregularVary('xmin',[.002*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(200)],5)
    rg41[0].irregularVary('Mh0',Mhlo,5)

    rg41[1].irregularVary("R",GetScaleLengths(100,Mh0=Mhhi,scatter=1.0e-10,multiple=4.1),5)
    rg41[1].irregularVary("vphiR",[220.0*(Mhhi[i]/1.0e12)**(1.0/3.0) for i in range(100)],5)
    rg41[1].irregularVary("diskScaleLength",GetScaleLengths(100,Mh0=Mhhi,scatter=1.0e-10,multiple=0.35),5)
    rg41[1].irregularVary("accScaleLength",GetScaleLengths(100,Mh0=Mhhi,scatter=1.0e-10,multiple=0.7),5)
    rg41[1].irregularVary("b",0)
    rg41[1].irregularVary("mu",[0.5*(Mhhi[i]/1.0e12)**(-1.0/3.0) for i in range(100)],5)
    rg41[1].irregularVary('whichAccretionHistory',[i+1500 for i in range(100)],5)
    rg41[1].irregularVary('Mh0',Mhhi,5)

    # Vary the inner powerlaw, i.e. beta0
    rg42=NewSetOfExperiments(rg01,"rg42",N=3)
    rg42[0].vary('innerPowerLaw',-0.5,-0.05,4,0)
    rg42[0].irregularVary('softening',-2)
    rg42[1].vary('innerPowerLaw',0.0,0.5, 4,0)
    rg42[2].vary('innerPowerLaw',.55,.95, 4,0)

    # Mass study but use lognormal accretion histories.
    rg43=NewSetOfExperiments(rg41,"rg43")
    [rg43[i].irregularVary('fscatter',0.3) for i in range(len(rg43))]
    rg43[0].irregularVary('whichAccretionHistory',[-2000+i for i in range(200)],5)
    rg43[1].irregularVary('whichAccretionHistory',[-1500+i for i in range(100)],5)

    # Mass study with lognormal accreiton histories, but use bursts uniform in time, not redshift.
    rg44=NewSetOfExperiments(rg43,"rg44")
    [rg44[i].irregularVary('dbg',2+2**14) for i in range(len(rg44))]

    rg45=NewSetOfExperiments(rg43,"rg45")
    [rg45[i].irregularVary('NChanges',1) for i in range(len(rg45))]

    rg46=NewSetOfExperiments(rg43,"rg46")
    [rg46[i].irregularVary('NChanges',1000) for i in range(len(rg46))]


    rg47=NewSetOfExperiments(rg43,"rg47")
    [rg47[i].irregularVary('fscatter',0.1) for i in range(len(rg47))]

    rg48=NewSetOfExperiments(rg47,"rg48")
    [rg48[i].irregularVary('NChanges',1) for i in range(len(rg48))]

    rg49=NewSetOfExperiments(rg47,"rg49")
    [rg49[i].irregularVary('NChanges',1000) for i in range(len(rg49))]


    rg50=NewSetOfExperiments(rg43,"rg50")
    [rg50[i].irregularVary('fscatter',0.5) for i in range(len(rg50))]

    rg51=NewSetOfExperiments(rg50,"rg51")
    [rg51[i].irregularVary('NChanges',1) for i in range(len(rg51))]

    rg52=NewSetOfExperiments(rg50,"rg52")
    [rg52[i].irregularVary('NChanges',1000) for i in range(len(rg52))]

    rg53=NewSetOfExperiments(rg02,"rg53")
    [rg53[i].irregularVary('alphaAccretionProfile',1.0) for i in range(len(rg53))]

    rg54=NewSetOfExperiments(rg01,"rg54",N=2)
    rg54[0].vary('alphaAccretionProfile',0,.3,3,0)
    rg54[1].vary('alphaAccretionProfile',.35,1.0,6,0)

    rg55=NewSetOfExperiments(rg20,"rg55")
    z0ScaleLengths = GetScaleLengths(400,Mh0=[1.0e12 for k in range(400)],scatter=0.3,multiple=0.7)
    [rg55[i].irregularVary('accScaleLength',z0ScaleLengths,3) for i in range(len(rg55))]
    [rg55[i].irregularVary('R',[z0SL*7 for z0SL in z0ScaleLengths],3) for i in range(len(rg55))]

    rg56=NewSetOfExperiments(rg01,"rg56",N=2)
    rg56[0].vary('yREC',.05,.052,2,0)
    rg56[1].vary('yREC',.056,.07,8,0)

    # Allllright. Investigate issues with acc history and B/T. smaller sample: only 101.
    rg57=NewSetOfExperiments(rg01,"rg57",N=8)
    # rg57a: a repeat of rg20b w/ smaller sample size
    rg57[0].irregularVary('deltaOmega',0.5) 
    rg57[0].vary('whichAccretionHistory',1000,1100,101,0,3)
    # rg57b: same as 57a but with a stricter merger cut
    rg57[1].irregularVary('deltaOmega',0.5) # same but with strict merger cut
    rg57[1].vary('whichAccretionHistory',1000,1100,101,0,3)
    rg57[1].irregularVary('invMassRatio',0.3)
    # rg57c:  an approximation to smooth, with fsc=0
    rg57[2].irregularVary('deltaOmega',0.5) 
    rg57[2].irregularVary('whichAccretionHistory',113)
    rg57[2].irregularVary('fscatter',0.0)
    # rg57d: # another smooth history, with much higher delta Omega smapling
    rg57[3].irregularVary('whichAccretionHistory',-1000) 
    rg57[3].irregularVary('fscatter',0.0)
    rg57[3].irregularVary('NChanges',100)
    # rg57e: stochastic sample w/ many variation points
    rg57[4].vary('whichAccretionHistory',-1100,-1000,101,0,3) 
    rg57[4].irregularVary('NChanges',20)
    rg57[4].irregularVary('fscatter',.3)
    # rg57f: stochastic sample with fewer variation points
    rg57[5].vary('whichAccretionHistory',-1100,-1000,101,0,3) 
    rg57[5].irregularVary('NChanges',4)
    rg57[5].irregularVary('fscatter',.3)
    # rg57g: like 57c but with shorter delta omega
    rg57[6].irregularVary('deltaOmega',0.05)
    rg57[6].irregularVary('whichAccretionHistory',113)
    rg57[6].irregularVary('fscatter',0.0)
    # rg57h: like 57d, but longer changes.
    rg57[7].irregularVary('whichAccretionHistory',-1000)
    rg57[7].irregularVary('fscatter',0.0)
    rg57[7].irregularVary('NChanges',3)



    # Another mass study, with an eye towards understanding th emain squence
    rg58=NewSetOfExperiments(rg41x[0],"rg58",N=2)
    Nsample = 200 
    Mhlo = [1.0e11 * 10**(i/float(Nsample)) for i in range(Nsample)]
    Mhhi = [1.0e12 * 10**((i+1.0)/float(Nsample)) for i in range(Nsample)]
    sclLo = GetScaleLengths(Nsample,Mh0=Mhlo,scatter=1.0e-10,multiple=4.1)
    sclHi = GetScaleLengths(Nsample,Mh0=Mhhi,scatter=1.0e-10,multiple=4.1)

    rg58[0].irregularVary("R",sclLo,5)
    rg58[0].irregularVary("vphiR",[220.0*(Mhlo[i]/1.0e12)**(1.0/3.0) for i in range(Nsample)],5)
    rg58[0].irregularVary("diskScaleLength",[sclLo[i]*.35/4.1 for i in range(Nsample)],5)
    rg58[0].irregularVary("accScaleLength",[sclLo[i]*.7/4.1 for i in range(Nsample)],5)
    rg58[0].irregularVary("b",0)
    rg58[0].irregularVary("mu",[0.5*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg58[0].irregularVary('whichAccretionHistory',[i+1000 for i in range(Nsample)],5)
    rg58[0].irregularVary('xmin',[.002*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg58[0].irregularVary('Mh0',Mhlo,5)

    rg58[1].irregularVary("R",sclHi,5)
    rg58[1].irregularVary("vphiR",[220.0*(Mhhi[i]/1.0e12)**(1.0/3.0) for i in range(Nsample)],5)
    rg58[1].irregularVary("diskScaleLength",[sclHi[i]*.35/4.1 for i in range(Nsample)],5)
    rg58[1].irregularVary("accScaleLength",[sclHi[i]*.7/4.1 for i in range(Nsample)],5)
    rg58[1].irregularVary("b",0)
    rg58[1].irregularVary("mu",[0.5*(Mhhi[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg58[1].irregularVary('whichAccretionHistory',[i+1500 for i in range(Nsample)],5)
    rg58[1].irregularVary('Mh0',Mhhi,5)

    # Random scale lengths (and stochastic accr history)
    rg59=NewSetOfExperiments(rg01,"rg59",N=2)
    [rg59[i].vary('whichAccretionHistory',1001,1200,200,0,3) for i in range(len(rg59))]
    scl = GetScaleLengths(200,Mh0=1.0e12,scatter=0.4,multiple=0.7)
    [rg59[i].irregularVary('accScaleLength',scl,3) for i in range(len(rg59))]
    [rg59[i].irregularVary('R',[scl[j]*5 for j in range(len(scl))],3) for i in range(len(rg59))]
    rg59[0].irregularVary('deltaOmega',.2)
    rg59[1].irregularVary('deltaOmega',.5)

    # Study the Main Sequence using smaller stepsizes
    rg60=NewSetOfExperiments(rg58,"rg60")
    [rg60[i].irregularVary('deltaOmega',0.1) for i in range(len(rg60))]
    [rg60[i].irregularVary('dbg',2+2**3) for i in range(len(rg60))]


    rg61=NewSetOfExperiments(rg01,"rg61",N=1)
    rg61[0].irregularVary('accScaleLength',3.0)
    rg61[0].irregularVary('R',30.0)

    rg62=NewSetOfExperiments(rg36,"rg62")
    [rg62[i].irregularVary('accScaleLength',3.0) for i in range(len(rg62))]
    [rg62[i].irregularVary('R',30.0) for i in range(len(rg62))]

    rg63=NewSetOfExperiments(rg20,"rg63")
    [rg63[i].irregularVary('accScaleLength',3.0) for i in range(len(rg63))]
    [rg63[i].irregularVary('R',30.0) for i in range(len(rg63))]

    rg64=NewSetOfExperiments(rg01,"rg64",N=1)
    rg64[0].irregularVary('accScaleLength',3.0)
    rg64[0].irregularVary('R',30.0)
    rg64[0].irregularVary('dbg',2+2**3)
    rg64[0].irregularVary('deltaOmega',0.1)

    rg65=NewSetOfExperiments(rg36,"rg65")
    [rg65[i].irregularVary('accScaleLength',3.0) for i in range(len(rg65))]
    [rg65[i].irregularVary('R',30.0) for i in range(len(rg65))]
    [rg65[i].irregularVary('dbg',2+2**3) for i in range(len(rg65))]
    [rg65[i].irregularVary('deltaOmega',0.1) for i in range(len(rg65))]


    rg66=NewSetOfExperiments(rg20,"rg66")
    [rg66[i].irregularVary('accScaleLength',3.0) for i in range(len(rg66))]
    [rg66[i].irregularVary('R',30.0) for i in range(len(rg66))]
    [rg66[i].irregularVary('dbg',2+2**3) for i in range(len(rg66))]
    [rg66[i].irregularVary('deltaOmega',0.1) for i in range(len(rg66))]


    rg67=NewSetOfExperiments(rg01,"rg67",N=1)
    rg67[0].irregularVary('dbg',2+2**3)
    rg67[0].irregularVary('deltaOmega',0.1)

    rg68=NewSetOfExperiments(rg36,"rg68")
    [rg68[i].irregularVary('dbg',2+2**3) for i in range(len(rg68))]
    [rg68[i].irregularVary('deltaOmega',0.1) for i in range(len(rg68))]


    rg69=NewSetOfExperiments(rg20,"rg69")
    [rg69[i].irregularVary('dbg',2+2**3) for i in range(len(rg69))]
    [rg69[i].irregularVary('deltaOmega',0.1) for i in range(len(rg69))]

    rg70=NewSetOfExperiments(rg69,"rg70")
    scl = GetScaleLengths(400,Mh0=1.0e12,scatter=0.4,multiple=0.7)
    [rg70[i].irregularVary('accScaleLength',scl,3) for i in range(len(rg70))]
    [rg70[i].irregularVary('R',[scl[j]*5 for j in range(len(scl))],3) for i in range(len(rg70))]

    rg71=NewSetOfExperiments(rg70,"rg71")
    [rg71[i].irregularVary('xmin',.0006) for i in range(len(rg71))]

    # Try the rg20bish thing but with no invMassRatio cut
    rg72=NewSetOfExperiments(rg69,"rg72")
    [rg72[i].irregularVary('invMassRatio',10000.0) for i in range(len(rg72))]

    # Another mass study, with an eye towards understanding th emain squence
    rg73=NewSetOfExperiments(rg41x[0],"rg73",N=2)
    Nsample = 200 
    Mhlo = [1.0e11 * 10**(i/float(Nsample)) for i in range(Nsample)]
    Mhhi = [1.0e12 * 10**((i+1.0)/float(Nsample)) for i in range(Nsample)]
    sclLo = GetScaleLengths(Nsample,Mh0=Mhlo,scatter=0.4,multiple=4.1,upper = 100,lower=7)
    sclHi = GetScaleLengths(Nsample,Mh0=Mhhi,scatter=0.4,multiple=4.1,upper = 100,lower=7)
    [rg73[i].irregularVary('dbg',2+2**3) for i in range(len(rg73))]
    [rg73[i].irregularVary('deltaOmega',0.1) for i in range(len(rg73))]

    rg73[0].irregularVary("R",sclLo,5)
    rg73[0].irregularVary("vphiR",[220.0*(Mhlo[i]/1.0e12)**(1.0/3.0) for i in range(Nsample)],5)
    rg73[0].irregularVary("accScaleLength",[sclLo[i]*.7/4.1 for i in range(Nsample)],5)
    rg73[0].irregularVary("mu",[0.5*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg73[0].irregularVary('whichAccretionHistory',[i+1000 for i in range(Nsample)],5)
    #rg73[0].irregularVary('xmin',[.002*(Mhlo[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg73[0].irregularVary('Mh0',Mhlo,5)

    rg73[1].irregularVary("R",sclHi,5)
    rg73[1].irregularVary("vphiR",[220.0*(Mhhi[i]/1.0e12)**(1.0/3.0) for i in range(Nsample)],5)
    rg73[1].irregularVary("accScaleLength",[sclHi[i]*.7/4.1 for i in range(Nsample)],5)
    rg73[1].irregularVary("mu",[0.5*(Mhhi[i]/1.0e12)**(-1.0/3.0) for i in range(Nsample)],5)
    rg73[1].irregularVary('whichAccretionHistory',[i+1500 for i in range(Nsample)],5)
    rg73[1].irregularVary('Mh0',Mhhi,5)

    rg74=NewSetOfExperiments(rg73,"rg74")
    #rg74[0].irregularVary('xmin',[.002*sclLo[i]/40.0 for i in range(Nsample)],5)
    rg74[0].irregularVary('xmin',.01) 

    # Accretion histories are not stochastic, but accr scales are random
    rg75 = NewSetOfExperiments(rg74,"rg75")
    [rg75[i].irregularVary('whichAccretionHistory',0) for i in range(len(rg75))]

    # Reduce fscatter artificially.
    rg76 = NewSetOfExperiments(rg74,"rg76")
    [rg76[i].irregularVary('fscatter',0.3) for i in range(len(rg76))]

    # Go back to not varying radial scale lengths. Now we have all four combinations of stochastic lambda, Mdot.
    rg77 = NewSetOfExperiments(rg74,"rg77")
    sclLoNV = GetScaleLengths(Nsample,Mh0=Mhlo,scatter=1.0e-10,multiple=4.1,upper=100,lower=7)
    sclHiNV = GetScaleLengths(Nsample,Mh0=Mhhi,scatter=1.0e-10,multiple=4.1,upper=100,lower=7)
    rg77[0].irregularVary("R",sclLoNV,5)
    rg77[0].irregularVary('accScaleLength',[sclLoNV[i]*.7/4.1 for i in range(Nsample)],5)
    rg77[1].irregularVary("R",sclHiNV,5)
    rg77[1].irregularVary('accScaleLength',[sclHiNV[i]*.7/4.1 for i in range(Nsample)],5)

    # Very similar to rg73, designed for stacking analysis.
    # xmin s.t. rmin = 100 pc
    # R which is at least 10 kpc
    # but still a wide range in accScaleLength (0.4 dex)
    rg78=NewSetOfExperiments(rg73,"rg78")
    sclLo = GetScaleLengths(Nsample,Mh0=Mhlo,scatter=0.4,multiple=4.1,upper = 100,lower=1)
    sclLoTr = [sclLo[i] if sclLo[i] > 10 else 10 for i in range(Nsample)]
    sclHi = GetScaleLengths(Nsample,Mh0=Mhhi,scatter=0.4,multiple=4.1,upper = 100,lower=1)
    sclHiTr = [sclHi[i] if sclHi[i] > 10 else 10 for i in range(Nsample)]
    rg78[0].irregularVary('xmin',[.1/(sclLoTr[i]) for i in range(Nsample)],5)
    rg78[0].irregularVary('R',sclLoTr,5)
    rg78[0].irregularVary("accScaleLength",[sclLo[i]*.7/4.1 for i in range(Nsample)],5)
    rg78[1].irregularVary('xmin',[.1/(sclHiTr[i]) for i in range(Nsample)],5)
    rg78[1].irregularVary('R',sclHiTr,5)
    rg78[1].irregularVary("accScaleLength",[sclHi[i]*.7/4.1 for i in range(Nsample)],5)


    rg79=NewSetOfExperiments(rg78,"rg79")
    [rg79[i].irregularVary('whichAccretionHistory',0) for i in range(len(rg79))]

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
