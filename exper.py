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



allModels={} # a global dictionary of all experiments created.

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
        self.p=[name,200,1.5,.01,4.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,1.0,2.0,50.0,int(1e9),1.0e-3,1.0,0,.5,2.0,2.0,0,0.0,1.5,1,2.0,.001,1.0e12,5.0,3.0,0,2.0,-1.0]
        self.p_orig=self.p[:] # store a copy of p, possibly necessary later on.
        self.pl=[self.p[:]] # define a 1-element list containing a copy of p.
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat','analyticQ','cosmologyOn','xmin','NActive','NPassive','vphiR','R','gasTemp','Qlim','fg0','phi0','zstart','tmax','stepmax','TOL','mu','b','innerPowerLaw','softening','diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive','fixedQ','kappaMetals','Mh0','minSigSt','ndecay','dbg','accScaleLength','zquench']
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
        
    def vary(self,name,mmin,mmax,n,log,cov=0):
        '''Set up to vary a particular parameter, name, over the
        range [min,max] with n steps spaced linearly if log=0 or
        logarithmically if log=1. The cov flag means that all variables
        with the same value of cov, e.g. accScaleLength and whichAccretionHistory,
        will be varied together rather than independently. '''
        if(n>1 and log==0):
            self.p[self.keys[name]]=[mmin + (mmax-mmin)*type(self.p[self.keys[name]])(i)/type(self.p[self.keys[name]])(n-1) for i in range(n)]
        elif(n>1 and log==1 and mmin!=0 and mmin!=0.0):
            rat = (float(mmax)/float(mmin))
            typ = type(self.p[self.keys[name]])
            self.p[self.keys[name]]=[typ(float(mmin) * (rat)**(float(i)/float(n-1))) for i in range(n)]
        elif(n==1):
            self.p[self.keys[name]]=mmin # a mechanism for changing the parameter from its default value.
	self.covariables[self.keys[name]] = cov

    def irregularVary(self,name,values,cov=0):
        ''' Similar to vary, except, rather than generating the
        array of parameter values to use for parameter "name",
        this list is specified explicitly in "values"'''
        keyNum = self.keys[name]
        typ = type(self.p[keyNum])
        if(type(values) == type([1])): # if values is a list..
            # replace the given index of p with the list given by values.
            self.p[keyNum]=[typ(value) for value in values]
        else:
            self.p[keyNum]=[typ(values)]
	self.covariables[self.keys[name]] = cov
    
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
                        	else:
	                            base='a'
        	                    app=str((i-(i%26))/26)
                	        a_p[self.keys['name']]+=chr(ord(base)+(i%26))+app
	                    # end loop over p's in pl2[i]
        	        # end loop over range of this varied parameter
		
			# collapse pl to be a 1d array of p's.
	                self.pl=[element for sub in pl2 for element in sub]
		
            else : #just a regular non-varying parameter
                # If the user has used vary() but with N=1, the new argument will
                # not be a list! In this case, we set this parameter in each element 
		# of pl to be the value the user picked, instead of the default value.
                for i in range(len(self.pl)):
                    if(self.pl[i][j]!=param):
                        self.pl[i][j]=param
                    else:
                        pass # nothing to do- just use default value


            
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
                    print "Parameters: "
                    print [binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]]
                    os.chdir(expDir)
                    procs.append(subprocess.Popen([binary]+tmpap[:1]+[repr(el) for el in tmpap[1:]],stdout=stdo,stderr=stde))
                    nPrinted=True

            # we've started a process off and running, but there are
            # probably more processes waiting to go. We do not want
            # to exceed the number of processors nproc the user
            # is willing to let run at once, so let's check what
            # our processes are up to:
            while True:
                nStillRunning=HowManyStillRunning(procs)
    
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


def GetScaleLengths(N,mean=0.045,scatter=0.5,Mh0=1.0e12,sd=100,lower=0.0,upper=1.0e3):
    ''' Return a vector of scale lengths in kpc such that the haloes will have
        spin parameters of mean on average, with a given scatter in dex. These 
        are also scaled by halo mass, since the Virial radius goes as M_h^(1/3).
        sd seeds the random number generator so that this function is deterministic.'''
    random.seed(sd)
    scLengths=[]
    while len(scLengths) < N:
	spinParameter = mean  * (10.0 ** random.gauss(0.0,scatter)) # lognormal w/ scatter in dex
        length =  (311.34 * (Mh0/1.0e12)**(1.0/3.0) * spinParameter/(2.0**.5))
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



if __name__ == "__main__":

    # http://docs.python.org/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(description='Analyze data and/or run experiments associated with a list of experiment names.')
    parser.add_argument('models', metavar='experiment', type=str, nargs='+',
                   help='an experiment to be run / analyzed')
    parser.add_argument('--nproc',type=int,help="maximum number of processors to use (default: 16)",default=16)
    parser.add_argument('--start',metavar='startingModel',type=int,
                   help='The number of the model in the experiment (as ordered by GeneratePl) at which we will start sending experiments to the processors to run. (default: 0)',default=0)
    parser.add_argument('--xgrid',type=bool,help="run on an xGrid (requires the user to submit the generated file to an xGrid (default: False)",default=False)
    args = parser.parse_args()
    
    modelList=[]
    for modelName in args.models:
        modelList.append(modelName)

    
#    allModels={} # empty dictionary into which we'll place all model definitions
    # keys are experiment names, e.g. "rn27", values are experiment instances.
    # There is no particular requirement that the name of the instance in this script
    # be the same as the instance's name variable, nor the value of the key,
    # but confusion is likely if that is not the case.
 

    rn27=experiment('rn27') # detailed view of MW-size halo
    rn27.vary('nx',200,200,1,0)
    rn27.vary('diskScaleLength',2,2,1,0)
    rn27.vary('dbg',2**10,2**10,1,0)
    rn27.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn27.vary('alphaMRI',0,0,1,0)
    rn27.vary('minSigSt',10,10,1,0)
    rn27.vary('whichAccretionHistory',100,1099,1000,0)

    rn28=experiment('rn28') # detailed view of a lower-mass halo
    rn28.vary('nx',200,200,1,0)
    rn28.vary('diskScaleLength',2,2,1,0)
    rn28.vary('dbg',2**10,2**10,1,0)
    rn28.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn28.vary('alphaMRI',0,0,1,0)
    rn28.vary('minSigSt',5,5,1,0)
    rn28.vary('whichAccretionHistory',1100,2099,1000,0)
    rn28.vary('Mh0',1.0e11,1.0e11,1,0)

    rn29=experiment('rn29') # More accretion!
    rn29.vary('nx',200,200,1,0)
    rn29.vary('diskScaleLength',2,2,1,0)
    rn29.vary('dbg',2**10,2**10,1,0)
    rn29.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn29.vary('alphaMRI',0,0,1,0)
    rn29.vary('minSigSt',10,10,1,0)
    rn29.vary('whichAccretionHistory',2100,3099,1000,0)
    rn29.vary('Mh0',3.0e12,3.0e12,1,0)
#    allModels['rn29']=rn29

    rn34=experiment('rn34') # variable metal diffusion coefficient; b-> lower zeta (.3)
    # c: zeta=.7, multiply variable coefficient by .1.
    rn34.vary('dbg',2**10+2**15,2**10+2**15,1,0)
    rn34.vary('minSigSt',10,10,1,0)
    rn34.vary('whichAccretionHistory',3100,4099,1000,0)
 #   allModels['rn34']=rn34

    rn35=experiment('rn35') # increased accr. history threshold.
    rn35.vary('dbg',2**10+2**16,2**10+2**16,1,0)
    rn35.vary('minSigSt',10,10,1,0)
    rn35.vary('whichAccretionHistory',4100,5099,1000,0)
#    allModels['rn35']=rn35

    rn36=experiment('rn36') # smaller eta
    rn36.vary('dbg',2**10,2**10,1,0)
    rn36.vary('minSigSt',10,10,1,0)
    rn36.vary('whichAccretionHistory',5100,6099,1000,0)
    rn36.vary('eta',.5,.5,1,0)
#    allModels['rn36']=rn36

    rn37=experiment('rn37') # larger eta
    rn37.vary('dbg',2**10,2**10,1,0)
    rn37.vary('minSigSt',10,10,1,0)
    rn37.vary('whichAccretionHistory',6100,7099,1000,0)
    rn37.vary('eta',4.5,4.5,1,0)
#    allModels['rn37']=rn37


    rn38=experiment('rn38') # lower fg
    rn38.vary('dbg',2**10,2**10,1,0)
    rn38.vary('minSigSt',10,10,1,0)
    rn38.vary('whichAccretionHistory',6100,7099,1000,0)
    rn38.vary('fg0',.1,.1,1,0)
#    allModels['rn38']=rn38

    rn39=experiment('rn39') # hotter stars!
    rn39.vary('dbg',2**10,2**10,1,0)
    rn39.vary('minSigSt',10,10,1,0)
    rn39.vary('whichAccretionHistory',6100,7099,1000,0)
    rn39.vary('phi0',2,2,1,0)
#    allModels['rn39']=rn39
    

    rn40=experiment('rn40') # higher fg
    rn40.vary('dbg',2**10,2**10,1,0)
    rn40.vary('minSigSt',10,10,1,0)
    rn40.vary('whichAccretionHistory',6100,7099,1000,0)
    rn40.vary('fg0',.9,.9,1,0)
#    allModels['rn40']=rn40

    rn41=experiment('rn41') # weaker feedback
    rn41.vary('dbg',2**10,2**10,1,0)
    rn41.vary('minSigSt',10,10,1,0)
    rn41.vary('whichAccretionHistory',7100,8099,1000,0)
    rn41.vary('fg0',.5,.5,1,0)
    rn41.vary('mu',.1,.1,1,0)

    rn42=experiment('rn42') # modified potential
    rn42.vary('dbg',2**10,2**10,1,0)
    rn42.vary('minSigSt',10,10,1,0)
    rn42.vary('whichAccretionHistory',8100,9099,1000,0)
    rn42.vary('innerPowerLaw',.5,.5,1,0)
    rn42.vary('softening',4,4,1,0)
    rn42.vary('b',10,10,1,0)

    rn46=experiment('rn46') # low SF efficiency & weak feedback
    rn46.vary('dbg',2**10,2**10,1,0)
    rn46.vary('minSigSt',10,10,1,0)
    rn46.vary('whichAccretionHistory',6100,7099,1000,0)
    rn46.vary('eta',1.5,1.5,1,0)
    rn46.vary('epsff',.005,.005,1,0)
    rn46.vary('mu',.1,.1,1,0)
    rn46.vary('NPassive',10,10,1,0)

    rn48=experiment('rn48') # low gas temperature
    rn48.vary('dbg',2**10,2**10,1,0)
    rn48.vary('minSigSt',10,10,1,0)
    rn48.vary('whichAccretionHistory',6100,7099,1000,0)
    rn48.vary('gasTemp',100,100,1,0)


    rn49=experiment('rn49') # const depletion time.
    rn49.vary('dbg',2**10+2**19,2**10+2**19,1,0)
    rn49.vary('minSigSt',10,10,1,0)
    rn49.vary('whichAccretionHistory',6100,7099,1000,0)


    rn50=experiment('rn50') # test new metallicity diffusion prescription
    rn50.vary('dbg',2**10,2**10,1,0)
    rn50.vary('minSigSt',10,10,1,0)
    rn50.vary('whichAccretionHistory',6100,7099,1000,0)


    # a new era, with several minor or confusing issues fixed, and using the same set accretion histories.
    rn51=experiment('rn51') # baseline.
    rn51.vary('dbg',2**10,2**10,1,0)
    rn51.vary('minSigSt',10,10,1,0)
    rn51.vary('whichAccretionHistory',1000,1399,400,0)


    rn52=experiment('rn52') # condition-dependent metal diffusion with fixed Z_boundary
    rn52.vary('dbg',2**10+2**15+2**1,2**10+2**15+2**1,1,0)
    rn52.vary('minSigSt',10,10,1,0)
    rn52.vary('whichAccretionHistory',1000,1399,400,0)

    rn53=experiment('rn53') # condition-dependent metal diffusion with fixed Z_boundary and zeta=.7
    rn53.vary('dbg',2**10+2**15+2**1+2**13,2**10+2**15+2**1+2**13,1,0)
    rn53.vary('minSigSt',10,10,1,0)
    rn53.vary('whichAccretionHistory',1000,1399,400,0)


    rn54=experiment('rn54') # condition-dependent metal diffusion with zero-outflow BC and zeta=.3
    rn54.vary('dbg',2**10+2**15+2**14,2**10+2**15+2**14,1,0)
    rn54.vary('minSigSt',10,10,1,0)
    rn54.vary('whichAccretionHistory',1000,1399,400,0)


    rn55=experiment('rn55') # IBC = -.01*x(0)
    rn55.vary('dbg',2**10+2**17,2**10+2**17,1,0)
    rn55.vary('minSigSt',10,10,1,0)
    rn55.vary('whichAccretionHistory',1000,1399,400,0)


    rn56=experiment('rn56') # IBC = -x(0)
    rn56.vary('dbg',2**10+2**18,2**10+2**18,1,0)
    rn56.vary('minSigSt',10,10,1,0)
    rn56.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn56']=rn56

    rn57=experiment('rn57') # const depletion time
    rn57.vary('dbg',2**10+2**19,2**10+2**19,1,0)
    rn57.vary('minSigSt',10,10,1,0)
    rn57.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn57']=rn57

    rn58=experiment('rn58') # F=0 only
    rn58.vary('dbg',2**10+2**12+2**9,2**10+2**12+2**9,1,0)
    rn58.vary('minSigSt',10,10,1,0)
    rn58.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn58']=rn58

    rn59=experiment('rn59') # low gas temperature
    rn59.vary('dbg',2**10,2**10,1,0)
    rn59.vary('gasTemp',100,100,1,0)
    rn59.vary('minSigSt',10,10,1,0)
    rn59.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn59']=rn59

    rn60=experiment('rn60') # rotation curve
    rn60.vary('dbg',2**10,2**10,1,0)
    rn60.vary('minSigSt',10,10,1,0)
    rn60.vary('whichAccretionHistory',1000,1399,400,0)
    rn60.vary('innerPowerLaw',.5,.5,1,0)
    rn60.vary('softening',4,4,1,0)
    rn60.vary('b',5,5,1,0)
#    allModels['rn60']=rn60

    rn61=experiment('rn61') # try everything...
    rn61.vary('minSigSt',10,10,1,0)
    rn61.vary('eta',.5,2.5,3,0)   # 3
    rn61.vary('epsff',.005,.02,3,1) # 3
    rn61.vary('gasTemp',100,7000,2,0) # 2
    rn61.vary('fg0',.5,.5,1,0)        # 1
    rn61.vary('phi0',1,3.16,2,1)   # 2
    rn61.vary('mu',1,1,1,1)     # 1
    rn61.vary('innerPowerLaw',.5,.5,1,0) 
    rn61.vary('softening',4,4,1,0)
    rn61.vary('b',0,5,2,0)            # 2
    rn61.vary('diskScaleLength',2,3.5,2,1) # should include 2.  --- 2
    rn61.vary('fixedQ',2.0,3.0,2,1) # should include 2 --- 2
    rn61.irregularVary('dbg',[2**10,2**10+2**19,2**10+2**17,2**10+2**15]) # 4 -- base, const depl time, IBC, metal diffusion
#    allModels['rn61']=rn61

    rn62=experiment('rn62') #test new code with non-outer-boundary-accretion infrastructure, but keep that switch off.
    rn62.vary('dbg',2**10,2**10,1,0)
    rn62.vary('minSigSt',10,10,1,0)
    rn62.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn62']=rn62

    rn63=experiment('rn63') # test new code with non-outer-boundary-accretion
    rn63.vary('dbg',2**10+2**8,2**10+2**8,1,0)
    rn63.vary('minSigSt',10,10,1,0)
    rn63.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn63']=rn63

    rn64=experiment('rn64') # test new code with non-outer-boundary-accretion; this time with acc scale length 3x the stellar scale length
    rn64.vary('dbg',2**10+2**8,2**10+2**8,1,0)
    rn64.vary('minSigSt',10,10,1,0)
    rn64.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn64']=rn64

    # a series of very-small-number experiments to really see the effects of varying individual parameters.
    rn65 = experiment('rn65')
    rn65.vary('epsff',.003,.03,10,1)
    rn65.vary('dbg',2**10,2**10,1,0)
    rn65.vary('minSigSt',10,10,1,0)
#    allModels['rn65']=rn65

    rn66 = experiment('rn66')
    rn66.vary('fg0',.1,.9,9,0)
    rn66.vary('dbg',2**10,2**10,1,0)
    rn66.vary('minSigSt',10,10,1,0)
#    allModels['rn66']=rn66

    rn67 = experiment('rn67')
    rn67.vary('fixedQ',1.5,3.5,9,1)
    rn67.vary('dbg',2**10,2**10,1,0)
    rn67.vary('minSigSt',10,10,1,0)
#    allModels['rn67']=rn67

    rn68 = experiment('rn68')
    rn68.vary('Qlim',1.5,3.5,9,1)
    rn68.vary('dbg',2**10,2**10,1,0)
    rn68.vary('minSigSt',10,10,1,0)
#    allModels['rn68']=rn68

    rn6768 = experiment('rn6768')
    rn6768.vary('Qlim',1.5,3.5,9,1)
    rn6768.vary('fixedQ',1.5,3.5,9,1)
    rn6768.vary('dbg',2**10,2**10,1,0)
    rn6768.vary('minSigSt',10,10,1,0)
#    allModels['rn6768']=rn6768

    rn69 = experiment('rn69')
    rn69.vary('eta',.5,4.5,9,1)
    rn69.vary('dbg',2**10,2**10,1,0)
    rn69.vary('minSigSt',10,10,1,0)
#    allModels['rn69']=rn69

    rn70 = experiment('rn70')
    rn70.irregularVary('dbg',[2**10+2**19,2**10+2**17])
    rn70.vary('minSigSt',10,10,1,0)
#    allModels['rn70']=rn70  

    rn90 = experiment('rn90')
    rn90.vary('minSigSt',10,10,1,0)
    rn90.vary('dbg',2**10+2**8,2**10+2**8,1,0)
    rn90.vary('accScaleLength',1,10,10,0)
#    allModels['rn90']=rn90

    rn91 = experiment('rn91')
    rn91.vary('minSigSt',10,10,1,0)
    rn91.vary('dbg',2**10+2**8+4,2**10+2**8+4,1,0)
    rn91.vary('accScaleLength',2,2,1,0)
    rn91.vary('whichAccretionHistory',100,1399,400,0)
#    allModels['rn91']=rn91

    rn92 = experiment('rn92')
    rn92.vary('minSigSt',10,10,1,0)
    rn92.vary('dbg',2**10+2**8+4,2**10+2**8+4,1,0)
    rn92.vary('accScaleLength',10,10,1,0)
    rn92.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn92']=rn92

    rn93 = experiment('rn93')
    rn93.vary('minSigSt',10,10,1,0)
    rn93.vary('dbg',2**10+2**8+4,2**10+2**8+4,1,0)
    rn93.vary('accScaleLength',300,300,1,0)
    rn93.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn93']=rn93

    rn94 = experiment('rn94')
    rn94.vary('minSigSt',10,10,1,0)
    rn94.vary('dbg',2**10+4,2**10+4,1,0)
    rn94.vary('accScaleLength',10,10,1,0)
    rn94.vary('whichAccretionHistory',1000,1399,400,0)
#    allModels['rn94']=rn94


    # Unrelaxed initial conditions.
    # 104: simple exponentials
    # 105: include a correction for very large scalelengths in which virtually no mass was put in the disk initially
    # 106: rerun with a low initial value of sigma (cold disks). Seems to fail pretty badly.
    # The above had 400 runs each. Decrease to 200 to speed up the process in the subsequent runs:
    rn104 = experiment('rn106')
    rn104.vary('nx',200,200,1,0)
    rn104.vary("TOL",1.0e-3,1.0e-3,1,0)
    rn104.vary('dbg',2**10+2**8+2**4+2**5,2**10+2**8+2**4+2**5,1,0)
    scLengths = GetScaleLengths(200,Mh0=1.0e12)
    rn104.irregularVary('accScaleLength',scLengths,1)
    rn104.irregularVary('diskScaleLength',scLengths,1)
    rn104.vary('whichAccretionHistory',1000,1199,200,0,1)

    # start w/ sig0= 50 km/s
    # rn108: back up to 100 km/s.
    # rn109: add in dbg.opt(12), which will set tau=0 instead of the mass flux in stable regions.
    # rn110... and back down to 8km/s
    rn107=experiment('rn110')
    rn107.vary('nx',200,200,1,0)
    rn107.vary("TOL",1.0e-3,1.0e-3,1,0)
    rn107.irregularVary('dbg',[2**10+2**8+2**4+2**5+2**12])
    rn107.vary('accScaleLength',.1,200,30,1,1)
    rn107.vary('diskScaleLength',.1,200,30,1,1)
    

    # 115: simple variation with halo mass, with no monotonicity in Mdot(!)
    # 116: 115 seems to work, so now let's try it with variable accretion histories.
    rn115=experiment('rn116')
    rn115.vary('nx',200,200,1,0)
    rn115.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn115.irregularVary('dbg',[2**10+2**8+2**5+2**12]) # without the Mdot>=0 fix. i.e. dbg.opt(4)
    rn115.vary('Mh0',1.0e10,1.0e12,21,1)
    rn115.vary('whichAccretionHistory',301,320,20,0)

    rn119=experiment("rp1")
    rn119.irregularVary('dbg',[2**10+2**8+2**5+2**12])
    rn119.vary('Mh0',1.0e11,1.0e11,1,0)
    scl = GetScaleLengths(100,Mh0=1.0e11, sd=112, upper=10.0)
    rn119.irregularVary('accScaleLength',scl,1)
    rn119.irregularVary('diskScaleLength',scl,1)
    rn119.vary('zquench',1.0,1.0,1,0)
    rn119.vary('fg0',.5,.9,2,0)

    # rp2: zquench=-1 (with failed runs with zquench=1)
    # rp2a: bug fixed, only run for zquench=1
    # rp2b: h'okay.
    rp2=experiment("rp2b")
    rp2.irregularVary('dbg',[2**10+2**8+2**5+2**12])
    rp2.vary('Mh0',1.0e10,1.0e12,15,1)
    rp2.vary('accScaleLength',2.0,2.0,1,0)
    rp2.vary('diskScaleLength',2.0,2.0,1,0)
    rp2.vary('zquench',-1.0,-1.0,1,0)

    # rp8 varies scalelengths & kappa but not diffusion prescription
    rp8=experiment("rp10")
    rp8.irregularVary('dbg',[2**10+2**8+2**5+2**12+2**15, 2**10+2**8+2**5+2**12])
    scl=GetScaleLengths(1,Mh0=1.0e12,sd=1.0e-10)[0]
    #scls=[scl/3,scl,3*scl]
    scls=[scl]
    rp8.irregularVary('accScaleLength',scls,2)
    rp8.irregularVary('diskScaleLength',scls,2)
    rp8.vary('kappaMetals',.000001,1000000,5,1)

    rp10a=experiment("rp10b") #10a: time-dependent diffusion. 10b - time-independent.
    rp10a.irregularVary('dbg',[2**10+2**8+2**5+2**12])
    rp10a.irregularVary('accScaleLength',scls,2)
    rp10a.irregularVary('diskScaleLength',scls,2)
    rp10a.vary('kappaMetals',1.0e-6,1.0,7,1)


#    nacc=90
#    nmh0=61
#    mh0Min = 1.0e9
#    mh0Max = 1.0e12
#    rn33 = HaloMassStudy("rn33",allModels,nacc,nmh0,mh0Min,mh0Max)
#    rn33.vary('nx',200,200,1,0)
#    rn33.vary('diskScaleLength',2,2,1,0)
#    rn33.irregularVary('dbg',[2**10])

    # rn73: R \propto 10 kpc, up to 1e13 MSol = Mh0
    # rn83: R \propto 20 kpc, up to 1e13 MSol = Mh0
    # rn84: R \propto 20 kpc, up to 1e12 MSol = Mh0
    # rn85: R \propto 20 kpc, smaller initial stellar disk.
    # rn86: exponential accretion
    # 87 - exponential accretion w/ bugfix
    # 97 - exponential accretion w/ scatter in spin parameter!
    # 99 - variable spin, 100, w/ primordial gas infall.
    # 101 - fixed spin
    # 102 - test a fix for gaps in gas column density.
    # 103 - back to varying spin parameter. (102 seems to work well)
    # rn111 - try exp initialization with smaller scatter in lambda.

    # 113: latest reasonable spread in lambda, w/ upper cutoff in scale length.
    # 114: constant lambda.
    # 117: abandon monotonicity, return to realistic spread in lambda.
    # 118: zquench = .5
    # rp3: increase baryon content by factor of 5 to be ~ consistent w/ Behroozi 2012b
    # rp4: give low-mass galaxies higher gas fractions.
    # rp5: fix lambda to be .045
    # rp6: include dbg 16, which ~eliminates the major merger restriction on accretion histories
    # rp7: back to normal restrictions, try out irregular rotation curves.
    # rp9: increased SF efficiency
    # rp11: back to regular SF efficiency. Try eps_in=1.0 (dbg.opt(6))
    # rp12: dbg.opt(3), i.e. Neistein 2010, turn off eps_in=1.0 for the moment. It seems this will also require dbg.opt(16), i.e. keep the major mergers.
    # rp13: for comparison to rp12. Keep dbg.opt(3), turn off dbg.opt(16). Also, from here on in, 100% of runs will be used, not 30ish %, so I'm reducing nacc from 90 to 32.
    # rp14: try hybrid SF law (dbg.opt(7)). Turn off dbg.opt(3), i.e. should be comparable to rp5.
    # rp15: again with the hybrid SF law. This time include dbg.opt(6) which implements CAFG efficienciency, not 1.0
    # rp16: keep all dbg options of rp15, but do the experiment of having a wide variance in accretion scale lengths but no variance in initial scale lengths.
    # rp17: same as rp15, but now do uniform col accretion, rather than exponential (dbg.opt(0)). Note that this requires, to be self-consistent in terms of angular momentum conservation of the accretion, that the accretion scale length be 3x larger.

    theCompositeName="rp17"



    experiments=[]
    nacc = 32
    nmh0 = 41
    compositeNames=[]
    for i in range(nmh0): # 100 different M_h values, each with 100 different accr. histories.
        theName = theCompositeName+'_'+str(i).zfill(3)+'_'
        compositeNames.append(theName)
        experiments.append(experiment(theName))
        experiments[i].vary('nx',200,200,1,0)
        experiments[i].irregularVary('dbg',[2**10+2**8+2**5+2**12+2**7+2**6+2**0])
	experiments[i].vary('TOL',1.0e-3,1.0e-3,1,0)
        experiments[i].vary('alphaMRI',0,0,1,0)
        Mh0 = 1.0e10 * (100.0)**(float(i)/float(nmh0-1))
        experiments[i].irregularVary('minSigSt',[4.0+(float(i)/float(nmh0-1))*(10.0-4.0)])
        experiments[i].irregularVary('R',[30.0])
        scLengths  =GetScaleLengths(nacc,mean=0.045, scatter=.5, Mh0=Mh0, sd=4+20*nacc, lower=1.0,upper=10.0)
        scLengthsNV=GetScaleLengths(nacc,mean=0.045, scatter=1.0e-14, Mh0=Mh0, sd = 4+20*nacc,lower=1.0,upper=10.0)
	experiments[i].irregularVary('accScaleLength',scLengthsNV*3,1)
        experiments[i].irregularVary('diskScaleLength',scLengthsNV,1)
        #experiments[i].irregularVary('b',scLengths,1)
        #experiments[i].irregularVary('innerPowerLaw',[.5])
        #experiments[i].irregularVary('softening',[4])
        #experiments[i].vary('epsff',.02,.03,1,0)
        experiments[i].irregularVary('Mh0',[Mh0])
        experiments[i].irregularVary('fg0',.4*(Mh0/1.0e12)**(-.13))
        experiments[i].vary('whichAccretionHistory',4+nacc*(i+20),4+nacc*(i+21)-1,nacc,0,1)
        #  experiments[i].vary('zquench',0.5,0.5,1,0)
        #  allModels[theName]=experiments[i]

    compositeFlag = False    
    if (theCompositeName in modelList):
        compositeFlag = True
        for ex in experiments:
            modelList.append(ex.expName)
        del modelList[modelList.index(theCompositeName)]

    successTables=[]
    # run all the models, and record which ones succeed.
    for model in modelList:
        if(not args.xgrid): #local run
          allModels[model].localRun(args.nproc,args.start)
        else: # write a file to run on the xgrid
          allModels[model].write('runExperiment_'+model+'.txt')
#        successTables.append(allModels[model].ExamineExperiment())

    if(compositeFlag):
	os.mkdir('/Users/jforbes/gidget/analysis/'+theCompositeName)
        for dirname in compositeNames:
            files=glob.glob('/Users/jforbes/gidget/analysis/'+dirname+'/*')
            for aFile in files:
	        os.symlink(aFile,'/Users/jforbes/gidget/analysis/'+theCompositeName+'/'+aFile[32+len(dirname):])
           

#    PrintSuccessTables(successTables)


