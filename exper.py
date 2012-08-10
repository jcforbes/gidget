import copy
import os, argparse
import sys
import shutil
import subprocess
import time
import math
import numpy as np

# This is a somewhat straightforward script to allow you to run the GIDGET code with
# a wide variety of systematically varied parameters. The functions here are described
# in some detail, but if you skip to the end of the script, there are plenty of 
# examples of what you can run from the command line, e.g.
#  $ python exper.py 'ex1'
# will execute the example used in the README file.


class covary:
    def __init__(self):
        self.names=[]
        # listOfVariants contains all the necessary information to covary two parameters together.
        self.listOfVariants=[] # [ [p_ind1^A, p_ind1^B, p_ind1^C, ...], [p_ind2^A, p_ind2^B,...], ...]
    def addVariable(self,name,theList):
        if(len(self.listOfVariants) == 0 and len(listOfVariants[0]) == len(theList)):
            self.names.append(name)
            self.listOfVariants.append(theList)
        else:
            raise ValueError('Length of given list incompatible with first list')
    def getIthVar(self,i):
        return self.names[i],self.listOfVariants[i]
    def getNVars(self):
        return len(self.names)
    def getNVals(self):
        return len(self.listOfVariants[0])
            

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
    def __init__(self,name):
        # fiducial model
        self.p=[name,200,1.5,.01,4.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,2.0,2.0,50.0,int(1e9),1.0e-3,1.0,0,.5,2.0,2.0,0,0.0,1.5,1,2.0,.001,1.0e12,5.0,3.0,0]
        self.p_orig=self.p[:] # store a copy of p, possibly necessary later on.
        self.pl=[self.p[:]] # define a 1-element list containing a copy of p.
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat','analyticQ','cosmologyOn','xmin','NActive','NPassive','vphiR','R','gasTemp','Qlim','fg0','phi0','zstart','tmax','stepmax','TOL','mu','b','innerPowerLaw','softening','diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive','fixedQ','kappaMetals','Mh0','minSigSt','ndecay','dbg']
        self.keys={}
        ctr=0
        for n in self.names:
            self.keys[n]=ctr
            ctr=ctr+1
        self.expName=name
	self.covariables=[]

        # store the location of various expected subdirectories in the gidget distribution.
        self.base=os.getcwd() # Assume we are in the base directory - alter this to /path/to/gidget/directory if necessary
        self.src=self.base+'/src'
        self.analysis=self.base+'/analysis'
        self.bin=self.base+'/bin'
        self.out=self.base+'/output'
        self.xgrid=self.base+'/xgrid'
        
    def vary(self,name,mmin,mmax,n,log):
        '''Set up to vary a particular parameter, name, over the
        range [min,max] with n steps spaced linearly if log=0 or
        logarithmically if log=1 '''
        if(n>1 and log==0):
            self.p[self.keys[name]]=[mmin + (mmax-mmin)*type(self.p[self.keys[name]])(i)/type(self.p[self.keys[name]])(n-1) for i in range(n)]
        elif(n>1 and log==1 and mmin!=0 and mmin!=0.0):
            rat = (float(mmax)/float(mmin))
            typ = type(self.p[self.keys[name]])
            self.p[self.keys[name]]=[typ(float(mmin) * (rat)**(float(i)/float(n-1))) for i in range(n)]
        elif(n==1):
            self.p[self.keys[name]]=mmin # a mechanism for changing the parameter from its default value.

    def irregularVary(self,name,values):
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
                # make len(param) copies of the current pl
                pl2=[]
                for i in range(len(param)):
                    f=copy.deepcopy(self.pl)
                    pl2.append(f) # pl2 ~ [pl,pl,pl]~[[p,p,p],[p,p,p],[p,p,p]]
                    # in each copy, set the jth parameter in p to param[i]
                    for a_p in pl2[i]: # each element is a list of parameters for
                        a_p[j]=param[i]
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


        # OK. The strategy is to take the existing pl ~ [[p,p,p..],[p,p,p..],...]
        # and create cov.getNVals() copies of it stored in pl2 ~ [pl, pl, ...]
        # For jth pl in pl2, change all the parameters listed in cov to the jth value of their list.
        # Repeat for every set of covarying variables
        
#        # For every set of covarying variables:
#	for i in range(len(self.covariables)):
#            pl2=[]
#            cov = self.covariables[i]
#            # for each variable which is covarying
#            for varj in range(cov.getNVals()):
#                preCovPl = copy.deepcopy(self.pl)
#                pl2.append(preCovPl)
#                for k in range(cov.getNVars()):
#                name,listForThisVar = cov.GetIthVar(varj)
#                self.keys[name]

            
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

#        for 
#        for flatInd in range(len(successTableKeys.size)):
#            theKey = successTableKeys.flat[flatInd]
#	    message=""
#            for i in range(len(theKey)):
#                message+=(theLists[theKey[i]]+' ')
#
#        for a_p in self.pl: # for each run
#            for success in range(4): # for each possible success code
#                successTableKey = [] 
#                for i in range(len(varied)): # for each quantity varied
#                    j = varied[i]
#                    if(successTabl
#
#        for varA in range(len(varied)):
#	    for varB in range(len(varied)):
#		if(varA!=varB) :
#		    pass
	            

    def fast(self):
        '''Set a few parameters to still-reasonable parameters
        which will speed up the simulation '''
        self.vary('nx',200,200,1,0)
        self.vary('minSigSt',10,10,1,0)
        self.vary('diskScaleLength',2.0,2.0,1,0) # not actually faster, but I often forget it.

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
                nStillRunning=0
                # count how many of our processes are still running
                for proc in procs:
                    if(proc.poll()==None): # still running
                        nStillRunning+=1
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
            nStillRunning=0
            # count how many processes are still running
            for proc in procs:
                if(proc.poll()==None): # still running
                    nStillRunning+=1
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


# Take the results of all experiments named when this script was called
# and sum them up such that it's obvious which runs succeeded and which failed.
# A schematic example output for a situation where 5 different experiments were
# run, each of which varied the same 2 variables, 2 values for 1, 3 for the other:
#
#   [ [ 00000,  00000 ], [ 00010, 00010], [00011, 00011] ]
#
# indicates that, as far as whether the runs crash or succeed,
# the 2-value variable doesn't make a difference, but for increasing
# values of the 3-value variable, fewer runs succeed. The numbers which
# appear here are explained in the successCode function above.
def PrintSuccessTables(successTables):
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

    
    allModels={} # empty dictionary into which we'll place all model definitions
 
#    expName=sys.argv[1]
#    a=experiment(expName)


    rn1=experiment('rn1') # originally rk69
    rn1.fast()
    rn1.vary('dbg',0,31,32,0)
    rn1.vary('ndecay',-.02,.1,4,0) #-.02, .02, .06, .1 
    rn1.vary('whichAccretionHistory',526,526,1,0)
    allModels['rn1']=rn1
    
    rn2=experiment('rn2') #originally rk70k1c, then rk71
    rn2.fast()
    rn2.vary('dbg',0,31,32,0)
    rn2.vary('ndecay',-.02,.1,4,0)
    rn2.vary('whichAccretionHistory',636,636,1,0)
    allModels['rn2']=rn2    

    rn3=experiment('rn3') # originally rk70b5c, then rk72
    rn3.fast()
    rn3.vary('dbg',0,31,32,0)
    rn3.vary('ndecay',-.02,.1,4,0)
    rn3.vary('whichAccretionHistory',731,731,1,0)
    allModels['rn3']=rn3

    rn4=experiment('rn4') # originally rk70b6c, then rk73
    rn4.fast()
    rn4.vary('dbg',0,31,32,0)
    rn4.vary('ndecay',-.02,0.1,4,0)
    rn4.vary('whichAccretionHistory',757,757,1,0)
    allModels['rn4']=rn4

    rn5=experiment('rn5') # originally rk70f6c, then rk73 (again?!)
    rn5.fast()
    rn5.vary('dbg',0,31,32,0)
    rn5.vary('ndecay',-.02,0.1,4,0)
    rn5.vary('whichAccretionHistory',761,761,1,0)
    allModels['rn5']=rn5

    rn6=experiment('rn6') # originally rk70o6c
    rn6.fast()
    rn6.vary('dbg',0,31,32,0)
    rn6.vary('ndecay',-.02,0.1,4,0)
    rn6.vary('whichAccretionHistory',770,770,1,0)
    allModels['rn6']=rn6


    # similar to rk70, except with only ndecay=-1, alpha=.1, and dbg=16
    rn7=experiment('rn7') 
    rn7.vary('nx',200,200,1,0)
    rn7.vary('minSigSt',8,8,1,0)
    rn7.vary('diskScaleLength',2,2,1,0)
    rn7.vary('dbg',16,16,1,0)
    rn7.vary('ndecay',-1,-1,1,0)
    rn7.vary('alphaMRI',0.1,0.1,1,0)
    rn7.vary('whichAccretionHistory',600,799,200,0)
    rn7.vary('Mh0',1.0e10,1.0e12,3,1)
    allModels['rn7']=rn7

    # similar to rk70, except with only ndecay=-1, alpha=.1, and dbg=16
    rn8=experiment('rn8') 
    rn8.vary('nx',200,200,1,0)
    rn8.vary('minSigSt',8,8,1,0)
    rn8.vary('diskScaleLength',2,2,1,0)
    rn8.vary('dbg',16,16,1,0)
    rn8.vary('ndecay',-1,-1,1,0)
    rn8.vary('alphaMRI',0.1,0.1,1,0)
    rn8.vary('whichAccretionHistory',100,599,500,0)
    rn8.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn8']=rn8

    # similar to rk70, except with only ndecay=-1, alpha=.1, and dbg=16
    rn9=experiment('rn9') 
    rn9.vary('nx',200,200,1,0)
    rn9.vary('minSigSt',8,8,1,0)
    rn9.vary('diskScaleLength',2,2,1,0)
    rn9.irregularVary('dbg',[16,32,48,80,112])
    rn9.vary('ndecay',1,8,4,1)
    rn9.vary('alphaMRI',0.1,0.1,1,0)
    rn9.vary('whichAccretionHistory',100,139,40,0)
    rn9.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn9']=rn9

    # similar to rk70, except with only ndecay=-1, alpha=.1, and dbg=16
    rn10=experiment('rn10') 
    rn10.vary('nx',200,200,1,0)
    rn10.vary('minSigSt',8,8,1,0)
    rn10.vary('diskScaleLength',2,2,1,0)
    rn10.irregularVary('dbg',[16,32,48,80,112])
    rn10.vary('ndecay',1,8,4,1)
    rn10.vary('alphaMRI',0.1,0.1,1,0)
    rn10.vary('whichAccretionHistory',100,139,40,0)
    rn10.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn10']=rn10

    rn11=experiment('rn11') # 1st cosmological run in a while. Using ysmoothing w/ ndecay = 8
    rn11.vary('nx',200,200,1,0)
    rn11.vary('minSigSt',8,8,1,0)
    rn11.vary('diskScaleLength',2,2,1,0)
    rn11.vary('dbg',16,16,1,0)
    rn11.vary('alphaMRI',0.1,0.1,1,0)
    rn11.vary('whichAccretionHistory',4,1003,1000,0)
    rn11.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn11']=rn11

    rn12=experiment('rn12') # ndecay = 2
    rn12.vary('nx',200,200,1,0)
    rn12.vary('minSigSt',8,8,1,0)
    rn12.vary('diskScaleLength',2,2,1,0)
    rn12.vary('dbg',16,16,1,0) # smooth y, ndecay=nsmooth
    rn12.vary('ndecay',2,2,1,0)
    rn12.vary('alphaMRI',.1,.1,1,0)
    rn12.vary('whichAccretionHistory',4,1003,1000,0)
    rn12.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn12']=rn12

    rn13=experiment('rn13') # ndecay = 8, hires
    rn13.vary('nx',1000,1000,1,0)
    rn13.vary('minSigSt',8,8,1,0)
    rn13.vary('diskScaleLength',2,2,1,0)
    rn13.vary('dbg',144,144,1,0)
    rn13.vary('alphaMRI',0.1,0.1,1,0)
    rn13.vary('whichAccretionHistory',4,1003,1000,0)
    rn13.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn13']=rn13

    rn14=experiment('rn14') # just do rn12, with lower masses!
    rn14.vary('nx',200,200,1,0)
    rn14.vary('minSigSt',4,4,1,0)
    rn14.vary('diskScaleLength',2,2,1,0)
    rn14.vary('dbg',16,16,1,0) # smooth y, ndecay=nsmooth
    rn14.vary('ndecay',2,2,1,0)
    rn14.vary('alphaMRI',.1,.1,1,0)
    rn14.vary('whichAccretionHistory',4,1003,1000,0)
    rn14.vary('Mh0',1.0e9,1.0e9,1,0)
    allModels['rn14']=rn14

    rn15=experiment('rn15') # just do rn12, with lower masses!
    rn15.vary('nx',200,200,1,0)
    rn15.vary('minSigSt',4,4,1,0)
    rn15.vary('diskScaleLength',2,2,1,0)
    rn15.vary('dbg',16,16,1,0) # smooth y, ndecay=nsmooth
    rn15.vary('ndecay',2,2,1,0)
    rn15.vary('alphaMRI',.1,.1,1,0)
    rn15.vary('whichAccretionHistory',4,1003,1000,0)
    rn15.vary('Mh0',1.0e10,1.0e10,1,0)
    allModels['rn15']=rn15

    rn16=experiment('rn16') # just do rn12, with lower masses!
    rn16.vary('nx',200,200,1,0)
    rn16.vary('minSigSt',4,4,1,0)
    rn16.vary('diskScaleLength',2,2,1,0)
    rn16.vary('dbg',16,16,1,0) # smooth y, ndecay=nsmooth
    rn16.vary('ndecay',2,2,1,0)
    rn16.vary('alphaMRI',.1,.1,1,0)
    rn16.vary('whichAccretionHistory',4,1003,1000,0)
    rn16.vary('Mh0',1.0e11,1.0e11,1,0)
    allModels['rn16']=rn16

    rn17=experiment('rn17') # just do rn12, with lower masses!
    rn17.vary('nx',200,200,1,0)
    rn17.vary('minSigSt',4,4,1,0)
    rn17.vary('diskScaleLength',2,2,1,0)
    rn17.vary('dbg',16,16,1,0) # smooth y, ndecay=nsmooth
    rn17.vary('ndecay',2,2,1,0)
    rn17.vary('alphaMRI',0,0,1,0)
    rn17.vary('whichAccretionHistory',4,1003,1000,0)
    rn17.vary('Mh0',1.0e9,1.0e9,1,0)
    allModels['rn17']=rn17

    rn18=experiment('rn18')
    rn18.vary('nx',200,200,1,0)
    rn18.vary('minSigSt',10,10,1,0)
    rn18.vary('diskScaleLength',2,2,1,0)
    rn18.vary('dbg',16,16,1,0)
    rn18.vary('ndecay',2,2,1,0)
    rn18.vary('alphaMRI',.1,.1,1,0)
    rn18.irregularVary('whichAccretionHistory',[-100,-316,-1000]) #try out oscillation acc history
    rn18.vary('Mh0',1.0e12,1.0e12,1,0) # shoudln't matter
    allModels['rn18']=rn18

    rn19=experiment('rn19')
    rn19.vary('nx',200,200,1,0)
    rn19.vary('minSigSt',10,10,1,0)
    rn19.vary('diskScaleLength',2,2,1,0)
    rn19.vary('epsff',0,.01,2,0)
    rn19.irregularVary('dbg',[16, 2**8+16, 2**9, 2**9+2**8, 2**9+2**10, 2**9+16, 2**10 ])
    rn19.vary('ndecay',2,2,1,0)
    rn19.vary('alphaMRI',.1,.1,1,0)
    rn19.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn19.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn19']=rn19

    rn21=experiment('rn21') # similar to rn19.
    rn21.vary('nx',200,200,1,0)
    rn21.vary('minSigSt',10,10,1,0)
    rn21.vary('diskScaleLength',2,2,1,0)
    rn21.vary('epsff',0,.01,2,0)
    rn21.irregularVary('dbg',[2**10,2**10+2**11])
    rn21.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn21.vary('alphaMRI',.1,.1,1,0)
    rn21.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn21.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn21']=rn21

    rn22=experiment('rn22') # similar to rn19.
    rn22.vary('nx',200,200,1,0)
    rn22.vary('minSigSt',10,10,1,0)
    rn22.vary('diskScaleLength',2,2,1,0)
    rn22.vary('epsff',0,.01,2,0)
    rn22.irregularVary('dbg',[2**10,2**10+2**11])
    rn22.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn22.vary('alphaMRI',0.0,0.0,1,0)
    rn22.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn22.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn22']=rn22

    rn23=experiment('rn23') # similar to rn21.
    rn23.vary('nx',200,200,1,0)
    rn23.vary('minSigSt',10,10,1,0)
    rn23.vary('diskScaleLength',2,2,1,0)
    rn23.vary('epsff',0,.01,2,0)
    rn23.irregularVary('dbg',[2**10,2**10+2**11])
    rn23.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn23.vary('alphaMRI',0.1,0.1,1,0)
    rn23.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn23.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn23']=rn23

    rn24=experiment('rn24') # similar to rn21.
    rn24.vary('nx',200,200,1,0)
    rn24.vary('minSigSt',10,10,1,0)
    rn24.vary('diskScaleLength',2,2,1,0)
    rn24.vary('epsff',0,.01,2,0)
    rn24.irregularVary('dbg',[2**10,2**10+2**11])
    rn24.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn24.vary('alphaMRI',0.1,0.1,1,0)
    rn24.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn24.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn24']=rn24

    rn25=experiment('rn25') # similar to rn21.
    rn25.vary('nx',200,200,1,0)
    rn25.vary('minSigSt',10,10,1,0)
    rn25.vary('diskScaleLength',2,2,1,0)
    rn25.vary('epsff',0,.01,2,0)
    rn25.irregularVary('dbg',[2**10,2**10+2**11])
    rn25.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn25.vary('alphaMRI',0.0,0.1,2,0)
    rn25.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn25.vary('Mh0',1.0e12,1.0e12,1,0)
    allModels['rn25']=rn25

    rn26=experiment('rn26') # same as rn25 but with Mdot_- =0, not tau=0
    rn26.vary('nx',200,200,1,0)
    rn26.vary('minSigSt',10,10,1,0)
    rn26.vary('diskScaleLength',2,2,1,0)
    rn26.vary('epsff',0.01,.01,1,0)
    #  9: only set F=0, not tau.
    # 10: never set y=0, i.e. smooth y.
    # 11: d/dt equations not explicitly mass-conserving
    # 12: turn off if F<0; else turn off if Q>Qf
    # 13: 1-cell padding for KTO
#    rn26.irregularVary('dbg',[2**9+2**10, 2**10, 2**10+2**11, 2**10+2**12, 2**10+2**12+2**13])
#    rn26.irregularVary('dbg',[2**10+2**14, 2**9+2**10, 2**10, 2**9+2**10+2**12, 2**10+2**12, 2**10+2**12+2**13])
    rn26.irregularVary('dbg',[2**10, 2**10+2**12, 2**10+2**12+2**9])
    rn26.vary('ndecay',2,2,1,0) # this line shouldn't matter
    rn26.vary('alphaMRI',0.0,0.05,2,0)
    rn26.irregularVary('whichAccretionHistory',[-1000,-3000,-6000,-9000])
    rn26.vary('Mh0',1.0e12,1.0e12,1,0)
    rn26.irregularVary('TOL',[1.0e-4,1.0e-3])
    allModels['rn26']=rn26



    rn27=experiment('rn27a') # detailed view of MW-size halo
    rn27.vary('nx',200,200,1,0)
    rn27.vary('diskScaleLength',2,2,1,0)
    rn27.vary('dbg',2**10,2**10,1,0)
    rn27.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn27.vary('alphaMRI',0,0,1,0)
    rn27.vary('minSigSt',10,10,1,0)
    rn27.vary('whichAccretionHistory',100,1099,1000,0)
    allModels['rn27']=rn27

    rn28=experiment('rn28a') # detailed view of a lower-mass halo
    rn28.vary('nx',200,200,1,0)
    rn28.vary('diskScaleLength',2,2,1,0)
    rn28.vary('dbg',2**10,2**10,1,0)
    rn28.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn28.vary('alphaMRI',0,0,1,0)
    rn28.vary('minSigSt',5,5,1,0)
    rn28.vary('whichAccretionHistory',1100,2099,1000,0)
    rn28.vary('Mh0',1.0e11,1.0e11,1,0)
    allModels['rn28']=rn28

    rn29=experiment('rn29a') # More accretion!
    rn29.vary('nx',200,200,1,0)
    rn29.vary('diskScaleLength',2,2,1,0)
    rn29.vary('dbg',2**10,2**10,1,0)
    rn29.vary('TOL',1.0e-3,1.0e-3,1,0)
    rn29.vary('alphaMRI',0,0,1,0)
    rn29.vary('minSigSt',10,10,1,0)
    rn29.vary('whichAccretionHistory',2100,3099,1000,0)
    rn29.vary('Mh0',3.0e12,3.0e12,1,0)
    allModels['rn29']=rn29

    rn34=experiment('rn34a') # variable metal diffusion coefficient
    rn34.vary('dbg',2**10+2**15,2**10+2**15,1,0)
    rn34.vary('minSigSt',10,10,1,0)
    rn34.vary('whichAccretionHistory',3100,4099,1000,0)
    allModels['rn34']=rn34

    rn35=experiment('rn35a') # increased accr. history threshold.
    rn35.vary('dbg',2**10+2**16,2**10+2**16,1,0)
    rn35.vary('minSigSt',10,10,1,0)
    rn35.vary('whichAccretionHistory',4100,5099,1000,0)
    allModels['rn35']=rn35

    rn36=experiment('rn36a') # smaller eta
    rn36.vary('dbg',2**10,2**10,1,0)
    rn36.vary('minSigSt',10,10,1,0)
    rn36.vary('whichAccretionHistory',5100,6099,1000,0)
    rn36.vary('eta',.5,.5,1,0)
    allModels['rn36']=rn36
    

    experiments=[]
    nacc = 90
    nmh0 = 61
    compositeNames=[]
    for i in range(nmh0): # 100 different M_h values, each with 100 different accr. histories.
        theName = 'rn33_'+str(i).zfill(3)+'_'
        compositeNames.append(theName)
        experiments.append(experiment(theName))
        experiments[i].vary('nx',200,200,1,0)
        experiments[i].vary('diskScaleLength',2,2,1,0)
        experiments[i].irregularVary('dbg',[2**10])
	experiments[i].irregularVary('TOL',[1.0e-3])
        experiments[i].vary('alphaMRI',0,0,1,0)
        Mh0 = 1.0e9 * (1000.0)**(float(i)/float(nmh0-1))
        experiments[i].irregularVary('minSigSt',[4.0+(float(i)/float(nmh0-1))*(10.0-4.0)])
#	experiments[i].irregularVary('R',[5.0+(float(i)/float(nmh0-1))*(20.0-5.0)])
        experiments[i].irregularVary('R',[20.0 * (Mh0/1.0e12)**(1.0/3.0)])
        experiments[i].irregularVary('Mh0',[Mh0])
        experiments[i].vary('whichAccretionHistory',4+nacc*i,4+nacc*(i+1)-1,nacc,0)
        allModels[theName]=experiments[i]
        
    if ('rn33' in modelList):
        for ex in experiments:
            modelList.append(ex.expName)
        del modelList[modelList.index('rn33')]

    successTables=[]
    # run all the models, and record which ones succeed.
    for model in modelList:
        if(not args.xgrid): #local run
          allModels[model].localRun(args.nproc,args.start)
        else: # write a file to run on the xgrid
          allModels[model].write('runExperiment_'+model+'.txt')
#        successTables.append(allModels[model].ExamineExperiment())

    if ('rn33' in modelList):
        os.mkdir('analysis/rn33')
        for dirname in compositeNames:
            files=glob.glob('./analysis/'+dirname+'/*')
            for aFile in files:
	        os.symlink(aFile,'./analysis/rn33/')
           

#    PrintSuccessTables(successTables)


