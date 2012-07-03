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

def successCode(filename):
    ''' Given a file assumed to contain the standard error output of
    a run of the gidget code, return 0 for success (file is empty)
    1 for a time step below floor result, 2 for a problem with Qst,
    or 3 for a miscellaneous error code.  '''
    if (not os.path.exists(filename)):
      return -1

    if os.stat(filename).st_size == 0:
      return 0

    with open(filename,'r') as f:
        contents = f.read()
        if 'floor' in contents:
            return 1
        if 'Qst' in contents:
            return 2
    return 3
    

class experiment:
    def __init__(self,name):
        # fiducial model
        self.p=[name,500,1.5,.01,4.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,2.0,2.0,50.0,int(1e9),1.0e-4,1.0,0,.5,2.0,-1.0,0,0.0,1.5,1,2.0,.001,1.0e12,1.0,3.0,0]
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
	successTable = np.zeros(tuple(Nvaried[:]))
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
	            


    def localRun(self,nproc,startAt):
        ''' Run the specified experiment on this machine,
        using no more than nproc processors, and starting
        at the "startAt"th run in the experiment '''
        self.generatePl()
        procs=[]
        ctr=0
        binary=self.bin+'/gidget'
        expDir=self.analysis+'/'+self.expName #directory for output files
        if(os.path.exists(expDir)):
            print "This directory already contains output. CANCELLING this run!"
            print
            print "This is fine if you've already run the desired models and are"
            print "just interested in looking at the results of ExamineExperiment."
            print "Otherwise, just delete the appropriate directory and run this"
            print "script again."
            return
        else:
            os.mkdir(expDir)

        for a_p in self.pl[startAt:]:
            ctr+=1
            tmpap=a_p[:]
            with open(expDir+'/'+a_p[self.keys['name']]+'_stdo.txt','w') as stdo:
                with open(expDir+'/'+a_p[self.keys['name']]+'_stde.txt','w') as stde:
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
            if(nStillRunning == nPrev):
                # do nothing except wait a little bit
                time.sleep(120)
            else:
                nPrev=nStillRunning
                if(nPrev == 0):
                    break # we're done!
                print "Still waiting for ",nPrev, " processes to finish; I'll check every two minutes for changes."
        print "Local run complete!"



if __name__ == "__main__":

    # http://docs.python.org/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(description='Analyze data and/or run experiments associated with a list of experiment names.')
    parser.add_argument('models', metavar='model', type=str, nargs='+',
                   help='a model to be run / analyzed')
    parser.add_argument('--nproc',type=int,help="maximum number of processors to use (default: 16)",default=16)
    parser.add_argument('--start',metavar='startingModel',type=int,
                   help='The number of the model in the experiment (as ordered by GeneratePl) at which we will start sending experiments to the processors to run. (default: 0)',default=0)
    args = parser.parse_args()
    
    modelList=[]
    for modelName in args.models:
        modelList.append(modelName)

    
    allModels={} # empty dictionary into which we'll place all model definitions
 
#    expName=sys.argv[1]
#    a=experiment(expName)

    rk1=experiment('rk1') #  test high-res y computation.
    rk1.vary('nx',200,200,1,0)
    rk1.vary('diskScaleLength',2.0,2.0,1,0)
    rk1.vary('whichAccretionHistory',4,103,100,0)
    rk1.vary('R',20.0,20.0,1,0)
    rk1.vary('Mh0',1.0e12,1.0e12,1,0)
    rk1.vary('fg0',.5,.5,1,0)
    allModels['rk1']=rk1

    rk2=experiment('rk2')# - find out in what parameter space we're allowed to work in terms of rotation curves
    rk2.vary('nx',200,200,1,0)
    rk2.vary('diskScaleLength',2.0,2.0,1,0)
    rk2.vary('whichAccretionHistory',4,4,1,0)
    rk2.vary('R',20.0,20.0,1,0)
    rk2.vary('Mh0',1.0e12,1.0e12,1,0)
    rk2.vary('b',0.0,5.0,10,0)
    rk2.vary('softening',1.0,4.0,3,1)
    rk2.vary('innerPowerLaw',.1,1,10,0)
    allModels['rk2']=rk2


    rk59=experiment("rk59") # Vary debug and ndecay parameters on one particular run that fails by Qst error.
    rk59.vary('whichAccretionHistory',23,23,1,0)
    rk59.vary('dbg',0,32,33,0)
    rk59.vary('ndecay',-1,7,5,0)
    allModels['rk59']=rk59

    rk60=experiment("rk60") # fails by dt falling below floor
    rk60.vary('whichAccretionHistory',31,31,1,0)
    rk60.vary('dbg',0,32,33,0)
    rk60.vary('ndecay',-1,7,5,0)
    allModels['rk60']=rk60

    rk61=experiment('rk61')  # no failures yet
    rk61.vary('whichAccretionHistory',4,4,1,0)
    rk61.vary('dbg',0,32,33,0)
    rk61.vary('ndecay',-1,7,5,0)
    allModels['rk61']=rk61

    rk62=experiment('rk62') # originally rk58b3d
    rk62.vary('whichAccretionHistory',83,83,1,0)
    rk62.vary('dbg',0,32,33,0)
    rk62.vary('ndecay',-1,7,5,0)
    allModels['rk62']=rk62

    rk63=experiment('rk63') # originally rk58g1d
    rk63.vary('whichAccretionHistory',36,36,1,0)
    rk63.vary('dbg',0,32,33,0)
    rk63.vary('ndecay',-1,7,5,0)
    allModels['rk63']=rk63

    rk64=experiment('rk64') # originally rk58o1d
    rk64.vary('whichAccretionHistory',44,44,1,0)
    rk64.vary('dbg',0,32,33,0)
    rk64.vary('ndecay',-1,7,5,0)
    allModels['rk64']=rk64

    # expand all the vary-ing into the appropriate number of 
    # parameter lists for individual runs.
##    a.generatePl()  
##    a.write('runExperiment_'+expName+'.txt') # to use with an xgrid
#    a.localRun(16,0) # run (in serial) on four processors.
#    a.ExamineExperiment()

    successTables=[]   
    for model in modelList:
        allModels[model].localRun(args.nproc,args.start)
        successTables.append(allModels[model].ExamineExperiment())

#    print len(successTables)
#    print type(successTables[0])

    sumOfSuccessTables = np.zeros(successTables[0].shape)
    for i in range(len(successTables)):
        table = successTables[i]
        sumOfSuccessTables += (table * (10**i))

    print sumOfSuccessTables

