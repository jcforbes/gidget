import copy
import os
import sys
import shutil
import subprocess
import time
import math

# This is a somewhat straightforward script to allow you to run the GIDGET code with
# a wide variety of systematically varied parameters. The functions here are described
# in some detail, but if you skip to the end of the script, there are plenty of 
# examples of what you can run from the command line, e.g.
#  $ python exper.py 'ex1'
# will execute the example used in the README file.
class experiment:
    def __init__(self,name):
        # fiducial model
        self.p=[name,500,1.5,.01,2.0,1,1,.01,1,10,220.0,20.0,7000.0,2.5,.5,2.0,2.0,50.0,int(1e9),1e-4,1.0,0,-1.0,0,0.0,1.5,1,2.0,.001]
        self.pl=[self.p[:]]
        # store some keys and the position to which they correspond in the p array
        self.names=['name','nx','eta','epsff','tauHeat','analyticQ','cosmologyOn','xmin','NActive','NPassive','vphiR','R','gasTemp','Qlim','fg0','phi0','zstart','tmax','stepmax','TOL','mu','b','diskScaleLength','whichAccretionHistory','alphaMRI','thickness','migratePassive','fixedQ','kappaMetals']
        self.keys={}
        ctr=0
        for n in self.names:
            self.keys[n]=ctr
            ctr=ctr+1
        self.expName=name
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
            self.p[self.keys[name]]=[type(self.p[self.keys[name]])(float(mmin) * (float(mmax)/float(mmin))**(float(i)/float(n-1))) for i in range(n)]
        elif(n==1):
            self.p[self.keys[name]]=mmin # a mechanism for changing the parameter from its default value.
            
    def generatePl(self):
        '''vary() will change a certain element of self.p into a list
        of values for that variable. This method assumes that some
        (perhaps many) variables have been varied in this manner, and
        expands p into a list of well-defined parameter lists, each
        one of which corresponds to a run of the program. This list can
        then be sent to the xgrid or run on your machine. See other
        methods for details on how to do that.'''
        for j in range(len(self.p)):
            param = self.p[j]
            if(type(param)==list):
                # make len(param) copies of the current pl
                pl2=[]
                for i in range(len(param)):
                    f=copy.deepcopy(self.pl)
                    pl2.append(f) # pl2 ~ [pl,pl,pl]~[[p,p,p],[p,p,p],[p,p,p]]
                    # in each copy, set the jth parameter in p to param[i]
                    for a_p in pl2[i]:
                        a_p[j]=param[i]
                        # in each copy, append to the name a...z corresponding to which copy is currently being edited
                        if(i<=25):
                            base='a'
                            app=''
                        else:
                            base='a'
                            app=str((i-(i%26))/26)
                        a_p[self.keys['name']]+=chr(ord(base)+(i%26))+app
                    # end loop over p's in pl2
                # end loop over range of this varied parameter
                self.pl=[element for sub in pl2 for
 element in sub] # collapse pl to be a 1d array of p's.
            else : #just a regular non-varying parameter
                # If the user has used vary() but with N=1, the new argument will
                # not be a list!
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


    def localRun(self,nproc):
        procs=[]
        ctr=0
        for a_p in self.pl:
            ctr+=1
            tmpap=a_p[:]
            binary=self.bin+'/gidget'
            expDir=self.analysis+'/'+self.expName #directory for output files
            if(os.path.exists(expDir)):
                pass # no need to do anything
            else:
                os.mkdir(expDir)
            with open(expDir+'/'+a_p[self.keys['name']]+'_stdo.txt','w') as stdo:
                with open(expDir+'/'+a_p[self.keys['name']]+'_stde.txt','w') as stde:
                    print "Sending run #",ctr,"/",len(self.pl)," , ",tmpap[0]," to a local core."
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
                print "Still waiting for ",nPrev, " processes to finish."


if __name__ == "__main__":
    expName=sys.argv[1]
    a=experiment(expName)

    # Examples of what one might be interested in co-varying


    # ri29- do only 1000 runs, i.e. similar to ri28, but low resolution. ri30 is high res, 
    # ri31- low res w/ fixed mdot normalization, ri33, higher res, smaller sample, ri34- ri31 w/ poly interp.
    ### vary the mass
#    a.vary('epsff',.01,.01,1,0) # not varying
#    a.vary('nx',200,200,1,0) # not varying, lowres
#    a.vary('Qlim',2.5,2.5,1,0) # not varying
#    a.vary('mu',1.0,1.0,1,0) # not varying
#    a.vary('phi0',2.0,2.0,1,0) # not varying
#    a.vary('whichAccretionHistory',4,3003,3000,0)
#    a.vary('thickness',1.5,1.5,1,0) # not varying
#    a.vary('fixedQ',2.0,2.0,1,0) # not varying
#    a.vary('kappaMetals',1.0e-3,1.0e-3,1,0) # not varying


#    ri32 - production run with t_ff replaced by its exact value (1/sqrt(G rho) -> sqrt( 3 pi / (32 G rho))
#    a.vary('epsff',.01,.01,1,0) # not varying
#    a.vary('Qlim',2.5,2.5,1,0) #not varying
#    a.vary('mu',1.0,1.0,1,0) # not varying
#    a.vary('phi0',2.0,2.0,2,0) # not varying
#    a.vary('whichAccretionHistory',0,2,3,0)
#    a.vary('migratePassive',0,1,2,0) 
#    a.vary('thickness',1.5,1.5,1,0) # not varying
#    a.vary('fixedQ',2.0,2.0,1,0) #not varying
#    a.vary('kappaMetals',1.0e-4,1.0e-2,9,1)


#   ri35 - back to simplicity: vary fg0 and phi0
#    a.vary('fg0',.1,.9,9,0)
#    a.vary('phi0',1.0,2.0,2,0)

#    ri36- set a fixed turnover radius / bulge radius, vary the accretion histories, see if a feature appears
    a.vary('nx',200,200,1,0)
    a.vary('b',1.5,1.5,1,0)
    a.vary('whichAccretionHistory',4,3003,3000,0)

    # expand all the vary-ing into the appropriate number of 
    # parameter lists for individual runs.
    a.generatePl()  
#    a.write('runExperiment_'+expName+'.txt') # to use with an xgrid
    a.localRun(16) # run (in serial) on four processors.
