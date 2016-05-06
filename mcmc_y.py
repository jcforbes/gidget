import shutil
import pdb
import sys
import copy
import emcee
import pickle
import triangle
from mpi4py import MPI 
from emcee.utils import MPIPool
import numpy as np
import exper
import readoutput
import os, glob
import matplotlib
import verticalProfile
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.stats
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


chaindirrel = 'mcmcIndFromMax25'
analysisdir = os.environ['GIDGETDIR']+'/analysis/'
chaindir = analysisdir+chaindirrel
bolshoidir = '~/bolshoi/'

proccounter=0
runnumber = 0

class bolshoireader:
    def __init__(self,fn, minmass, maxmass):
        self.filename = fn
        trees=[]
        with open(bolshoidir+fn,'r') as f:
            line = f.readline()
            while line!='':
                sline = line.split()
                treenum = int(sline[0])
                mh0 = float(sline[1])
                weightedmdot = float(sline[2])
                if minmass<mh0 and mh0<maxmass:
                    trees.append([treenum, mh0, weightedmdot])
                line=f.readline()
        self.trees = sorted(trees, key=lambda tree: tree[2])
        self.keys = [tree[2] for tree in self.trees]
    def copynearest(self, dirname, filename, value):
        ''' copy a file with the closest value of tree[2]'''
        #ind = np.searchsorted(self.keys, value)
        ind = int(value)
        if ind == len(self.keys):
            ind=ind-1
        assert ind<len(self.keys) and ind>-1
        #rffile = bolshoidir+ repr(self.trees[ind][0])+'_rf.txt'
        #shutil.copy(rffile, dirname)
        # This chunk of code searches the file rf_data.txt for 
        #  the appropriate tree number from the register as identified by ind above.
        # Once found, the subsequent lines (containing the values of the random factors
        #  for the main progenitor accretion history of that tree) are copied to a file
        #  in the working directory of the gidget run.
        with open(bolshoidir+'rf_data.txt','r') as f:
            line = f.readline()
            while line!='':
                if 'tree '+repr(self.trees[ind][0]) in line:
                    break # we found the tree we're looking for
                line = f.readline()
            if line=='':
                print "Failed to find line tree "+repr(self.trees[ind][0])
            line = f.readline()
            collectlines=[]
            while not 'tree' in line and line!='':
                collectlines.append(line)
                line = f.readline()
            with open(dirname+'/'+filename,'w') as ff:
                assert len(collectlines)>0
                for line in collectlines:
                    ff.write(line)

globalBolshoiReader = bolshoireader('rf_registry.txt',3.0e11,3.0e12)
bolshoiSize =  len(globalBolshoiReader.keys)



def lnbetadensity(theta, a,b):
    if theta<0 or theta>1:
        return -np.inf
    return (a-1.0)*np.log(theta) + (b-1.0)*np.log(1.0-theta)
def lngammadensity(theta, a,b):
    if theta<0:
        return -np.inf
    return (a-1.0)*np.log(theta) - b*theta
def lnlognormaldensity(theta, mean, var):
    if theta<=0:
        return -np.inf
    return -np.log(theta) - 0.5*(np.log(theta) - mean)**2.0/var
def lnnormaldensity(theta, mean, var):
    return -0.5*(theta-mean)**2.0/var

def lnloguniformdensity(theta, a, b):
    if theta<a or theta>b:
        return -np.inf
    return -np.log(theta)

def lnuniformdensity(theta, a, b):
    if theta<a or theta>b:
        return -np.inf
    return 0.0 

def samplefrombetadensity(a,b):
    assert a>0 and b>0
    return np.random.beta(a,b)
def samplefromgammadensity(a,b):
    assert a>0 and b>0
    return np.random.gamma(a,1.0/b) # numpy's definition of the gamma uses e^{-x/scale}/scale^a
def samplefromlognormaldensity(mean,var):
    return np.exp(samplefromnormaldensity(mean,var))
def samplefromnormaldensity(mean,var):
    return np.random.normal(mean,np.sqrt(var))
def samplefromloguniformdensity(a,b):
    return a*(b/a)**np.random.uniform()
def samplefromuniformdensity(a,b):
    return np.random.uniform(a,b)

# define the base experiment we want.
def emceeparameterspacetogidgetexperiment(emceeparams,name=None):
    global proccounter
    global runnumber


    # unpack emceeparams
    eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc, systematicVphiError, bolshoiWeight = emceeparams

    # if any of the following assertions fail, you should probably adjust / check the prior
    assert fg0>=0
    assert fg0<=1.0
    assert munorm>=0
    assert accscalelength>=0
    assert 0<=fcool and fcool<=1

    # create experiment
    basename = chaindirrel+'_'
    
    if name is None:
        name = basename+str(runnumber).zfill(3)+'_'+str(rank).zfill(5)+'_'+str(proccounter).zfill(5)
    thisexper = exper.experiment(copy.copy(name))

    # we need to put the random factors for the accretion history into a file for gidget to read.
    directory = chaindir[0:-len(chaindirrel)]+'/'+name+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    globalBolshoiReader.copynearest(directory, name+'_inputRandomFactors.txt', bolshoiWeight)
    #np.savetxt(directory+name+'_inputrandomfactors.txt', [x0,x1,0.0]  ) # fscatter is accounted for in gidget proper.
    yscatter = 0.5 # check!
    #np.savetxt(directory+name+'_inputrandomfactorsy.txt', [y0*yscatter,y1*yscatter]  )

    proccounter+=1

    # set up a range of masses and computational domains
    thisexper.irregularVary('Mh0', mh0)
    #thisexper.irregularVary('r', 50*accscalelength/.05 *(mh0/1.0e12)**(1.0/3.0) )
    thisexper.irregularVary('R', 80.0)
    thisexper.irregularVary('xmin', .005)
    thisexper.irregularVary('TOL', 1.0e-3)
    thisexper.irregularVary('NChanges', 301)

    # set up some common parameters.
    thisexper.irregularVary('dbg', 2**4 + 2**1 + 2**0)
    thisexper.irregularVary('alphaMRI', 0.005)
    thisexper.irregularVary('zstart',3.99)
    thisexper.irregularVary('zrelax',4.0)
    thisexper.irregularVary('NPassive',2)
    thisexper.irregularVary('Noutputs',12)


    thisexper.irregularVary('yREC', .054*.54)
    thisexper.irregularVary('eta',eta)
    thisexper.irregularVary('epsff',epsff)
    thisexper.irregularVary('fg0',fg0)
    thisexper.irregularVary('muNorm',munorm)
    thisexper.irregularVary('muColScaling',mucolscaling)
    thisexper.irregularVary('muFgScaling',mufgscaling)
    thisexper.irregularVary('fixedQ',fixedq)
    thisexper.irregularVary('Qlim',qlim)
    thisexper.irregularVary('accScaleLength',accscalelength)
    thisexper.irregularVary('fscatter',0.45)
    thisexper.irregularVary('accNorm',epsilonacc)
    #thisexper.irregularVary('accalphaz',accalphaz)
    #thisexper.irregularVary('accalphamh',accalphamh)
    #thisexper.irregularVary('accceiling',accceiling)
    thisexper.irregularVary('fcool',fcool)
    thisexper.irregularVary('ZIGM',zigm)
    thisexper.irregularVary('concentrationRandomFactor',conrf)

    thisexper.irregularVary('whichAccretionHistory', -344)


    # r0 and v0 are irrelevant to running the physical model and only matter in the comparison to data.
    return thisexper, name

# define the base experiment we want.
def emceeparameterspacetogidgetexperimentrerun(emceeparams,name=None):
    global proccounter
    global runnumber


    # unpack emceeparams
    # eta, epsff, fg0, munorm, mumhscaling, fixedq, accscalelength, xirec, fcool, kappametals, zigm, mh0, alphamri, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsscale, conrf, muhgscaling = emceeparams
    eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc, systematicVphiError, bolshoiWeight = emceeparams


    # create experiment
    basename = chaindirrel+'-rerun-ppd_'
    
    if name is None:
        name = basename+str(rank).zfill(5)+'_'+str(proccounter).zfill(5)


    # if any of the following assertions fail, you should probably adjust / check the prior
    try:
        assert fg0>=0
        assert fg0<=1.0
        assert munorm>=0
        assert accscalelength>=0
        assert 0<=fcool and fcool<=1
    except:
        print "assertions failed. ", emceeparams, name
        assert False


    thisexper = exper.experiment(copy.copy(name))

    # we need to put the random factors for the accretion history into a file for gidget to read.
    directory = chaindir[0:-len(chaindirrel)]+'/'+name+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    #np.savetxt(directory+name+'_inputrandomfactors.txt', [x0,x1, 0.0]  )
    globalBolshoiReader.copynearest(directory, name+'_inputRandomFactors.txt', bolshoiWeight)
    yscatter=0.5
    #np.savetxt(directory+name+'_inputrandomfactorsy.txt', [y0*yscatter,y1*yscatter]  )

    proccounter+=1

    # set up a range of masses and computational domains
    thisexper.irregularVary('Mh0', mh0)
    #thisexper.irregularVary('r', 50*accscalelength/.05 *(mh0/1.0e12)**(1.0/3.0) )
    thisexper.irregularVary('R', 80.0)
    thisexper.irregularVary('xmin', .005)
    thisexper.irregularVary('TOL', 1.0e-3)
    thisexper.irregularVary('NChanges', 301)

    # set up some common parameters.
    thisexper.irregularVary('dbg', 2**4+2**1+2**0)
    thisexper.irregularVary('alphaMRI', 0.005)
    thisexper.irregularVary('zstart',3.99)
    thisexper.irregularVary('zrelax',4.0)
    thisexper.irregularVary('NPassive',20)
    thisexper.irregularVary('Noutputs',120)
    thisexper.irregularVary('fscatter',1.0)
    thisexper.irregularVary('yREC', .054*.54)

    thisexper.irregularVary('eta',eta)
    thisexper.irregularVary('epsff',epsff)
    thisexper.irregularVary('fg0',fg0)
    thisexper.irregularVary('muNorm',munorm)
    thisexper.irregularVary('muColScaling',mucolscaling)
    thisexper.irregularVary('muFgScaling',mufgscaling)
    thisexper.irregularVary('fixedQ',fixedq)
    thisexper.irregularVary('accScaleLength',accscalelength)
    thisexper.irregularVary('Qlim',qlim)
    thisexper.irregularVary('fcool',fcool)
    thisexper.irregularVary('ZIGM',zigm)
    thisexper.irregularVary('concentrationRandomFactor',conrf)
    thisexper.irregularVary('accNorm',epsilonacc)

    thisexper.irregularVary('whichAccretionHistory', -344)
    return thisexper, name


def lnprior(emceeparams):
    eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc, systematicvphierror, bolshoiWeight = emceeparams
    accum = 0.0
    accum += lnlognormaldensity(eta, np.log(1.5), np.log(3)**2.0)
    accum += lnlognormaldensity(epsff, np.log(0.01), np.log(3)**2.0)
    accum += lnbetadensity(fg0, 1.0, 1.0)
    accum += lnlognormaldensity(munorm, np.log10(13.0), np.log10(3.0)**2.0)
    accum += lnnormaldensity(mucolscaling, -1.15, 0.5**2.0)
    accum += lnnormaldensity(mufgscaling, 0.16, 0.5**2.0)
    accum += lnlognormaldensity(fixedq, np.log(2.0), np.log(3.0)**2.0)
    accum += lnlognormaldensity(qlim, np.log(2.0), np.log(3.0)**2.0)
    accum += lnlognormaldensity(accscalelength, np.log(.05), np.log(10)**2.0)
    accum += lnbetadensity(fcool, 1.0, 1.0 )
    accum += lnlognormaldensity(zigm, np.log(.0002), np.log(10.0)**2.0 )
    accum += lnlognormaldensity(mh0, np.log(1.0e12), np.log(2.0)**2.0 )
    accum += lnnormaldensity(conrf, 0, 0.5**2.0) # dex.
    accum += lnloguniformdensity(systematicvphierror, 0.5, 30.0)
    accum += lnuniformdensity(bolshoiWeight,0, bolshoiSize) 
    corr = 0.5
    infl = 1.0
    accum += np.log(scipy.stats.multivariate_normal.pdf( [r0mcmc, v0mcmc], [ 8.1, 246], [[0.6*0.6*infl, 0.6*24 * corr*infl], [0.6*24*corr*infl, 24*24*infl]] ))
    accum += lnbetadensity(epsilonacc, 1.0, 1.0)
    if not np.isfinite(accum):
        return -np.inf
    return accum

def samplefromprior():
    # eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, x0, x1, obsscale, conrf, mufgscaling, y0,y1, zigm, r0mcmc, v0mcmc, epsilonacc = emceeparams
    r=np.sqrt(np.random.uniform())
    phi=np.random.uniform()*2.0*np.pi
    corr = 0.5
    infl = 1.0
    rr,vv = np.random.multivariate_normal( [ 8.1, 246], [[0.6*0.6*infl, 0.6*24 * corr*infl], [0.6*24*corr*infl, 24*24*infl]] )

    return [ \
            samplefromlognormaldensity( np.log(1.5), np.log(3.0)**2.0), # eta
            samplefromlognormaldensity( np.log(0.01), np.log(3.0)**2.0), # epsff
            samplefrombetadensity(1.0,1.0), # fg0
            samplefromnormaldensity(np.log10(13.0), np.log10(3.0)**2.0), # munorm
            samplefromnormaldensity(-1.15, .5**2.0), # mucolscaling
            samplefromlognormaldensity( np.log(2.0), np.log(3.0)**2.0), # fixedq
            samplefromlognormaldensity( np.log(2.0), np.log(3.0)**2.0), # qlim
            samplefromlognormaldensity(np.log(.05),np.log(10.0)**2.0), # accscalelength
            samplefrombetadensity( 1.0, 1.0 ), #fcool
            samplefromlognormaldensity( np.log(1.0e12), np.log(2.0)**2.0), # mh0
            samplefromnormaldensity(0,0.5),  # concentrationrandomfactor
            samplefromnormaldensity(0.16, 0.5**2.0), # mufgscaling
            samplefromlognormaldensity(np.log10(0.0002), np.log10(10.0)**2.0), # zigm
            rr, # r0mcmc
            vv, # v0mcmc
            samplefrombetadensity(1.0,1.0),  # epsilonacc
            samplefromloguniformdensity(.5,30.0),  # systematicVphiError
            samplefromuniformdensity(0,bolshoiSize) ] # bolshoiWeight


def constlikelihood(emceeparams, modelname=None):
    ''' just return a constant so that posterior=prior. for testing! '''
    assert False
    return 0

def lnlikelihood(emceeparams, modelname=None):
    # set up the experiment

    time0 = time.time()
    eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc, systematicVphiError, bolshoiWeight = emceeparams

    ## if we're being given a model that's already been run, don't run it again.
    if modelname is None:
        expertorun, name = emceeparameterspacetogidgetexperiment(emceeparams)

        # run the experiment.
        expertorun.localRun(1,0,maxTime=3600*3)
    else:
        name = modelname

    # read the results of the model
    output = readoutput.Experiment(name)
    # ... but only keep the radial functions to which we will compare real data.
    output.read(keepOnly=['vPhi','colst','colH2','colHI'])


    if len(output.models)==0:
        print "warning: model did not return sensible results, setting likelihood to zero"
        return -np.inf

    model0 = output.models[0]
    zs = model0.var['z'].sensible()

    accum = 0


    # in this next section we take the rotation curve from Bhattacharjee (2014) and 
    # compare it to our model. this involves several steps. 
    # next we interpolate the gidget model's circular velocity curve onto the radii from b14.
    # at each radius, we take the data to be drawn from a normal distribution with b14's quoted errors.
    # we assume each measurement is independent, i.e. the likelihood of each point is simply multiplied.
    # to begin to account for possible systematics in the data, we allow the quoted errors to be
    # scaled by an additional nuissance parameter obsscale.

    #rotcurve = np.loadtxt(chaindir+"/../../Bhattacharjee2014.txt", skiprows=15)
    rotcurve = np.loadtxt(chaindir[0:-len(chaindirrel)]+"/../Bhattacharjee2014.txt", skiprows=15)
    #whichrotcurve = convertpairtolabel(l0,l1)
    whichrotcurve = 1
    if whichrotcurve==0:
        rc = rotcurve[0:51, 2:5]
    if whichrotcurve==1:
        rc = rotcurve[51:101, 2:5]
    if whichrotcurve==2:
        rc = rotcurve[101:151, 2:5]

    rmodel = model0.var['r'].sensible(timeIndex=-1)




    rdatah2 = np.array( [.5, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, \
                7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75])
    colh2data =np.array( [26.7, 6.548974943052395, 4.498861047835993, 2.1526195899772205, 2.7220956719817764, \
            4.316628701594535, 4.316628701594535, 11.241457858769932, 12.220956719817767, \
            10.19362186788155, 10.125284738041003, 10.444191343963555, 7.551252847380409, \
            4.088838268792713, 4.7038724373576315, 3.6332574031890665, 3.5193621867881575, \
            2.5626423690205016, 3.041002277904326, 2.927107061503417, 1.9476082004555835, \
            2.175398633257405, 1.6059225512528492, 1.287015945330296, 2.1070615034168583, \
            1.150341685649206, 1.7425968109339394, 1.6514806378132114, 2.835990888382689] )
    colh2error = np.array( [20.3+26.7, 9.09185803757829, 6.628392484342381, 4.1022964509394555, 4.77035490605428, \
            6.816283924843425, 7.0459290187891455, 14.749478079331944, 16.29436325678497, \
            14.373695198329855, 14.436325678496871, 14.937369519832988, 11.784968684759917, \
            7.839248434237998, 8.298538622129438, 7.3382045929018815, 8.653444676409187, \
            6.9624217118997915, 7.171189979123174, 7.171189979123174, 6.064718162839251, \
            6.127348643006265, 4.937369519832988, 4.54070981210856, 5.041753653444678, \
            3.810020876826723, 4.039665970772443, 3.6638830897703585, 4.979123173277662 ] ) - colh2data


    rdatahi = np.array([.75+i*.5 for i in range(39)])
    colhidata = np.array([ 2.175398633257405, 1.6514806378132114, 1.5375854214123024, 1.5831435079726646, \
            1.7198177676537618, 2.1070615034168583, 2.357630979498861, 2.9498861047836016, \
            3.3826879271070602, 3.610478359908882, 3.7243735763097945, 3.7927107061503413, \
            3.5193621867881575, 3.2915717539863323, 3.0865603644646917, 3.1093394077448764, \
            4.56719817767654, 3.906605922551254, 3.5876993166287043, 3.5193621867881575, \
            3.6332574031890665, 3.4738041002277917, 2.995444191343964, 2.927107061503417, \
            2.6082004555808673, 2.3120728929384953, 1.9476082004555835, 1.6287015945330303, \
            1.5148063781321177, 1.3097949886104807, 1.3781321184510276, 1.4009111617312051, \
            1.1958997722095681, 1.3097949886104807, 1.1731207289293835, 1.1958997722095681, \
            1.3781321184510276, 1.264236902050115, 1.287015945330296 ])
    colhierror = colhidata* 0.2



    #rdatah2 = np.arange(11)+0.5  # quoted bin centers
    vcmodelath2data = model0.var['vPhi'].atR( rdatah2, rmodel, -1)
    vcmodelathidata = model0.var['vPhi'].atR( rdatahi, rmodel, -1)
    v0assumed = 220 # double check this value
    r0assumed = 8.5 # double check this value
    rdatah2 = r0mcmc / ( (v0mcmc/ vcmodelath2data) * ( ( v0assumed/v0mcmc) * (1.0 * (r0assumed/rdatah2) - 1.0) + 1.0) ) # adjust for model's circular velocity and this mcmc iteration's values of v0 and r0.
    rdatahi = r0mcmc / ( (v0mcmc/ vcmodelathidata) * ( ( v0assumed/v0mcmc) * (1.0 * (r0assumed/rdatahi) - 1.0) + 1.0) ) # adjust for model's circular velocity and this mcmc iteration's values of v0 and r0.
    #colh2data = [26.7, 6.0, 3.5, 3.5, 4.6, 3.9, 2.8, 1.9, 0.9, 0.5, 0.3]
    #colh2error= [20.3, 5.4, 1.8, 1.6, 2.0, 2.0, 1.4, 1.0, 0.7, 0.4, 0.2]
    nr = np.shape(rc)[0]
    rdata = rc[:,0]
    vcmodel = model0.var['vPhi'].atR(rdata, rmodel, -1)
    colh2model = model0.var['colH2'].atR(rdatah2, rmodel, -1)
    colhimodel = model0.var['colHI'].atR(rdatahi, rmodel, -1)
    minr = np.min(rmodel)
    maxr = np.max(rmodel)

    for i in range(len(rdatah2)):
        r = rdatah2[i]
        if r<minr or r>maxr:
            pass
        else:
            accum += -0.5*(colh2model[i] - colh2data[i])**2.0/colh2error[i]**2.0

        r = rdatahi[i]
        if r<minr or r>maxr:
            pass
        else:
            accum += -0.5*(colhimodel[i] - colhidata[i])**2.0/colhierror[i]**2.0
    
    for i in range(nr):
        r = rc[i,0]
        if r<minr or r>maxr:
            # if the observational data lie outside the computational domain don't count them.
            # note that this would be bad if the computational domain depended on the model parameters!
            # but we don't do that.
            pass 
            
        else:
            r0assumed = 8.3
            v0assumed = 244
            vc = (r0assumed/r0mcmc)*( rc[i,1]  + ( v0mcmc - v0assumed)*rc[i,0]/r0assumed)
            dvc = np.sqrt( np.power(np.array(rc[i,2]),2.0) + systematicVphiError**2.0 )
            vcm = vcmodel[i]
            accum += - 0.5*(vc-vcm)**2.0/dvc**2.0

    # next up we compare to some results enumerated in licquia and newman (2014) arxiv 1407.1078
    #r0 = 8.3 #samplefromnormaldensity(8.33, 0.35**2.0) # first adopt a value of r0, the distance of the sun from the galactic center
    r0 = r0mcmc

    
    # this is the boxy & rix (2013) value of the solar neighborhood stellar surface density
    #rinterp = [8.0,8.3,8.5]
    rinterp = [r0mcmc, 20.0]
    sigstmodel = model0.var['colst'].atR(rinterp,rmodel,-1)[0]
    accum += - 0.5*(sigstmodel-38.0)**2.0/(4.0)**2.0

    rscale = model0.var['scaleLength'].sensible(timeIndex=-1)
#    accum +=  - 0.5*(rscale-2.15)**2.0 / (0.14*obsscale)**2.0

    #sfr
    sfr = model0.var['sfr'].sensible(timeIndex=-1)
    accum +=  - 0.5*(sfr-1.65)**2.0 / (0.19)**2.0

    # total stellar mass
    mstar = model0.var['mstar'].sensible(timeIndex=-1)
    accum +=  -0.5*((6.08e10 - mstar)/1.14e10)**2.0

#    # bulge:total ratio
    bt = model0.var['BT'].sensible(timeIndex=-1)
    mean = 0.150 + (0.028 - 0.019)/2.0
    #accum += -0.5*((mean-bt)/0.028/obsscale)**2.0

    maxcolstindex = np.argmax(model0.var['colst'].sensible(timeIndex=-1))
    time1 = time.time()
    
    
    print "with params ",emceeparams," we get bt=",bt," sfr=",sfr,' rscale=',rscale,' mstar=',np.log10(mstar)," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. the maximum of colst is at ",maxcolstindex,", the number of time outputs is ",model0.nt,len(zs),zs[-1]

    return accum        





#def lnlikelihood(emceeparams, modelname=None):
#    # set up the experiment
#    #eta, epsff, fg0, munorm, mumhscaling, fixedq, accscalelength, xirec, fcool, kappametals, zigm, mh0, alphamri, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsscale, conrf, muhgscaling= emceeparams
#    eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, x0, x1, obsscale, conrf, mufgscaling, y0,y1, kappametals, zigm, xirec, r0mcmc, v0mcmc, epsilonacc = emceeparams
#    time0 = time.clock()
#
#    ## if we're being given a model that's already been run, don't run it again.
#    if modelname is None:
#        expertorun, name = emceeparameterspacetogidgetexperiment(emceeparams)
#
#        # run the experiment.
#        expertorun.localRun(1,0,maxTime=3600)
#    else:
#        name = modelname
#
#    # read the results of the model
#    output = readoutput.Experiment(name)
#    # ... but only keep the radial functions to which we will compare real data.
#    output.read(keeponly=['vphi','colst','colh2'])
#
#
#    if len(output.models)==0:
#        print "warning: model did not return sensible results, setting likelihood to zero"
#        return -np.inf
#
#    model0 = output.models[0]
#    zs = model0.var['z'].sensible()
#
#    accum = 0
#
#
#    # in this next section we take the rotation curve from bhattacharjee (2014) and 
#    # compare it to our model. this involves several steps. first, we decide which 
#    # model from b14 to use - this is controlled by l0 and l1, two "nuissance parameters" in our model.
#    # next we interpolate the gidget model's circular velocity curve onto the radii from b14.
#    # at each radius, we take the data to be drawn from a normal distribution with b14's quoted errors.
#    # we assume each measurement is independent, i.e. the likelihood of each point is simply multiplied.
#    # to begin to account for possible systematics in the data, we allow the quoted errors to be
#    # scaled by an additional nuissance parameter obsscale.
#
#    #rotcurve = np.loadtxt(chaindir+"/../../Bhattacharjee2014.txt", skiprows=15)
#    rotcurve = np.loadtxt(chaindir[0:-len(chaindirrel)]+"/../Bhattacharjee2014.txt", skiprows=15)
#    #whichrotcurve = convertpairtolabel(l0,l1)
#    whichrotcurve = 1
#    if whichrotcurve==0:
#        rc = rotcurve[0:51, 2:5]
#    if whichrotcurve==1:
#        rc = rotcurve[51:101, 2:5]
#    if whichrotcurve==2:
#        rc = rotcurve[101:151, 2:5]
#
#
#    # similar data for hi available from: http://adsabs.harvard.edu/abs/2003pasj...55..191n
#    # hii from http://adsabs.harvard.edu/cgi-bin/bib_query?arxiv:astro-ph/0301598
#    # h2 from http://adsabs.harvard.edu/abs/2006pasj...58..847n
#
#
#    rdatah2 = np.arange(11)+0.5
#    colh2data = [26.7, 6.0, 3.5, 3.5, 4.6, 3.9, 2.8, 1.9, 0.9, 0.5, 0.3]
#    colh2error= [20.3, 5.4, 1.8, 1.6, 2.0, 2.0, 1.4, 1.0, 0.7, 0.4, 0.2]
#    nr = np.shape(rc)[0]
#    rmodel = model0.var['r'].sensible(timeIndex=-1)
#    rdata = rc[:,0]
#    vcmodel = model0.var['vphi'].atr(rdata, rmodel, -1)
#    colh2model = model0.var['colh2'].atr(rdatah2, rmodel, -1)
#    minr = np.min(rmodel)
#    maxr = np.max(rmodel)
#
#    for i in range(len(rdatah2)):
#        r = rdatah2[i]
#        if r<minr or r>maxr:
#            pass
#        else:
#            accum += -0.5*(colh2model[i] - colh2data[i])**2.0/colh2error[i]**2.0
#    
#    for i in range(nr):
#        r = rc[i,0]
#        if r<minr or r>maxr:
#            # if the observational data lie outside the computational domain don't count them.
#            # note that this would be bad if the computational domain depended on the model parameters!
#            # but we don't do that.
#            pass 
#            
#        else:
#            vc = rc[i,1]
#            dvc = rc[i,2] * obsscale
#            vcm = vcmodel[i]
#            accum += - 0.5*(vc-vcm)**2.0/dvc**2.0
#
#    # next up we compare to some results enumerated in licquia and newman (2014) arxiv 1407.1078
#    r0 = 8.3 #samplefromnormaldensity(8.33, 0.35**2.0) # first adopt a value of r0, the distance of the sun from the galactic center
#
#    
#    # this is the boxy & rix (2013) value of the solar neighborhood stellar surface density
#    rinterp = [8.0,8.3,8.5]
#    sigstmodel = model0.var['colst'].atr(rinterp,rmodel,-1)[1]
#    accum += - 0.5*(sigstmodel-38.0)**2.0/(4.0)**2.0
#
#    rscale = model0.var['scalelength'].sensible(timeindex=-1)
##    accum +=  - 0.5*(rscale-2.15)**2.0 / (0.14*obsscale)**2.0
#
#    #sfr
#    sfr = model0.var['sfr'].sensible(timeindex=-1)
#    accum +=  - 0.5*(sfr-1.65)**2.0 / (0.19)**2.0
#
#    # total stellar mass
#    mstar = model0.var['mstar'].sensible(timeindex=-1)
#    accum +=  -0.5*((6.08e10 - mstar)/1.14e10)**2.0
#
##    # bulge:total ratio
#    bt = model0.var['bt'].sensible(timeindex=-1)
##    mean = 0.150 + (0.028 - 0.019)/2.0
##    accum += -0.5*((mean-bt)/0.028/obsscale)**2.0
#
#    maxcolstindex = np.argmax(model0.var['colst'].sensible(timeindex=-1))
#    time1 = time.clock()
#    
#    
#    print "with params ",emceeparams," we get bt=",bt," sfr=",sfr,' rscale=',rscale,' mstar=',np.log10(mstar)," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. the maximum of colst is at ",maxcolstindex,", the number of time outputs is ",model0.nt,len(zs),zs[-1]
#
#    return accum        

    

def lnprob(emceeparams):
    pr = lnprior(emceeparams)
    if np.isfinite(pr):
        return lnlikelihood(emceeparams) + pr
        #return constlikelihood(emceeparams) + pr
    return pr


def dic(restart,burnin,nspace):
    ''' compute deviance information criterion'''


    epsffs = restart['chain'][:,burnin::nspace, 1].flatten()
    etas = restart['chain'][:,burnin::nspace, 0].flatten()


    chainkeys=[]
    keyindices={}
    for i in range(len(epsffs)):
        key = ("%.5e" % epsffs[i]) +'_'+ ("%.5e" % etas[i])
        chainkeys.append(key)
        keyindices[key] = i


    names = glob.glob(chaindir+'-rerun-ppd_*')
    dbar = 0.0
    counter = 0
    for name in names:
        shortname = name[ name.find('-rerun-ppd_'): ]
        ## most of these don't matter since they're not explicitly
        ##  used in the likelihood.
        #emceeparams = [0]*25 
        ## the exception is the scaling of the observational errorbars.
        ##  to do this we need to identify the place in the chain 
        ##  where we got this model.
        thisoutput = readoutput.Experiment(shortname)
        thisoutput.read(paramsonly=True)
        if len(thisoutput.models)==0:
            print 'warning: skipping ',shortname
        else:
            model= thisoutput.models[0]
            key = ("%.5e" % model.p['epsff']) +'_'+ ("%.5e" % model.p['eta'])
            index = keyindices[key] # find our index in the chain
            emceeparams = []
            # assemble the emcee parameters.
            for k in range(np.shape(restart['chain'])[2]):
                emceeparams.append(restart['chain'][: , burnin::nspace, k].flatten()[index])
            dbar += -2.0 * lnlikelihood(emceeparams, modelname=shortname)
            counter += 1
        print "current dbar: ", dbar/counter
        print "counter: ", counter
    dbar = dbar/counter
    print "counter = ",counter
    assert np.isfinite(dbar)

    npa = restart['chain'][:,burnin::nspace, :] # nwalkers x niter x nparam
    pdb.set_trace()
    thetabar = np.mean(npa, axis=0)[0]
    assert len(thetabar)==25
    dthetabar = -2.0 * lnlikelihood(thetabar)

    pd = dbar - dthetabar
    return pd+dbar, pd



def rerunposteriorpredictive():
    ''' rerun the posterior predictive distribution. this can be used to e.g. increase the resolution
        spatially or in terms of the age of stellar populations, or vary some parameter systematically.
        the mandatory argument func is a user-provided function that specifies how a model with known
        parameters should be modified and (re) run.'''
    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    output = readoutput.Experiment(chaindirrel+'-ppd') # read in the posterior predictive distribution.
    output.read(paramsonly=True,keepstars=False)
    emcee_params = []
    print "output.models: ",len(output.models)
    # for each model, take the parameters we have read in and construct the corresponding emcee parameters.
    # eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling,  kappametals, zigm, xirec, r0mcmc, v0mcmc, epsilonacc = emceeparams
    for model in output.models:
        eta = model.p['eta']
        epsff = model.p['epsff']
        fg0 = model.p['fg0']
        munorm = model.p['muNorm']
        mucolscaling = model.p['muColScaling']
        fixedq = model.p['fixedq']
        qlim = model.p['Qlim']
        accscalelength = model.p['accScaleLength']
        fcool = model.p['fcool']
        mh0 = model.p['Mh0']
        mufgscaling = model.p['muFgScaling']
        zigm = model.p['ZIGM'] 
        r0mcmc = 8
        v0mcmc = 220
        epsilonacc = model.p['accNorm']

        # we have everything except obsscale, but actually that doesn't matter,
        # since it only affects the model in post-processing, i.e. in comparing to the data,
        # not the running of the model itself. so..... we good!
        thelist = [ eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc ]
        try:
            assert eta>0 and epsff>0 and fg0>0 and fg0<=1 and fixedq>0 and munorm>=0 and fcool>=0 and fcool<=1 and mh0>0
        except:
            print 'unexpected ppd params: ',thelist
        emcee_params.append( copy.deepcopy(thelist) )
    # ok, from here on out, we just need to emulate parts of the run() function to trick emcee into running a single iteration of the algorithm with this ic.


    ndim =  18

    restart = {}
    restart['currentPosition'] = emcee_params
    restart['chain'] = None
    restart['state'] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0



    nwalkers = len(emcee_params) # need one walker per sample from posterior predictive distribution
    print "starting up the ensemble sampler!"
    sampler = emcee.EnsembleSampler(nwalkers, ndim, fakeprob, pool=pool)
    #pos, prob, state = sampler.run_mcmc(restart['currentPosition'], n, rstate0=restart['state'], lnprob0=restart['prob'])
    print "take a step with the ensemble sampler"

    # take a single step with the ensemble sampler.
    print np.shape(restart['currentPosition']), np.shape(np.random.uniform(0,1,nwalkers))
    sampler._get_lnprob(pos = restart['currentPosition'])

    #result = sampler.sample(restart['currentPosition'], iterations=1, lnprob0=None, rstate0=None)

    #pos, prob, state = result
    print "close the pool"

    pool.close()


def fakeprob(params):
    # get the gidget experiment object we want to run.
    print "running with params ",params
    expertorun, name = emceeparameterspacetogidgetexperimentrerun(params)
    print "experiment name: ",name

    # run the experiment.
    expertorun.localRun(1,0,maxTime=3600*6)
    return np.random.uniform(0,1) # value doesn't matter -- all we care about is re-running the experiment.


def getposteriorpredictive(restart, burnin=0, nspace=10):
    ''' we have to find the set of models over some period of time in the chain (after burnin)
        that represents the posterior predictive distribution of models. this is not the same as
        just taking a set of all models run after you think the sampler has converged because
        some (most!) of those models are not accepted! it's also a bit non-trivial because
        when the new model isn't accepted, you need to include the identical model again.'''
    #allruns = glob.glob(chaindir+'*')
    frac=1.0
    output = readoutput.Experiment(chaindirrel)
    output.read(paramsOnly=True)
    modeldict = {}
    for model in output.models:
        key = ("%.5e" % model.p['epsff']) +'_'+ ("%.5e" % model.p['eta'])
        modeldict[key] = model

    # the following block of code ends up working but not being effective because for some reason i changed the number
    # of outputs in between the two runs, so when i can't find a model in 06, it's gone forever.

    # these are all the accscalelengths in the posterior distribution.
    # this is a somewhat subjective decision: the user needs to have picked a
    #  burnin time, and spacing to cut down on autocorrelation, and possibly also
    #  a fraction of models to take on board (if e.g. you don't want the full sample).
    epsffs = restart['chain'][:,burnin::nspace, 1].flatten()
    etas = restart['chain'][:,burnin::nspace, 0].flatten()
    assert len(epsffs)==len(etas)
    print "copying over models as ",len(epsffs)," samples from the ppd."
    for i in range(len(epsffs)):
        if np.random.uniform()<frac:
            # we have decided to take the model, identified by the key below, and copy it from
            # the mcmc chain into our posterior predictive distribution.
            key = ("%.5e" % epsffs[i]) +'_'+ ("%.5e" % etas[i])
            print "using key ", key
            # to do so, we find the model in our dictionary of models
            if key in modeldict.keys():
                model = modeldict[key]
            else:
                print "didn't find the key in modeldict! diagnostics:"
                print "i: ",i
                print "key: ", key
                print "len(keys())", len(modeldict.keys())
                print "examples: ", modeldict.keys()[0]
                print "examples: ", modeldict.keys()[1]
                print "examples: ", modeldict.keys()[2]
                print "chain shape: ", np.shape(restart['chain'])
                print "len(epsffs): ",len(epsffs)
                raise keyerror


            # assign it a new name
            destname = chaindirrel+'-ppd_'+str(i).zfill(5)
            # and copy the original run to a new directory.
            shutil.copytree( model.dirname, analysisdir+'/'+destname )
            print "copied over key ", key
            # for every file in the copied folder, replace its old prefix with its new name.
            for filename in os.listdir(analysisdir+'/'+destname):
                filenamedest = filename[len(model.name):] # strip off the model name, leaving only e.g. _evolution.dat
                filenamedest = destname+filenamedest # replace with the current model name.
                os.rename(analysisdir+'/'+destname+'/'+filename, analysisdir+'/'+destname+'/'+filenamedest)
#            except keyerror:
#                print "warning: skipping selected model because it's not in "+chaindirrel
#                print "failed to find key ",key," in this list of keys "
#                #print sorted(modeldict.keys())


def run(n, p00=None, nwalkers=500):


    fn = chaindirrel+'.pickle'
    ndim = 18

    if p00 is not None:
        p0 = [p00*(1.0+0.01*np.random.randn( ndim )) for i in range(nwalkers)]
    else:
        p0 = [samplefromprior() for i in range(nwalkers)]

    restart = {}
    restart['currentPosition'] = p0
    restart['chain'] = None
    restart['state'] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0

    # read in our past progress unless we've been given a new starting location.
    if p00 is None:
        updaterestart(fn,restart)

    if restart['chain'] is not None:
        # this may save some time if you change something and forget to delete the .pickle file.
        restartedshape = np.shape(restart['chain'])
        print restartedshape, nwalkers, ndim
        assert restartedshape[0] == nwalkers
        assert restartedshape[2] == ndim

    global runnumber
    runnumber = restart['mcmcRunCounter']

    restart['iterationCounter'] += n
    restart['mcmcRunCounter'] += 1

    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
    #pos, prob, state = sampler.run_mcmc(restart['currentPosition'], n, rstate0=restart['state'], lnprob0=restart['prob'])

    for result in sampler.sample(restart['currentPosition'], iterations=n, lnprob0=restart['prob'], rstate0=restart['state']):

        pos, prob, state = result

        restart['acor'] = sampler.acor[:] # autocorr length for each param (ndim)
        restart['accept'] = sampler.acceptance_fraction[:]  # acceptance frac for each walker.
        restart['currentPosition'] = pos # same shape as p0: nwalkers x ndim
        restart['state'] = state # random number generator state
        restart['prob'] = prob # nwalkers x __
        if restart['chain'] is None:
            restart['chain'] = np.expand_dims(sampler.chain[:,0,:],1) # nwalkers x niterations x ndim
            restart['allProbs'] = np.expand_dims(prob,1)  # nwalkers x niterations
        else:
            print np.shape(restart['chain']), np.shape(sampler.chain[:,-1,:]), np.shape(sampler.chain)
            print restart['mcmcRunCounter'], restart['iterationCounter']
            #restart['chain'] = np.concatenate((restart['chain'], sampler.chain[:,-1,:]), axis=1)
            print "dbg1: ",np.shape(restart['chain']), np.shape(np.zeros((nwalkers, 1, ndim))), np.shape(np.expand_dims(pos,1))
            restart['chain'] = np.concatenate((restart['chain'], np.expand_dims(pos, 1)),axis=1)
            restart['allProbs'] = np.concatenate((restart['allProbs'], np.expand_dims(prob, 1)),axis=1)

        
        saverestart(fn,restart)

    pool.close()


def traceplots(chain, fn, burnin=0):
    ndim = np.shape(chain)[2]
    sq = np.sqrt(float(ndim))
    nr = int(np.ceil(sq))
    nc=nr
    while nr*(nc-1)>ndim:
        nc=nc-1

    # plot the trace for every parameter, for a subset of the walkers.
    fig,ax = plt.subplots(nrows=nr,ncols=nc, figsize=(nc*4,nr*4))
    for dim in range(ndim):
        i = np.mod(dim, nr)
        j = ( dim -i )/nr
        for walker in range(np.shape(chain)[0]):
            if np.random.uniform(0,1) < 1.0e-1: # print every tenth-ish
                ax[i,j].plot(chain[walker,burnin:,dim],alpha=.3,ls='--')

    plt.tight_layout()
    plt.savefig(fn+'.png')
    plt.close(fig)
def ksprob(arr1, arr2):
    from scipy.stats import ks_2samp
    return ks_2samp(arr1,arr2)[1]

def probsplots(allprobs, fn, chain, burnin=0):
    nwalker = np.shape(allprobs)[0]
    iters = np.shape(allprobs)[1]
    ndim = np.shape(chain)[2]

    # plot the trace of the probabilities for every walker.
    fig,ax = plt.subplots()
    for walker in range(nwalker):
        ax.plot(allprobs[walker,burnin:],alpha=.3,ls='--')
    plt.savefig(fn+'_probs.png')
    plt.close(fig)
    print "saved "+fn+'_probs.png'

    fig,ax=plt.subplots()
    kss=[]
    kssthetas = np.zeros( (ndim, iters-1) )
    # chain ~ walker x iter x dim
    for iter in range(iters-1):
        kss.append( ksprob(allprobs[:,iter], allprobs[:,iters-1] ) )
        for k in range(ndim):
            kssthetas[k,iter] = ksprob( chain[:, iter, k], chain[:,iters-1, k])
    ax.scatter(range(iters-1)[-25:], kss[-25:], s=30, marker='+')
    colors = ['k','r','green','orange','lightblue','grey','purple','pink','yellow','blue','lightgreen','darkgreen']*5
    for k in range(ndim):
        ax.plot(range(iters-1)[-25:], kssthetas[k,-25:], color=colors[k])
    ax.set_xlabel('iter')
    ax.set_ylabel('ks prob vs last iteration')
    ax.set_ylim(-0.01,1.02)
    plt.savefig(fn+'_ks.png')
    plt.close(fig)

    



    changes = np.zeros(nwalker)
    for walker in range(nwalker):
        for iter in range(iters-1):
            if allprobs[walker,iter]!=allprobs[walker,iter+1]:
                changes[walker]+=1.0
    changes = changes/float(iters-1.0)
    print "long-term acceptance fraction stats: "
    stats(changes)

    acor = np.zeros(nwalker)
    for walker in range(nwalker):
        acor[walker] = emcee.autocorr.integrated_time(allprobs[walker,burnin:],  window=min([50,iters/2]))
    print "acor stats: "
    stats(acor)

def stats(array):
    print "length: ",len(array)
    print "mean: ",np.mean(array)
    print "gauss percentiles: ",np.percentile(array,[2.5,16,50,84,97.5])
    print "5-spaced percentiles: ",np.percentile(array,range(96)[1::5])
    print "min,max: ",np.min(array), np.max(array)
    print "********************************"

def printrestart(restart):
    ''' print quick summary info about the current state of the sampler. '''
    print "restart info: "
    print " current shape of chain: (nwalkers x niterations x ndim) ",np.shape(restart['chain'])
    print " autocorrelation lengths for each parameter: ",restart['acor']
    stats(restart['acor'])
    print " acceptance rate for each walker: ",restart['accept']
    stats(restart['accept'])

def trianglemovie(restart):
    shp = np.shape(restart['chain'])
    bn = 'triangle_frame_'
    for i in range(shp[1]):
        triangleplot( restart, bn+str(i).zfill(3)+'.png', burnin=i, nspace=10000 )

def triangleplot(restart,fn,burnin=0, nspace=10):
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnin)*shp[2]
    prior = np.array([samplefromprior() for i in range(prs)])
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    #eta, epsff, fg0, munorm, mucolscaling, fixedq, qlim, accscalelength, fcool, mh0, conrf, mufgscaling, zigm, r0mcmc, v0mcmc, epsilonacc, systematicVphiError, bolshoiWeight = emceeparams
    labels = [r"$\eta$",r"$\epsilon_\mathrm{ff}$",r"$f_{g,0}$",r"$\mu_0$", \
            r"$\mu_{\Sigma}$",r"$Q_f$",r"$Q_\mathrm{lim}$",r"$r_\mathrm{acc}/r_\mathrm{vir}$", \
            r"$f_\mathrm{cool}$", r"$M_{h,0}$",  r'c offset (dex)', r'$\mu_{f_g}$', \
            r'$Z_\mathrm{IGM}$', r'$R_0$', r'$v_0$', r'$\epsilon_\mathrm{acc}$', \
            r'$\sigma_{v_\phi,\mathrm{sys}}$', r'Accretion History']

    logs = [0,1,0,1, \
            0,0,0,1, \
            0,1,0,0, \
            1,0,0,0, \
            1,0]

    ### this code selects a fraction frac of the samples after burnin iterations.
    ### the selection is random, meaning consecutive samples from the same walker
    ### will be common. in this case we had a kwarg frac instead of nspace.
    #samplered = restart['chain'][:,burnin:,:].reshape((-1,ndim)) # (nwalkers*(niterations-burnin)) x ndim

    #ns = np.shape(samplered)[0]
    #selection = np.empty((int(ns*frac),) , dtype='int')
    #np.fix(np.random.uniform(0,ns-1,int(ns*frac)), y=selection)
    #samplered = samplered[selection, :]

    ### instead we take a sample not at random, but evenly spaced:
    samplered = copy.deepcopy(restart['chain'][:,burnin::nspace,:].reshape((-1,ndim)))


    extents=[]
    for i in range(np.shape(samplered)[1]):
        if logs[i]==1:
            samplered[:,i] = np.log10(np.clip(samplered[:,i],1.0e-5,np.inf))
            prior[:,i] = np.log10(prior[:,i])
            labels[i] = r'$\log_{10}$ '+labels[i]

        mi = np.min(samplered[:,i])
        ma = np.max(samplered[:,i])
        if(mi==ma):
            extents.append([mi-1,ma+1])
        else:
            extents.append([mi,ma])

    trifig = triangle.corner(samplered, labels=labels, extents=extents)
    #trifigprior = triangle.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, extents=extents)
    trifig.savefig(fn)
    plt.close(trifig)

    # print a summary of the current variables!
    print "**********************************"
    print "name, mean, pr mean, std dev, prior std dev"
    for i in range(np.shape(samplered)[1]):
        print labels[i],"   ",np.mean(samplered[:,i]),"   ", np.mean(prior[:,i]),"   ",np.std(samplered[:,i]), "    ",np.std(prior[:,i])
    print "*********************************"
    print "correlation coefficients"
    for i in range(np.shape(samplered)[1]):
        for j in range(np.shape(samplered)[1]):
            print labels[i], " x ",labels[j]," : ",np.corrcoef(samplered[:,i],samplered[:,j])[0][1]

    # now we're going to follow the formalism of p 284 of bda3 to assess convergence as a fn of iteration.
    wholechain = copy.deepcopy(restart['chain'][:,:,:]) # ~ walker x iter x dim
    m = np.shape(wholechain)[0]*2 # number of walkers times 2 (since we're splitting each chain in 2)
    print "     ",labels
    for finaliter in range(16,np.shape(wholechain)[1],4)[:-1]:
        if finaliter%4 == 0:
            n = finaliter/4  # throw away the first 1/4, then split the second half in two
            psi = np.zeros( (m,n,ndim) ) # j=1...m, i=1...n, k=1...ndim
            psi[0:m/2,: ,:] = wholechain[:,2*n:3*n,:]
            psi[m/2:m,: ,:] = wholechain[:,3*n:4*n,:]
            #psibarj = np.zeros( (m,ndim) ) # average over i=1...n for each j=1..m
            psibarj = np.mean(psi,axis=1) 
            assert np.shape(psibarj)[0] == m and np.shape(psibarj)[1] == ndim
            #psibar = np.zeros( ndim ) # average of psibarj over j=1..m for each variable k
            psibar = np.mean( psibarj, axis=0)
            assert np.shape(psibar)[0] == ndim
            #sj = np.zeros( (m,ndim) ) # variance of psi_ij about psibarj
            sj = np.var( psi, axis=1, ddof=1 )
            assert np.shape(sj)[0] == m and np.shape(sj)[1] ==ndim
            # w = np.zeros(ndim)
            w = np.mean(sj, axis=0)
            b = np.zeros(ndim)
            for k in range(ndim):
                b[k] = float(n)/float(m-1) * np.sum(np.power( psibarj[:,k] - psibar[k], 2.0) )
            vark = float(n-1)/float(n) * w + b/float(n)
            rhat = np.sqrt( vark/w )

            if np.shape(wholechain)[1]>500 and finaliter%100==0:
                print finaliter, "    ", rhat
            elif np.shape(wholechain)[1]<=500:
                print finaliter, "    ", rhat
            else:
                pass



def saverestart(fn,restart):
    with open(fn,'wb') as f:
        print "checkpoint to "+fn
        pickle.dump(restart,f,2)

def updaterestart(fn,restart):
    if os.path.isfile(fn):
        with open(fn,'rb') as f:
            tmp_dict = pickle.load(f)
            restart.update(tmp_dict)

def maximumposteriorprob(x0):
    def tominimize(x):
        xx = copy.deepcopy(x)
        xx[11] = xx[11]*1.0e12
        prob = lnprob(xx)
        print "found ln prob = ",prob
        print "from x = ",x
        if not np.isfinite(prob) or prob<-1.0e30:
            prob=-1.0e30
        return -prob
    result = scipy.optimize.minimize(tominimize, x0)
    print result
    return result

if __name__=="__main__":


    


    ##x0 = [1.5, 0.01, .99, .5, -2./3., 2, .03, .01, .1, 1.0, .002, 1.0, 0.01, .45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    #x0 = sampleFromPrior()
    #x0[11] /= 1.0e12
    #maximumPosteriorProb(x0)

    #eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, conRF, muFgScaling, kappaMetals, ZIGM, xiREC, R0MCMC, V0MCMC, epsilonAcc = emceeParams

    xmax = [ 0.2, 0.01, 0.7, 1.4, # eta, epsff, fg0, muNorm
            -0.43, 1.0, 1.0, .025, # muColScaling, fixedQ, Qlim, accScaleLength
            0.3, 4.5e12,  # fcool, Mh0, 
            0.4, .23,  # conRF, muFgScaling,
            .002, # ZIGM, 
            8, 220, .5, # R0, V0, epsAcc 
            3.0, bolshoiSize/2] # systematicVphiError, bolshoiWeight


#    xmax = [  1.89572200e-01,   1.04076859e-02,   7.38152246e-01,   1.37673100e+00,
#             -4.37140916e-01,   5.58972092e-01,   2.49790419e-02,   2.93190618e-01,
#              4.57221629e+12,   2.33209376e-01,  -1.63343636e-01,  -3.21647816e-01,
#              9.05444175e-02,   1.84943961e-01,   1.11707726e+01,   3.98292819e-01,
#              2.33958073e-01,   7.97679715e-01,   1.83307648e-01,   2.26457747e-01,
#              3.51223250e-01,   3.06985476e-01,   2.26597455e+00,   1.76156619e-03,
#              6.22481329e-04]


#    xmax = [  1.38110984e-01,   1.56054436e-02,   6.77659749e-01,   7.51281733e-01,
#      -1.61484088e-01,   1.52195432e+00,   2.00000000e+00,   2.40937635e-02,
#       1.00530879e-01,   5.65076067e+12,   2.13798472e-01,  -5.05746147e-01,
#       4.36293557e+00,   3.85208913e-01,   2.86150739e-01,   3.54157890e-01,
#       2.19773960e-01,   2.01738761e-03,   
#       8.73335819e+00,   2.68739699e+02,   1.03621353e-01]
#
#
#    xmax = [  1.47974343e-01,   1.25535599e-02,   6.76028939e-01,   8.05514278e-01,
#             -1.10342067e-01,   1.36792855e+00,   1.09758830e+00,   2.90302105e-02,
#              1.90124109e-01,   6.18665396e+12,   2.91009898e-01,  -5.59652735e-01,
#              4.13462858e+00,   3.01332231e-01,   1.44540920e-01,   1.06956484e-01,
#              3.13000089e-01,   1.78235790e-03,  
#              9.06850006e+00,   2.54572880e+02,   1.16864456e-01]
#
#    xmax = [  1.34085076e-01,   1.27581051e-02,   6.72328130e-01,   4.99389544e-01,
#             -8.98714744e-02,   1.75081080e+00,   1.32953768e+00,   4.17966204e-02,
#              2.90496695e-01,   5.87203706e+12,   1.98650145e-01,  -5.10407939e-01,
#              6.60228652e+00,   2.96848022e-01,   1.78408267e-01,   9.12929719e-02,
#              2.24720233e-01,   7.26599888e-04,   9.48338953e+00,   2.75070157e+02,
#              6.60213847e-02]

   
    #run(400, nwalkers=1024,p00=xmax) 
    run(400, nwalkers=1024) 
    #rerunPosteriorPredictive()

    if True:
        # Load in the resulting chain:
        restart={}
        updaterestart(chaindirrel+'.pickle', restart)
        printrestart(restart)
        #dic(restart, 210, 100)
        triangleplot(restart,chaindirrel+'_triangle.png',burnin=67, nspace=300)
        #triangleMovie(restart)
        traceplots(restart['chain'], chaindirrel+'_trace', burnin=0)
        probsplots(restart['allProbs'], chaindirrel+'_allprob', restart['chain'], burnin=0)

        # Find the maximum among all models sampled so far.
        allprobs = restart['allProbs'].flatten()
        index = np.argmax(allprobs)
        print "probs stats: ", np.max(allprobs), np.percentile(allprobs,[5,50,90,95,97.5,99])
        print allprobs[index]
        indices = np.unravel_index(index, np.shape(restart['allProbs']))
        xmax = restart['chain'][indices[0],indices[1],:]
        print "Favorite coordinates: ", xmax

        #getposteriorpredictive( restart, burnin=67, nspace=100)




        # Now use xmax as p00 in a new chain
        #shutil.move(chainDirRel+'.pickle', chainDirRel+'_initial.pickle')

        #run(20, nwalkers=1024, p00=xmax) # now do the real chain initialized from a hopefully-reasonable location.


        ##pdb.set_trace()

        #xmax[11]  /= 1.0e12
        
        #print "Starting a new maximization procedure at this location!"
        #maximumPosteriorProb(xmax)


