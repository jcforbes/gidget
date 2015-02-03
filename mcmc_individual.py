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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


chainDirRel = 'mcmcIndFromMax13'
analysisDir = os.environ['GIDGETDIR']+'/analysis/'
chainDir = analysisDir+chainDirRel

procCounter=0
runNumber = 0

def lnBetaDensity(theta, a,b):
    if theta<0 or theta>1:
        return -np.inf
    return (a-1.0)*np.log(theta) + (b-1.0)*np.log(1.0-theta)
def lnGammaDensity(theta, a,b):
    if theta<0:
        return -np.inf
    return (a-1.0)*np.log(theta) - b*theta
def lnLogNormalDensity(theta, mean, var):
    if theta<=0:
        return -np.inf
    return -np.log(theta) - 0.5*(np.log(theta) - mean)**2.0/var
def lnNormalDensity(theta, mean, var):
    return -0.5*(theta-mean)**2.0/var

def lnLogUniformDensity(theta, a, b):
    if theta<a or theta>b:
        return -np.inf
    return -np.log(theta)

def lnUniformDensity(theta, a, b):
    if theta<a or theta>b:
        return -np.inf
    return -np.log(1.0e-5) # for historical reasons.

def sampleFromBetaDensity(a,b):
    assert a>0 and b>0
    return np.random.beta(a,b)
def sampleFromGammaDensity(a,b):
    assert a>0 and b>0
    return np.random.gamma(a,1.0/b) # numpy's definition of the gamma uses e^{-x/scale}/scale^a
def sampleFromLogNormalDensity(mean,var):
    return np.exp(sampleFromNormalDensity(mean,var))
def sampleFromNormalDensity(mean,var):
    return np.random.normal(mean,np.sqrt(var))
def sampleFromLogUniformDensity(a,b):
    return a*(b/a)**np.random.uniform()
def sampleFromUniformDensity(a,b):
    return np.random.uniform(a,b)

# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperiment(emceeParams,name=None):
    global procCounter
    global runNumber


    # unpack emceeParams
    # eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling = emceeParams

    # if any of the following assertions fail, you should probably adjust / check the prior
    assert fg0>=0
    assert fg0<=1.0
    assert muNorm>=0
    assert accScaleLength>=0
    assert 0<=fcool and fcool<=1

    # Create experiment
    basename = chainDirRel+'_'
    
    if name is None:
        name = basename+str(runNumber).zfill(3)+'_'+str(rank).zfill(5)+'_'+str(procCounter).zfill(5)
    thisExper = exper.experiment(copy.copy(name))

    # We need to put the random factors for the accretion history into a file for gidget to read.
    directory = chainDir[0:-len(chainDirRel)]+'/'+name+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    np.savetxt(directory+name+'_inputRandomFactors.txt', [x0,x1,x2,x3]  )

    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.irregularVary('Mh0', Mh0)
    #thisExper.irregularVary('R', 50*accScaleLength/.05 *(Mh0/1.0e12)**(1.0/3.0) )
    thisExper.irregularVary('R', 20.0)
    thisExper.irregularVary('xmin', .01)
    thisExper.irregularVary('TOL', 1.0e-3)
    thisExper.irregularVary('NChanges', 4)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4+2**1+2**0)
    thisExper.irregularVary('alphaMRI', 0.005)
    thisExper.irregularVary('zstart',4.99)
    thisExper.irregularVary('zrelax',5.0)
    thisExper.irregularVary('NPassive',2)
    thisExper.irregularVary('Noutputs',12)

    #thisExper.irregularVary('eta',1.5)
    #thisExper.irregularVary('epsff',0.01)
    #thisExper.irregularVary('fixedQ',2.0)
    #thisExper.irregularVary('xiREC',0.0)
    #thisExper.irregularVary('kappaMetals',1.0)
    #thisExper.irregularVary('ZIGM',.02*.01)

    thisExper.irregularVary('eta',eta)
    thisExper.irregularVary('epsff',epsff)
    thisExper.irregularVary('fg0',fg0)
    thisExper.irregularVary('muNorm',muNorm)
    thisExper.irregularVary('muMhScaling',muMhScaling)
    thisExper.irregularVary('muHgScaling',muHgScaling)
    thisExper.irregularVary('fixedQ',fixedQ)
    thisExper.irregularVary('accScaleLength',accScaleLength)
    #thisExper.irregularVary('xiREC',xiREC)
    thisExper.irregularVary('fscatter',fscatter)
    #thisExper.irregularVary('accNorm',accNorm)
    #thisExper.irregularVary('accAlphaZ',accAlphaZ)
    #thisExper.irregularVary('accAlphaMh',accAlphaMh)
    #thisExper.irregularVary('accCeiling',accCeiling)
    thisExper.irregularVary('fcool',fcool)
    #thisExper.irregularVary('kappaMetals',kappaMetals)
    #thisExper.irregularVary('ZIGM',ZIGM)
    thisExper.irregularVary('concentrationRandomFactor',conRF)

    thisExper.irregularVary('whichAccretionHistory', -344)
    return thisExper, name

# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperimentRerun(emceeParams,name=None):
    global procCounter
    global runNumber


    # unpack emceeParams
    # eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling = emceeParams


    # Create experiment
    basename = chainDirRel+'-rerun-ppd_'
    
    if name is None:
        name = basename+str(rank).zfill(5)+'_'+str(procCounter).zfill(5)


    # if any of the following assertions fail, you should probably adjust / check the prior
    try:
        assert fg0>=0
        assert fg0<=1.0
        assert muNorm>=0
        assert accScaleLength>=0
        assert 0<=fcool and fcool<=1
    except:
        print "Assertions failed. ", emceeParams, name
        assert False


    thisExper = exper.experiment(copy.copy(name))

    # We need to put the random factors for the accretion history into a file for gidget to read.
    directory = chainDir[0:-len(chainDirRel)]+'/'+name+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    np.savetxt(directory+name+'_inputRandomFactors.txt', [x0,x1,x2,x3]  )

    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.irregularVary('Mh0', Mh0)
    #thisExper.irregularVary('R', 50*accScaleLength/.05 *(Mh0/1.0e12)**(1.0/3.0) )
    thisExper.irregularVary('R', 20.0)
    thisExper.irregularVary('xmin', .01)
    thisExper.irregularVary('TOL', 1.0e-3)
    thisExper.irregularVary('NChanges', 4)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4+2**1+2**0)
    thisExper.irregularVary('alphaMRI', 0.005)
    thisExper.irregularVary('zstart',4.99)
    thisExper.irregularVary('zrelax',5.0)
    thisExper.irregularVary('NPassive',20)
    thisExper.irregularVary('Noutputs',120)

    #thisExper.irregularVary('eta',1.5)
    #thisExper.irregularVary('epsff',0.01)
    #thisExper.irregularVary('fixedQ',2.0)
    #thisExper.irregularVary('xiREC',0.0)
    #thisExper.irregularVary('kappaMetals',1.0)
    #thisExper.irregularVary('ZIGM',.02*.01)

    thisExper.irregularVary('eta',eta)
    thisExper.irregularVary('epsff',epsff)
    thisExper.irregularVary('fg0',fg0)
    thisExper.irregularVary('muNorm',muNorm)
    thisExper.irregularVary('muMhScaling',muMhScaling)
    thisExper.irregularVary('muHgScaling',muHgScaling)
    thisExper.irregularVary('fixedQ',fixedQ)
    thisExper.irregularVary('accScaleLength',accScaleLength)
    #thisExper.irregularVary('xiREC',xiREC)
    thisExper.irregularVary('fscatter',fscatter)
    #thisExper.irregularVary('accNorm',accNorm)
    #thisExper.irregularVary('accAlphaZ',accAlphaZ)
    #thisExper.irregularVary('accAlphaMh',accAlphaMh)
    #thisExper.irregularVary('accCeiling',accCeiling)
    thisExper.irregularVary('fcool',fcool)
    #thisExper.irregularVary('kappaMetals',kappaMetals)
    #thisExper.irregularVary('ZIGM',ZIGM)
    thisExper.irregularVary('concentrationRandomFactor',conRF)

    thisExper.irregularVary('whichAccretionHistory', -344)
    return thisExper, name


def lnprior(emceeParams):
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling = emceeParams
    #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams
    #fg0, muNorm, muScaling, accScaleLength, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool = emceeParams
    accum = 0.0
    accum += lnLogNormalDensity(eta, np.log(1.5), np.log(3)**2.0)
    accum += lnLogNormalDensity(epsff, np.log(0.01), np.log(3)**2.0)
    accum += lnBetaDensity(fg0, 1.0, 1.0)
    accum += lnGammaDensity(muNorm, 1, 1)
    accum += lnNormalDensity(muMhScaling, -.5, 3.0**2.0)
    accum += lnNormalDensity(muHgScaling, 0.0, 2.0**2.0)
    #accum += lnGammaDensity(fixedQ-1.0, 2, 2)
    accum += lnLogNormalDensity(fixedQ, np.log(2.0), np.log(3.0)**2.0)
    #accum += lnBetaDensity(accScaleLength, .5, 9.5)
    accum += lnLogNormalDensity(accScaleLength, np.log(.05), np.log(10)**2.0)
    #accum += lnBetaDensity(xiREC, .1, .9)
    #accum += lnGammaDensity(accNorm, 0.3, 1)
    #accum += lnNormalDensity(accAlphaZ, 0.38, 0.5)
    #accum += lnNormalDensity(accAlphaMh, -.25, 0.5) 
    #accum += lnBetaDensity(accCeiling, 1.0, 1.0)
    accum += lnBetaDensity(fcool, 1.0, 1.0 )
    #accum += lnLogNormalDensity(kappaMetals, 0, np.log(3)**2.0 )
    #accum += lnBetaDensity(ZIGM, 2, 998)
    #accum += lnLogNormalDensity(ZIGM, np.log(.0002), np.log(10.0)**2.0 )
    accum += lnLogNormalDensity(Mh0, np.log(1.0e12), np.log(3)**2.0 )
    #accum += lnBetaDensity(alphaMRI, 1, 100)
    for x in [x0,x1,x2,x3]:
        accum += lnNormalDensity(x, 0, 1)
    accum += lnUniformDensity(fscatter, 0.0, 0.5)
    accum += lnLogNormalDensity(obsScale, 0, np.log(2.0)**2.0)
    accum += lnNormalDensity(conRF, 0, 0.5**2.0) # dex.
    if not np.isfinite(accum):
        return -np.inf
    return accum

def sampleFromPrior():
    fscatter = sampleFromUniformDensity(0.0, 0.5)
    r=np.sqrt(np.random.uniform())
    phi=np.random.uniform()*2.0*np.pi
    return [ \
            sampleFromLogNormalDensity( np.log(1.5), np.log(3.0)**2.0), # eta
            sampleFromLogNormalDensity( np.log(0.01), np.log(3.0)**2.0), # epsff
            sampleFromBetaDensity(1.0,1.0), # fg0
            sampleFromGammaDensity(1.0, 1.0), # muNorm
            sampleFromNormalDensity(-.5, 3.0**2.0), # muMhScaling
            sampleFromLogNormalDensity( np.log(2.0), np.log(3.0)**2.0), # fixedQ
            sampleFromLogNormalDensity(np.log(.05),np.log(10.0)**2.0), # accScaleLength
            sampleFromBetaDensity( 1.0, 1.0 ), #fcool
            sampleFromLogNormalDensity( np.log(1.0e12), np.log(3.0)**2.0), # Mh0
            fscatter,
            sampleFromNormalDensity(0,1), #x0
            sampleFromNormalDensity(0,1), #x1
            sampleFromNormalDensity(0,1), #x2
            sampleFromNormalDensity(0,1), #x3
            sampleFromLogNormalDensity( 0, np.log(2.0)**2.0),  # obsScale
            sampleFromNormalDensity(0,0.5),  # concentrationRandomFactor
            sampleFromNormalDensity(0.0, 2.0**2.0) ] # muHgScaling


def lnlikelihood(emceeParams):
    # Set up the experiment
    #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling= emceeParams
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling= emceeParams
    time0 = time.clock()

    experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)

    # Run the experiment.
    experToRun.localRun(1,0,maxTime=3600)

    output = readoutput.Experiment(name)
    output.read(keepOnly=['vPhi','colst'])


    if len(output.models)==0:
        print "WARNING: Model did not return sensible results, setting likelihood to zero"
        return -np.inf

    model0 = output.models[0]
    zs = model0.var['z'].sensible()

    accum = 0


    # In this next section we take the rotation curve from Bhattacharjee (2014) and 
    # compare it to our model. This involves several steps. First, we decide which 
    # model from B14 to use - this is controlled by l0 and l1, two "nuissance parameters" in our model.
    # Next we interpolate the GIDGET model's circular velocity curve onto the radii from B14.
    # At each radius, we take the data to be drawn from a Normal distribution with B14's quoted errors.
    # We assume each measurement is independent, i.e. the likelihood of each point is simply multiplied.
    # To begin to account for possible systematics in the data, we allow the quoted errors to be
    # scaled by an additional nuissance parameter obsScale.

    #rotcurve = np.loadtxt(chainDir+"/../../Bhattacharjee2014.txt", skiprows=15)
    rotcurve = np.loadtxt(chainDir[0:-len(chainDirRel)]+"/../Bhattacharjee2014.txt", skiprows=15)
    #whichRotCurve = convertPairToLabel(l0,l1)
    whichRotCurve = 1
    if whichRotCurve==0:
        rc = rotcurve[0:51, 2:5]
    if whichRotCurve==1:
        rc = rotcurve[51:101, 2:5]
    if whichRotCurve==2:
        rc = rotcurve[101:151, 2:5]

    nr = np.shape(rc)[0]
    rModel = model0.var['r'].sensible(timeIndex=-1)
    rData = rc[:,0]
    vcModel = model0.var['vPhi'].atR(rData, rModel, -1)
    minR = np.min(rModel)
    maxR = np.max(rModel)
    
    for i in range(nr):
        r = rc[i,0]
        if r<minR or r>maxR:
            # if the observational data lie outside the computational domain don't count them.
            # Note that this would be bad if the computational domain depended on the model parameters!
            # But we don't do that.
            pass 
            
        else:
            vc = rc[i,1]
            dvc = rc[i,2] * obsScale
            vcM = vcModel[i]
            accum += - 0.5*(vc-vcM)**2.0/dvc**2.0

    # Next up we compare to some results enumerated in Licquia and Newman (2014) arXiv 1407.1078
    r0 = 8.3 #sampleFromNormalDensity(8.33, 0.35**2.0) # first adopt a value of R0, the distance of the sun from the galactic center

    
    # This is the Boxy & Rix (2013) value of the solar neighborhood stellar surface density
    rInterp = [8.0,8.3,8.5]
    SigStModel = model0.var['colst'].atR(rInterp,rModel,-1)[1]
    accum += - 0.5*(SigStModel-38.0)**2.0/(4.0)**2.0

    rScale = model0.var['scaleLength'].sensible(timeIndex=-1)
#    accum +=  - 0.5*(rScale-2.15)**2.0 / (0.14*obsScale)**2.0

    #SFR
    sfr = model0.var['sfr'].sensible(timeIndex=-1)
    accum +=  - 0.5*(sfr-1.65)**2.0 / (0.19)**2.0

    # Total stellar mass
    mstar = model0.var['mstar'].sensible(timeIndex=-1)
    accum +=  -0.5*((6.08e10 - mstar)/1.14e10)**2.0

#    # Bulge:Total Ratio
    BT = model0.var['BT'].sensible(timeIndex=-1)
#    mean = 0.150 + (0.028 - 0.019)/2.0
#    accum += -0.5*((mean-BT)/0.028/obsScale)**2.0

    maxColStIndex = np.argmax(model0.var['colst'].sensible(timeIndex=-1))
    time1 = time.clock()
    
    
    print "With params ",emceeParams," we get BT=",BT," sfr=",sfr,' rScale=',rScale,' mstar=',np.log10(mstar)," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. The maximum of ColSt is at ",maxColStIndex,", the number of time outputs is ",model0.nt,len(zs),zs[-1]

    return accum        

    

def lnProb(emceeParams):
    pr = lnprior(emceeParams)
    if np.isfinite(pr):
        return lnlikelihood(emceeParams) + pr
    return pr



def rerunPosteriorPredictive():
    ''' Rerun the posterior predictive distribution. This can be used to e.g. increase the resolution
        spatially or in terms of the age of stellar populations, or vary some parameter systematically.
        The mandatory argument func is a user-provided function that specifies how a model with known
        parameters should be modified and (re) run.'''
    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    output = readoutput.Experiment(chainDirRel+'-ppd') # read in the posterior predictive distribution.
    output.read(paramsOnly=True,keepStars=False)
    emcee_params = []
    print "output.models: ",len(output.models)
    # For each model, take the parameters we have read in and construct the corresponding emcee parameters.
    for model in output.models:
        #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling = emceeParams
        eta = model.p['eta']
        epsff = model.p['epsff']
        fg0 = model.p['fg0']
        muNorm = model.p['muNorm']
        muMhScaling = model.p['muMhScaling']
        fixedQ = model.p['fixedQ']
        accScaleLength = model.p['accScaleLength']
        fcool = model.p['fcool']
        Mh0 = model.p['Mh0']
        fscatter = model.p['fscatter']
        x0 = model.p['x0']
        x1 = model.p['x1']
        x2 = model.p['x2']
        x3 = model.p['x3']
        obsScale = 1.0 # doesn't matter.. see below
        conRF = model.p['concentrationRandomFactor']
        muHgScaling = model.p['muHgScaling']
        # We have everything except obsScale, but actually that doesn't matter,
        # since it only affects the model in post-processing, i.e. in comparing to the data,
        # not the running of the model itself. So..... we good!
        theList = [ eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling]
        try:
            assert eta>0 and epsff>0 and fg0>0 and fg0<=1 and fixedQ>0 and muNorm>=0 and fcool>=0 and fcool<=1 and Mh0>0
        except:
            print 'Unexpected ppd params: ',theList
        emcee_params.append( copy.deepcopy(theList) )
    # OK, from here on out, we just need to emulate parts of the run() function to trick emcee into running a single iteration of the algorithm with this IC.


    ndim =  17

    restart = {}
    restart['currentPosition'] = emcee_params
    restart['chain'] = None
    restart['state'] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0



    nwalkers = len(emcee_params) # Need one walker per sample from posterior predictive distribution
    print "Starting up the ensemble sampler!"
    sampler = emcee.EnsembleSampler(nwalkers, ndim, fakeProb, pool=pool)
    #pos, prob, state = sampler.run_mcmc(restart['currentPosition'], N, rstate0=restart['state'], lnprob0=restart['prob'])
    print "Take a step with the ensemble sampler"

    # Take a single step with the ensemble sampler.
    print np.shape(restart['currentPosition']), np.shape(np.random.uniform(0,1,nwalkers))
    sampler._get_lnprob(pos = restart['currentPosition'])

    #result = sampler.sample(restart['currentPosition'], iterations=1, lnprob0=None, rstate0=None)

    #pos, prob, state = result
    print "Close the pool"

    pool.close()


def fakeProb(params):
    # Get the gidget experiment object we want to run.
    print "Running with params ",params
    experToRun, name = emceeParameterSpaceToGidgetExperimentRerun(params)
    print "Experiment name: ",name

    # Run the experiment.
    experToRun.localRun(1,0,maxTime=3600*3)
    return np.random.uniform(0,1) # value doesn't matter -- all we care about is re-running the experiment.


def getPosteriorPredictive(restart, burnIn=0, nspace=10):
    ''' We have to find the set of models over some period of time in the chain (after burnIn)
        that represents the posterior predictive distribution of models. This is NOT the same as
        just taking a set of all models run after you think the sampler has converged because
        some (most!) of those models are not accepted! It's also a bit non-trivial because
        when the new model isn't accepted, you need to include the identical model again.'''
    #allRuns = glob.glob(chainDir+'*')
    frac=1.0
    output = readoutput.Experiment(chainDirRel)
    output.read(paramsOnly=True)
    modelDict = {}
    for model in output.models:
        key = ("%.5e" % model.p['epsff']) +'_'+ ("%.5e" % model.p['eta'])
        modelDict[key] = model

    # The following block of code ends up working but not being effective because for some reason I changed the number
    # of outputs in between the two runs, so when I can't find a model in 06, it's gone forever.

    # These are all the accScaleLengths in the posterior distribution.
    # This is a somewhat subjective decision: the user needs to have picked a
    #  burnIn time, and spacing to cut down on autocorrelation, and possibly also
    #  a fraction of models to take on board (if e.g. you don't want the full sample).
    epsffs = restart['chain'][:,burnIn::nspace, 1].flatten()
    etas = restart['chain'][:,burnIn::nspace, 0].flatten()
    assert len(epsffs)==len(etas)
    for i in range(len(epsffs)):
        if np.random.uniform()<frac:
            # We have decided to take the model, identified by the key below, and copy it from
            # the MCMC chain into our Posterior Predictive Distribution.
            key = ("%.5e" % epsffs[i]) +'_'+ ("%.5e" % etas[i])
            try:
                # To do so, we find the model in our dictionary of models
                model = modelDict[key]
                # Assign it a new name
                destname = chainDirRel+'-ppd_'+str(i).zfill(5)
                # and copy the original run to a new directory.
                shutil.copytree( model.dirname, analysisDir+'/'+destname )
                # For every file in the copied folder, replace its old prefix with its new name.
                for filename in os.listdir(analysisDir+'/'+destname):
                    filenameDest = filename[len(model.name):] # strip off the model name, leaving only e.g. _evolution.dat
                    filenameDest = destname+filenameDest # replace with the current model name.
                    os.rename(analysisDir+'/'+destname+'/'+filename, analysisDir+'/'+destname+'/'+filenameDest)
            except KeyError:
                print "WARNING: skipping selected model because it's not in mcmcIndFromMax06."
                #print "failed to find key ",key," in this list of keys ",
                #print sorted(modelDict.keys())


def run(N, p00=None, nwalkers=500):
    fn = chainDirRel+'.pickle'
    ndim =  17

    if p00 is not None:
        p0 = [p00*(1.0+0.001*np.random.randn( ndim )) for i in range(nwalkers)]
    else:
        p0 = [sampleFromPrior() for i in range(nwalkers)]

    restart = {}
    restart['currentPosition'] = p0
    restart['chain'] = None
    restart['state'] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0

    # Read in our past progress UNLESS we've been given a new starting location.
    if p00 is None:
        updateRestart(fn,restart)

    if restart['chain'] is not None:
        # This may save some time if you change something and forget to delete the .pickle file.
        restartedShape = np.shape(restart['chain'])
        print restartedShape, nwalkers, ndim
        assert restartedShape[0] == nwalkers
        assert restartedShape[2] == ndim

    global runNumber
    runNumber = restart['mcmcRunCounter']

    restart['iterationCounter'] += N
    restart['mcmcRunCounter'] += 1

    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnProb, pool=pool)
    #pos, prob, state = sampler.run_mcmc(restart['currentPosition'], N, rstate0=restart['state'], lnprob0=restart['prob'])

    for result in sampler.sample(restart['currentPosition'], iterations=N, lnprob0=restart['prob'], rstate0=restart['state']):

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

        
        saveRestart(fn,restart)

    pool.close()


def tracePlots(chain, fn, burnIn=0):
    ndim = np.shape(chain)[2]
    sq = np.sqrt(float(ndim))
    nr = int(np.ceil(sq))
    nc=nr
    while nr*(nc-1)>ndim:
        nc=nc-1

    # Plot the trace for every parameter, for a subset of the walkers.
    fig,ax = plt.subplots(nrows=nr,ncols=nc, figsize=(nc*4,nr*4))
    for dim in range(ndim):
        i = np.mod(dim, nr)
        j = ( dim -i )/nr
        for walker in range(np.shape(chain)[0]):
            if np.random.uniform(0,1) < 1.0e-1: # print every tenth-ish
                ax[i,j].plot(chain[walker,burnIn:,dim],alpha=.3,ls='--')

    plt.tight_layout()
    plt.savefig(fn+'.png')
    plt.close(fig)

def probsPlots(allProbs, fn, burnIn=0):
    ndim = np.shape(allProbs)[0]
    iters = np.shape(allProbs)[1]

    # Plot the trace of the probabilities for every walker.
    fig,ax = plt.subplots()
    for walker in range(ndim):
        ax.plot(allProbs[walker,burnIn:],alpha=.3,ls='--')
    plt.savefig(fn+'_probs.png')
    plt.close(fig)
    print "Saved "+fn+'_probs.png'

    changes = np.zeros(ndim)
    for walker in range(ndim):
        for iter in range(iters-1):
            if allProbs[walker,iter]!=allProbs[walker,iter+1]:
                changes[walker]+=1.0
    changes = changes/float(iters-1.0)
    print "Long-term acceptance fraction stats: "
    stats(changes)

    acor = np.zeros(ndim)
    for walker in range(ndim):
        acor[walker] = emcee.autocorr.integrated_time(allProbs[walker,burnIn:],  window=min([50,iters/2]))
    print "acor stats: "
    stats(acor)

def stats(array):
    print "Length: ",len(array)
    print "Mean: ",np.mean(array)
    print "Gauss percentiles: ",np.percentile(array,[2.5,16,50,84,97.5])
    print "5-spaced percentiles: ",np.percentile(array,range(96)[1::5])
    print "Min,Max: ",np.min(array), np.max(array)
    print "********************************"

def printRestart(restart):
    ''' Print quick summary info about the current state of the sampler. '''
    print "restart info: "
    print " current shape of chain: (nwalkers x niterations x ndim) ",np.shape(restart['chain'])
    print " autocorrelation lengths for each parameter: ",restart['acor']
    stats(restart['acor'])
    print " acceptance rate for each walker: ",restart['accept']
    stats(restart['accept'])

def trianglePlot(restart,fn,burnIn=0, nspace=10):
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnIn)*shp[2]
    prior = np.array([sampleFromPrior() for i in range(prs)])
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    labels = [r"$\eta$",r"$\epsilon_\mathrm{ff}$",r"$f_{g,0}$",r"$\mu_0$",r"$\mu_{M_h}$",r"Q",r"$r_\mathrm{acc}/r_\mathrm{vir}$", \
            r"$f_\mathrm{cool}$", r"$M_{h,0}$", r"$\sigma$ (dex)", \
            r"$x_0$",r"$x_1$",r"$x_2$",r"$x_3$",r'Error scaling',r'c offset', r'$\mu_{h_g}$']

    logs = [1,1,0,0,0,0,1,0,1,0,0,0,0,0,1,0,0]

    ### This code selects a fraction frac of the samples after burnIn iterations.
    ### The selection is random, meaning consecutive samples from the same walker
    ### will be common. In this case we had a kwarg frac instead of nspace.
    #sampleRed = restart['chain'][:,burnIn:,:].reshape((-1,ndim)) # (nwalkers*(niterations-burnIn)) x ndim

    #ns = np.shape(sampleRed)[0]
    #selection = np.empty((int(ns*frac),) , dtype='int')
    #np.fix(np.random.uniform(0,ns-1,int(ns*frac)), y=selection)
    #sampleRed = sampleRed[selection, :]

    ### Instead we take a sample not at random, but evenly spaced:
    sampleRed = restart['chain'][:,burnIn::nspace,:].reshape((-1,ndim))


    extents=[]
    for i in range(np.shape(sampleRed)[1]):
        if logs[i]==1:
            sampleRed[:,i] = np.log10(np.clip(sampleRed[:,i],1.0e-5,np.inf))
            prior[:,i] = np.log10(prior[:,i])
            labels[i] = r'$\log_{10}$ '+labels[i]

        mi = np.min(sampleRed[:,i])
        ma = np.max(sampleRed[:,i])
        if(mi==ma):
            extents.append([mi-1,ma+1])
        else:
            extents.append([mi,ma])

    trifig = triangle.corner(sampleRed, labels=labels, extents=extents)
    #trifigPrior = triangle.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, extents=extents)
    trifig.savefig(fn)

def saveRestart(fn,restart):
    with open(fn,'wb') as f:
        print "checkpoint to "+fn
        pickle.dump(restart,f,2)

def updateRestart(fn,restart):
    if os.path.isfile(fn):
        with open(fn,'rb') as f:
            tmp_dict = pickle.load(f)
            restart.update(tmp_dict)

def maximumPosteriorProb(x0):
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, muHg= emceeParams
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, x4, obsScale= emceeParams
    def toMinimize(x):
        xx = copy.deepcopy(x)
        xx[11] = xx[11]*1.0e12
        prob = lnProb(xx)
        print "Found ln prob = ",prob
        print "From x = ",x
        if not np.isfinite(prob) or prob<-1.0e30:
            prob=-1.0e30
        return -prob
    result = scipy.optimize.minimize(toMinimize, x0)
    print result
    return result

if __name__=="__main__":

    ##x0 = [1.5, 0.01, .99, .5, -2./3., 2, .03, .01, .1, 1.0, .002, 1.0, 0.01, .45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0]
    #x0 = sampleFromPrior()
    #x0[11] /= 1.0e12
    #maximumPosteriorProb(x0)

    #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling= emceeParams
    #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muHgScaling= emceeParams


    xmax = [  2.23796988e-01,   3.57956297e-03,   4.63899560e-01,   5.35669200e-01,
             -2.28085601e-01,   1.50434451e+00,   1.74728368e-02,   1.45029007e-01,
              2.60707088e+12,   4.88308818e-01,  -7.50932976e-01,  -9.15802005e-01,
              2.84786722e-01,   9.11660610e-02,   6.88968259e+00,   4.91919834e-01,
             -1.85137413e+00]

    
    #run(400, nwalkers=2048) 
    rerunPosteriorPredictive()
    if False: 
        # Load in the resulting chain:
        restart={}
        updateRestart(chainDirRel+'.pickle', restart)
        printRestart(restart)
        trianglePlot(restart,chainDirRel+'_triangle.png',burnIn=255, nspace=40)
        tracePlots(restart['chain'], chainDirRel+'_trace', burnIn=0)
        probsPlots(restart['allProbs'], chainDirRel+'_allProb',burnIn=0)

        # Find the maximum among all models sampled so far.
        allProbs = restart['allProbs'].flatten()
        index = np.argmax(allProbs)
        print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])
        print allProbs[index]
        indices = np.unravel_index(index, np.shape(restart['allProbs']))
        xmax = restart['chain'][indices[0],indices[1],:]
        print "Favorite coordinates: ", xmax

    #getPosteriorPredictive( restart, burnIn=290, nspace=100)




    # Now use xmax as p00 in a new chain
    #shutil.move(chainDirRel+'.pickle', chainDirRel+'_initial.pickle')

    #run(20, nwalkers=1024, p00=xmax) # now do the real chain initialized from a hopefully-reasonable location.


    ##pdb.set_trace()

    #xmax[11]  /= 1.0e12
    
    #print "Starting a new maximization procedure at this location!"
    #maximumPosteriorProb(xmax)

