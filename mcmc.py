import sys
import copy
import emcee
import pickle
from mpi4py import MPI
from emcee.utils import MPIPool
import numpy as np
import exper
import readoutput
import os, glob
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


chainDirRel = 'mcmcChain04'
chainDir = '/Users/jforbes/gidget/analysis/'+chainDirRel

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

# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperiment(emceeParams):
    global procCounter
    global runNumber


    # unpack emceeParams
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams
    fg0, muNorm, muScaling, accScaleLength, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool = emceeParams
    # Create experiment
    basename = chainDirRel+'_'
    
    name = basename+str(runNumber).zfill(3)+'_'+str(rank).zfill(5)+'_'+str(procCounter).zfill(5)
    thisExper = exper.experiment(copy.copy(name))
    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.vary('Mh0', 1.0e10, 5.0e11, 4, 1, 4)
    thisExper.vary('R', 10, 50, 4, 1, 4)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4+2**1+2**0)
    thisExper.irregularVary('alphaMRI', 0)
    thisExper.irregularVary('zstart',4.99)
    thisExper.irregularVary('zrelax',5.0)
    thisExper.irregularVary('Noutputs',20)

    thisExper.irregularVary('eta',1.5)
    thisExper.irregularVary('epsff',0.01)
    thisExper.irregularVary('fixedQ',2.0)
    thisExper.irregularVary('xiREC',0.0)
    thisExper.irregularVary('kappaMetals',1.0)
    thisExper.irregularVary('ZIGM',.02*.01)

    #thisExper.irregularVary('eta',eta)
    #thisExper.irregularVary('epsff',epsff)
    thisExper.irregularVary('fg0',fg0)
    thisExper.irregularVary('muNorm',muNorm)
    thisExper.irregularVary('muScaling',muScaling)
    #thisExper.irregularVary('fixedQ',fixedQ)
    thisExper.irregularVary('accScaleLength',accScaleLength)
    #thisExper.irregularVary('xiREC',xiREC)
    thisExper.irregularVary('accNorm',accNorm)
    thisExper.irregularVary('accAlphaZ',accAlphaZ)
    thisExper.irregularVary('accAlphaMh',accAlphaMh)
    thisExper.irregularVary('accCeiling',accCeiling)
    thisExper.irregularVary('fcool',fcool)
    #thisExper.irregularVary('kappaMetals',kappaMetals)
    #thisExper.irregularVary('ZIGM',ZIGM)

    return thisExper, name


def lnprior(emceeParams):
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams
    fg0, muNorm, muScaling, accScaleLength, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool = emceeParams
    accum = 0.0
#    accum += lnLogNormalDensity(eta, np.log(1.5), np.log(2)**2.0)
#    accum += lnLogNormalDensity(epsff, np.log(0.01), np.log(2)**2.0)
    accum += lnBetaDensity(fg0, 1.0, 0.1)
    accum += lnGammaDensity(muNorm, 1, 1)
    accum += lnNormalDensity(muScaling, -.5, 3.0)
#    accum += lnGammaDensity(fixedQ-1.0, 2, 2)
    accum += lnBetaDensity(accScaleLength, .05, .95)
#    accum += lnBetaDensity(xiREC, .1, .9)
    accum += lnGammaDensity(accNorm, 0.3, 1)
    accum += lnNormalDensity(accAlphaZ, 0.38, 0.5)
    accum += lnNormalDensity(accAlphaMh, -.25, 0.5) 
    accum += lnBetaDensity(accCeiling, 1.0, 1.0)
    accum += lnBetaDensity(fcool, 1.0, 1.0 )
#    accum += lnLogNormalDensity(kappaMetals, 0, np.log(2)**2.0 )
#    accum += lnBetaDensity(ZIGM, 2, 998)
    if not np.isfinite(accum):
        return -np.inf
    return accum

def lnlikelihood(emceeParams):
    # Set up the experiment
    experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)

    # Run the experiment.
    experToRun.localRun(1,0)

    output = readoutput.Experiment(name)
    output.read()

    if(len(output.models) < 4):
        print "WARNING: setting likelihood to zero because ",len(output.models)," of the 4 models produced sensible results"
        return -np.inf

    model0 = output.models[0]
    zs = model0.var['z'].sensible()

    accum = 0
    for model in output.models:
        for ti in range(len(zs)):
            Mh = model.var['Mh'].sensible(timeIndex=ti)
            lo, hi, mid = efficiency(Mh, zs[ti])
            eff = model.var['mstar'].sensible(timeIndex=ti)/Mh
            logdist = np.abs(np.log(eff/mid)/np.log(hi/mid))
            accum += -0.5 * logdist*logdist

    return accum        

    

def lnProb(emceeParams):
    return lnlikelihood(emceeParams) + lnprior(emceeParams)


# Formula from Moster, assuming no correlation bewteen parameters.
def efficiency(thisMh, z):
    def Moster(Mh, mparams):
        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
        logM1z = M10 + M11*z/(z+1.0)
        Nz = N10 + N11*z/(z+1.0)
        betaz = beta10 + beta11*z/(z+1.0)
        gammaz = gamma10 + gamma11*z/(z+1.0)
        M1 = np.power(10.0, logM1z)
        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
        return eff
    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329])
    unc =     np.array([0.236, 0.353, 0.0058, 0.0069, 0.153, 0.225, 0.059, 0.173])
    Mhs = np.array([thisMh])
    eff = Moster(Mhs, central)
    for i in range(len(unc)):
        theseParams = copy.copy(central)
        theseParams[i] = theseParams[i]+unc[i]
        eff = np.vstack([eff, Moster(Mhs, theseParams)])
        theseParams = copy.copy(central)
        theseParams[i] = theseParams[i]-unc[i]
        eff = np.vstack([eff, Moster(Mhs, theseParams)])
    effM = np.min(eff, axis=0)
    effP = np.max(eff, axis=0)
    return effM, effP, Moster(Mhs,central)

def run(N):
    fn = chainDirRel+'.pickle'
    nwalkers = 50
    ndim =  9 # 15
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams
    #p00 = np.array([1.5, 0.01, .5, .5, -2./3., 2, .05, .01, .30959, .38, -.25, .7, .5, 1, .002])
    p00 = np.array([ .5, .5, -2./3., .05, .30959, .38, -.25, .7, .01 ])



    p0 = [p00*(1.0+0.2*np.random.randn( ndim )) for i in range(nwalkers)]

    restart = {}
    restart['currentPosition'] = p0
    restart['chain'] = None
    restart['state'] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0


    updateRestart(fn,restart)

    global runNumber
    runNumber = restart['mcmcRunCounter']

    restart['iterationCounter'] += N
    restart['mcmcRunCounter'] += 1

    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)


    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnProb, pool=pool)
    pos, prob, state = sampler.run_mcmc(restart['currentPosition'], N, rstate0=restart['state'], lnprob0=restart['prob'])

    restart['acor'] = sampler.acor[:] # autocorr length for each param (ndim)
    restart['accept'] = sampler.acceptance_fraction[:]  # acceptance frac for each walker.
    restart['currentPosition'] = pos # same shape as p0: nwalkers x ndim
    restart['state'] = state # random number generator state
    restart['prob'] = prob # nwalkers x dim
    if restart['chain'] is None:
        restart['chain'] = sampler.chain # nwalkers x niterations x ndim
    else:
        restart['chain'] = np.concatenate((restart['chain'], sampler.chain), axis=1)
    
    saveRestart(fn,restart)

    pool.close()


def tracePlots(chain, fn):
    ndim = np.shape(chain)[2]
    sq = np.sqrt(float(ndim))
    nr = int(np.ceil(sq))
    fig,ax = plt.subplots(nrows=nr,ncols=nr)

    for dim in range(ndim):
        i = np.mod(dim, nr)
        j = ( dim -i )/nr
        for walker in range(np.shape(chain)[0]):
            ax[i,j].plot(chain[walker,:,dim])

    plt.savefig(fn+'.png')


def saveRestart(fn,restart):
    with open(fn,'wb') as f:
        print "checkpoint to "+fn
        pickle.dump(restart,f,2)

def updateRestart(fn,restart):
    if os.path.isfile(fn):
        with open(fn,'rb') as f:
            tmp_dict = pickle.load(f)
            restart.update(tmp_dict)


run(200)

restart={}
updateRestart(chainDirRel+'.pickle', restart)
tracePlots(restart['chain'], chainDirRel+'_trace')









