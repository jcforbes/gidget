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


chainDirRel = 'mcmcChain06'
chainDir = '/pfs/jforbes/gidget/analysis/'+chainDirRel

procCounter=0
runNumber = 0

def lnBetaDensity(theta, a,b):
    if theta<=0 or theta>=1:
        return -np.inf
    return (a-1.0)*np.log(theta) + (b-1.0)*np.log(1.0-theta)
def lnGammaDensity(theta, a,b):
    if theta<=0:
        return -np.inf
    return (a-1.0)*np.log(theta) - b*theta
def lnLogNormalDensity(theta, mean, var):
    if theta<=0:
        return -np.inf
    return -np.log(theta) - 0.5*(np.log(theta) - mean)**2.0/var
def lnNormalDensity(theta, mean, var):
    return -0.5*(theta-mean)**2.0/var

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


# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperiment(emceeParams):
    global procCounter
    global runNumber


    # unpack emceeParams
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams
    fg0, muNorm, muScaling, accScaleLength, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool = emceeParams

    # if any of the following assertions fail, you should probably adjust / check the prior
    assert fg0>=0
    assert fg0<=1.0
    assert muNorm>=0
    assert accScaleLength>=0
    assert accCeiling<=1 and accCeiling>=0
    assert 0<=fcool and fcool<=1

    # Create experiment
    basename = chainDirRel+'_'
    
    name = basename+str(runNumber).zfill(3)+'_'+str(rank).zfill(5)+'_'+str(procCounter).zfill(5)
    thisExper = exper.experiment(copy.copy(name))
    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.vary('Mh0', 1.0e10, 5.0e11, 3, 1, 4)
    thisExper.vary('R', 10*accScaleLength/.05, 50*accScaleLength/.05, 3, 1, 4)

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
    #accum += lnBetaDensity(accScaleLength, .5, 9.5)
    accum += lnLogNormalDensity(accScaleLength, np.log(.05), np.log(2)**2.0)
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

def sampleFromPrior():
    return [sampleFromBetaDensity(1.0,0.1), # fg0
    sampleFromGammaDensity(1.0, 1.0), # muNorm
    sampleFromNormalDensity(-.5, 3.0), # muScaling
    sampleFromLogNormalDensity(np.log(.05),np.log(2)**2.0), # accScaleLength
    sampleFromGammaDensity(0.3, 1), # accNorm
    sampleFromNormalDensity(0.38, 0.5), # accAlphaZ
    sampleFromNormalDensity(-0.25, 0.5), # accAlphaMh
    sampleFromBetaDensity( 1.0, 1.0 ), # accCeiling
    sampleFromBetaDensity( 1.0, 1.0 )] #fcool


def lnlikelihood(emceeParams):
    # Set up the experiment

    experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)

    # Run the experiment.
    print "Evaluating likelihood for params ",emceeParams
    experToRun.localRun(1,0,maxTime = 3000)

    output = readoutput.Experiment(name)
    output.read()

    if(len(output.models) < 3):
        print "WARNING: setting likelihood to zero because ",len(output.models)," of the 3 models produced sensible results"
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
            accum += -0.5 * logdist*logdist - np.log(np.log(hi/mid))



    return accum        

    

def lnProb(emceeParams):
    pr = lnprior(emceeParams)
    if np.isfinite(pr):
        return lnlikelihood(emceeParams) + pr
    return pr


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
    nwalkers = 500
    ndim =  9 # 15
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool, kappaMetals, ZIGM = emceeParams

    #p00 = np.array([ .9, .1, -1., .08, .50959, .38, -.25, .7, .01 ])
    #p0 = [p00*(1.0+0.2*np.random.randn( ndim )) for i in range(nwalkers)]

    p0 = [sampleFromPrior() for i in range(nwalkers)]

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
    #pos, prob, state = sampler.run_mcmc(restart['currentPosition'], N, rstate0=restart['state'], lnprob0=restart['prob'])

    counter = 0
    for result in sampler.sample(restart['currentPosition'], iterations=N, lnprob0=restart['prob'], rstate0=restart['state']):
        print "Beginning iteration number ",counter," of ",N

        pos, prob, state = result

        restart['acor'] = sampler.acor[:] # autocorr length for each param (ndim)
        restart['accept'] = sampler.acceptance_fraction[:]  # acceptance frac for each walker.
        restart['currentPosition'] = pos # same shape as p0: nwalkers x ndim
        restart['state'] = state # random number generator state
        restart['prob'] = prob # nwalkers x dim
        if restart['chain'] is None:
            restart['chain'] = sampler.chain # nwalkers x niterations x ndim
        else:
            print np.shape(restart['chain']), np.shape(sampler.chain[:,-1,:]), np.shape(sampler.chain)
            print restart['mcmcRunCounter'], restart['iterationCounter']
            #restart['chain'] = np.concatenate((restart['chain'], sampler.chain[:,-1,:]), axis=1)
            print "dbg1: ",np.shape(restart['chain']), np.shape(np.zeros((nwalkers, 1, ndim))), np.shape(np.expand_dims(pos,1))
            restart['chain'] = np.concatenate((restart['chain'], np.expand_dims(pos, 1)),axis=1)

        
        saveRestart(fn,restart)
        counter+=1

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


def printRestart(restart):
    ''' Print quick summary info about the current state of the sampler. '''
    print "restart info: "
    print " current shape of chain: (nwalkers x niterations x ndim) ",np.shape(restart['chain'])
    print " autocorrelation lengths for each parameter: ",restart['acor']
    print " acceptance rate for each walker: ",restart['accept']

def trianglePlot(restart,fn,burnIn=0):
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    trifig = triangle.corner(restart['chain'][:,burnIn:,:].reshape((-1,ndim)), \
             labels=[r'$f_{g,0}$',r'$\mu_0$',r'$\alpha_\mu$', r'$r_\mathrm{acc}/r_\mathrm{vir}$', \
             r'$\epsilon_0$',r'$\alpha_z$',r'$\alpha_{M_h}$',r'$\epsilon_\mathrm{max}$',r'$f_\mathrm{cool}$'])
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

if __name__=="__main__":
    run(20)

    restart={}
    updateRestart(chainDirRel+'.pickle', restart)
    printRestart(restart)
    #trianglePlot(restart,chainDirRel+'_triangle.png',burnIn=50)
    tracePlots(restart['chain'], chainDirRel+'_trace')









