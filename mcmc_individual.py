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


#chainDirRel = 'mcmcChain10'
chainDirRel = 'mcmcIndFromMax03'
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


# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperiment(emceeParams,name=None):
    global procCounter
    global runNumber


    # unpack emceeParams
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
    #fg0, muNorm, muScaling, accScaleLength, accNorm, accAlphaZ, accAlphaMh, accCeiling, fcool = emceeParams

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
    np.savetxt(directory+name+'_inputRandomFactors.txt', [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9]  )

    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.irregularVary('Mh0', Mh0)
    #thisExper.irregularVary('R', 50*accScaleLength/.05 *(Mh0/1.0e12)**(1.0/3.0) )
    thisExper.irregularVary('R', 30.0)
    thisExper.irregularVary('xmin', .01)
    thisExper.irregularVary('TOL', 1.0e-3)
    thisExper.irregularVary('NChanges', 10)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4+2**1+2**0)
    thisExper.irregularVary('alphaMRI', alphaMRI)
    thisExper.irregularVary('zstart',4.99)
    thisExper.irregularVary('zrelax',5.0)
    thisExper.irregularVary('NPassive',2)
    thisExper.irregularVary('Noutputs',10)

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
    thisExper.irregularVary('xiREC',xiREC)
    thisExper.irregularVary('fscatter',fscatter)
    #thisExper.irregularVary('accNorm',accNorm)
    #thisExper.irregularVary('accAlphaZ',accAlphaZ)
    #thisExper.irregularVary('accAlphaMh',accAlphaMh)
    #thisExper.irregularVary('accCeiling',accCeiling)
    thisExper.irregularVary('fcool',fcool)
    thisExper.irregularVary('kappaMetals',kappaMetals)
    thisExper.irregularVary('ZIGM',ZIGM)
    thisExper.irregularVary('concentrationRandomFactor',conRF)

    thisExper.irregularVary('whichAccretionHistory', -344)
    return thisExper, name


def lnprior(emceeParams):
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
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
    accum += lnBetaDensity(xiREC, .1, .9)
    #accum += lnGammaDensity(accNorm, 0.3, 1)
    #accum += lnNormalDensity(accAlphaZ, 0.38, 0.5)
    #accum += lnNormalDensity(accAlphaMh, -.25, 0.5) 
    #accum += lnBetaDensity(accCeiling, 1.0, 1.0)
    accum += lnBetaDensity(fcool, 1.0, 1.0 )
    accum += lnLogNormalDensity(kappaMetals, 0, np.log(3)**2.0 )
    #accum += lnBetaDensity(ZIGM, 2, 998)
    accum += lnLogNormalDensity(ZIGM, np.log(.0002), np.log(10.0)**2.0 )
    accum += lnLogNormalDensity(Mh0, np.log(1.0e12), np.log(3)**2.0 )
    accum += lnBetaDensity(alphaMRI, 1, 100)
    for x in [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9]:
        accum += lnNormalDensity(x, 0, 1)
    accum += lnLogUniformDensity(fscatter, .001, 10)
    accum += lnLogNormalDensity(obsScale, 0, np.log(2.0)**2.0)
    accum += lnNormalDensity(conRF, 0, 0.5**2.0) # dex.
    if not np.isfinite(accum):
        return -np.inf
    return accum

def sampleFromPrior():
    fscatter = sampleFromLogUniformDensity(.001, 10)
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
            sampleFromBetaDensity(0.1, 0.9), # xiREC
            sampleFromBetaDensity( 1.0, 1.0 ), #fcool
            sampleFromLogNormalDensity( 0, np.log(2.0)**2.0), # kappaMetals
            sampleFromLogNormalDensity( np.log(.0002), np.log(10.0)**2.0), # ZIGM
            sampleFromLogNormalDensity( np.log(1.0e12), np.log(3.0)**2.0), # Mh0
            sampleFromBetaDensity(1, 100), # alphaMRI
            fscatter,
            sampleFromNormalDensity(0,1), #x0
            sampleFromNormalDensity(0,1), #x1
            sampleFromNormalDensity(0,1), #x2
            sampleFromNormalDensity(0,1), #x3
            sampleFromNormalDensity(0,1), #x4
            sampleFromNormalDensity(0,1), #x5
            sampleFromNormalDensity(0,1), #x6
            sampleFromNormalDensity(0,1), #x7
            sampleFromNormalDensity(0,1), #x8
            sampleFromNormalDensity(0,1), #x9
            sampleFromLogNormalDensity( 0, np.log(2.0)**2.0),  # obsScale
            sampleFromNormalDensity(0,0.5),  # concentrationRandomFactor
            sampleFromNormalDensity(0.0, 2.0**2.0) ] # muHgScaling


def lnlikelihood(emceeParams):
    # Set up the experiment
    eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling= emceeParams

    experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)

    # Run the experiment.
    time0 = time.clock()
    experToRun.localRun(1,0,maxTime=3600)
    time1 = time.clock()

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
    accum += - 0.5*(SigStModel-38.0)**2.0/4.0**2.0

    rScale = model0.var['scaleLength'].sensible(timeIndex=-1)
    accum +=  - 0.5*(rScale-2.15)**2.0 / 0.14**2.0

    #SFR
    sfr = model0.var['sfr'].sensible(timeIndex=-1)
    accum +=  - 0.5*(sfr-1.65)**2.0 / 0.19**2.0

    # Total stellar mass
    mstar = model0.var['mstar'].sensible(timeIndex=-1)
    accum +=  -0.5*((6.08e10 - mstar)/1.14e10)**2.0

    # Bulge:Total Ratio
    BT = model0.var['BT'].sensible(timeIndex=-1)
    mean = 0.150 + (0.028 - 0.019)/2.0
    accum += -0.5*((mean-BT)/0.028)**2.0

    maxColStIndex = np.argmax(model0.var['colst'](timeIndex=-1))
    
    
    print "With params ",emceeParams," we get BT=",BT," sfr=",sfr,' rScale=',rScale,' mstar=',mstar," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. The maximum of ColSt is at ",maxColStIndex

    return accum        

    

def lnProb(emceeParams):
    pr = lnprior(emceeParams)
    if np.isfinite(pr):
        return lnlikelihood(emceeParams) + pr
    return pr



def run(N, p00=None, nwalkers=500):
    fn = chainDirRel+'.pickle'
    ndim =  27

    if p00 is not None:
        p0 = [p00*(1.0+0.2*np.random.randn( ndim )) for i in range(nwalkers)]
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


def tracePlots(chain, fn):
    ndim = np.shape(chain)[2]
    sq = np.sqrt(float(ndim))
    nr = int(np.ceil(sq))
    fig,ax = plt.subplots(nrows=nr,ncols=nr)

    for dim in range(ndim):
        i = np.mod(dim, nr)
        j = ( dim -i )/nr
        for walker in range(np.shape(chain)[0]):
            if np.random.uniform(0,1) < 1.0e-2: # print every hundredth-ish
                ax[i,j].plot(chain[walker,:,dim],alpha=.3,ls='--')

    plt.savefig(fn+'.png')


def printRestart(restart):
    ''' Print quick summary info about the current state of the sampler. '''
    print "restart info: "
    print " current shape of chain: (nwalkers x niterations x ndim) ",np.shape(restart['chain'])
    print " autocorrelation lengths for each parameter: ",restart['acor']
    print " acceptance rate for each walker: ",restart['accept']

def trianglePlot(restart,fn,burnIn=0):
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnIn)*shp[2]
    prior = [sampleFromPrior() for i in range(prs)]
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 = emceeParams
    labels = [r"$\eta$",r"$\epsilon_\mathrm{ff}$",r"$f_{g,0}$",r"$\mu_0$",r"$\mu_{M_h}$",r"Q",r"$r_\mathrm{acc}/r_\mathrm{vir}$", \
            r"$\xi$", r"$f_\mathrm{cool}$", r"$\kappa_Z$", r"$Z_\mathrm{IGM}$", r"$M_{h,0}$", r"$\alpha_\mathrm{MRI}$", r"$\sigma$", \
            r"$x_0$",r"$x_1$",r"$x_2$",r"$x_3$",r"$x_4$",r"$x_5$",r"$x_6$",r"$x_7$",r"$x_8$",r"$x_9$",r'Error scaling',r'c offset', r'$\mu_{h_g}$']

    sampleRed = restart['chain'][:,burnIn:,:].reshape((-1,ndim))

    extents=[]
    for i in range(np.shape(sampleRed)[1]):
        mi = np.min(sampleRed[:,i])
        ma = np.max(sampleRed[:,i])
        if(mi==ma):
            extents.append([mi-1,ma+1])
        else:
            extents.append([mi,ma])

    trifig = triangle.corner(sampleRed, labels=labels, extents=extents)
    trifigPrior = triangle.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, extents=extents)
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
    #eta, epsff, fg0, muNorm, muScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale= emceeParams
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

    xmax = [  1.65513740e+00,   1.48901222e-02,   4.16536362e-01,   1.47648169e+00,
             -1.89256167e-02,   2.61430452e+00,   6.91524490e-03,   3.19210057e-02,
              5.43699341e-01,   2.18900245e+00,   3.26680422e-05,   1.23144143e+12,
              7.41721728e-03,   2.97664659e-03,   1.00254762e+00,  -7.73690253e-01,
             -4.04853927e-01,  -1.84838373e+00,  -4.47083465e-01,   1.36789494e-01,
             -8.46508322e-01,   4.17355813e-01,  -2.43255337e-01,   1.53598702e+00,
              4.47703558e+00,   5.92558508e-01] # Manually add this from an initial run of ~5 iterations over 5000 walkers starting from samples of the prior.

    #run(100, nwalkers=1024, p00=None) 
     
    # Load in the resulting chain:
    restart={}
    updateRestart(chainDirRel+'.pickle', restart)
    printRestart(restart)
    trianglePlot(restart,chainDirRel+'_triangle.png',burnIn=60)
    tracePlots(restart['chain'], chainDirRel+'_trace')

    # Find the maximum among all models sampled so far.
    allProbs = restart['allProbs'].flatten()
    index = np.argmax(allProbs)
    print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])
    print allProbs[index]
    indices = np.unravel_index(index, np.shape(restart['allProbs']))
    xmax = restart['chain'][indices[0],indices[1],:]
    print "Favorite coordinates: ", xmax

    # Now use xmax as p00 in a new chain
    #shutil.move(chainDirRel+'.pickle', chainDirRel+'_initial.pickle')

    #run(20, nwalkers=1024, p00=xmax) # now do the real chain initialized from a hopefully-reasonable location.


    ##pdb.set_trace()

    #xmax[11]  /= 1.0e12
    
    #print "Starting a new maximization procedure at this location!"
    #maximumPosteriorProb(xmax)

