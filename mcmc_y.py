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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.stats
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


chainDirRel = 'mcmcIndFromMax23'
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
    eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling,  y0,y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams

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
    np.savetxt(directory+name+'_inputRandomFactors.txt', [x0,x1,0.0]  ) # fscatter is accounted for in gidget proper.
    yscatter = 0.5 # Check!
    np.savetxt(directory+name+'_inputRandomFactorsY.txt', [y0*yscatter,y1*yscatter]  )

    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.irregularVary('Mh0', Mh0)
    #thisExper.irregularVary('R', 50*accScaleLength/.05 *(Mh0/1.0e12)**(1.0/3.0) )
    thisExper.irregularVary('R', 30.0)
    thisExper.irregularVary('xmin', .02)
    thisExper.irregularVary('TOL', 1.0e-3)
    thisExper.irregularVary('NChanges', 3)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4 + 2**1 + 2**0)
    thisExper.irregularVary('alphaMRI', 0.005)
    thisExper.irregularVary('zstart',3.99)
    thisExper.irregularVary('zrelax',4.0)
    thisExper.irregularVary('NPassive',2)
    thisExper.irregularVary('Noutputs',12)


    thisExper.irregularVary('yREC', .054*.54)
    thisExper.irregularVary('eta',eta)
    thisExper.irregularVary('epsff',epsff)
    thisExper.irregularVary('fg0',fg0)
    thisExper.irregularVary('muNorm',muNorm)
    thisExper.irregularVary('muColScaling',muColScaling)
    thisExper.irregularVary('muFgScaling',muFgScaling)
    thisExper.irregularVary('fixedQ',fixedQ)
    thisExper.irregularVary('Qlim',Qlim)
    thisExper.irregularVary('accScaleLength',accScaleLength)
    thisExper.irregularVary('fscatter',0.45)
    thisExper.irregularVary('accNorm',epsilonAcc)
    #thisExper.irregularVary('accAlphaZ',accAlphaZ)
    #thisExper.irregularVary('accAlphaMh',accAlphaMh)
    #thisExper.irregularVary('accCeiling',accCeiling)
    thisExper.irregularVary('fcool',fcool)
    thisExper.irregularVary('ZIGM',ZIGM)
    thisExper.irregularVary('concentrationRandomFactor',conRF)

    thisExper.irregularVary('whichAccretionHistory', -344)


    # R0 and V0 are irrelevant to running the physical model and only matter in the comparison to data.
    return thisExper, name

# Define the base experiment we want.
def emceeParameterSpaceToGidgetExperimentRerun(emceeParams,name=None):
    global procCounter
    global runNumber


    # unpack emceeParams
    # eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling = emceeParams
    eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0, y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams


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
    np.savetxt(directory+name+'_inputRandomFactors.txt', [x0,x1, 0.0]  )
    yscatter=0.5
    np.savetxt(directory+name+'_inputRandomFactorsY.txt', [y0*yscatter,y1*yscatter]  )

    procCounter+=1

    # Set up a range of masses and computational domains
    thisExper.irregularVary('Mh0', Mh0)
    #thisExper.irregularVary('R', 50*accScaleLength/.05 *(Mh0/1.0e12)**(1.0/3.0) )
    thisExper.irregularVary('R', 30.0)
    thisExper.irregularVary('xmin', .01)
    thisExper.irregularVary('TOL', 1.0e-3)
    thisExper.irregularVary('NChanges', 3)

    # Set up some common parameters.
    thisExper.irregularVary('dbg', 2**4+2**1+2**0)
    thisExper.irregularVary('alphaMRI', 0.005)
    thisExper.irregularVary('zstart',3.99)
    thisExper.irregularVary('zrelax',4.0)
    thisExper.irregularVary('NPassive',20)
    thisExper.irregularVary('Noutputs',120)
    thisExper.irregularVary('fscatter',0.45)
    thisExper.irregularVary('yREC', .054*.54)

    thisExper.irregularVary('eta',eta)
    thisExper.irregularVary('epsff',epsff)
    thisExper.irregularVary('fg0',fg0)
    thisExper.irregularVary('muNorm',muNorm)
    thisExper.irregularVary('muColScaling',muColScaling)
    thisExper.irregularVary('muFgScaling',muFgScaling)
    thisExper.irregularVary('fixedQ',fixedQ)
    thisExper.irregularVary('accScaleLength',accScaleLength)
    thisExper.irregularVary('Qlim',Qlim)
    thisExper.irregularVary('fcool',fcool)
    thisExper.irregularVary('ZIGM',ZIGM)
    thisExper.irregularVary('concentrationRandomFactor',conRF)
    thisExper.irregularVary('accNorm',epsilonAcc)

    thisExper.irregularVary('whichAccretionHistory', -344)
    return thisExper, name


def lnprior(emceeParams):
    eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams
    accum = 0.0
    accum += lnLogNormalDensity(eta, np.log(1.5), np.log(3)**2.0)
    accum += lnLogNormalDensity(epsff, np.log(0.01), np.log(3)**2.0)
    accum += lnBetaDensity(fg0, 1.0, 1.0)
    accum += lnLogNormalDensity(muNorm, np.log10(13.0), np.log10(3.0)**2.0)
    accum += lnNormalDensity(muColScaling, -1.15, 0.5**2.0)
    accum += lnNormalDensity(muFgScaling, 0.16, 0.5**2.0)
    accum += lnLogNormalDensity(fixedQ, np.log(2.0), np.log(3.0)**2.0)
    accum += lnLogNormalDensity(Qlim, np.log(2.0), np.log(3.0)**2.0)
    accum += lnLogNormalDensity(accScaleLength, np.log(.05), np.log(10)**2.0)
    accum += lnBetaDensity(fcool, 1.0, 1.0 )
    accum += lnLogNormalDensity(ZIGM, np.log(.0002), np.log(10.0)**2.0 )
    accum += lnLogNormalDensity(Mh0, np.log(1.0e12), np.log(3.0)**2.0 )
    for x in [x0,x1]:
        accum += lnNormalDensity(x, 0, 1)
    for y in [y0,y1]:
        accum += lnNormalDensity(y, 0, 1)
    accum += lnLogNormalDensity(obsScale, 0, np.log(2.0)**2.0)
    accum += lnNormalDensity(conRF, 0, 0.5**2.0) # dex.
    corr = 0.5
    infl = 1.0
    accum += np.log(scipy.stats.multivariate_normal.pdf( [R0MCMC, V0MCMC], [ 8.1, 246], [[0.6*0.6*infl, 0.6*24 * corr*infl], [0.6*24*corr*infl, 24*24*infl]] ))
    accum += lnBetaDensity(epsilonAcc, 1.0, 1.0)
    if not np.isfinite(accum):
        return -np.inf
    return accum

def sampleFromPrior():
    # eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams
    r=np.sqrt(np.random.uniform())
    phi=np.random.uniform()*2.0*np.pi
    corr = 0.5
    infl = 1.0
    rr,vv = np.random.multivariate_normal( [ 8.1, 246], [[0.6*0.6*infl, 0.6*24 * corr*infl], [0.6*24*corr*infl, 24*24*infl]] )

    return [ \
            sampleFromLogNormalDensity( np.log(1.5), np.log(3.0)**2.0), # eta
            sampleFromLogNormalDensity( np.log(0.01), np.log(3.0)**2.0), # epsff
            sampleFromBetaDensity(1.0,1.0), # fg0
            sampleFromNormalDensity(np.log10(13.0), np.log10(3.0)**2.0), # muNorm
            sampleFromNormalDensity(-1.15, .5**2.0), # muColScaling
            sampleFromLogNormalDensity( np.log(2.0), np.log(3.0)**2.0), # fixedQ
            sampleFromLogNormalDensity( np.log(2.0), np.log(3.0)**2.0), # Qlim
            sampleFromLogNormalDensity(np.log(.05),np.log(10.0)**2.0), # accScaleLength
            sampleFromBetaDensity( 1.0, 1.0 ), #fcool
            sampleFromLogNormalDensity( np.log(1.0e12), np.log(3.0)**2.0), # Mh0
            sampleFromNormalDensity(0,1), #x0
            sampleFromNormalDensity(0,1), #x1
            sampleFromLogNormalDensity( 0, np.log(2.0)**2.0),  # obsScale
            sampleFromNormalDensity(0,0.5),  # concentrationRandomFactor
            sampleFromNormalDensity(0.16, 0.5**2.0), # muFgScaling
            sampleFromNormalDensity(0,1), #y0
            sampleFromNormalDensity(0,1), #y1
            sampleFromLogNormalDensity(np.log10(0.0002), np.log10(10.0)**2.0), # ZIGM
            rr, # R0MCMC
            vv, # V0MCMC
            sampleFromBetaDensity(1.0,1.0) ] # epsilonAcc


def constlikelihood(emceeParams, modelName=None):
    ''' Just return a constant so that posterior=prior. For testing! '''
    assert False
    return 0

def lnlikelihood(emceeParams, modelName=None):
    # Set up the experiment

    time0 = time.time()
    eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams

    ## If we're being given a model that's already been run, don't run it again.
    if modelName is None:
        experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)

        # Run the experiment.
        experToRun.localRun(1,0,maxTime=3600*3)
    else:
        name = modelName

    # Read the results of the model
    output = readoutput.Experiment(name)
    # ... but only keep the radial functions to which we will compare real data.
    output.read(keepOnly=['vPhi','colst','colH2','colHI'])


    if len(output.models)==0:
        print "WARNING: Model did not return sensible results, setting likelihood to zero"
        return -np.inf

    model0 = output.models[0]
    zs = model0.var['z'].sensible()

    accum = 0


    # In this next section we take the rotation curve from Bhattacharjee (2014) and 
    # compare it to our model. This involves several steps. 
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

    rModel = model0.var['r'].sensible(timeIndex=-1)




    rDataH2 = np.array( [.5, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, \
                7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25, 13.75, 14.25, 14.75])
    colH2Data =np.array( [26.7, 6.548974943052395, 4.498861047835993, 2.1526195899772205, 2.7220956719817764, \
            4.316628701594535, 4.316628701594535, 11.241457858769932, 12.220956719817767, \
            10.19362186788155, 10.125284738041003, 10.444191343963555, 7.551252847380409, \
            4.088838268792713, 4.7038724373576315, 3.6332574031890665, 3.5193621867881575, \
            2.5626423690205016, 3.041002277904326, 2.927107061503417, 1.9476082004555835, \
            2.175398633257405, 1.6059225512528492, 1.287015945330296, 2.1070615034168583, \
            1.150341685649206, 1.7425968109339394, 1.6514806378132114, 2.835990888382689] )
    colH2Error = np.array( [20.3+26.7, 9.09185803757829, 6.628392484342381, 4.1022964509394555, 4.77035490605428, \
            6.816283924843425, 7.0459290187891455, 14.749478079331944, 16.29436325678497, \
            14.373695198329855, 14.436325678496871, 14.937369519832988, 11.784968684759917, \
            7.839248434237998, 8.298538622129438, 7.3382045929018815, 8.653444676409187, \
            6.9624217118997915, 7.171189979123174, 7.171189979123174, 6.064718162839251, \
            6.127348643006265, 4.937369519832988, 4.54070981210856, 5.041753653444678, \
            3.810020876826723, 4.039665970772443, 3.6638830897703585, 4.979123173277662 ] ) - colH2Data


    rDataHI = np.array([.75+i*.5 for i in range(39)])
    colHIData = np.array([ 2.175398633257405, 1.6514806378132114, 1.5375854214123024, 1.5831435079726646, \
            1.7198177676537618, 2.1070615034168583, 2.357630979498861, 2.9498861047836016, \
            3.3826879271070602, 3.610478359908882, 3.7243735763097945, 3.7927107061503413, \
            3.5193621867881575, 3.2915717539863323, 3.0865603644646917, 3.1093394077448764, \
            4.56719817767654, 3.906605922551254, 3.5876993166287043, 3.5193621867881575, \
            3.6332574031890665, 3.4738041002277917, 2.995444191343964, 2.927107061503417, \
            2.6082004555808673, 2.3120728929384953, 1.9476082004555835, 1.6287015945330303, \
            1.5148063781321177, 1.3097949886104807, 1.3781321184510276, 1.4009111617312051, \
            1.1958997722095681, 1.3097949886104807, 1.1731207289293835, 1.1958997722095681, \
            1.3781321184510276, 1.264236902050115, 1.287015945330296 ])
    colHIError = colHIData* 0.2



    #rDataH2 = np.arange(11)+0.5  # quoted bin centers
    vcModelAtH2Data = model0.var['vPhi'].atR( rDataH2, rModel, -1)
    vcModelAtHIData = model0.var['vPhi'].atR( rDataHI, rModel, -1)
    V0Assumed = 220 # double check this value
    R0Assumed = 8.5 # double check this value
    rDataH2 = R0MCMC / ( (V0MCMC/ vcModelAtH2Data) * ( ( V0Assumed/V0MCMC) * (1.0 * (R0Assumed/rDataH2) - 1.0) + 1.0) ) # Adjust for model's circular velocity and this MCMC iteration's values of V0 and R0.
    rDataHI = R0MCMC / ( (V0MCMC/ vcModelAtHIData) * ( ( V0Assumed/V0MCMC) * (1.0 * (R0Assumed/rDataHI) - 1.0) + 1.0) ) # Adjust for model's circular velocity and this MCMC iteration's values of V0 and R0.
    #colH2Data = [26.7, 6.0, 3.5, 3.5, 4.6, 3.9, 2.8, 1.9, 0.9, 0.5, 0.3]
    #colH2Error= [20.3, 5.4, 1.8, 1.6, 2.0, 2.0, 1.4, 1.0, 0.7, 0.4, 0.2]
    nr = np.shape(rc)[0]
    rData = rc[:,0]
    vcModel = model0.var['vPhi'].atR(rData, rModel, -1)
    colH2Model = model0.var['colH2'].atR(rDataH2, rModel, -1)
    colHIModel = model0.var['colHI'].atR(rDataHI, rModel, -1)
    minR = np.min(rModel)
    maxR = np.max(rModel)

    for i in range(len(rDataH2)):
        r = rDataH2[i]
        if r<minR or r>maxR:
            pass
        else:
            accum += -0.5*(colH2Model[i] - colH2Data[i])**2.0/colH2Error[i]**2.0

        r = rDataHI[i]
        if r<minR or r>maxR:
            pass
        else:
            accum += -0.5*(colHIModel[i] - colHIData[i])**2.0/colHIError[i]**2.0
    
    for i in range(nr):
        r = rc[i,0]
        if r<minR or r>maxR:
            # if the observational data lie outside the computational domain don't count them.
            # Note that this would be bad if the computational domain depended on the model parameters!
            # But we don't do that.
            pass 
            
        else:
            R0Assumed = 8.3
            V0Assumed = 244
            vc = (R0Assumed/R0MCMC)*( rc[i,1]  + ( V0MCMC - V0Assumed)*rc[i,0]/R0Assumed)
            dvc = rc[i,2] * obsScale
            vcM = vcModel[i]
            accum += - 0.5*(vc-vcM)**2.0/dvc**2.0

    # Next up we compare to some results enumerated in Licquia and Newman (2014) arXiv 1407.1078
    #r0 = 8.3 #sampleFromNormalDensity(8.33, 0.35**2.0) # first adopt a value of R0, the distance of the sun from the galactic center
    r0 = R0MCMC

    
    # This is the Boxy & Rix (2013) value of the solar neighborhood stellar surface density
    #rInterp = [8.0,8.3,8.5]
    rInterp = [R0MCMC, 20.0]
    SigStModel = model0.var['colst'].atR(rInterp,rModel,-1)[0]
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
    mean = 0.150 + (0.028 - 0.019)/2.0
    #accum += -0.5*((mean-BT)/0.028/obsScale)**2.0

    maxColStIndex = np.argmax(model0.var['colst'].sensible(timeIndex=-1))
    time1 = time.time()
    
    
    print "With params ",emceeParams," we get BT=",BT," sfr=",sfr,' rScale=',rScale,' mstar=',np.log10(mstar)," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. The maximum of ColSt is at ",maxColStIndex,", the number of time outputs is ",model0.nt,len(zs),zs[-1]

    return accum        





#def lnlikelihood(emceeParams, modelName=None):
#    # Set up the experiment
#    #eta, epsff, fg0, muNorm, muMhScaling, fixedQ, accScaleLength, xiREC, fcool, kappaMetals, ZIGM, Mh0, alphaMRI, fscatter, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, obsScale, conRF, muHgScaling= emceeParams
#    eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, kappaMetals, ZIGM, xiREC, R0MCMC, V0MCMC, epsilonAcc = emceeParams
#    time0 = time.clock()
#
#    ## If we're being given a model that's already been run, don't run it again.
#    if modelName is None:
#        experToRun, name = emceeParameterSpaceToGidgetExperiment(emceeParams)
#
#        # Run the experiment.
#        experToRun.localRun(1,0,maxTime=3600)
#    else:
#        name = modelName
#
#    # Read the results of the model
#    output = readoutput.Experiment(name)
#    # ... but only keep the radial functions to which we will compare real data.
#    output.read(keepOnly=['vPhi','colst','colH2'])
#
#
#    if len(output.models)==0:
#        print "WARNING: Model did not return sensible results, setting likelihood to zero"
#        return -np.inf
#
#    model0 = output.models[0]
#    zs = model0.var['z'].sensible()
#
#    accum = 0
#
#
#    # In this next section we take the rotation curve from Bhattacharjee (2014) and 
#    # compare it to our model. This involves several steps. First, we decide which 
#    # model from B14 to use - this is controlled by l0 and l1, two "nuissance parameters" in our model.
#    # Next we interpolate the GIDGET model's circular velocity curve onto the radii from B14.
#    # At each radius, we take the data to be drawn from a Normal distribution with B14's quoted errors.
#    # We assume each measurement is independent, i.e. the likelihood of each point is simply multiplied.
#    # To begin to account for possible systematics in the data, we allow the quoted errors to be
#    # scaled by an additional nuissance parameter obsScale.
#
#    #rotcurve = np.loadtxt(chainDir+"/../../Bhattacharjee2014.txt", skiprows=15)
#    rotcurve = np.loadtxt(chainDir[0:-len(chainDirRel)]+"/../Bhattacharjee2014.txt", skiprows=15)
#    #whichRotCurve = convertPairToLabel(l0,l1)
#    whichRotCurve = 1
#    if whichRotCurve==0:
#        rc = rotcurve[0:51, 2:5]
#    if whichRotCurve==1:
#        rc = rotcurve[51:101, 2:5]
#    if whichRotCurve==2:
#        rc = rotcurve[101:151, 2:5]
#
#
#    # Similar data for HI available from: http://adsabs.harvard.edu/abs/2003PASJ...55..191N
#    # HII from http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/0301598
#    # H2 from http://adsabs.harvard.edu/abs/2006PASJ...58..847N
#
#
#    rDataH2 = np.arange(11)+0.5
#    colH2Data = [26.7, 6.0, 3.5, 3.5, 4.6, 3.9, 2.8, 1.9, 0.9, 0.5, 0.3]
#    colH2Error= [20.3, 5.4, 1.8, 1.6, 2.0, 2.0, 1.4, 1.0, 0.7, 0.4, 0.2]
#    nr = np.shape(rc)[0]
#    rModel = model0.var['r'].sensible(timeIndex=-1)
#    rData = rc[:,0]
#    vcModel = model0.var['vPhi'].atR(rData, rModel, -1)
#    colH2Model = model0.var['colH2'].atR(rDataH2, rModel, -1)
#    minR = np.min(rModel)
#    maxR = np.max(rModel)
#
#    for i in range(len(rDataH2)):
#        r = rDataH2[i]
#        if r<minR or r>maxR:
#            pass
#        else:
#            accum += -0.5*(colH2Model[i] - colH2Data[i])**2.0/colH2Error[i]**2.0
#    
#    for i in range(nr):
#        r = rc[i,0]
#        if r<minR or r>maxR:
#            # if the observational data lie outside the computational domain don't count them.
#            # Note that this would be bad if the computational domain depended on the model parameters!
#            # But we don't do that.
#            pass 
#            
#        else:
#            vc = rc[i,1]
#            dvc = rc[i,2] * obsScale
#            vcM = vcModel[i]
#            accum += - 0.5*(vc-vcM)**2.0/dvc**2.0
#
#    # Next up we compare to some results enumerated in Licquia and Newman (2014) arXiv 1407.1078
#    r0 = 8.3 #sampleFromNormalDensity(8.33, 0.35**2.0) # first adopt a value of R0, the distance of the sun from the galactic center
#
#    
#    # This is the Boxy & Rix (2013) value of the solar neighborhood stellar surface density
#    rInterp = [8.0,8.3,8.5]
#    SigStModel = model0.var['colst'].atR(rInterp,rModel,-1)[1]
#    accum += - 0.5*(SigStModel-38.0)**2.0/(4.0)**2.0
#
#    rScale = model0.var['scaleLength'].sensible(timeIndex=-1)
##    accum +=  - 0.5*(rScale-2.15)**2.0 / (0.14*obsScale)**2.0
#
#    #SFR
#    sfr = model0.var['sfr'].sensible(timeIndex=-1)
#    accum +=  - 0.5*(sfr-1.65)**2.0 / (0.19)**2.0
#
#    # Total stellar mass
#    mstar = model0.var['mstar'].sensible(timeIndex=-1)
#    accum +=  -0.5*((6.08e10 - mstar)/1.14e10)**2.0
#
##    # Bulge:Total Ratio
#    BT = model0.var['BT'].sensible(timeIndex=-1)
##    mean = 0.150 + (0.028 - 0.019)/2.0
##    accum += -0.5*((mean-BT)/0.028/obsScale)**2.0
#
#    maxColStIndex = np.argmax(model0.var['colst'].sensible(timeIndex=-1))
#    time1 = time.clock()
#    
#    
#    print "With params ",emceeParams," we get BT=",BT," sfr=",sfr,' rScale=',rScale,' mstar=',np.log10(mstar)," and total lnlikelihood = ",accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. The maximum of ColSt is at ",maxColStIndex,", the number of time outputs is ",model0.nt,len(zs),zs[-1]
#
#    return accum        

    

def lnProb(emceeParams):
    pr = lnprior(emceeParams)
    if np.isfinite(pr):
        return lnlikelihood(emceeParams) + pr
        #return constlikelihood(emceeParams) + pr
    return pr


def dic(restart,burnIn,nspace):
    ''' Compute deviance information criterion'''


    epsffs = restart['chain'][:,burnIn::nspace, 1].flatten()
    etas = restart['chain'][:,burnIn::nspace, 0].flatten()


    chainKeys=[]
    keyIndices={}
    for i in range(len(epsffs)):
        key = ("%.5e" % epsffs[i]) +'_'+ ("%.5e" % etas[i])
        chainKeys.append(key)
        keyIndices[key] = i


    names = glob.glob(chainDir+'-rerun-ppd_*')
    dbar = 0.0
    counter = 0
    for name in names:
        shortName = name[ name.find('-rerun-ppd_'): ]
        ## Most of these don't matter since they're not explicitly
        ##  used in the likelihood.
        #emceeParams = [0]*25 
        ## The exception is the scaling of the observational errorbars.
        ##  To do this we need to identify the place in the chain 
        ##  where we got this model.
        thisOutput = readoutput.Experiment(shortName)
        thisOutput.read(paramsOnly=True)
        if len(thisOutput.models)==0:
            print 'WARNING: Skipping ',shortName
        else:
            model= thisOutput.models[0]
            key = ("%.5e" % model.p['epsff']) +'_'+ ("%.5e" % model.p['eta'])
            index = keyIndices[key] # find our index in the chain
            emceeParams = []
            # Assemble the emcee parameters.
            for k in range(np.shape(restart['chain'])[2]):
                emceeParams.append(restart['chain'][: , burnIn::nspace, k].flatten()[index])
            dbar += -2.0 * lnlikelihood(emceeParams, modelName=shortName)
            counter += 1
        print "Current dbar: ", dbar/counter
        print "Counter: ", counter
    dbar = dbar/counter
    print "counter = ",counter
    assert np.isfinite(dbar)

    npa = restart['chain'][:,burnIn::nspace, :] # nwalkers x niter x nparam
    pdb.set_trace()
    thetaBar = np.mean(npa, axis=0)[0]
    assert len(thetaBar)==25
    dThetaBar = -2.0 * lnlikelihood(thetaBar)

    pD = dbar - dThetaBar
    return pD+dbar, pD



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
    # eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling,  y0, y1, kappaMetals, ZIGM, xiREC, R0MCMC, V0MCMC, epsilonAcc = emceeParams
    for model in output.models:
        eta = model.p['eta']
        epsff = model.p['epsff']
        fg0 = model.p['fg0']
        muNorm = model.p['muNorm']
        muColScaling = model.p['muColScaling']
        fixedQ = model.p['fixedQ']
        Qlim = model.p['Qlim']
        accScaleLength = model.p['accScaleLength']
        fcool = model.p['fcool']
        Mh0 = model.p['Mh0']
        x0 = model.p['x0']
        x1 = model.p['x1']
        obsScale = 1.0 # doesn't matter.. see below
        conRF = model.p['concentrationRandomFactor']
        muFgScaling = model.p['muFgScaling']
        y0 = model.p['y0']
        y1 = model.p['y1']
        ZIGM = model.p['ZIGM']
        R0MCMC = 8
        V0MCMC = 220
        epsilonAcc = model.p['accNorm']

        # We have everything except obsScale, but actually that doesn't matter,
        # since it only affects the model in post-processing, i.e. in comparing to the data,
        # not the running of the model itself. So..... we good!
        theList = [ eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0, y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc ]
        try:
            assert eta>0 and epsff>0 and fg0>0 and fg0<=1 and fixedQ>0 and muNorm>=0 and fcool>=0 and fcool<=1 and Mh0>0
        except:
            print 'Unexpected ppd params: ',theList
        emcee_params.append( copy.deepcopy(theList) )
    # OK, from here on out, we just need to emulate parts of the run() function to trick emcee into running a single iteration of the algorithm with this IC.


    ndim =  21

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
    experToRun.localRun(1,0,maxTime=3600*6)
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
    print "Copying over models as ",len(epsffs)," samples from the ppd."
    for i in range(len(epsffs)):
        if np.random.uniform()<frac:
            # We have decided to take the model, identified by the key below, and copy it from
            # the MCMC chain into our Posterior Predictive Distribution.
            key = ("%.5e" % epsffs[i]) +'_'+ ("%.5e" % etas[i])
            print "Using key ", key
            # To do so, we find the model in our dictionary of models
            if key in modelDict.keys():
                model = modelDict[key]
            else:
                print "Didn't find the key in modelDict! Diagnostics:"
                print "i: ",i
                print "key: ", key
                print "len(keys())", len(modelDict.keys())
                print "Examples: ", modelDict.keys()[0]
                print "Examples: ", modelDict.keys()[1]
                print "Examples: ", modelDict.keys()[2]
                print "Chain shape: ", np.shape(restart['chain'])
                print "len(epsffs): ",len(epsffs)
                raise KeyError


            # Assign it a new name
            destname = chainDirRel+'-ppd_'+str(i).zfill(5)
            # and copy the original run to a new directory.
            shutil.copytree( model.dirname, analysisDir+'/'+destname )
            print "Copied over key ", key
            # For every file in the copied folder, replace its old prefix with its new name.
            for filename in os.listdir(analysisDir+'/'+destname):
                filenameDest = filename[len(model.name):] # strip off the model name, leaving only e.g. _evolution.dat
                filenameDest = destname+filenameDest # replace with the current model name.
                os.rename(analysisDir+'/'+destname+'/'+filename, analysisDir+'/'+destname+'/'+filenameDest)
#            except KeyError:
#                print "WARNING: skipping selected model because it's not in "+chainDirRel
#                print "failed to find key ",key," in this list of keys "
#                #print sorted(modelDict.keys())


def run(N, p00=None, nwalkers=500):
    fn = chainDirRel+'.pickle'
    ndim = 21

    if p00 is not None:
        p0 = [p00*(1.0+0.01*np.random.randn( ndim )) for i in range(nwalkers)]
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
def ksProb(arr1, arr2):
    from scipy.stats import ks_2samp
    return ks_2samp(arr1,arr2)[1]

def probsPlots(allProbs, fn, chain, burnIn=0):
    nwalker = np.shape(allProbs)[0]
    iters = np.shape(allProbs)[1]
    ndim = np.shape(chain)[2]

    # Plot the trace of the probabilities for every walker.
    fig,ax = plt.subplots()
    for walker in range(nwalker):
        ax.plot(allProbs[walker,burnIn:],alpha=.3,ls='--')
    plt.savefig(fn+'_probs.png')
    plt.close(fig)
    print "Saved "+fn+'_probs.png'

    fig,ax=plt.subplots()
    kss=[]
    kssThetas = np.zeros( (ndim, iters-1) )
    # chain ~ walker x iter x dim
    for iter in range(iters-1):
        kss.append( ksProb(allProbs[:,iter], allProbs[:,iters-1] ) )
        for k in range(ndim):
            kssThetas[k,iter] = ksProb( chain[:, iter, k], chain[:,iters-1, k])
    ax.scatter(range(iters-1)[-25:], kss[-25:], s=30, marker='+')
    colors = ['k','r','green','orange','lightblue','grey','purple','pink','yellow','blue','lightgreen','darkgreen']*5
    for k in range(ndim):
        ax.plot(range(iters-1)[-25:], kssThetas[k,-25:], color=colors[k])
    ax.set_xlabel('iter')
    ax.set_ylabel('KS prob vs last iteration')
    ax.set_ylim(-0.01,1.02)
    plt.savefig(fn+'_ks.png')
    plt.close(fig)

    



    changes = np.zeros(nwalker)
    for walker in range(nwalker):
        for iter in range(iters-1):
            if allProbs[walker,iter]!=allProbs[walker,iter+1]:
                changes[walker]+=1.0
    changes = changes/float(iters-1.0)
    print "Long-term acceptance fraction stats: "
    stats(changes)

    acor = np.zeros(nwalker)
    for walker in range(nwalker):
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

def triangleMovie(restart):
    shp = np.shape(restart['chain'])
    bn = 'triangle_frame_'
    for i in range(shp[1]):
        trianglePlot( restart, bn+str(i).zfill(3)+'.png', burnIn=i, nspace=10000 )

def trianglePlot(restart,fn,burnIn=0, nspace=10):
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnIn)*shp[2]
    prior = np.array([sampleFromPrior() for i in range(prs)])
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    #eta, epsff, fg0, muNorm, muColScaling, fixedQ, accScaleLength, fcool, Mh0, fscatter, x0, x1, x2, x3, obsScale, conRF, muFgScaling, y0,y1,y2,y3, yscatter, kappaMetals, ZIGM, xiREC = emceeParams
    #eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, ZIGM, R0MCMC, V0MCMC, epsilonAcc = emceeParams
    labels = [r"$\eta$",r"$\epsilon_\mathrm{ff}$",r"$f_{g,0}$",r"$\mu_0$",r"$\mu_{\Sigma}$",r"$Q_f$",r"$Q_\mathrm{lim}$",r"$r_\mathrm{acc}/r_\mathrm{vir}$", \
            r"$f_\mathrm{cool}$", r"$M_{h,0}$",  \
            r"$x_0$",r"$x_1$",r'Error scaling',r'c offset (dex)', r'$\mu_{f_g}$', \
            r'$y_0$', r'$y_1$',  \
            r'$Z_\mathrm{IGM}$', r'$R_0$', r'$V_0$', r'$\epsilon_\mathrm{acc}$']

    logs = [1,1,0,0,0,0,0,1,0,1,0,0,1,0,0, 0,0,1,0,0,0]

    ### This code selects a fraction frac of the samples after burnIn iterations.
    ### The selection is random, meaning consecutive samples from the same walker
    ### will be common. In this case we had a kwarg frac instead of nspace.
    #sampleRed = restart['chain'][:,burnIn:,:].reshape((-1,ndim)) # (nwalkers*(niterations-burnIn)) x ndim

    #ns = np.shape(sampleRed)[0]
    #selection = np.empty((int(ns*frac),) , dtype='int')
    #np.fix(np.random.uniform(0,ns-1,int(ns*frac)), y=selection)
    #sampleRed = sampleRed[selection, :]

    ### Instead we take a sample not at random, but evenly spaced:
    sampleRed = copy.deepcopy(restart['chain'][:,burnIn::nspace,:].reshape((-1,ndim)))


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
    plt.close(trifig)

    # Print a summary of the current variables!
    print "**********************************"
    print "Name, mean, pr mean, std dev, prior std dev"
    for i in range(np.shape(sampleRed)[1]):
        print labels[i],"   ",np.mean(sampleRed[:,i]),"   ", np.mean(prior[:,i]),"   ",np.std(sampleRed[:,i]), "    ",np.std(prior[:,i])
    print "*********************************"
    print "Correlation coefficients"
    for i in range(np.shape(sampleRed)[1]):
        for j in range(np.shape(sampleRed)[1]):
            print labels[i], " x ",labels[j]," : ",np.corrcoef(sampleRed[:,i],sampleRed[:,j])[0][1]

    # Now we're going to follow the formalism of p 284 of BDA3 to assess convergence as a fn of iteration.
    wholeChain = copy.deepcopy(restart['chain'][:,:,:]) # ~ walker x iter x dim
    m = np.shape(wholeChain)[0]*2 # number of walkers times 2 (since we're splitting each chain in 2)
    print "     ",labels
    for finalIter in range(16,np.shape(wholeChain)[1],4)[:-1]:
        if finalIter%4 == 0:
            n = finalIter/4  # throw away the first 1/4, then split the second half in two
            psi = np.zeros( (m,n,ndim) ) # j=1...m, i=1...n, k=1...ndim
            psi[0:m/2,: ,:] = wholeChain[:,2*n:3*n,:]
            psi[m/2:m,: ,:] = wholeChain[:,3*n:4*n,:]
            #psibarj = np.zeros( (m,ndim) ) # average over i=1...n for each j=1..m
            psibarj = np.mean(psi,axis=1) 
            assert np.shape(psibarj)[0] == m and np.shape(psibarj)[1] == ndim
            #psibar = np.zeros( ndim ) # average of psibarj over j=1..m for each variable k
            psibar = np.mean( psibarj, axis=0)
            assert np.shape(psibar)[0] == ndim
            #sj = np.zeros( (m,ndim) ) # variance of psi_ij about psibarj
            sj = np.var( psi, axis=1, ddof=1 )
            assert np.shape(sj)[0] == m and np.shape(sj)[1] ==ndim
            # W = np.zeros(ndim)
            W = np.mean(sj, axis=0)
            B = np.zeros(ndim)
            for k in range(ndim):
                B[k] = float(n)/float(m-1) * np.sum(np.power( psibarj[:,k] - psibar[k], 2.0) )
            vark = float(n-1)/float(n) * W + B/float(n)
            Rhat = np.sqrt( vark/W )

            if np.shape(wholeChain)[1]>500 and finalIter%100==0:
                print finalIter, "    ", Rhat
            elif np.shape(wholeChain)[1]<=500:
                print finalIter, "    ", Rhat
            else:
                pass



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

    #eta, epsff, fg0, muNorm, muColScaling, fixedQ, Qlim, accScaleLength, fcool, Mh0, x0, x1, obsScale, conRF, muFgScaling, y0,y1, kappaMetals, ZIGM, xiREC, R0MCMC, V0MCMC, epsilonAcc = emceeParams

    xmax = [ 0.2, 0.01, 0.7, 1.4, # eta, epsff, fg0, muNorm
            -0.43, 1.0, 1.0, .025, # muColScaling, fixedQ, Qlim, accScaleLength
            0.3, 4.5e12, 0.23, -0.32, # fcool, Mh0, x0, x1
            3.0, 0.4, .23, .8, # obsScale, conRF, muFgScaling, y0
            .2, .002, # y1, ZIGM, 
            8, 220, .5 ] # R0, V0, epsAcc 


#    xmax = [  1.89572200e-01,   1.04076859e-02,   7.38152246e-01,   1.37673100e+00,
#             -4.37140916e-01,   5.58972092e-01,   2.49790419e-02,   2.93190618e-01,
#              4.57221629e+12,   2.33209376e-01,  -1.63343636e-01,  -3.21647816e-01,
#              9.05444175e-02,   1.84943961e-01,   1.11707726e+01,   3.98292819e-01,
#              2.33958073e-01,   7.97679715e-01,   1.83307648e-01,   2.26457747e-01,
#              3.51223250e-01,   3.06985476e-01,   2.26597455e+00,   1.76156619e-03,
#              6.22481329e-04]


    xmax = [  1.38110984e-01,   1.56054436e-02,   6.77659749e-01,   7.51281733e-01,
      -1.61484088e-01,   1.52195432e+00,   2.00000000e+00,   2.40937635e-02,
       1.00530879e-01,   5.65076067e+12,   2.13798472e-01,  -5.05746147e-01,
       4.36293557e+00,   3.85208913e-01,   2.86150739e-01,   3.54157890e-01,
       2.19773960e-01,   2.01738761e-03,   
       8.73335819e+00,   2.68739699e+02,   1.03621353e-01]


    xmax = [  1.47974343e-01,   1.25535599e-02,   6.76028939e-01,   8.05514278e-01,
             -1.10342067e-01,   1.36792855e+00,   1.09758830e+00,   2.90302105e-02,
              1.90124109e-01,   6.18665396e+12,   2.91009898e-01,  -5.59652735e-01,
              4.13462858e+00,   3.01332231e-01,   1.44540920e-01,   1.06956484e-01,
              3.13000089e-01,   1.78235790e-03,  
              9.06850006e+00,   2.54572880e+02,   1.16864456e-01]

    xmax = [  1.34085076e-01,   1.27581051e-02,   6.72328130e-01,   4.99389544e-01,
             -8.98714744e-02,   1.75081080e+00,   1.32953768e+00,   4.17966204e-02,
              2.90496695e-01,   5.87203706e+12,   1.98650145e-01,  -5.10407939e-01,
              6.60228652e+00,   2.96848022e-01,   1.78408267e-01,   9.12929719e-02,
              2.24720233e-01,   7.26599888e-04,   9.48338953e+00,   2.75070157e+02,
              6.60213847e-02]

   
    #run(400, nwalkers=1024,p00=xmax) 
    run(400, nwalkers=1024) 
    #rerunPosteriorPredictive()

    if True:
        # Load in the resulting chain:
        restart={}
        updateRestart(chainDirRel+'.pickle', restart)
        printRestart(restart)
        #dic(restart, 210, 100)
        trianglePlot(restart,chainDirRel+'_triangle.png',burnIn=70, nspace=300)
        #triangleMovie(restart)
        tracePlots(restart['chain'], chainDirRel+'_trace', burnIn=0)
        probsPlots(restart['allProbs'], chainDirRel+'_allProb', restart['chain'], burnIn=0)

        # Find the maximum among all models sampled so far.
        allProbs = restart['allProbs'].flatten()
        index = np.argmax(allProbs)
        print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])
        print allProbs[index]
        indices = np.unravel_index(index, np.shape(restart['allProbs']))
        xmax = restart['chain'][indices[0],indices[1],:]
        print "Favorite coordinates: ", xmax

        #getPosteriorPredictive( restart, burnIn=142, nspace=100)




        # Now use xmax as p00 in a new chain
        #shutil.move(chainDirRel+'.pickle', chainDirRel+'_initial.pickle')

        #run(20, nwalkers=1024, p00=xmax) # now do the real chain initialized from a hopefully-reasonable location.


        ##pdb.set_trace()

        #xmax[11]  /= 1.0e12
        
        #print "Starting a new maximization procedure at this location!"
        #maximumPosteriorProb(xmax)


