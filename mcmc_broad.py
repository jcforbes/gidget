import shutil
import pdb
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
import scipy.optimize
import scipy.stats
import time

import bolshoireader
import analyticDistributions

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


chaindirrel = 'broadPostSim153a'
haloMass = 1.0e13
analysisdir = os.environ['GIDGETDIR']+'/analysis/'
chaindir = analysisdir+chaindirrel
bolshoidir = os.environ['GIDGETDIR']+'/../bolshoi/' 

proccounter=0
runnumber = 0

#class bolshoireader:
#    def __init__(self,fn, minmass, maxmass):
#        self.filename = fn
#        trees=[]
#        with open(bolshoidir+fn,'r') as f:
#            line = f.readline()
#            while line!='':
#                sline = line.split()
#                treenum = int(sline[0])
#                mh0 = float(sline[1])
#                weightedmdot = float(sline[2])
#                if minmass<mh0 and mh0<maxmass:
#                    trees.append([treenum, mh0, weightedmdot])
#                line=f.readline()
#        self.trees = sorted(trees, key=lambda tree: tree[2])
#        self.keys = [tree[2] for tree in self.trees]
#    def copynearest(self, dirname, filename, value):
#        ''' copy a file with the closest value of tree[2]'''
#        #ind = np.searchsorted(self.keys, value)
#        ind = int(value)
#        if ind == len(self.keys):
#            ind=ind-1
#        assert ind<len(self.keys) and ind>-1
#        #rffile = bolshoidir+ repr(self.trees[ind][0])+'_rf.txt'
#        #shutil.copy(rffile, dirname)
#        # This chunk of code searches the file rf_data.txt for 
#        #  the appropriate tree number from the register as identified by ind above.
#        # Once found, the subsequent lines (containing the values of the random factors
#        #  for the main progenitor accretion history of that tree) are copied to a file
#        #  in the working directory of the gidget run.
#        with open(bolshoidir+'rf_data.txt','r') as f:
#            line = f.readline()
#            while line!='':
#                if line=='tree '+repr(self.trees[ind][0]):
#                    break # we found the tree we're looking for
#                line = f.readline()
#            line = f.readline()
#            collectlines=[]
#            while not 'tree' in line and line!='':
#                collectlines.append(line)
#                line = f.readline()
#            with open(dirname+'/'+filename,'w') as ff:
#                for line in collectlines:
#                    ff.write(line)

globalBolshoiReader = bolshoireader.bolshoireader('rf_registry4.txt',3.0e11,3.0e14, bolshoidir)
bolshoiSize =  len(globalBolshoiReader.keys)

globalPrior = analyticDistributions.jointDistribution(\
        [ analyticDistributions.simpleDistribution( 'loguniform', [haloMass*0.99, haloMass*1.01], 'Mh0'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(0.1), np.log(4.0)**2.0], 'alphaR' ), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(2.0), np.log(2.0)**2.0], 'alphaRSt0' ), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(2.0), np.log(2.0)**2.0], 'alphaRGa0' ), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.5), np.log(2.0)**2.0], 'chifg0' ), \
        analyticDistributions.simpleDistribution( 'normal', [-0.5, 0.5**2.0], 'muColScaling'), \
        analyticDistributions.simpleDistribution( 'normal', [0.3, 0.2**2.0], 'muFgScaling'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(0.1), np.log(10.0)**2.0], 'muNorm'), \
        analyticDistributions.simpleDistribution( 'normal', [-.5, 1.0**2.0], 'muMhScaling'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.0), np.log(10.0)**2.0], 'ZIGMfac'), \
        analyticDistributions.simpleDistribution( 'beta', [1,1], 'zmix' ), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.5), np.log(2.0)**2.0], 'eta'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.5), np.log(2.0)**2.0], 'Qf'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(0.03), np.log(2.0)**2.0], 'alphaMRI'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.0e-3), np.log(10.0)**2.0], 'epsquench'), \
        analyticDistributions.simpleDistribution( 'beta', [1.0,1.0], 'accCeiling'), \
        analyticDistributions.simpleDistribution( 'normal', [0.0, 0.3**2.0], 'conRF'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(0.1), np.log(2.0)**2.0], 'kZ'), \
        analyticDistributions.simpleDistribution( 'beta', [1,2], 'xiREC'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.0e-2), np.log(2.0)**2.0], 'epsff'), \
        analyticDistributions.simpleDistribution( 'normal', [0.0, 0.4**2], 'scaleAdjust'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.0e12), np.log(4.0)**2.0], 'mquench'), \
        analyticDistributions.simpleDistribution( 'lognormal', [np.log(1.0), np.log(4.0)**2.0], 'enInjFac'), \
        analyticDistributions.simpleDistribution( 'normal', [0.3, 0.3**2], 'chiZslope') ] )



# define the base experiment we want.
def emceeparameterspacetogidgetexperiment(emceeparams,name=None):
    global proccounter
    global runnumber


    # unpack emceeparams
    #accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams
    #Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC = emceeparams
    Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, initialSlope, mquench, enInjFac, chiZslope = emceeparams

    # create experiment
    basename = chaindirrel+'_'
    
    if name is None:
        name = basename+str(runnumber).zfill(3)+'_'+str(rank).zfill(5)+'_'+str(proccounter).zfill(5)
    thisexper = exper.experiment(copy.copy(name))


    # we need to put the random factors for the accretion history into a file for gidget to read.
    directory = chaindir[0:-len(chaindirrel)]+'/'+name+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    proccounter+=1

    #### Paramters to vary:
    # Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC

    thisexper.irregularVary('kappaMetals',kZ)
    thisexper.irregularVary('xiREC',xiREC)
    #Mhz0 = np.array(Mhz0)
    Mhz4 = Mhz0/(.75*33.3 * np.power(Mhz0/1.5e13,0.359)) # very roughly... We know for Mh0 = 3e11, Mhz4 =4e10, and same for Mh0=2e13->Mhz4=6e11. Using those 4 numbers and assuming a powerlaw, we get this relation.
    def Moster(Mh, mparams):
        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
        zti = 4.0
        logM1z = M10 + M11*zti/(zti+1.0)
        Nz = N10 + N11*zti/(zti+1.0)
        betaz = beta10 + beta11*zti/(zti+1.0)
        gammaz = gamma10 + gamma11*zti/(zti+1.0)
        M1 = np.power(10.0, logM1z)
        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
        return eff
    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376 + initialSlope, -0.826, 0.608, 0.329])
    eff = Moster(Mhz4,central)
    mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!
    f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
    tau4 = 12.27/(12.27+1.60) # fractional lookback time at z=4
    fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
    reff4 = 5.28*np.power(mst/1.0e10, 0.4)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    reff411 = 5.28*np.power(1.0e11/1.0e10, 0.4)*np.power(1.0+4.0,-0.6) # kpc (eq B3) at z=4
    weight1 = np.exp(-Mhz0/1.0e12)  #the ad-hoc reductions in fg and fcool make things bad for high-mass galaxies. weight1 is for tiny galaxies.
    weight2 = 1-weight1
    rat4 = 1.0/(1.0/fgz4 - 1.0)
    rat4 *= fg0mult # double the gas content
    fgUsed =1.0/(1.0/rat4 + 1.0)*(1.0*weight1+1.0*weight2) 
    thisexper.irregularVary( 'fg0', fgUsed, 5)
    thisexper.irregularVary( 'fg0mult', fg0mult)
    thisexper.irregularVary( 'Mh0', Mhz0, 5)
    fcools = 1.0/(1.0-fgUsed) * mst/(0.17*Mhz4) * (1.0*weight1+1.0*weight2)  # The factor of 0.3 is a fudge factor to keep the galaxies near Moster for longer, and maybe shift them to be a bit more in line with what we would expect for e.g. the SF MS.
    thisexper.irregularVary('fcool', fcools, 5)
    thisexper.irregularVary('NChanges', 1001)

    # This little section is pretty important - it sets the size of the computational domain.
    width = 0.4/100000.0
    asls = raccRvir*np.power(10.0, np.random.normal()*width) # * np.power(Mhz0/1.0e12,alpharmh)
    if asls<0.005:
        asls=0.005
    thisexper.irregularVary('accScaleLength', asls, 5)
    #Rs = np.power(reff4/reff411, 1.0)*50* asls/0.042 
    Rs = reff4 * 100
    rmin = reff4/3
    xmin = rmin/Rs
    #print "Mh0 is ", np.log10(Mhz0)
    #print "mst is ", np.log10(mst)
    #print "reff4 is ", reff4
    #print "accScaleLength is ", asls
    #print "Setting outer radius: ", Rs
    #print "Setting inner radius: ", rmin 
    thisexper.irregularVary( 'R', Rs , 5)
    thisexper.irregularVary('xmin',rmin/Rs,5)

    bolweights =  np.random.random()
    thisexper.irregularVary('bolshoiWeight', bolweights ,5)
    thisexper.irregularVary('Noutputs',400)
    thisexper.irregularVary('zstart',3.98)
    thisexper.irregularVary('zrelax',4.0)
    thisexper.irregularVary('muNorm', muNorm, 5)
    thisexper.irregularVary('muFgScaling', muFgScaling)
    thisexper.irregularVary('muColScaling', muColScaling)
    thisexper.irregularVary('muMhScaling', muMhScaling)
    thisexper.irregularVary('fscatter', 1.0)
    thisexper.irregularVary('accCeiling',accCeiling)
    thisexper.irregularVary('NPassive',4)
    thisexper.irregularVary('eta',eta)
    thisexper.irregularVary('yREC',0.03)
    thisexper.irregularVary( 'nx', 256 ) # FFT o'clock.
    ZLee = np.power(10.0, 5.65 + chiZslope*np.log10(mst/1.0e10) + 0.3*np.log10(1.0e10) - 8.7 ) # Z in Zsun
    thisexper.irregularVary('ZIGM', ZIGMfac*0.05*ZLee*0.02, 5)
    thisexper.irregularVary('ZIGMfac', ZIGMfac)
    thisexper.irregularVary('chiZslope', chiZslope)
    thisexper.irregularVary('deltaBeta', initialSlope)
    thisexper.irregularVary('concentrationRandomFactor', conRF) ## units of dex! # 0.7

    thisexper.irregularVary('whichAccretionHistory', -112)
    thisexper.irregularVary('ksuppress', 10.0 )
    thisexper.irregularVary('kpower', 2.0)
    #thisexper.irregularVary('dbg',2**1 + 2**0 + 2**4 + 2**13 + 2**7 + 2**5 + 2**16 + 2**2 + 2**10 )  
    thisexper.irregularVary('dbg', 2**0 + 2**1 + 2**2 + 2**4 + 2**5 + 2**6 + 2**7 + 2**10 + 2**16 )
    thisexper.irregularVary('fscatter', 1.0)
    thisexper.irregularVary('MQuench', 2.0e12)
    thisexper.irregularVary('epsquench', epsquench)
    thisexper.irregularVary('muQuench', 0.0)
    thisexper.irregularVary('gaScaleReduction', rgasRed )
    thisexper.irregularVary('stScaleReduction', rstarRed )
    thisexper.irregularVary('ZMix', zmix)
    thisexper.irregularVary('alphaMRI', alphaMRI)
    thisexper.irregularVary('Qlim', Qf+0.05)
    thisexper.irregularVary('fixedQ', Qf)
    thisexper.irregularVary('TOL', 6.0e-4)
    thisexper.irregularVary('minSigSt', 10.0)
    thisexper.irregularVary('energyInjectionFactor', enInjFac)

    return thisexper, name




def samplefromposterior():
    mu = np.array([ 0.147451569889, 2.52941811659, 2.59034734186, 2.41120695741, -0.124283831858, -0.0523679435879, 9.35680382698, -0.974822093888, 1.89905286619, 0.511551421578, 2.0337488747, 2.02929369251, 0.0642244458824, 0.00988965146683, 0.509787819545, 0.279394476293, 1.74417214913, 0.342311450585, 0.0107366934282, 0.414472066814, 1.50430105103e+12, 1.21653967757, -0.028389697679, 0.497607252288])
    sig = np.array([ 0.15602724955, 0.125382594594, 0.134442070724, 0.13119698371, 0.581097952545, 0.7667336685, 0.440834644703, 0.385460062581, 0.205840416096, 0.109821737135, 0.10631070297, 0.134348461, 0.133969583672, 0.389980322595, 0.100853519509, 0.106265874887, 0.177096457663, 0.0850124838557, 0.11420926966, 0.162960229799, 0.187855261279, 0.159803314605, 0.172088242166, 0.109903788948])*1.1

    sig[0] = sig[0]*1.3 # need a bit of a wider range on alphar, since formally the posterior includes a wider range via the inclusion of a parameter that allows alphar to depend on halo mass.

    logs = [1,1,1,1,0,0,1,0,1,0,1,1,1,1,0,0,1,0,1,0,1,1,0,0]
    for i in range(len(mu)):
        if logs[i]==1:
            mu[i] = np.log10(mu[i])

    cov = np.diag(sig*sig)
    out = np.random.multivariate_normal( mu, cov)

    for i in range(len(mu)):
        if logs[i]==1:
            out[i] = np.power(10.0, out[i])

    ### The numbers above refer to physical parameters from \eta through to alpharmh and fh.
    ### For our purposes here, the thing from which we're drawing includes mass and not the last two parameters, which are done in post-processing
    ### Therefore the coordinates we actually want are mass + the physical parameters we're varying, and not alpharmh and fh.
    return [ analyticDistributions.samplefromloguniformdensity(5.0e9, 2.0e13) ] + list(out[:-2])




def constlikelihood(emceeparams, modelname=None):
    ''' just return a constant so that posterior=prior. for testing! '''
    assert False
    return 0

def lnlikelihoodRehash(emceeparams, modelname=None):
    # set up the experiment
    #Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC = emceeparams

    time0 = time.time()
    runExists = os.path.isfile( analysisdir+'/'+name+'_sampleInfo.txt' )
    if not runExists:
        return 0 # run doesn't exist - no need to include it in our sample obv
    old_sample_info = np.loadtxt( analysisdir+'/'+name+'_sampleInfo.txt' )
    if len(old_sample_info) != 19:
        return 0 # something strange has happened.

    output = readoutput.Experiment(name)
    output.read(keepOnly=['vPhi','col','colst','ageSt','Z'])




    ## if we're being given a model that's already been run, don't run it again.
    if modelname is None:
        expertorun, name = emceeparameterspacetogidgetexperiment(emceeparams)

        # run the experiment.
        expertorun.localRun(1,0,maxTime=3600*1.2)
    else:
        name = modelname

    # read the results of the model
    output = readoutput.Experiment(name)
    # ... but only keep the radial functions to which we will compare real data.
    output.read(keepOnly=['vPhi'])


    successfullyRun=1
    if len(output.models)==0:
        print "warning: model did not return sensible results, setting likelihood to zero"
        successfullyRun=0
        #return -np.inf

    np.savetxt( analysisdir+'/'+name+'_sampleInfo.txt', list(emceeparams)+[successfullyRun] )

    return 0.0

def lnlikelihood(emceeparams, modelname=None):
    # set up the experiment
    Mh0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, scaleAdjust, mquench, enInjFac, chiZslope = emceeparams
    #accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams

    time0 = time.time()

    ## if we're being given a model that's already been run, don't run it again.
    if modelname is None:
        expertorun, name = emceeparameterspacetogidgetexperiment(emceeparams)

        # run the experiment.
        expertorun.localRun(1,0,maxTime=3600*2)
    else:
        name = modelname

    # read the results of the model
    output = readoutput.Experiment(name)
    # ... but only keep the radial functions to which we will compare real data.
    radialVars = ['vPhi', 'col', 'colst', 'colsfr', 'Z', 'ageRadial']
    output.read(keepOnly=radialVars, keepStars=True)


    successfullyRun=1
    if len(output.models)==0:
        print "warning: model did not return sensible results, setting likelihood to zero"
        successfullyRun=0
        #return -np.inf
    variables = ['Mh', 'mstar', 'sSFR', 'sfZ', 'stZ', 'gasToStellarRatioH2', 'gasToStellarRatioHI', 'halfMassStars', 'vPhi22', 'c82', 'Sigma1', 'specificJStars', 'metallicityGradientR90', 'sfsig', 'mdotBulgeG', 'fractionGI', 'tdep', 'tDepH2', 'broeilsHI', 'mStellarHalo']
    nz=4
    nRadii=20
    toFit = np.zeros( nRadii*len(radialVars) + len(variables)*nz )
    if successfullyRun==1:
        model = output.models[0]
        zinds = [ readoutput.Nearest( model.var['z'].sensible(), z)[0] for z in [0,1,2,3]]
        for j in range(len(variables)):
            for k in range(nz):
                toFit[ j*nz+k] = model.var[variables[j]].sensible(timeIndex=zinds[k])
        for l in range(len(radialVars)):
            rNew = np.linspace(0.1, 3.0, nRadii)
            rVec = model.var['rx'].sensible(timeIndex=zinds[0])
            toFit[len(variables)*nz + l*nRadii : len(variables)*nz + (l+1)*nRadii ] = model.var[radialVars[l]].atR(rNew, rVec, zinds[0], sensible=True)


    # 200 + 1 + 24
    outputList = list(emceeparams)+[successfullyRun]+list(toFit)
    np.savetxt( analysisdir+'/'+name+'_sampleInfo.txt', outputList )


    # shutil.rmtree( analysisdir+'/'+name+'/' )
    # remove the large files in the output, since we no longer need them
    try:
        os.remove( analysisdir+'/'+name+'/'+name+'_evolution.dat')
        os.remove( analysisdir+'/'+name+'/'+name+'_radial.dat')
        os.remove( analysisdir+'/'+name+'/'+name+'_stars.dat')
    except:
        print "WARNING: did not successfully remove evolution.dat radial.dat or stars.dat from ",analysisdir+'/'+name+'/'

    return 0.0

    model0 = output.models[0]
    zs = model0.var['z'].sensible()
    ts = model0.var['t'].sensible()

    accum = np.zeros((len(zs), 5)) # contribution to likelihood from each redshift & each of the 5 relations.


    def Moster(Mh, mparams):
        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
        zti = 4.0
        logM1z = M10 + M11*zti/(zti+1.0)
        Nz = N10 + N11*zti/(zti+1.0)
        betaz = beta10 + beta11*zti/(zti+1.0)
        gammaz = gamma10 + gamma11*zti/(zti+1.0)
        M1 = np.power(10.0, logM1z)
        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
        return eff
    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329])
    #eff = Moster(Mhz4,central)
    #mst = eff*Mhz4 # mstar according to the moster relation. ## at z=4 !!



    print "For the record, the redshifts we're analyzing here are ",zs
    for ti in range(len(zs[1:-1])):
        # make a list of the mstar's at this redshift.
        mstar = np.array([model.var['mstar'].sensible(timeIndex=ti) for model in output.models])

        ZHayward = -8.69 + 9.09*np.power(1.0+zs[ti],-0.017) - 0.0864*np.power(np.log10(mstar) - 11.07*np.power(1.0+zs[ti],0.094),2.0)
        ZHayward = np.power(10.0, ZHayward) * 0.02
        reff = 5.28 * np.power(mstar/1.0e10, 0.25) * np.power(1.0+zs[ti],-0.6) # kpc (eq B3) at z=4
        f0 = 1.0/(1.0 + np.power(mstar/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
        tau4 = (12.27-ts[ti])/(12.27+1.60) # fractional lookback time at z=4
        fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
        obsgasratio = fgz4/(1-fgz4)
        dataVphi = 147.0*np.power(mstar/1.0e10, 0.23) # independent of redshift, km/s.
        for i, model in enumerate( output.models ):
            modelZ = model.var['sfZ'].cgs(timeIndex=ti) 
            # Compare metallicities
            accum[ti,0] += -0.5* ( np.log10(modelZ) - np.log10(ZHayward[i]) )**2.0 / np.log10(2.0)**2.0
            
            modelReff = model.var['halfMassStars'].sensible(timeIndex=ti)
            # Compare half mass radii
            accum[ti,1] += -0.5* ( np.log10(modelReff) - np.log10(reff[i]) )**2.0 / np.log10(2.0)**2.0

            ## fgas = Mgas/(Mstar+Mgas)
            ## fgas*(Mstar/Mgas+1) = 1
            ## Mgas/Mstar = fgas / (1 - fgas)
            modelfg = model.var['fg'].sensible(timeIndex=ti)
            modelgasratio = modelfg/(1-modelfg)
            # Compare gas:stellar mass ratio
            accum[ti,2] += -0.5* ( np.log10(modelgasratio) - np.log10(obsgasratio[i]) )**2.0 / np.log10(2.0**2.0)
            modelVphi = model.var['vPhiOuter'].sensible(timeIndex=ti)
            # Compare TF relation (mstar vs vphi)
            accum[ti,3] += -0.5* ( np.log10(modelVphi) - np.log10(dataVphi[i]) )**2.0 / np.log10(2.0)**2.0

            modelMh = model.var['Mh'].sensible(timeIndex=ti)
            modelMst = mstar[i]
            dataMst = Moster(modelMh,central) * modelMh
            accum[ti,4] += -0.5 * (np.log10(dataMst) - np.log10(modelMst))**2.0 / np.log10(2.0)**2.0


    time1 = time.time()
    
    totaccum = np.sum(accum)
    
    print "with params ",emceeparams," we get total lnlikelihood = ",totaccum, accum, " requring a model runtime of ",(time1-time0)/60.0,"minutes. The number of time outputs is ",model0.nt,len(zs),zs[-1]

    return totaccum        





    

def lnprob(emceeparams):
    pr = globalPrior.lndensity(emceeparams)
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
        '''
    pool = MPIPool(comm=comm, loadbalance=True)
    if not pool.is_master():
        pool.wait()
        sys.exit(0)

    output = readoutput.Experiment(chaindirrel+'-ppd') # read in the posterior predictive distribution.
    output.read(paramsonly=True,keepstars=False)
    emcee_params = []
    print "output.models: ",len(output.models)
    # for each model, take the parameters we have read in and construct the corresponding emcee parameters.
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams
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
    output.read(paramsonly=True)
    modeldict = {}
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams
    for model in output.models:
        key = ("%.5e" % model.p['accScaleLength']) +'_'+ ("%.5e" % model.p['eta'])
        modeldict[key] = model

    # the following block of code ends up working but not being effective because for some reason i changed the number
    # of outputs in between the two runs, so when i can't find a model in 06, it's gone forever.

    # these are all the accscalelengths in the posterior distribution.
    # this is a somewhat subjective decision: the user needs to have picked a
    #  burnin time, and spacing to cut down on autocorrelation, and possibly also
    #  a fraction of models to take on board (if e.g. you don't want the full sample).
    accScaleLength = restart['chain'][:,burnin::nspace, 0].flatten()
    eta = restart['chain'][:,burnin::nspace, 6].flatten()
    assert len(accScaleLength)==len(etas)
    print "copying over models as ",len(eta)," samples from the ppd."
    for i in range(len(eta)):
        if np.random.uniform()<frac:
            # we have decided to take the model, identified by the key below, and copy it from
            # the mcmc chain into our posterior predictive distribution.
            key = ("%.5e" % accScaleLength[i]) +'_'+ ("%.5e" % eta[i])
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


def run(n, p00=None, nwalkers=500, prior=True, pikfilename=None):
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams


    fn = chaindirrel+'.pickle'
    ndim = 24

    if p00 is not None:
        p0 = [p00*(1.0+0.01*np.random.randn( ndim )) for i in range(nwalkers)]
    else:
        if prior:
            p0 = [globalPrior.sample() for i in range(nwalkers)]
        else:
            assert pikfilename is not None

            arr = pickle.load( open(pikfilename,'r') )
            assert nwalkers>len(arr[:,0])
            print "Simulating cases for ",pikfilename, "with ",nwalkers,"walkers and ",len(arr[:,0]),"samples"
            p0 = []
            for i in range(nwalkers):
                j = i % len(arr[:,0])
                p0.append( [haloMass]+list(arr[j,:-1]) ) 

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

    if False:
        pool = MPIPool(comm=comm, loadbalance=True)
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
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

    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

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
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
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
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import corner
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnin)*shp[2]
    prior = np.array([globalPrior.sample() for i in range(prs)])
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams
    labels = [r"$r_\mathrm{acc}/r_\mathrm{vir}$", r"$\mu_0$", r"$\mu_M$", r'$\mu_{f_g}$', r"$\mu_{\Sigma}$", r'$\epsilon_\mathrm{acc}$', r"$\eta$",r"$Q_f$",r"$Q_\mathrm{lim}$", \
            r'c offset (dex)',  \
            r'$\kappa_0$', r'$\kappa_\mathrm{scaling}$']

    logs = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]

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

    trifig = corner.corner(samplered, labels=labels, extents=extents)
    #trifigprior = corner.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, extents=extents)
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


    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalization, kappaMassScaling = emceeparams

    #xmax = [ 0.05, 1.0, -2.0/3.0, 0.23, -0.43, 0.5, 
    #         0.5, 2.0, 3.0, .3, 0.2, 2.0/3.0]


#    xmax = [ 0.2, 0.01, 0.7, 1.4, # eta, epsff, fg0, muNorm
#            -0.43, 1.0, 1.0, .025, # muColScaling, fixedQ, Qlim, accScaleLength
#            0.3, 4.5e12,  # fcool, Mh0, 
#            0.4, .23,  # conRF, muFgScaling,
#            .002, # ZIGM, 
#            8, 220, .5, # R0, V0, epsAcc 
#            3.0, bolshoiSize/2] # systematicVphiError, bolshoiWeight

   
    #run(400, nwalkers=1024,p00=xmax) 
    #run(2, nwalkers=131072)
    run(1, nwalkers=131072, prior=False, pikfilename='py/fakemcmc153c_filteredposterior.pickle')
    #rerunPosteriorPredictive()

    if False:
        # Load in the resulting chain:
        restart={}
        updaterestart(chaindirrel+'.pickle', restart)
        printrestart(restart)
        #dic(restart, 210, 100)
        triangleplot(restart,chaindirrel+'_triangle.png',burnin=48, nspace=300)
        #triangleMovie(restart)
        traceplots(restart['chain'], chaindirrel+'_trace', burnin=0)
        probsplots(restart['allProbs'], chaindirrel+'_allProb', restart['chain'], burnin=0)

        # Find the maximum among all models sampled so far.
        allProbs = restart['allProbs'].flatten()
        index = np.argmax(allProbs)
        print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])
        #print [ "{:0.10f}".format(f)  for f in allProbs[index] ]
        indices = np.unravel_index(index, np.shape(restart['allProbs']))
        xmax = restart['chain'][indices[0],indices[1],:]
        print "Favorite coordinates: ", xmax
        print [ "{:0.15f}".format(f)  for f in xmax ]

        #getPosteriorPredictive( restart, burnIn=190, nspace=100)




        # Now use xmax as p00 in a new chain
        #shutil.move(chainDirRel+'.pickle', chainDirRel+'_initial.pickle')

        #run(20, nwalkers=1024, p00=xmax) # now do the real chain initialized from a hopefully-reasonable location.


        ##pdb.set_trace()

        #xmax[11]  /= 1.0e12
        
        #print "Starting a new maximization procedure at this location!"
        #maximumPosteriorProb(xmax)


