import numpy as np
import pickle
import matplotlib.pyplot as plt
import pdb
import observationalData
import emcee
from sklearn import svm, linear_model

# Useful information about the setup of the linear models and their fits....
logVars = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0]
logPreds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1]
logRadials = [1,1,1,1,1,1]
nRadii=20
nz=4
nvars = len(logVars)

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
    assert var>0
    return np.exp(samplefromnormaldensity(mean,var))
def samplefromnormaldensity(mean,var):
    assert var>0
    return np.random.normal(mean,np.sqrt(var))
def samplefromloguniformdensity(a,b):
    assert b>a
    assert a>0
    return a*(b/a)**np.random.uniform()
def samplefromuniformdensity(a,b):
    assert b>a
    return np.random.uniform(a,b)











class lassoRidge:
    def __init__(self, sklLasso, sklRidge, expon):
        self.sklLasso = sklLasso
        self.sklRidge = sklRidge
        self.expon = expon
    def predict(self, X):
        X_lasso = X[:, :-1000]    
        X_ridge = X[:, -1000: ]
        if self.expon:
            X_ridge = np.exp(X_ridge)
        Y_lasso = self.sklLasso.predict( X_lasso )
        Y_ridge = self.sklRidge.predict( X_ridge )
        return Y_lasso + Y_ridge
# Temporarily comment out so that pickle will successfully load list of models
#    def fractionalVariances(self, X_test, Y_test):
#        '''  Report the fraction of the variance which is explained by ridge, lasso, and neither. Also just tell us the total var.'''
#        overallVariance = np.var(Y_test)
#        X_lasso = X[:, :-1000]    
#        X_ridge = X[:, -1000: ]
#        if self.expon:
#            X_ridge = np.exp(X_ridge)
#        Y_lasso = self.sklLasso.predict( X_lasso )
#        Y_ridge = self.sklRidge.predict( X_ridge )
#        lassoVariance = np.var(Y_lasso)
#        ridgeVariance = np.var(Y_ridge)
#        unexplainedVariance = np.var(Y_test - Y_ridge - Y_lasso)
#
#        return lassoVariance/overallVariance, ridgeVariance/overallVariance, unexplainedVariance/overallVariance, overallVariance

def runEmcee():
    models = pickle.load( open( 'validatedModels.pickle', 'r' ) )
    import bolshoireader
    bolshoidir = '/Users/jforbes/bolshoi/'
    globalBolshoiReader = bolshoireader.bolshoireader('rf_registry4.txt',3.0e11,3.0e14, bolshoidir)
    #bolshoiSize =  len(globalBolshoiReader.keys)
    MhGrid = np.power(10.0, np.linspace(10,13,1000))

    def lnlikelihood(emceeparams):
        # First transform the emceeparams into the same format used by 'X' in the fit of the linear models
        # emceeparams is the set of 18 parameters that we fit with Lasso, minus mass, so we have..
        lnlik = 0.0
        for k,Mh in enumerate(MhGrid):
            X1 = np.array([Mh]+list(emceeparams)).reshape((1,len(emceeparams)+1))
            for i in range(len(logVars)):
                # Take the log of the input variable if it makes sense to do so
                if logVars[i] == 1:
                    X1[:,i] = np.log10(X1[:,i])
            # Add in squares and cross-terms for every pair of input parameters
            nxv = np.shape(X1)[1]
            for i in range(nxv):
                for j in range(nxv):
                    # don't repeat combinations, only look at squares instead of cross-terms, or cross-terms with mass alone
                    if j>=i and (i==j or i==0):
                        to_add = X1[:,i]*X1[:,j]
                        X1 = np.vstack( [ X1.T, to_add ] ).T
            # OK, we've now successfully constructed X1, i.e. the part that's fit with Lasso. Now let's get X2:
            X2 = globalBolshoiReader.returnnearest( '', '', np.random.random(), np.log10(Mh))
            while len(X2)!=1000:
                print "X2 wasn't the right length! ", len(X2)
                X2 = globalBolshoiReader.returnnearest( '', '', np.random.random(), np.log10(Mh))
            X2 = np.array(X2).reshape((1,1000))
            X = np.vstack( [ X1.T, X2.T ] ).T # mc params + random factors

            # Then use the linear models to make predicitons for every quantity we will compare to observations!
            logMh0 = models[0].predict(X)
            logMh1 = models[1].predict(X)
            logMh2 = models[2].predict(X)
            logMh3 = models[3].predict(X)
            logMst0 = models[4].predict(X)
            logMst1 = models[5].predict(X)
            logMst2 = models[6].predict(X)
            logMst3 = models[7].predict(X)
            logsSFRz0  = models[8].predict(X)
            logsSFRz1 = models[9].predict(X)
            logsSFRz2 = models[10].predict(X)
            logsSFRz3 = models[11].predict(X)
            logsfZz0 = models[12].predict(X)
            logsfZz1 = models[13].predict(X)
            logsfZz2 = models[14].predict(X)
            logsfZz3 = models[15].predict(X)
            logstZz0 = models[16].predict(X)
            logstZz1 = models[17].predict(X)
            logstZz2 = models[18].predict(X)
            logstZz3 = models[19].predict(X)
            loggasToStellarRatioH2z0 = models[20].predict(X)
            loggasToStellarRatioH2z1 = models[21].predict(X)
            loggasToStellarRatioH2z2  = models[22].predict(X)
            loggasToStellarRatioH2z3  = models[23].predict(X)
            loggasToStellarRatioHIz0  = models[24].predict(X)
            loghalfMassStarsz0 = models[28].predict(X)
            loghalfMassStarsz1  = models[29].predict(X)
            loghalfMassStarsz2 = models[30].predict(X)
            loghalfMassStarsz3 = models[31].predict(X)
            logvPhi22z0  = models[32].predict(X)
            logvPhi22z1 = models[33].predict(X)
            logvPhi22z2  = models[34].predict(X)
            logvPhi22z3 = models[35].predict(X)
            logc82z0  = models[36].predict(X)
            logSigma1z0  = models[40].predict(X)
            logSigma1z1  = models[41].predict(X)
            logSigma1z2  = models[42].predict(X)
            logSigma1z3  = models[43].predict(X)
            logspecificJStarsz0  = models[44].predict(X)
            logspecificJStarsz1  = models[45].predict(X)
            logspecificJStarsz2  = models[46].predict(X)
            logspecificJStarsz3  = models[47].predict(X)
            logbroeilsHIz0  = models[72].predict(X)
            logmStellarHaloz0  = models[76].predict(X)
            logmStellarHaloz1  = models[77].predict(X)
            logmStellarHaloz2  = models[78].predict(X)
            logmStellarHaloz3 = models[79].predict(X)
            #  logPreds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1]
            # corresponding to
            # Mhz0 mstarz0 sSFRz0 sfZz0 stZz0 gasToStellarRatioH2z0 gasToStellarRatioHIz0 halfMassStarsz0 vPhi22z0 c82z0 Sigma1z0 specificJStarsz0 metallicityGradientR90z0 maxsigz0 mdotBulgeGz0 fractionGIz0 tdepz0 tDepH2z0 broeilsHIz0 mStellarHaloz0 
            # i.e. only fractionGIz0 and metallicityGradientR90z0 are not logarithmic

            # Label points that fit Moster SMHM
            distancesSMHM0 = observationalData.datasets['Moster10z0'].distance(10.0**logMh0,10.0**logMst0)
            distancesSMHM1 = observationalData.datasets['Moster10z1'].distance(10.0**logMh1,10.0**logMst1)
            distancesSMHM2 = observationalData.datasets['Moster10z2'].distance(10.0**logMh2,10.0**logMst2)
            distancesSMHM3 = observationalData.datasets['Moster10z3'].distance(10.0**logMh3,10.0**logMst3)

            # Label points that fit Brinchmann04 
            distancesMS0 = observationalData.datasets['Brinchmann04Specz0'].distance(10.0**logMst0, 10.0**logsSFRz0)
            distancesMS1 = observationalData.datasets['whitaker14Specz1'].distance(10.0**logMst1, 10.0**logsSFRz1)
            distancesMS2 = observationalData.datasets['whitaker14Specz2'].distance(10.0**logMst2, 10.0**logsSFRz2)
            distancesMS3 = observationalData.datasets['lilly13Specz3'].distance(10.0**logMst3, 10.0**logsSFRz3)

            # Label points that fit stellar MZR 
            distancesSTMZR0 = np.min(np.abs( [observationalData.datasets['kirby13'].distance(10.0**logMst0,  10.0**logstZz0 ), observationalData.datasets['gallazi05'].interior(10.0**logMst0, 10.0**logstZz0) ] ))

            # Label points that fit gas MZR 
            distancesGMZR0 = np.min( np.abs([observationalData.datasets['Tremonti04'].distance(10.0**logMst0, 10.0**logsfZz0 ), observationalData.datasets['Lee06'].distance(10.0**logMst0, 10.0**logsfZz0)] ) )
            distancesGMZR1 = observationalData.datasets['Genzel15Z1'].distance(10.0**logMst1, 10.0**logsfZz1)
            distancesGMZR2 = observationalData.datasets['Genzel15Z2'].distance(10.0**logMst2, 10.0**logsfZz2)
            distancesGMZR3 = observationalData.datasets['Genzel15Z3'].distance(10.0**logMst3, 10.0**logsfZz3)

            # Label points that fit gas fractions
            distancesHMOL0 = np.min(np.abs([ observationalData.datasets['genzel15COz0'].distance(10.0**logMst0, 10.0**loggasToStellarRatioH2z0), observationalData.datasets['genzel15Dustz0'].interior(10.0**logMst0, 10.0**loggasToStellarRatioH2z0)] ))
            distancesHMOL1 = np.min(np.abs([ observationalData.datasets['genzel15COz1'].distance(10.0**logMst1, 10.0**loggasToStellarRatioH2z1), observationalData.datasets['genzel15Dustz1'].interior(10.0**logMst1, 10.0**loggasToStellarRatioH2z1)] ))
            distancesHMOL2 = np.min(np.abs([ observationalData.datasets['genzel15COz2'].distance(10.0**logMst2, 10.0**loggasToStellarRatioH2z2), observationalData.datasets['genzel15Dustz2'].interior(10.0**logMst2, 10.0**loggasToStellarRatioH2z2)] ))
            distancesHMOL3 = np.min(np.abs([ observationalData.datasets['genzel15COz3'].distance(10.0**logMst3, 10.0**loggasToStellarRatioH2z3), observationalData.datasets['genzel15Dustz3'].interior(10.0**logMst3, 10.0**loggasToStellarRatioH2z3)] ))

            # Label points that fit gas fractions
            distancesHI = np.min( np.abs( [observationalData.datasets['papastergis12'].distance(10.0**logMst0, 10.0**loggasToStellarRatioHIz0), observationalData.datasets['peeples11'].interior(10.0**logMst0, 10.0**loggasToStellarRatioHIz0) ] ) )

            # Label points that fit structural relations
            distancesSIG10 = np.min(np.abs([ observationalData.datasets['Fang13'].distance(10.0**logMst0, 10.0**logSigma1z0 ), observationalData.datasets['Barro15HS'].interior(10.0**logMst0, 10.0**logSigma1z0) , observationalData.datasets['Barro15HQ'].interior(10.0**logMst0, 10.0**logSigma1z0)]))
            distancesSIG11 = np.min(np.abs([  observationalData.datasets['Barro151S'].interior(10.0**logMst1, 10.0**logSigma1z1) , observationalData.datasets['Barro151Q'].interior(10.0**logMst1, 10.0**logSigma1z1)]))
            distancesSIG12 = np.min(np.abs([  observationalData.datasets['Barro152S'].interior(10.0**logMst2, 10.0**logSigma1z2) , observationalData.datasets['Barro152Q'].interior(10.0**logMst2, 10.0**logSigma1z2)]))
            distancesSIG13 = np.min(np.abs([  observationalData.datasets['Barro153S'].interior(10.0**logMst3, 10.0**logSigma1z3) , observationalData.datasets['Barro153Q'].interior(10.0**logMst3, 10.0**logSigma1z3)]))

            # Label points that fit galaxy sizes
            distancesMRE0 = np.min( np.abs( [ observationalData.datasets['vdW14ETG0'].distance(10.0**logMst0, 10.0**loghalfMassStarsz0), observationalData.datasets['vdW14LTG0'].interior(10.0**logMst0, 10.0**loghalfMassStarsz0)]) )
            distancesMRE1 = np.min( np.abs( [ observationalData.datasets['vdW14ETG1'].distance(10.0**logMst1, 10.0**loghalfMassStarsz1), observationalData.datasets['vdW14LTG1'].interior(10.0**logMst1, 10.0**loghalfMassStarsz1)]) )
            distancesMRE2 = np.min( np.abs( [ observationalData.datasets['vdW14ETG2'].distance(10.0**logMst2, 10.0**loghalfMassStarsz2), observationalData.datasets['vdW14LTG2'].interior(10.0**logMst2, 10.0**loghalfMassStarsz2)]) )
            distancesMRE3 = np.min( np.abs( [ observationalData.datasets['vdW14ETG3'].distance(10.0**logMst3, 10.0**loghalfMassStarsz3), observationalData.datasets['vdW14LTG3'].interior(10.0**logMst3, 10.0**loghalfMassStarsz3)]) )

            # Label points that fit Broeils
            distancesHIR = observationalData.datasets['broeils97'].distance(10.0**(loggasToStellarRatioHIz0+logMst0),  10.0**logbroeilsHIz0  )

            distancesC82 = observationalData.datasets['Dutton09All'].distance(10.0**logMst0, 10.0**logc82z0) 




            lnlik -= distancesSMHM0**2.0 / 2.0 
            lnlik -= distancesSMHM1**2.0 / 2.0
            lnlik -= distancesSMHM2**2.0 / 2.0
            lnlik -= distancesSMHM3**2.0 / 2.0

            lnlik -= distancesMS0**2.0 / 2.0
            lnlik -= distancesMS1**2.0 / 2.0
            lnlik -= distancesMS2**2.0 / 2.0
            lnlik -= distancesMS3**2.0 / 2.0

            lnlik -= distancesSTMZR0 **2.0 / 2.0

            lnlik -= distancesGMZR0**2.0 / 2.0
            lnlik -= distancesGMZR1**2.0 / 2.0
            lnlik -= distancesGMZR2 **2.0 / 2.0
            lnlik -= distancesGMZR3 **2.0 / 2.0

            lnlik -= distancesHMOL0**2.0 / 2.0
            lnlik -= distancesHMOL1**2.0 / 2.0
            lnlik -= distancesHMOL2**2.0 / 2.0
            lnlik -= distancesHMOL3**2.0 / 2.0

            lnlik -= distancesHI**2.0 / 2.0

            lnlik -= distancesSIG10**2.0 / 2.0
            lnlik -= distancesSIG11**2.0 / 2.0
            lnlik -= distancesSIG12**2.0 / 2.0
            lnlik -= distancesSIG13**2.0 / 2.0

            lnlik -= distancesMRE0**2.0 / 2.0
            lnlik -= distancesMRE1**2.0 / 2.0
            lnlik -= distancesMRE2**2.0 / 2.0
            lnlik -= distancesMRE3**2.0 / 2.0

            lnlik -= distancesHIR**2.0 / 2.0

            lnlik -= distancesC82**2.0 / 2.0
        return lnlik


    def lnprior(emceeparams):
        raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC = emceeparams
        #accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalizat     ion, kappaMassScaling = emceeparams
        accum = 0.0
    
        #varRedFac = 10000.0
        varRedFac = 1.0
    
        accum += lnlognormaldensity( raccRvir, np.log(0.141), np.log(3.0)**2.0 )
        accum += lnlognormaldensity( rstarRed, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( rgasRed, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( fg0mult, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
        accum += lnnormaldensity( muColScaling, 0.6, 1.0**2.0 / varRedFac )
        accum += lnnormaldensity( muFgScaling, -0.1, 1.0**2.0 / varRedFac )
        accum += lnlognormaldensity( muNorm, np.log(1.0), np.log(10.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( ZIGMfac, np.log(1.0), np.log(3.0)**2.0/ varRedFac )
        accum += lnbetadensity( zmix, 1.0, 1.0 ) # not accurate
        accum += lnlognormaldensity( eta, np.log(1.5), np.log(2.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( Qf, np.log(1.5), np.log(2.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( alphaMRI, np.log(0.05), np.log(2.0)**2.0/ varRedFac )
        accum += lnlognormaldensity( epsquench, np.log(1.0e-3), np.log(10.0)**2.0/ varRedFac )
        accum += lnbetadensity( accCeiling, 1.0, 1.0 ) # not accurate
        accum += lnnormaldensity( conRF, 0.3, 0.3**2.0 / varRedFac )
        accum += lnlognormaldensity( kZ, np.log(1.0), np.log(3.0)**2.0/ varRedFac )
        accum += lnbetadensity( xiREC, 1.0, 2.0 ) # not accurate
    
    
        if not np.isfinite(accum):
            return -np.inf
        return accum


    def samplefromprior():
        # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormaliza     tion, kappaMassScaling = emceeparams
        # Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsqu     ench, accCeiling, conRF, kZ, xiREC = emceeparams

        #varRedFac = 10000.0
        varRedFac = 1.0

        return [ \
            samplefromlognormaldensity( np.log(0.141), np.log(3.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
            samplefromnormaldensity( 0.6, 1.0**2.0 / varRedFac ),
            samplefromnormaldensity( -0.1, 1.0**2.0 / varRedFac ),
            samplefromlognormaldensity( np.log(1.0), np.log(10.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(1.0), np.log(3.0)**2.0/ varRedFac ),
            samplefrombetadensity( 1.0 * varRedFac, 1.0 * varRedFac),
            samplefromlognormaldensity( np.log(1.5), np.log(2.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(1.5), np.log(2.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(0.05), np.log(2.0)**2.0/ varRedFac),
            samplefromlognormaldensity( np.log(1.0e-3), np.log(10.0)**2.0/ varRedFac),
            samplefrombetadensity( 1.0 * varRedFac, 1.0 ),
            samplefromnormaldensity( 0.3, 0.3**2.0 / varRedFac ),
            samplefromlognormaldensity( np.log(1.0), np.log(3.0)**2.0/ varRedFac),
            samplefrombetadensity( 1.0, 2.0*varRedFac ) ]

    def lnprob(emceeparams):
        pr = lnprior(emceeparams)
        if np.isfinite(pr):
            return lnlikelihood(emceeparams) + pr
            #return constlikelihood(emceeparams) + pr
        return pr
    
    ndim, nwalkers = 17, 64 
    nsteps = 300
    p0 = [ samplefromprior() for w in range(nwalkers) ]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
    sampler.run_mcmc(p0, nsteps) 
    chain = sampler.chain # nwalkers x nsteps x ndim
    samples = chain[:,150:,:].reshape((-1,ndim))

    fig,ax = plt.subplots((4,5))
    for q in range(ndim):
        for r in range(nwalkers):
            ax.flatten()[q].plot( range(nsteps),chain[r,:,q], c=(np.random.random(), np.random.random(), np.random.random()) )
    plt.savefig('mcmc_fake_trace.png')
    plt.close(fig)

    import corner
    fig = corner.corner(samples, labels=[ 'raccRvir', 'rstarRed', 'rgasRed', 'fg0mult', 'muColScaling', 'muFgScaling', 'muNorm', 'ZIGMfac', 'zmix', 'eta', 'Qf', 'alphaMRI', 'epsquench', 'accCeiling', 'conRF', 'kZ', 'xiREC']) 
    fig.savefig("mcmc_fake_corner.png") 

def appendLabels():
    arr = np.loadtxt('broad_to_fail.txt')
    Mh = arr[:,19]
    mstar = arr[:,20]
    sSFR = arr[:,21]
    sfZ = arr[:,22]
    stZ = arr[:,23]
    gasToStellarRatioH2 = arr[:,24]
    gasToStellarRatioHI = arr[:,25]
    halfMassStars = arr[:,26]
    vPhi22 = arr[:,27] 
    c82  = arr[:,28] 
    Sigma1  = arr[:,29] 
    specificJStars  = arr[:,30] 
    metallicityGradientR90  = arr[:,31] 
    maxsig  = arr[:,32] 
    mdotBulgeG  = arr[:,33] 
    fractionGI  = arr[:,34] 
    tdep = arr[:,35] 
    tDepH2  = arr[:,36] 
    broeilsHI = arr[:,37] 
    MHI = gasToStellarRatioHI*mstar

    # Label points that fit Moster SMHM
    successesSMHM = observationalData.datasets['Moster10z0'].interior(Mh,mstar)
    labels = np.zeros(len(Mh))
    labels[successesSMHM] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit Brinchmann04 
    successesMS = observationalData.datasets['Brinchmann04Specz0'].interior(mstar, sSFR)
    labels = np.zeros(len(Mh))
    labels[successesMS] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit stellar MZR 
    successesSTMZR = np.logical_or( observationalData.datasets['kirby13'].interior(mstar, stZ), observationalData.datasets['gallazi05'].interior(mstar, stZ) )
    labels = np.zeros(len(Mh))
    labels[successesSTMZR] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit gas MZR 
    successesGMZR = np.logical_or( observationalData.datasets['Tremonti04'].interior(mstar, sfZ), observationalData.datasets['Lee06'].interior(mstar, sfZ) )
    labels = np.zeros(len(Mh))
    labels[successesGMZR] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit gas fractions
    successesHMOL = np.logical_or( observationalData.datasets['genzel15COz0'].interior(mstar, gasToStellarRatioH2), observationalData.datasets['genzel15Dustz0'].interior(mstar, gasToStellarRatioH2) )
    labels = np.zeros(len(Mh))
    labels[successesHMOL] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit gas fractions
    successesHI = np.logical_or( observationalData.datasets['papastergis12'].interior(mstar, gasToStellarRatioHI), observationalData.datasets['peeples11'].interior(mstar, gasToStellarRatioHI) )
    labels = np.zeros(len(Mh))
    labels[successesHI] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit structural relations
    successesSIG1 = np.logical_or( observationalData.datasets['Fang13'].interior(mstar, Sigma1), observationalData.datasets['Barro15HS'].interior(mstar, Sigma1) )
    labels = np.zeros(len(Mh))
    labels[successesSIG1] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit galaxy sizes
    successesMRE = np.logical_or( observationalData.datasets['vdW14ETG0'].interior(mstar, halfMassStars), observationalData.datasets['vdW14LTG0'].interior(mstar, halfMassStars) )
    labels = np.zeros(len(Mh))
    labels[successesMRE] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T

    # Label points that fit Broeils
    successesHIR = observationalData.datasets['broeils97'].interior(MHI, broeilsHI)
    labels = np.zeros(len(Mh))
    labels[successesHIR] = 1
    ne = np.shape(arr)[1]
    arr = np.vstack([arr.T, labels]).T


    arr = np.vstack([arr.T, MHI]).T # add in a column giving us the HI mass explicitly

    totalSuccess = np.logical_and( successesMRE, np.logical_and( successesHIR, np.logical_and( successesSIG1, np.logical_and( successesHI, np.logical_and( successesHMOL, np.logical_and( successesGMZR, np.logical_and( successesSTMZR, np.logical_and( successesMS, successesSMHM ))))))))
    labels = np.zeros(len(Mh))
    labels[totalSuccess] = 1
    arr = np.vstack([arr.T, labels]).T

    with open('broad_to_fail.txt', 'r') as f:
        header = f.readline()
    header = header[2:-1] + 'smhm ms stmzr gmzr hmol hi sigmacen mre broeils mhi totsuccess '
    np.savetxt('broad_to_fail_labeled.txt', arr, header=header)


def ksTestAnalysis():
    import scipy.stats
    with open('broad_to_fail_labeled.txt', 'r') as f:
        header = f.readline()[2:]
    arr = np.loadtxt('broad_to_fail_labeled.txt')
    def oneComparison(successColumn, massColumn, paramColumn, lowerMassRange, upperMassRange):
        assert lowerMassRange >=0 and lowerMassRange < upperMassRange and upperMassRange<=1.0
        successes =  arr[:,successColumn]==1
        masses = arr[:,massColumn]
        minx = np.min(masses[successes])
        maxx = np.max(masses[successes])
        massRange = [minx*(maxx/minx)**lowerMassRange, minx*(maxx/minx)**upperMassRange]
        selection = np.logical_and(successes, np.logical_and( masses>=massRange[0], masses<=massRange[1] ))
        modelsRun = arr[:,20]>1.0 # if mstar has a non-zero value the model was actually run. selection should be guaranteed to be a subset of modelsRun
        ksstat, pval = scipy.stats.ks_2samp( arr[:, paramColumn][modelsRun], arr[:,paramColumn][selection] ) # compare the full sample of the variable to the subset picked out by its match to a certain observational dataset AND a certain mass range
        return pval
    logs = [ 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0]  
    inputVariables = range(1,18) # the first 18 columns are the inpute variables. First one is halo mass at z=0
    successColumns = [38, 39, 40, 41, 42, 43, 44, 45, 46, 48]
    #massColumns = [19, 20, 20, 20, 20, 20, 20, 20, 47, 19 ] # the columns of the x-variable in the observational relation corresponding to the success column, i.e. the column telling us whether this model falls in the acceptable range of this observation
    massColumns = [19]*10 ### guarantees a more even number of groups
    pvalMatrix = np.zeros(( len(inputVariables)*3, len(successColumns) ))
    for i in range(np.shape(pvalMatrix)[0]):
        for j in range(np.shape(pvalMatrix)[1]):
            pvalMatrix[i,j] = oneComparison( successColumns[j], massColumns[j], inputVariables[(i - i%3)/3], (i%3)*(1.0/3.0), (i%3 +1)*(1.0/3.0) )
    headerList = header.split()
    tableString = '    '
    for sc in successColumns:
        tableString += headerList[sc] + ' '
    print tableString
    for i in range(np.shape(pvalMatrix)[0]):
        # set up a string for this line of the table
        tableString = headerList[ inputVariables[(i-i%3)/3] ]+ ' '+"{0:.2g}".format( (i%3)*(1.0/3.0))+'-'+ "{0:.2g}".format( (i%3+1)*(1.0/3.0))+' '
        for j in range(np.shape(pvalMatrix)[1]):
            tableString += "{0:.2g}".format(pvalMatrix[i,j])
            if pvalMatrix[i,j]<0.05/(np.shape(pvalMatrix)[0]*np.shape(pvalMatrix)[1]):
                tableString += '** '
            elif pvalMatrix[i,j]<0.05/3.0:
                tableString += '* '
            else:
                tableString += ' '
        print tableString

def filterArrayRows(array, valid):
    new_arr = []
    for i in range(len(valid)):
        if valid[i]:
            new_arr.append(array[i,:])
    return np.array(new_arr)

def readData(trainFrac=0.5, validateFrac=0.4):
    ''' trainFrac denotes what fraction of the sample should be used for training 
        validateFrac denotes what fraction of the sample should be used for validation
        1 - trainFrac - validateFrac is the out-of-sample test.'''
    arr = np.loadtxt('broad03_to_lasso.txt')
    with open('broad03_to_lasso.txt', 'r') as f:
        header = f.readline()[2:]
    labels = header.split()[nvars+1:nvars+1+len(logPreds)*nz+nRadii*len(logRadials)]
    X1 = arr[:,:nvars] # the 18 input variables we varied
    X2 = arr[:, nvars+1+len(logPreds)*nz+nRadii*len(logRadials):] # the random factors that define the accretion history
    for i in range(len(logVars)):
        # Take the log of the input variable if it makes sense to do so
        if logVars[i] == 1:
            X1[:,i] = np.log10(X1[:,i])
    # Add in squares and cross-terms for every pair of input parameters
    nxv = np.shape(X1)[1]
    for i in range(nxv):
        for j in range(nxv):
            # don't repeat combinations, only look at squares instead of cross-terms, or cross-terms with mass alone
            if j>=i and (i==j or i==0):
                to_add = X1[:,i]*X1[:,j]
                X1 = np.vstack( [ X1.T, to_add ] ).T

    X = np.vstack( [ X1.T, X2.T ] ).T # mc params + random factors

    #nxvOrig = np.shape(X)[1]
    Ys = arr[:,nvars+1:nvars+1+len(logPreds)*nz+nRadii*len(logRadials)]
    # Mh mstar sSFR sfZ stZ gasToStellarRatioH2 gasToStellarRatioHI halfMassStars vPhi22 c82 Sigma1 specificJStars metallicityGradientR90 maxsig mdotBulgeG fractionGI tdep tDepH2 broeilsHI, mStellarHalo
    valid = Ys[:,0]>0
    X =  filterArrayRows( X, valid )
    Ys = filterArrayRows( Ys, valid)
    ns = np.shape(arr)[0]
    nv = np.shape(X)[1]
    trainingSubset = int(ns*trainFrac) # 
    validateSubset = int(ns*(validateFrac+trainFrac))
    X_train = X[:trainingSubset, :]
    print "Training set shape: ",np.shape(X_train)
    Ys_train = Ys[:trainingSubset, :]
    X_validate = X[trainingSubset:validateSubset, :]
    Ys_validate = Ys[trainingSubset:validateSubset, :]
    X_test = X[validateSubset: , :]
    Ys_test = Ys[validateSubset: , :]
#    for i in range(len(logVars)):
#        to_replace = X_train[:,i]
#        if logVars[i]==1:
#            to_replace = np.log10(to_replace)
#        minx = np.min(to_replace)
#        maxx = np.max(to_replace)
#        X_train[:,i] = (to_replace-minx)/(maxx-minx)
#        if logVars[i]==1:
#            X_test[:,i] = (np.log10(X_test[:,i]) - minx)/(maxx-minx)
#        else:
#            X_test[:,i] = (X_test[:,i] - minx)/(maxx-minx)
    # Add terms \propto x1*x2
    #for i in range(nxv):
    #    for j in range(nxv):
    #        if j>=i and i==0 :
    #            if (j>nvars and j%3==0) or j<=nvars:
    #                to_add = X_train[:,i]*X_train[:,j]
    #                X_train = np.vstack([X_train.T, to_add]).T
    #                to_add = X_test[:,i]*X_test[:,j]
    #                X_test = np.vstack([X_test.T, to_add]).T

    for j in range(len(logPreds)):
        if logPreds[j] == 1:
            print "Taking log for ", labels[j*nz+0]
            for k in range(nz):
                # mStellarHalo is sometimes identically zero which is bad for taking the log!
                if j==19:
                    Ys_train[:,j*nz+k] += 1
                    Ys_validate[:,j*nz+k] += 1
                    Ys_test[:,j*nz+k] += 1
                Ys_train[:,j*nz+k] = np.log10(Ys_train[:,j*nz+k])
                Ys_validate[:,j*nz+k] = np.log10(Ys_validate[:,j*nz+k]) 
                Ys_test[:,j*nz+k] = np.log10(Ys_test[:,j*nz+k])
    for i in range(len(logRadials)):
        if logRadials[i]==1:
            for j in range( len(logPreds)*nz + i*nRadii, len(logPreds)*nz + (i+1)*nRadii ):
                Ys_train[:,j] = np.log10(Ys_train[:,j])
                Ys_validate[:,j] = np.log10(Ys_validate[:,j])
                Ys_test[:,j] = np.log10(Ys_test[:,j])

    return X_train, X_validate, X_test,  Ys_train, Ys_validate, Ys_test, labels


def learnLassoRidge(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, alphaLasso, alphaRidge, exponentiateRidge=False):
    errors_train = np.zeros(len(labels))
    errors_validate = np.zeros(len(labels))
    errors_test = np.zeros(len(labels))

    ridge_coefs = np.zeros( (len(labels), 1000) )
    lasso_coefs = np.zeros( (len(labels), np.shape(X_train)[1]-1000) )
    theModel = np.empty( len(labels), dtype=object)
    for k in range(len(labels)):
        print "Fitting k=",k,labels[k]
        invalid = Ys_train[:,k]!=Ys_train[:,k]
        if np.any(invalid):
            pdb.set_trace()
        if np.any(np.logical_not(np.isfinite(X_train))):
            pdb.set_trace()

        #regOLS = linear_model.LinearRegression()
        regLasso = linear_model.Lasso(alpha=alphaLasso, selection='random', tol=1.0e-5, max_iter=1e6, normalize=True)
        X_train_inter = X_train[:, :-1000] #np.vstack([ X_train[:, :-1000].T, X_train[:, -1]]).T
        X_validate_inter = X_validate[:, :-1000] # np.vstack([ X_validate[:, 0:18].T, X_validate[:, -1]]).T
        X_test_inter = X_test[:, :-1000]
        try:
            regLasso.fit(  X_train_inter, Ys_train[:,k])
        
            # Fit the physics params w/ Lasso. Now compute residuals!
            Y_train_inter = regLasso.predict( X_train_inter) - Ys_train[:,k]

            regRidge = linear_model.Ridge(alpha=alphaRidge, tol=1.0e-5, max_iter=1e6, normalize=True)
            X_train_ridge = X_train[:, -1000:]
            if exponentiateRidge:
                X_train_ridge = np.exp(X_train_ridge)
            regRidge.fit(X_train_ridge, Y_train_inter)

            final_residuals = regRidge.predict(X_train_ridge) - Y_train_inter
            errors_train[k] = np.std(final_residuals)

            
            Y_validate_inter = regLasso.predict( X_validate_inter) - Ys_validate[:,k]
            Y_test_inter = regLasso.predict( X_test_inter) - Ys_test[:,k]
            X_validate_ridge = X_validate[:, -1000:]
            X_test_ridge = X_test[:, -1000:]
            if exponentiateRidge:
                X_validate_ridge = np.exp(X_validate_ridge)
                X_test_ridge = np.exp(X_test_ridge)
            validation_residuals = regRidge.predict(X_validate_ridge) - Y_validate_inter
            test_residuals = regRidge.predict(X_test_ridge) - Y_test_inter
            errors_validate[k] = np.std(validation_residuals)
            errors_test[k] = np.std(test_residuals)
            ridge_coefs[k, : ] = regRidge.coef_[:]
            lasso_coefs[k, : ] = regLasso.coef_[:]
            theModel[k] = lassoRidge( regLasso, regRidge, exponentiateRidge)
        except:
            print "Failed for ", labels[k]
            errors_validate[k] = 1.0e10
            errors_train[k] = 1.0e10


    return errors_train, errors_validate, errors_test, labels, lasso_coefs, ridge_coefs, test_residuals, theModel


def learnLasso(X_train, X_test, Ys_train, Ys_test, alpha=0.1):



#    Ymins = []
#    Ymaxs = []
#    for j in range(len(logPreds)):
#        to_replace = Ys_train[:,j]
#        if logPreds[j]==1:
#            to_replace = np.log10(to_replace)
#        miny = np.min(to_replace)
#        maxy = np.max(to_replace)
#        Ys_train[:,j] = (to_replace-miny)/(maxy-miny)
#        Ymins.append(miny)
#        Ymaxs.append(maxy)
#        if logPreds[j]==1:
#            Ys_test[:,j] = (np.log10(Ys_test[:,j]) - miny)/(maxy-miny)
#        else:
#            Ys_test[:,i] = (Ys_test[:,i] - miny)/(maxy-miny)
    errors_train = np.zeros(len(labels))
    errors_test = np.zeros(len(labels))
    for k in range(len(labels)):
        regOLS = linear_model.LinearRegression()
        regLasso = linear_model.Lasso(alpha=alpha, selection='random', tol=1.0e-5, max_iter=1e6, normalize=True)
        regRidge = linear_model.Ridge(alpha=alpha, tol=1.0e-5, max_iter=1e6, normalize=True)
        invalid = Ys_train[:,k]!=Ys_train[:,k]
        if np.any(invalid):
            pdb.set_trace()
        regLasso.fit(X_train, Ys_train[:,k])
        regRidge.fit(X_train, Ys_train[:,k])
        regOLS.fit(X_train, Ys_train[:,k])
        Y_pred_train = regLasso.predict(X_train)
        Y_pred_test = regLasso.predict(X_test)
        Y_pred_train_ols = regOLS.predict(X_train)
        Y_pred_test_ols = regOLS.predict(X_test)
        Y_pred_train_ridge = regRidge.predict(X_train)
        Y_pred_test_ridge = regRidge.predict(X_test)
        errors_train[k] = np.std(Y_pred_train - Ys_train[:,k])
        errors_test[k] = np.std(Y_pred_test - Ys_test[:,k])
        print labels[k], '--', np.std(Y_pred_train - Ys_train[:,k]), np.std(Y_pred_test - Ys_test[:,k]), '--', np.std(Y_pred_train_ols - Ys_train[:,k]), np.std(Y_pred_test_ols - Ys_test[:,k]), '--', np.std(Y_pred_train_ridge - Ys_train[:,k]), np.std(Y_pred_test_ridge - Ys_test[:,k])
        print "L1 norm for lasso, ridge fit: ", np.sum(np.abs(regLasso.coef_)), np.sum(np.abs(regRidge.coef_))
        print "Number of coeffs used in lasso/ridge fit: ", np.sum( np.ones(len(regLasso.coef_))[np.abs(regLasso.coef_)>1.0e-4] ),  np.sum( np.ones(len(regRidge.coef_))[np.abs(regRidge.coef_)>1.0e-4] ), "of", len(regLasso.coef_)
        print "------"
        #except:
        #    print "Failed for k=",k, labels[k]
        #    print " "
    return errors_train, errors_test, labels
    

def learnSuccesses(C, gamma):
    arr = np.loadtxt('broad02_to_lasso.txt')
    print "Read in array of shape: ", np.shape(arr)
    #X = np.hstack( [arr[:,:18], arr[:,19:]] )
    X1 = arr[:,:18]
    X2 = arr[:, 19+len(logPreds)+nRadii*len(logRadials):]
    X = np.vstack( [ X1.T, X2.T ] ).T # mc params + random factors
    #X = arr[:,:18]
    # Take an output quantity, and check whether it's zero or non-zero. zero means we weren't able to read it in, so the run failed.
    Y = np.clip(arr[:,19],0, 1)
    Mh = arr[:,19]
    mstar = arr[:,20]
    #Y = observationalData.datasets['Moster10z0'].interior(Mh,mstar)
    ns = np.shape(arr)[0]
    nv = np.shape(X)[1]
    # Split the sample into a training set and a test set
    trainingSubset = ns*1/2 # 
    X_train = X[:trainingSubset, :]
    Y_train = Y[:trainingSubset]
    X_test = X[trainingSubset:, :]
    Y_test = Y[trainingSubset:]

    print "Proportion of successes in training, test sets: ", float(np.sum(Y_train))/len(Y_train), float(np.sum(Y_test))/len(Y_test)

    # Normalize the data. Rescale parameters to lie between 0 and 1. If the parameter should naturally be logarithmic, take the logarithm first. Subject the test set to the same transformation (taking no additional info about the test set into account!!!)
    logs = logVars + [0]*(nv-18)
    for i in range(len(logs)):
        to_replace = X_train[:,i]
        if logs[i]==1:
            to_replace = np.log10(to_replace)
        minx = np.min(to_replace)
        maxx = np.max(to_replace)
        X_train[:,i] = (to_replace-minx)/(maxx-minx)
        if logs[i]==1:
            X_test[:,i] = (np.log10(X_test[:,i]) - minx)/(maxx-minx)
        else:
            X_test[:,i] = (X_test[:,i] - minx)/(maxx-minx)

    clf = svm.SVC(kernel='rbf',gamma=gamma, C=C, class_weight='balanced')
    clf.fit(X_train, Y_train)

    # Test the result:
    y_predict_training = clf.predict( X_train )
    n_successes_training = np.sum(np.ones(len(y_predict_training))[ y_predict_training==Y_train ])
    in_sample_error = 1.0 - float(n_successes_training)/float(len(y_predict_training))
    print "In the training set:"
    print n_successes_training, " successes out of ", len(y_predict_training)
    print float(n_successes_training)/float(len(y_predict_training))
    y_predict = clf.predict( X_test )
    n_successes = np.sum(np.ones(len(y_predict))[ y_predict==Y_test ])
    out_sample_error = 1.0 - float(n_successes)/float(len(y_predict))
    print "In the test set:"
    print n_successes, " successes out of ", len(y_predict)
    print float(n_successes)/float(len(y_predict))
    return in_sample_error, out_sample_error

def validateSVM():
    Cs = [0.1, 1.0, 10.0, 100.0, 1000, 10000]
    ls=['-','--','-.',':']*5
    cs=['k','r','b','green', 'orange', 'purple', 'lightblue', 'grey']*5
    gammas = [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
    fig,ax = plt.subplots()
    for i, C in enumerate(Cs):
        errors_in = []
        errors_out = []
        for j, gamma in enumerate(gammas):
            ei, eo = learnSuccesses(C,gamma)
            errors_in.append( ei )
            errors_out.append( eo )
        ax.plot(gammas, errors_in, lw=2, ls=ls[i], c=cs[i])
        ax.plot(gammas, errors_out, lw=4, ls=ls[i], label='C='+str(C), c=cs[i])
    ax.set_xlabel(r'$\gamma$')
    ax.set_xscale('log')
    ax.set_ylabel(r'Validation Error')
    ax.legend(loc=0, frameon=False)
    plt.savefig('svm_validate_successes.png')
    plt.close(fig)

def validateLassoRidge():
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.15, validateFrac=0.4)

    alphasLasso = np.power(10.0, np.linspace(-4, -1, 5))
    alphasRidge = np.power(10.0, np.linspace(-4, -1, 6))

    errors_train = np.zeros((len(labels), len(alphasLasso), len(alphasRidge), 2))
    errors_validate = np.zeros((len(labels), len(alphasLasso), len(alphasRidge), 2))
    errors_test = np.zeros((len(labels), len(alphasLasso), len(alphasRidge), 2))
    residuals_test = np.zeros((len(labels), len(alphasLasso), len(alphasRidge), 2))
    model_array = np.empty((len(labels), len(alphasLasso), len(alphasRidge), 2), dtype=object)
    #ridge_coefs = np.zeros( (len(labels), np.shape(X_train)[1]-19) )
    #lasso_coefs = np.zeros( (len(labels), 19) )
    lasso_coefs = np.zeros((len(labels), np.shape(X_train)[1]-1000, len(alphasLasso), len(alphasRidge), 2))
    ridge_coefs = np.zeros((len(labels), 1000, len(alphasLasso), len(alphasRidge), 2))
    print "ridge_coefs shape: ", np.shape(ridge_coefs)
    print "lasso_coefs shape: ", np.shape(lasso_coefs)
    minIs = np.zeros( len(labels) )
    minJs = np.zeros( len(labels) )
    minQs = np.zeros( len(labels) )
    minValErrors = np.zeros( len(labels) ) + 1.0e30
    for i, aL in enumerate(alphasLasso):
        for j, aR in enumerate(alphasRidge):
            for q, tf in enumerate([True,False]):
                errors_train_this, errors_validate_this, errors_test_this, labels_this, lasso_coef_this, ridge_coef_this, test_residuals_this, theModel = learnLassoRidge(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, aL, aR, exponentiateRidge=tf)
                errors_train[:,i,j,q] = errors_train_this[:]
                errors_validate[:,i,j,q] = errors_validate_this[:]
                errors_test[:,i,j,q] = errors_test_this[:]
                lasso_coefs[:,:, i,j,q] = lasso_coef_this[:,:]
                ridge_coefs[:,:, i,j,q] = ridge_coef_this[:,:]
                model_array[:,i,j,q] = theModel[:]
    for i in range(len(alphasLasso)):
        for j in range(len(alphasRidge)):
            for q in range(2):
                for k in range(len(labels)):
                    if errors_validate[k,i,j,q] < minValErrors[k]:
                        minValErrors[k] = errors_validate[k,i,j,q]
                        minIs[k] = i
                        minJs[k] = j
                        minQs[k] = q
    print "DependentVar, alphaLasso, alphaRidge, #Lasso, testErr, validationErr, trainingErr, avgRidgeCoef, stdRidgeCoef"
    models = []
    for k in range(len(labels)):
        # For each dependent variable, show number of lasso coeffs, minimum validation error, corresponding training error, mean/std of ridge coeffs.
        print labels[k], alphasLasso[minIs[k]], alphasRidge[minJs[k]], minQs[k], np.sum(np.ones(len(lasso_coefs[0,:,0,0,0]))[np.abs(lasso_coefs[k,:,minIs[k], minJs[k],0])>1.0e-4]), errors_test[k,minIs[k],minJs[k],minQs[k]], minValErrors[k], errors_train[k,minIs[k], minJs[k], minQs[k]], np.mean(ridge_coefs[k, :, minIs[k], minJs[k], minQs[k]]), np.std(ridge_coefs[k, :, minIs[k], minJs[k], minQs[k]]) #, lasso_coefs[k,:,minIs[k], minJs[k]]
        models.append( model_array[k, minIs[k], minJs[k], minQs[k]] )

    pickle.dump(models, open('validatedModels.pickle', 'w'))
    return models 


def validateLasso():
    alphas = np.power(10.0, np.linspace(-9,-1.5,30))
    npar = 19
    cs=[ (1.0 - float(i)/npar/2.0, 0.5, 0.5+float(i)/npar/2.0) for i in range(npar)]
    fig,ax = plt.subplots()
    errors_train = np.zeros((npar, len(alphas)))
    errors_test = np.zeros((npar, len(alphas)))
    labels=None
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData() 
    for i, a in enumerate(alphas):
        errors_train_this, errors_test_this, labels_this = learnLasso(X_train, X_test, X_validate, Ys_validate, Ys_train, Ys_test, labels, alpha=a)
        errors_train[:,i] = errors_train_this[:]
        errors_test[:,i] = errors_test_this[:]
        labels=labels_this
    for j in range(npar):
        ax.plot(alphas, errors_train[j,:], c=cs[j], ls='--', lw=2)
        ax.plot(alphas, errors_test[j,:], c=cs[j], ls='-', lw=2, label=labels[j])
    #ax.legend(loc=0)
    ax.set_xlabel(r'$\alpha$')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\sigma_{y-y_\mathrm{pred}}$')
    plt.savefig('lasso_validate.png')
    plt.close(fig)

if __name__=='__main__':
    #learn(1.0e4, 1.0e-5)
    #validateSVM()
    #appendLabels()
    #ksTestAnalysis()
    #learnLasso(alpha=1.0e-3)
    #validateLasso()
    #validateLassoRidge()

    runEmcee()



