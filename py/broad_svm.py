import numpy as np
import os
import copy
import pickle
import matplotlib.pyplot as plt
import pdb
import observationalData
import emcee
from sklearn import svm, linear_model, ensemble

#import pyqt_fit.nonparam_regression as smooth
#from pyqt_fit import npr_methods

# Useful information about the setup of the linear models and their fits....
logVars = [1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0]
logPreds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1]
logRadials = [1,1,1,1,1,1]
#logVars = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#logPreds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#logRadials = [0, 0, 0, 0, 0, 0]
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

class xtransform:
    def __init__(self, X_train):
        xps = X_train[:,:] # All physical parameters and 
        # now normalize ordinate variables so that we can sensibly evaluate distances
        self.minxps = np.min(xps, axis=0)
        self.maxxps = np.max(xps, axis=0)
    def transform(self, X_other):
        minxpsTiled = np.tile( self.minxps, (np.shape(X_other)[0],1))
        maxxpsTiled = np.tile( self.maxxps, (np.shape(X_other)[0],1))
        return (X_other - minxpsTiled)/(maxxpsTiled-minxpsTiled)
    def inverseTransform(self, x_transformed):
        minxpsTiled = np.tile( self.minxps, (np.shape(x_transformed)[0],1))
        maxxpsTiled = np.tile( self.maxxps, (np.shape(x_transformed)[0],1))
        return x_transformed*(maxxpsTiled-minxpsTiled) + minxpsTiled

def plotResiduals():
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.25, validateFrac=0.1, fn='broad05_to_lasso.txt')
    tra = xtransform(X_train)
    X_train = tra.transform(X_train)
    X_validate = tra.transform(X_validate)
    X_test = tra.transform(X_test)
    models = pickle.load( open( 'validatedRF_nolog.pickle', 'r' ) )
    for i, model in enumerate(models):
        Y_test_predict = model.predict(X_test)
        Y_test_actual = Ys_test[:,i]
        residuals = Y_test_predict-Y_test_actual

        fig,ax = plt.subplots(1,2)
        plt.hist( residuals, bins='auto' )
        ax[0].scatter(X_test[:,0], residuals)
        plt.savefig('residualsRF_'+labels[i]+'.png')
        plt.close(fig)

def ridgeCoeffsPlot():
    models = pickle.load( open( 'validatedModels.pickle', 'r' ) )

    fig,ax = plt.subplots()
    cs=[]
    for i in range(4):
        a,b,c = (np.random.random(), np.random.random(), np.random.random())
        cs.append( (a,b,c) )
        cs.append( (a,b,c) )
        cs.append( (a,b,c) )
        cs.append( (a,b,c) )
    for i, model in enumerate(models[:16]):
        ax.plot(model.sklRidge.coef_, c=cs[i])
    plt.savefig('ridge_coefs.png')
    plt.close(fig)


def fractionalVariancePlot(normalize=False):
    models = pickle.load( open( 'validatedRF_nolog.pickle', 'r' ) )
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.25, validateFrac=0.1, fn='broad05_to_lasso.txt')
    ## Transform the ordinates to (perhaps) improve behavior of the fitting algorithms
    tra = xtransform(X_train)
    X_train = tra.transform(X_train)
    X_validate = tra.transform(X_validate)
    X_test = tra.transform(X_test)

    Ytra = xtransform(Ys_train)
    Ys_train = Ytra.transform(Ys_train)
    Ys_validate = Ytra.transform(Ys_validate)
    Ys_test = Ytra.transform(Ys_test)


    lefts = []
    heights = []
    bottoms = []
    colors = []
    runningleft = -0.2
    for i, model in enumerate(models[:80]):
        #print "i= ",i," model type is ", type(model)
        if type(model) is type(None):
            pass
        else:
            #lass, ridg, unex, tot = fractionalVariances(model, X_test, Ys_test[:,i])
            tot = np.var(Ys_test[:,i])
            ridg= 0.0
            lass = np.var( model.predict(X_test) ) /tot
            unex = np.var( Ys_test[:,i] - model.predict(X_test) ) / tot

            actualTot = lass+ridg+unex
            if normalize:
                lass/=actualTot
                ridg/=actualTot
                unex/=actualTot
            if i%4 == 0:
                runningleft+=0.4
            else:
                runningleft+=0.2

            heights.append( lass )
            heights.append( ridg )
            heights.append( unex )
            bottoms.append( 0 )
            bottoms.append( lass )
            bottoms.append( (lass+ridg) )
            colors.append( 'blue' )
            colors.append( 'red' )
            colors.append( 'gray' )
            for j in range(3):
                lefts.append(runningleft)

    print "labels, lefts: ", len(labels), len(heights), len(lefts), len(models)
    fig,ax = plt.subplots(figsize=(24,6))
    ax.bar( lefts, heights, width=0.2, bottom=np.clip(bottoms, 0.01, np.inf), color=colors )
    for i, label in enumerate(labels[:80]):
        if i%4==0:
            ax.text(lefts[i*3], heights[i*3+2] + bottoms[i*3+2]+0.05, label)
    ax.xaxis.set_ticks([])
    #ax.set_yscale('log')
    #ax.set_ylim(0.0, 2.5)
    plt.axhline(1.0, c='k', ls='--')
    plt.savefig('variance_fractions_RF_nolog.png')
    plt.close(fig)


def fractionalVariances(model, X_test, Y_test):
    '''  Report the fraction of the variance which is explained by ridge, lasso, and neither. Also just tell us the total var.'''
    overallVariance = np.var(Y_test)
    X_lasso = X_test[:, :-1000]    
    X_ridge = X_test[:, -1000: ]
    if model.expon:
        X_ridge = np.exp(X_ridge)
    Y_lasso = model.sklLasso.predict( X_lasso )
    Y_ridge = model.sklRidge.predict( X_ridge )
    lassoVariance = np.var(Y_lasso)
    ridgeVariance = np.var(Y_ridge)
    unexplainedVariance = np.var(Y_test - Y_ridge - Y_lasso)

    return lassoVariance/overallVariance, ridgeVariance/overallVariance, unexplainedVariance/overallVariance, overallVariance



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

def lnlikelihood(emceeparams, models):
    # First transform the emceeparams into the same format used by 'X' in the fit of the linear models
    # emceeparams is the set of 18 parameters that we fit with Lasso, minus mass, so we have..
    nmh = 5
    MhGrid = np.power(10.0, np.linspace(10.5,14,nmh))
    lnlik = 0.0
    for k,Mh in enumerate(MhGrid):
        X1 = np.array([Mh]+list(emceeparams[:-1])).reshape((1,len(emceeparams[:-1])+1))
        for i in range(len(logVars)):
            # Take the log of the input variable if it makes sense to do so
            if logVars[i] == 1:
                X1[:,i] = np.log10(X1[:,i])
        Y_eval = np.array( [models[j].predictFill(X1)[0][j] for j in range(80)] ).reshape(1,80)
        lnlik += np.sum( globalLikelihood(Y_eval, fh=emceeparams[-1], returnlikelihood=True) )


    print "Returning lnlik = ", lnlik
    if not np.isfinite(lnlik):
        return -np.inf
    return lnlik
def samplefromprior( varRedFac=1.0):
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormaliza     tion, kappaMassScaling = emceeparams
    # Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsqu     ench, accCeiling, conRF, kZ, xiREC = emceeparams


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
        samplefrombetadensity( 1.0, 2.0*varRedFac ),
        samplefrombetadensity( 1.0 *varRedFac, 1.0*varRedFac) ]

def lnprob(emceeparams, models):
    pr = lnprior(emceeparams)
    if np.isfinite(pr):
        return lnlikelihood(emceeparams, models) + pr
        #return constlikelihood(emceeparams) + pr
    return pr

def lnprior(emceeparams):
    raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, fh = emceeparams
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
    accum += lnbetadensity( fh, 1.0, 1.0 )

    if not np.isfinite(accum):
        return -np.inf
    return accum


def saveRestart(fn,restart):
    with open(fn,'wb') as f:
        print "checkpoint to "+fn
        pickle.dump(restart,f,2)

def updateRestart(fn,restart):
    if os.path.isfile(fn):
        with open(fn,'rb') as f:
            tmp_dict = pickle.load(f)
            restart.update(tmp_dict)

### some functions to make plots of mcmc
def tracePlots(restart, fn, burnIn=0):
    chain = restart['chain']
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
            if np.random.uniform(0,1) <= 1.0: # print every one ----tenth-ish----
                ax[i,j].plot(chain[walker,burnIn:,dim],alpha=.3,ls='--')

    plt.tight_layout()
    plt.savefig(fn+'.png')
    plt.close(fig)

def probsPlots(restart, fn, burnIn=0):
    allProbs = restart['allProbs']
    
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

    try:
        acor = np.zeros(ndim)
        for walker in range(ndim):
            acor[walker] = emcee.autocorr.integrated_time(allProbs[walker,burnIn:] )#,  window=min([50,iters/2]))
        print "acor stats: "
        stats(acor)
    except:
        print "failed to compute acor."

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
    #print " autocorrelation lengths for each parameter: ",restart['acor']
    #stats(restart['acor'])
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

    import corner
    trifig = corner.corner(sampleRed, labels=labels, extents=extents)
    #trifigPrior = triangle.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, extents=extents)
    trifig.savefig(fn)



def runEmcee(mpi=False, continueRun=False):
    import sys
    if mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        from emcee.utils import MPIPool

    models = [ pickle.load( open( 'rfnt_'+str(k)+'_0.pickle', 'r' ) ) for k in range(80) ]
    #import bolshoireader
    #bolshoidir = '/Users/jforbes/bolshoi/'
    #globalBolshoiReader = bolshoireader.bolshoireader('rf_registry4.txt',3.0e11,3.0e14, bolshoidir)
    #globalBolshoiReader.storeAll()
    #bolshoiSize =  len(globalBolshoiReader.keys)



    if mpi:
        pool = MPIPool(comm=comm, loadbalance=True)
        if not pool.is_master():
            pool.wait()
            sys.exit()
    
    ndim, nwalkers = 18, 40 
    fn = 'fakemcmc_restart.pickle'
    restart = {}
    nsteps = 3000 # test run
    p0 = [ samplefromprior(varRedFac=1000.0) for w in range(nwalkers) ]

    restart['currentPosition'] = p0
    restart['chain' ] = None
    restart['state' ] = None
    restart['prob'] = None
    restart['iterationCounter'] = 0
    restart['mcmcRunCounter'] = 0

    if continueRun:
        updateRestart(fn, restart)

    restart['iterationCounter'] += nsteps
    restart['mcmcRunCounter'] += 1


    if mpi:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[models])
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[models])
    
    for result in sampler.sample(restart['currentPosition'], iterations=nsteps, lnprob0=restart['prob'], rstate0=restart['state']) :
        pos, prob, state = result

        #restart['acor'] = sampler.acor[:] # autocorr length for each param (ndim)
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




    #sampler.run_mcmc(p0, nsteps) 


    chain = sampler.chain # nwalkers x nsteps x ndim
    samples = chain[:,15:,:].reshape((-1,ndim))

    fig,ax = plt.subplots(4,5)
    for q in range(ndim):
        for r in range(nwalkers):
            ax.flatten()[q].plot( range(nsteps),chain[r,:,q], c=(np.random.random(), np.random.random(), np.random.random()) )
    plt.savefig('mcmc_fake_trace.png')
    plt.close(fig)

    import corner
    fig = corner.corner(samples, labels=[ 'raccRvir', 'rstarRed', 'rgasRed', 'fg0mult', 'muColScaling', 'muFgScaling', 'muNorm', 'ZIGMfac', 'zmix', 'eta', 'Qf', 'alphaMRI', 'epsquench', 'accCeiling', 'conRF', 'kZ', 'xiREC', 'fh']) 
    fig.savefig("mcmc_fake_corner.png") 

    if mpi:
        pool.close()


def singleRelationLikelihood(x,y,datasets):
    lnlik = np.zeros(len(x))
    distances = np.zeros((len(x),len(datasets)))
    sigmas = np.zeros((len(x),len(datasets)))
    finalResiduals = np.zeros(len(x))
    for i,ds in enumerate(datasets):
        distance, sigma = observationalData.datasets[ds].distance(x,y)
        distances[:,i] = distance
        sigmas[:,i] = sigma
    finite = np.isfinite(sigmas)
    mask = np.zeros(np.shape(sigmas))
    mask[finite] = 1 
    counts = np.sum(mask, axis=1)

    nods = counts==0
    oneds = counts==1
    nds = counts>1

    ## pick out the sigmas that are not infinity. These same pts will also have non-zero distances.
    lnlik[oneds] = np.log(1.0/np.sqrt(2.0*np.pi* np.power(np.min(sigmas[oneds,:],axis=1),2.0))) - np.power(np.sum(distances[oneds,:],axis=1),2.0)/2.0
    finalResiduals[oneds] = np.sum(distances[oneds,:],axis=1)

    lnlik[nds] = np.log( 1.0/counts[nds] * np.sum(  1.0/np.sqrt(2.0*np.pi* np.power(sigmas[nds,:],2.0))*np.exp( - np.power(distances[nds,:],2.0)/2.0), axis=1) )
    ndsnonzero = np.abs(distances[nds])>1.0e-5
    finalResiduals[nds] = np.mean(distances[nds,:],axis=1)  ### this is a little sketchy

    lnlik[np.logical_not(np.isfinite(lnlik))] = -1.0e20 # floor
    lnlik[lnlik<-1.0e20] = -1.0e20

    # because of the floor above, this trace should never be triggered...
    if np.any(np.logical_not(np.isfinite(lnlik))):
        pdb.set_trace()

    return lnlik, finalResiduals

def globalLikelihood(Ys_train, fh=0, returnlikelihood=True):
    logMh0 = Ys_train[:,0]
    logMh1 = Ys_train[:,1]
    logMh2 = Ys_train[:,2]
    logMh3 = Ys_train[:,3]
    logsSFRz0  = Ys_train[:,8]
    logsSFRz1 = Ys_train[:,9]
    logsSFRz2 = Ys_train[:,10]
    logsSFRz3 = Ys_train[:,11]
    logsfZz0 = Ys_train[:,12]
    logsfZz1 = Ys_train[:,13]
    logsfZz2 = Ys_train[:,14]
    logsfZz3 = Ys_train[:,15]
    logstZz0 = Ys_train[:,16]
    logstZz1 = Ys_train[:,17]
    logstZz2 = Ys_train[:,18]
    logstZz3 = Ys_train[:,19]
    loggasToStellarRatioH2z0 = Ys_train[:,20]
    loggasToStellarRatioH2z1 = Ys_train[:,21]
    loggasToStellarRatioH2z2  = Ys_train[:,22]
    loggasToStellarRatioH2z3  = Ys_train[:,23]
    loggasToStellarRatioHIz0  = Ys_train[:,24]
    loghalfMassStarsz0 = Ys_train[:,28] 
    loghalfMassStarsz1 = Ys_train[:,29] 
    loghalfMassStarsz2 = Ys_train[:,30] 
    loghalfMassStarsz3 =Ys_train[:,31] 
    logvPhi22z0  = Ys_train[:,32] 
    logvPhi22z1 = Ys_train[:,33] 
    logvPhi22z2  = Ys_train[:,34] 
    logvPhi22z3 = Ys_train[:,35] 
    logc82z0  = Ys_train[:,36] 
    logSigma1z0  = Ys_train[:,40] 
    logSigma1z1  = Ys_train[:,41] 
    logSigma1z2  = Ys_train[:,42] 
    logSigma1z3  = Ys_train[:,43] 
    logspecificJStarsz0  = Ys_train[:,44] 
    logspecificJStarsz1  = Ys_train[:,45] 
    logspecificJStarsz2  = Ys_train[:,46] 
    logspecificJStarsz3  = Ys_train[:,47] 
    logbroeilsHIz0  = Ys_train[:,72]
    logmStellarHaloz0  = Ys_train[:,76]
    logmStellarHaloz1  = Ys_train[:,77]
    logmStellarHaloz2  = Ys_train[:,78]
    logmStellarHaloz3 = Ys_train[:,79]
    #  logPreds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1]
    # corresponding to
    # Mhz0 mstarz0 sSFRz0 sfZz0 stZz0 gasToStellarRatioH2z0 gasToStellarRatioHIz0 halfMassStarsz0 vPhi22z0 c82z0 Sigma1z0 specificJStarsz0 metallicityGradientR90z0 maxsigz0 mdotBulgeGz0 fractionGIz0 tdepz0 tDepH2z0 broeilsHIz0 mStellarHaloz0 
    # i.e. only fractionGIz0 and metallicityGradientR90z0 are not logarithmic
    logMst0 = np.log10( 10.0**Ys_train[:,4] + fh* 10.0**logmStellarHaloz0 )
    logMst1 = np.log10( 10.0**Ys_train[:,5] + fh* 10.0**logmStellarHaloz1 )
    logMst2 = np.log10( 10.0**Ys_train[:,6] + fh* 10.0**logmStellarHaloz2 )
    logMst3 = np.log10( 10.0**Ys_train[:,7] + fh* 10.0**logmStellarHaloz3 )

    # Label points that fit Moster SMHM

    lnlik = np.zeros((len(Ys_train[:,0]),29) )

    srlInd=0
    if not returnlikelihood:
        srlInd=1
    lnlik[:,0] = singleRelationLikelihood(10.0**logMh0,10.0**logMst0,['Moster10z0'])[srlInd]
    lnlik[:,1] = singleRelationLikelihood(10.0**logMh1,10.0**logMst1,['Moster10z1'])[srlInd]
    lnlik[:,2] += singleRelationLikelihood(10.0**logMh2,10.0**logMst2,['Moster10z2'])[srlInd]
    lnlik[:,3] += singleRelationLikelihood(10.0**logMh3,10.0**logMst3,['Moster10z3'])[srlInd]

    lnlik[:,4] += singleRelationLikelihood(10.0**logMst0,10.0**logsSFRz0,['Brinchmann04Specz0'])[srlInd]
    lnlik[:,5] += singleRelationLikelihood(10.0**logMst1,10.0**logsSFRz1,['whitaker14Specz1'])[srlInd]
    lnlik[:,6] += singleRelationLikelihood(10.0**logMst2,10.0**logsSFRz2,['whitaker14Specz2'])[srlInd]
    lnlik[:,7] += singleRelationLikelihood(10.0**logMst3,10.0**logsSFRz3,['lilly13Specz3'])[srlInd]

    lnlik[:,8] += singleRelationLikelihood(10.0**logMst0,10.0**logstZz0,['gallazi05','kirby13'])[srlInd]

    lnlik[:,9] += singleRelationLikelihood(10.0**logMst0,10.0**logsfZz0,['Tremonti04','Lee06'])[srlInd]
    lnlik[:,10] += singleRelationLikelihood(10.0**logMst1,10.0**logsfZz1,['Genzel15Z1'])[srlInd]
    lnlik[:,11] += singleRelationLikelihood(10.0**logMst2,10.0**logsfZz2,['Genzel15Z2'])[srlInd]
    lnlik[:,12] += singleRelationLikelihood(10.0**logMst3,10.0**logsfZz3,['Genzel15Z3'])[srlInd]

    
    lnlik[:,13] += singleRelationLikelihood(10.0**logMst0,10.0**loggasToStellarRatioH2z0,['genzel15COz0','genzel15Dustz0'])[srlInd]
    lnlik[:,14] += singleRelationLikelihood(10.0**logMst1,10.0**loggasToStellarRatioH2z1,['genzel15COz1','genzel15Dustz1'])[srlInd]
    lnlik[:,15] += singleRelationLikelihood(10.0**logMst2,10.0**loggasToStellarRatioH2z2,['genzel15COz2','genzel15Dustz2'])[srlInd]
    lnlik[:,16] += singleRelationLikelihood(10.0**logMst3,10.0**loggasToStellarRatioH2z3,['genzel15COz3','genzel15Dustz3'])[srlInd]

    lnlik[:,17] += singleRelationLikelihood(10.0**logMst0,10.0**loggasToStellarRatioHIz0,['papastergis12','peeples11'])[srlInd]

    #lnlik += singleRelationLikelihood(10.0**logMst0,10.0**logSigma1z0,['Fang13','Barro15HS','Barro15HQ'])
    lnlik[:,18] += singleRelationLikelihood(10.0**logMst0,10.0**logSigma1z0,['Fang13'])[srlInd]
    lnlik[:,19] += singleRelationLikelihood(10.0**logMst1,10.0**logSigma1z1,['Barro151S','Barro151Q'])[srlInd]
    lnlik[:,20] += singleRelationLikelihood(10.0**logMst2,10.0**logSigma1z2,['Barro152S','Barro152Q'])[srlInd]
    lnlik[:,21] += singleRelationLikelihood(10.0**logMst3,10.0**logSigma1z3,['Barro153S','Barro153Q'])[srlInd]


    # Label points that fit galaxy sizes
    lnlik[:,22] += singleRelationLikelihood(10.0**logMst0,10.0**loghalfMassStarsz0,['vdW14ETG0','vdW14LTG0'])[srlInd]
    lnlik[:,23] += singleRelationLikelihood(10.0**logMst1,10.0**loghalfMassStarsz1,['vdW14ETG1','vdW14LTG1'])[srlInd]
    lnlik[:,24] += singleRelationLikelihood(10.0**logMst2,10.0**loghalfMassStarsz2,['vdW14ETG2','vdW14LTG2'])[srlInd]
    lnlik[:,25] += singleRelationLikelihood(10.0**logMst3,10.0**loghalfMassStarsz3,['vdW14ETG3','vdW14LTG3'])[srlInd]

    lnlik[:,26] += singleRelationLikelihood(10.0**(logMst0+loggasToStellarRatioHIz0),10.0**logbroeilsHIz0,['broeils97'])[srlInd]

    lnlik[:,27] += singleRelationLikelihood(10.0**logMst0,10.0**logc82z0,['Dutton09All'])[srlInd]

    lnlik[:,28] += singleRelationLikelihood(10.0**logMst0,10.0**logvPhi22z0,['miller11'])[srlInd]

    if returnlikelihood:
        return np.sum(lnlik,axis=1)
    else:
        return lnlik # in this case, these are the residuals



def nuclearSearch(Nbins = 7, Niter=10000, Ninits=50):
    X_train, X_validate, X_test,  Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=1.0, validateFrac=0.0, fn='broad06_to_lasso.txt')
    logMh0s= Ys_train[:,0]
    ## Set up Mh0 bins. These selections won't change throughout this process!
    MhEdges = np.linspace(np.min(logMh0s), np.max(logMh0s), Nbins+1)
    binSelections = []
    for i in range(Nbins):
        binSelections.append( np.logical_and(MhEdges[i]<logMh0s, logMh0s<MhEdges[i+1]) )

    xps = X_train[:,1:18] # The physical parameters excluding mass
    # now normalize ordinate variables so that we can sensibly evaluate distances
    minxps = np.min(xps, axis=0)
    maxxps = np.max(xps, axis=0)
    minxpsTiled = np.tile( minxps, (np.shape(xps)[0],1))
    maxxpsTiled = np.tile( maxxps, (np.shape(xps)[0],1))
    normxps = (xps - minxpsTiled)/(maxxpsTiled-minxpsTiled)
    
    globalLikelihoodArray = globalLikelihood(Ys_train[:,:])
    
    def run():
        sequenceOfXps = []
        sequenceOfNormXps = []
        sequenceOfLiks = []
        sequenceOfFhs = []
        normxp = np.random.random(17) #np.ones(17)/2.0
        fh = np.random.random()

        for i in range(Niter):
            if i%100==0:
                print "Running nuclear search iteration i=",i
            normxp =  np.clip( np.random.multivariate_normal( normxp, np.diag(np.ones(len(normxp)))*1.0e-3 ), 0,1)
            fh = np.clip( np.random.normal(fh,1.0e3), 0, 1 )

            distances = np.sum( np.abs( normxps - np.tile(normxp, (np.shape(normxps)[0],1)) ), axis=1)
            
            sampleInds = []

            for j in range( Nbins ):
                closestInd = np.argmin(distances[ binSelections[j] ]) 
                sampleInds.append( np.arange(len(distances))[binSelections[j]][closestInd] )

            lnlik = np.sum( globalLikelihood(Ys_train[sampleInds,:], fh) )
            #lnlik = np.sum(globalLikelihoodArray[sampleInds])

            if i==0:
                sequenceOfNormXps.append(normxp)
                sequenceOfXps.append(normxp*(maxxps-minxps)+minxps)
                sequenceOfLiks.append(lnlik)
                sequenceOfFhs.append(fh)
            else:
                alpha= lnlik-sequenceOfLiks[-1]
                if alpha>=0:
                    sequenceOfNormXps.append(normxp)
                    sequenceOfXps.append(normxp*(maxxps-minxps)+minxps)
                    sequenceOfLiks.append(lnlik)
                    sequenceOfFhs.append(fh)
                elif np.random.random() < np.exp(alpha):
                    sequenceOfNormXps.append(normxp)
                    sequenceOfXps.append(normxp*(maxxps-minxps)+minxps)
                    sequenceOfLiks.append(lnlik)
                    sequenceOfFhs.append(fh)
                else:
                    sequenceOfXps.append( sequenceOfXps[-1] )
                    sequenceOfNormXps.append( sequenceOfNormXps[-1] )
                    sequenceOfLiks.append( sequenceOfLiks[-1] )
                    normxp = sequenceOfNormXps[-1]
                    sequenceOfFhs.append( sequenceOfFhs[-1] )
                    fh = sequenceOfFhs[-1]
            if i==Niter-1:
                lastSetOfResiduals = globalLikelihood(Ys_train[sampleInds,:], fh, returnlikelihood=False)
        return sequenceOfLiks, sequenceOfNormXps, sequenceOfXps, sequenceOfFhs, lastSetOfResiduals


    likelihoods = np.zeros((Ninits, Niter))
    chains = np.zeros(( Ninits, Niter, 17))
    allFhs = np.zeros((Ninits,Niter))
    residuals = np.zeros((Ninits,Nbins,29))
    for j in range(Ninits):
        lis, nXps, Xps, sfh, lsor = run()
        likelihoods[j,:] = np.array(lis)
        chains[j,:,:] = np.array(Xps)
        residuals[j,:,:] = np.array(lsor)

    fig,ax = plt.subplots()
    for j in range(Ninits):
        ax.plot(range(Niter), likelihoods[j,:], c='gray')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('lnlikelihood')
    plt.savefig('nuclear_mcmc_lnlik.png')
    plt.close(fig)

    fig,ax = plt.subplots(3,6)
    for i in range(17):
        for j in range(Ninits):
            ax.flatten()[i].plot( range(Niter), chains[j,:,i]  )
    for j in range(Ninits):
        ax.flatten()[17].plot( range(Niter), allFhs[j,:] )
    plt.savefig('nuclear_mcmc_trace.png')
    plt.close(fig)

    likind = np.argmax( likelihoods[:,-1] )

    lkrange = [ np.min(likelihoods[:,-1]), np.max(likelihoods[:,-1]) ]
    reslabels=['SMHM0','SMHM1','SMHM2','SMHM3','MS0','MS1','MS2','MS3','MZR*','MZR0','MZR1','MZR2','MZR3','FG0','FG1','FG2','FG3','FGHI','Den0','Den1','Den2','Den3','MRe0','MRe1','MRe2','MRe3','HIR','MC82','TF']
    

    fig,ax=plt.subplots(5,6, figsize=(14,14))
    for j in range(Ninits):
        cscalar = (likelihoods[j,-1] - lkrange[0])/(lkrange[1]-lkrange[0])
        lw=1
        color = (0.8, 0.8, 0.8)
        if j==likind:
            lw=3
            color = (0.98-cscalar/2.0, .7, 0.5+cscalar/2.1 )
        for k in range(29):
            ax.flatten()[k].plot( range(Nbins), residuals[j,:,k], color=color, lw=lw )
            ax.flatten()[k].set_ylabel(reslabels[k])
            ax.flatten()[k].set_ylim(-6,6)
    plt.savefig('nuclear_residuals.png')
    plt.close(fig)


    print "Maximum likelihood: ", likelihoods[likind,-1], chains[likind,-1,:], allFhs[likind,-1]




def appendLikelihood():
    ''' Use something analogous to lnlik in the faceMCMC above, but append it to our data file! '''
    X_train, X_validate, X_test,  Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=1.0, validateFrac=0.0, fn='broad05_to_lasso.txt')
    
    lnlik = globalLikelihood(Ys_train)


    arr = np.vstack([ X_train.T, Ys_train.T, lnlik.T]).T
    header = "Mh0 raccRvir rstarRed rgasRed fg0mult muColScaling muFgScaling muNorm ZIGMfac zmix eta Qf alphaMRI epsquench accCeiling conRF kZ xiREC"
    header += " Mh02 raccRvir2 rstarRed2 rgasRed2 fg0mult2 muColScaling2 muFgScaling2 muNorm2 ZIGMfac2 zmix2 eta2 Qf2 alphaMRI2 epsquench2 accCeiling2 conRF2 kZ2 xiREC2"
    k0s = [0,250,500,750]
    kss = [15,30,60,120,240]
    for i, k0 in enumerate( k0s ):
        for j, ks in enumerate( kss ):
            header+=" k0-"+str(k0)+"ks-"+str(ks)
    for label in labels:
        header += " "+label
    header += " lnlik"
    #with open('broad05_to_lasso.txt', 'r') as f:
    #    header = f.readline()
    #header = header[2:-1] + 'lnlik'
    np.savetxt('broad05_to_lasso_lnlik.txt', arr, header=header)




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

def loglin(x,p):
    return np.log10(x) #### for debugging purposes
    # p is a parameter varying from 0 to 1, with zero corresponding to log, and 1 corresponding to linear.
    assert 0<=p and p<=1
    if p==0:
        return np.log(x)
    else:
        return (np.power(x,p)-1)/p

def readData(trainFrac=0.5, validateFrac=0.4, fn='broad05_to_lasso.txt', naccr=50, arr=None, kscale=10, logscale=0, includeCrossTerms=False):
    ''' trainFrac denotes what fraction of the sample should be used for training 
        validateFrac denotes what fraction of the sample should be used for validation
        1 - trainFrac - validateFrac is the out-of-sample test.'''
    if arr is None:
        arr = np.loadtxt(fn)
    with open(fn, 'r') as f:
        header = f.readline()[2:]
    labels = header.split()[nvars+1:nvars+1+len(logPreds)*nz+nRadii*len(logRadials)]
    X1 = arr[:,:nvars] # the 18 input variables we varied
    X2 = arr[:, nvars+1+len(logPreds)*nz+nRadii*len(logRadials):] # the random factors that define the accretion history
    #X1[:,0] = X1[:,0]/1.0e12
    for i in range(len(logVars)):
        # Take the log of the input variable if it makes sense to do so
        if logVars[i] == 1:
            #X1[:,i] = np.log10(X1[:,i])
            X1[:,i] = loglin(X1[:,i], logscale)
    # Add in squares and cross-terms for every pair of input parameters
    nxv = np.shape(X1)[1]
#### In this block add squares of various quantities (useful for linear models, but superfluous (I think!) for tree-based methods
    if includeCrossTerms:
        for i in range(nxv):
            for j in range(nxv):
                # don't repeat combinations, only look at squares instead of cross-terms, or cross-terms with mass alone <--- no longer applies
                # at this point we're just adding squares, e.g. M_h0^2, rgasRed^2, ...
                if j>=i and (i==0):
                    to_add = X1[:,i]*X1[:,j]
                    X1 = np.vstack( [ X1.T, to_add ] ).T
    #    ## add in even higher order mass terms
        #to_add = np.power(X1[:,0]-12, 2.0)
        #X1 = np.vstack( [ X1.T, to_add] ).T
        to_add = np.power(X1[:,0], 3.0)
        X1 = np.vstack( [ X1.T, to_add] ).T
        to_add = np.power(X1[:,0], 4.0)
        X1 = np.vstack( [ X1.T, to_add] ).T

##### In this block, add in various weights of the accretion history.
    #k0s = [0,250,500,750] # only sum up k's higher than these values. Causality argument: e.g. SFR at z=2 shouldn't be affected by accr. at z<2
    k0s = [0]
    #kss = [4, 8, 16, 32, 64, 128] # falloff rate of accretion random factor indices to 
    kss = [kscale]
    X3 = np.zeros( (np.shape(X2)[0], len(kss)*len(k0s)) )
    for i, k0 in enumerate( k0s ):
        for j, ks in enumerate( kss ):
            ind = i*len(kss) + j
            multArray = np.zeros( np.shape(X2)[1] )
            multArray[k0:] = np.exp( - (np.arange( np.shape(X2)[1] )[k0:] - k0)/float(ks) )
            multArray = np.tile(multArray, (np.shape(X2)[0], 1) )
            print np.shape(multArray), np.shape(X2) # should be the same!

            # X3[:, ind] = np.log10( np.sum( np.power(10.0, X2) * multArray, axis=1 ) )
            #X3[:, ind] = np.sum( np.power(10.0, X2) * multArray, axis=1 ) 
            X3[:, ind] = loglin( np.sum( np.power(10.0, X2) * multArray, axis=1 ) , logscale )

##### In this block, just take the accretion history and reduce it down to a more manageable number by simple averaging.
    X4 = np.zeros( (np.shape(X2)[0], naccr) )
    for i in range(naccr):
        ind0 = i*int(1000/naccr)
        ind1 = (i+1)*int(1000/naccr)
        multArray = np.zeros( np.shape(X2)[1] )
        multArray[ind0:ind1] = 1.0
        multArray = np.tile(multArray, (np.shape(X2)[0], 1) )

        # X4[:, ind] = np.log10( np.sum( np.power(10.0, X2) * multArray, axis=1 ) )
        # X4[:, i] = np.sum( np.power(10.0, X2) * multArray, axis=1 ) 
        X4[:, i] = loglin( np.sum( np.power(10.0, X2) * multArray, axis=1 ) , logscale )

    X = np.vstack( [ X1.T, X3.T, X4.T ] ).T # mc params + random factors

    #nxvOrig = np.shape(X)[1]
    Ys = arr[:,nvars+1:nvars+1+len(logPreds)*nz+nRadii*len(logRadials)]
    # Mh mstar sSFR sfZ stZ gasToStellarRatioH2 gasToStellarRatioHI halfMassStars vPhi22 c82 Sigma1 specificJStars metallicityGradientR90 maxsig mdotBulgeG fractionGI tdep tDepH2 broeilsHI, mStellarHalo
    valid = Ys[:,0]>0
    X =  filterArrayRows( X, valid )
    Ys = filterArrayRows( Ys, valid)
    ns = np.shape(X)[0]
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

    for j in range(len(logPreds)):
        if logPreds[j] == 1:
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

    print "Are all X data finite? ", np.all(np.isfinite(X_train)), np.all(np.isfinite(X_test))
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
            pdb.set_trace()


    return errors_train, errors_validate, errors_test, labels, lasso_coefs, ridge_coefs, theModel



def learnNPR(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, q=1):


    errors_train = np.zeros(len(labels))
    errors_validate= np.zeros(len(labels))
    errors_test = np.zeros(len(labels))
    theModel = np.empty( len(labels), dtype=object)
    #lasso_coefs = np.zeros( (len(labels), np.shape(X_train)[1]) )
    for k in range(len(labels)):
        if True:
            #kq = smooth.NonParamRegression(X_train.T, Ys_train[:,k].T, method=npr_methods.LocalPolynomialKernel(q=q))
            kq = smooth.NonParamRegression(X_train.T, Ys_train[:,k].T, method=npr_methods.SpatialAverage())
            #regLasso = linear_model.Lasso(alpha=alpha, selection='random', tol=1.0e-5, max_iter=1e6, normalize=True)
            invalid = Ys_train[:,k]!=Ys_train[:,k]
            if not np.all(np.isfinite(Ys_train[:,k])): 
                pdb.set_trace()
            #regLasso.fit(X_train, Ys_train[:,k])
            kq.fit()
            # pdb.set_trace()
            Y_pred_train = np.zeros(np.shape(Ys_train[:,k].T))
            kq(X_train.T, out=Y_pred_train)
            Y_pred_validate = kq(X_validate.T)
            Y_pred_test = kq(X_test.T)
            #lasso_coefs[k, : ] = regLasso.coef_[:]
            errors_train[k] = np.std(Y_pred_train - Ys_train[:,k])
            errors_validate[k] = np.std(Y_pred_validate - Ys_validate[:,k])
            errors_test[k] = np.std(Y_pred_test - Ys_test[:,k])
            #theModel[k] = copy.deepcopy( kq )
            theModel[k] = kq 
            print "Success for ",labels[k]
        else:
            print "Failed for k=",k, labels[k]
            print " "
    return errors_train, errors_validate, errors_test, labels, theModel



def learnRF(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, k=0, n_estimators=1000, max_depth=10, max_features=10, min_per_leaf=5):

#    errors_train = np.zeros(len(labels))
#    errors_validate= np.zeros(len(labels))
#    errors_test = np.zeros(len(labels))
#    theModel = np.empty( len(labels), dtype=object)
    theModel= None
#    feature_importances = np.zeros( (len(labels), np.shape(X_train)[1]) )
#    for k in range(len(labels)):
#        if True:
    regRF = ensemble.RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth, max_features=max_features, min_samples_leaf=min_per_leaf, n_jobs=3)
    #regRF = ensemble.ExtraTreesRegressor(n_estimators=n_estimators, max_depth=max_depth, max_features=max_features, min_samples_leaf=min_per_leaf, n_jobs=3)
    invalid = Ys_train[:,k]!=Ys_train[:,k]
    if np.any(invalid):
        pdb.set_trace()
    try:
        regRF.fit(X_train, Ys_train[:,k])
    except:
        pdb.set_trace()
    Y_pred_train = copy.deepcopy( regRF.predict(X_train) )
    Y_pred_validate = copy.deepcopy( regRF.predict(X_validate) )
    Y_pred_test = copy.deepcopy( regRF.predict(X_test) )
    feature_importances = copy.deepcopy( regRF.feature_importances_[:] )
    errors_train = np.std(Y_pred_train - Ys_train[:,k])
    errors_validate = np.std(Y_pred_validate - Ys_validate[:,k]) 
    errors_test = np.std(Y_pred_test - Ys_test[:,k])
    theModel = regRF
    #theModel[k] = copy.deepcopy( regRF )
    print "Finished learning ",labels[k], "for ntrees=",n_estimators
#        else:
#            print "Failed for k=",k, labels[k]
#            print " "
    return errors_train, errors_validate, errors_test, labels, feature_importances, theModel



def learnLasso(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, alpha=0.1, k=0, tol=1.0e-5, l1_ratio=0.5):

    lasso_coefs = np.zeros(  np.shape(X_train)[1]) 

    #regOLS = linear_model.LinearRegression()
    #regLasso = linear_model.ElasticNet(alpha=alpha, l1_ratio=l1_ratio, selection='random', tol=tol, max_iter=1e6, normalize=True, warm_start=True)
    regLasso = linear_model.Lasso(alpha=alpha, selection='random', tol=tol, max_iter=1e6, normalize=True, warm_start=True)
    #regRidge = linear_model.Ridge(alpha=alpha, tol=1.0e-5, max_iter=1e6, normalize=True) 
    invalid = Ys_train[:,k]!=Ys_train[:,k]
    if np.any(invalid):
        pdb.set_trace()
    regLasso.fit(X_train, Ys_train[:,k])
    #regRidge.fit(X_train, Ys_train[:,k])
    #regOLS.fit(X_train, Ys_train[:,k])
    Y_pred_train = regLasso.predict(X_train)
    Y_pred_validate = regLasso.predict(X_validate)
    Y_pred_test = regLasso.predict(X_test)
    lasso_coefs[ : ] = regLasso.coef_[:]
    errors_train = np.std(Y_pred_train - Ys_train[:,k])
    errors_validate = np.std(Y_pred_validate - Ys_validate[:,k])
    errors_test = np.std(Y_pred_test - Ys_test[:,k])
    theModel = copy.deepcopy( regLasso )
    #print labels[k], '--', np.std(Y_pred_train - Ys_train[:,k]), np.std(Y_pred_test - Ys_test[:,k]), '--', np.std(Y_pred_train_ols - Ys_train[:,k]), np.std(Y_pred_test_ols - Ys_test[:,k]), '--', np.std(Y_pred_train_ridge - Ys_train[:,k]), np.std(Y_pred_test_ridge - Ys_test[:,k])
    #print "L1 norm for lasso, ridge fit: ", np.sum(np.abs(regLasso.coef_)), np.sum(np.abs(regRidge.coef_))
    #print "Number of coeffs used in lasso/ridge fit: ", np.sum( np.ones(len(regLasso.coef_))[np.abs(regLasso.coef_)>1.0e-4] ),  np.sum( np.ones(len(regRidge.coef_))[np.abs(regRidge.coef_)>1.0e-4] ), "of", len(regLasso.coef_)
    return errors_train, errors_validate, errors_test, labels, lasso_coefs, theModel
        ####errors_train_this, errors_validate_this, errors_test_this, labels_this, lasso_coef_this, test_residuals_this, theModel = 
        ### learnLasso(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, aL )
    

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
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.8, validateFrac=0.1)

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





def validateNPR():
    ## Read in data
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.1, validateFrac=0.4)

    ## For now we're only going to fit the first 80 variables
    labels=labels[:80]

    ## Choose the order of polynomials we'll use in the nonparametric fit
    #qs = range(1,3)
    qs = [1]
    
    ## Transform the ordinates to (perhaps) improve behavior of the fitting algorithms
    tra = xtransform(X_train)
    X_train = tra.transform(X_train)
    X_validate = tra.transform(X_validate)
    X_test = tra.transform(X_test)


    ## Initialize arrays to store results from each choice of hyperparameters for the learning algorithm
    errors_train = np.zeros((len(labels), len(qs)))
    errors_validate = np.zeros((len(labels), len(qs)))
    errors_test = np.zeros((len(labels), len(qs)))
    residuals_test = np.zeros((len(labels), len(qs)))
    model_array = np.empty((len(labels), len(qs)), dtype=object)
    minIs = np.zeros( len(labels), dtype=int )
    minValErrors = np.zeros( len(labels) ) + 1.0e30

    for i, q in enumerate(qs):
        errors_train_this, errors_validate_this, errors_test_this, labels_this, theModel = learnNPR(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, q )
        errors_train[:,i] = errors_train_this[:]
        errors_validate[:,i] = errors_validate_this[:]
        errors_test[:,i] = errors_test_this[:]
        model_array[:,i] = theModel[:]
    for i in range(len(qs)):
        for k in range(len(labels)):
            if errors_validate[k,i] < minValErrors[k]:
                minValErrors[k] = errors_validate[k,i]
                minIs[k] = i
    print "DependentVar, q, testErr, validationErr, trainingErr, "
    models = []
    for k in range(len(labels)):
        # For each dependent variable, show number of lasso coeffs, minimum validation error, corresponding training error, mean/std of ridge coeffs.
        toprint = labels[k]
        toprint += ' ' + str(qs[minIs[k]])
        toprint += ' ' + str(errors_test[k,minIs[k]])
        toprint += ' ' + str(minValErrors[k])
        toprint += ' ' + str(errors_train[k,minIs[k]]) #,  lasso_coefs[k,:,minIs[k]]
        print toprint
        models.append( model_array[k, minIs[k]] )

    pickle.dump(models, open('validatedNPRModels.pickle', 'w'))
    return models 





def validateLasso():
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.85, validateFrac=0.1)

    labels = labels[:80] # only look at first 80 variables for the moment
    #alphasLasso = np.power(10.0, np.linspace(-8, 1, 12)) # 3 just to test something out
    alphasLasso = [1.0e-7]

    ## Transform the ordinates to (perhaps) improve behavior of the fitting algorithms
    tra = xtransform(X_train)
    X_train = tra.transform(X_train)
    X_validate = tra.transform(X_validate)
    X_test = tra.transform(X_test)



    errors_train = np.zeros((len(labels), len(alphasLasso)))
    errors_validate = np.zeros((len(labels), len(alphasLasso)))
    errors_test = np.zeros((len(labels), len(alphasLasso)))
    residuals_test = np.zeros((len(labels), len(alphasLasso)))
    model_array = np.empty((len(labels), len(alphasLasso)), dtype=object)
    #ridge_coefs = np.zeros( (len(labels), np.shape(X_train)[1]-19) )
    #lasso_coefs = np.zeros( (len(labels), 19) )
    lasso_coefs = np.zeros((len(labels), np.shape(X_train)[1], len(alphasLasso)))
    print "lasso_coefs shape: ", np.shape(lasso_coefs)
    minIs = np.zeros( len(labels), dtype=int )
    minValErrors = np.zeros( len(labels) ) + 1.0e30
    for i, aL in enumerate(alphasLasso):
        errors_train_this, errors_validate_this, errors_test_this, labels_this, lasso_coef_this, theModel = learnLasso(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, aL )
        errors_train[:,i] = errors_train_this[:]
        errors_validate[:,i] = errors_validate_this[:]
        errors_test[:,i] = errors_test_this[:]
        lasso_coefs[:,:, i] = lasso_coef_this[:,:]
        model_array[:,i] = theModel[:]
    for i in range(len(alphasLasso)):
        for k in range(len(labels)):
            if errors_validate[k,i] < minValErrors[k]:
                minValErrors[k] = errors_validate[k,i]
                minIs[k] = i
    print "DependentVar, alphaLasso, alphaRidge, #Lasso, testErr, validationErr, trainingErr, avgRidgeCoef, stdRidgeCoef"
    models = []
    for k in range(len(labels)):
        # For each dependent variable, show number of lasso coeffs, minimum validation error, corresponding training error, mean/std of ridge coeffs.
        print labels[k], alphasLasso[minIs[k]], np.sum(np.ones(len(lasso_coefs[0,:,0]))[np.abs(lasso_coefs[k,:,minIs[k]])>1.0e-4]), errors_test[k,minIs[k]], minValErrors[k], errors_train[k,minIs[k]] #,  lasso_coefs[k,:,minIs[k]]
        models.append( model_array[k, minIs[k]] )

    pickle.dump(models, open('validatedLassoModels_p85.pickle', 'w'))
    return models 


def searchLinearModels(maxiter = 1000):
    #ntrees, max_depth, min_per_leaf, naccr = 200, 1240, 3, 5 # initial guesses
    #hyperparams = [ ntrees, max_depth, min_per_leaf, naccr ]
    #hyperparambounds = [[3, 2000], [2, 30000], [1,200], [1,1000]]

    alpha, tol, kscale, naccr, logscale, l1_ratio = 1.0e-8, 1.0e-5, 10, 10, 0.1, 0.5
    hyperparams = [ alpha, tol, kscale, naccr, logscale, l1_ratio ]
    hyperparambounds = [[1.0e-15, 100.0], [1.0e-10, 1.0e-1], [1,1000], [1,1000], [0,1], [0,1]]

    history = np.zeros( (len(hyperparams)+3, maxiter) )
    standardErrors = np.zeros( (3,maxiter) )

    arr = np.loadtxt('broad567.txt')

    def trial(hyperparams):
        print "Beginning trial with hyperparams: ", hyperparams
        #readData(trainFrac=0.5, validateFrac=0.4, fn='broad05_to_lasso.txt', naccr=50, arr=None, kscale=10, logscale=0):
        X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.85, validateFrac=0, kscale=hyperparams[2], logscale=0, naccr=int(hyperparams[3]), arr=copy.deepcopy(arr), fn='broad567.txt')
        ### Select a random subset of the X_train to be the validation set for this iteration
        nsamples = np.shape(X_train)[0]
        validationIndices = np.random.uniform(0,nsamples-1, size=int(nsamples*0.1))
        vi = [int(v) for v in validationIndices]
        validationSample = np.array([k in vi for k in range(nsamples)])
        X_validate = X_train[validationSample, :]
        X_train = X_train[np.logical_not(validationSample), :]
        Ys_validate = Ys_train[validationSample,:]
        Ys_train = Ys_train[np.logical_not(validationSample),:]

        # Regularize the input and output vectors to give the learning algorithm an easier time
        Xtra = xtransform(X_train)
        X_train = Xtra.transform(X_train)
        X_validate = Xtra.transform(X_validate)
        X_test = Xtra.transform(X_test)

        Ytra = xtransform(Ys_train)
        Ys_train = Ytra.transform(Ys_train)
        Ys_validate = Ytra.transform(Ys_validate)
        Ys_test = Ytra.transform(Ys_test)


        # k=8 corresponds to z=0 sfr. Try really hard to get this right!
        ###errors_train, errors_validate, errors_test, labels, lasso_coefs, theModel
        errors_train_this, errors_validate_this, errors_test_this, labels_this, coefs, theModel = learnLasso(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, alpha=hyperparams[0], k=8, tol=hyperparams[1], l1_ratio=hyperparams[5] )

        return errors_train_this, errors_validate_this, errors_test_this, np.std(Ys_train[:,8]), np.std(Ys_validate[:,8]), np.std(Ys_test[:,8])

    # initialize history
    history[:-3,0] = hyperparams[:]
    history[-3,0] = 100
    history[-2,0] = 100
    history[-1,0] = 100
    nacceptances=0

    for i in np.arange(1,maxiter):
        proposal = np.array(hyperparams) * np.power(10.0, np.random.normal(0,1,size=len(hyperparams))*0.1) # log-normal w/ 10% variance
        # clean up the proposal: make sure it's in bounds and is composed of integers
        for k in range(len(hyperparams)):
            if proposal[k]>hyperparambounds[k][1]:
                proposal[k] = hyperparambounds[k][1]
            if proposal[k]<hyperparambounds[k][0]:
                proposal[k] = hyperparambounds[k][0]

        tre, ve, tee, sTr, sV, sTe = trial(proposal)
        standardErrors[0,i] = sTr
        standardErrors[1,i] = sV
        standardErrors[2,i] = sTe
        accept=False
        # if this trial is better (in the sense that its validation error is lower), accept it
        if ve < history[-2,i-1]:
            accept=True
        # even if this trial is worse, accept with a probability which gets drastically lower as ve gets worse
        #if np.power(history[-2,i-1]/ve, 8) > np.random.uniform():
        #    accept=True
        if np.power(10.0, (history[-2,i-1]-ve)/0.003) > np.random.uniform():
            accept=True

        if accept:
            hyperparams = copy.deepcopy(proposal)
            history[:-3,i] = copy.deepcopy(proposal)
            history[-3,i] = tre
            history[-2,i] = ve
            history[-1,i] = tee
            nacceptances += 1
        else:
            history[:,i] = copy.deepcopy(history[:,i-1])
        print "finished iteration ",i," of ",maxiter," acceptance fraction = ",float(nacceptances)/float(i)

        if i%10==0 or i==maxiter-1:
            # now plot the results
            fig,ax = plt.subplots(7,1, figsize=(8,12))
            for j in range(len(hyperparams)):
                ax[j].plot(history[j,1:i], c='k')
            color=[ 'k', 'r', 'b']
            labels=[ 'train', 'validate', 'test']
            for j in range(3):
                ax[6].plot(history[6+j,1:i], c=color[j], label=labels[j])
                ax[6].plot(standardErrors[j,1:i], ls='--', color=color[j])
            #hyperparams = [ alpha, tol, kscale, naccr, logscale, l1_ratio ]
            ax[0].set_ylabel('Penalty')
            ax[0].set_yscale('log')
            ax[1].set_yscale('log')
            ax[1].set_ylabel('Tolerance')
            ax[2].set_ylabel('Accretion Timescale')
            ax[3].set_ylabel('Number of Accr Divisions')
            ax[4].set_ylabel('Log-Linear interpolation')
            ax[5].set_ylabel('L1 Ratio')
            ax[6].set_ylabel('Errors')
            ax[6].set_xlabel('Number of Iterations')
            plt.legend()
            plt.savefig('searchNetParams_forcelog_forcelasso.png')
            plt.close(fig)


def printTable(arr, columnnames, rownames, fn='table.tex'):
    assert len(np.shape(arr))==2
    nrows, ncols = np.shape(arr)
    #if not (nrows==len(rownames) and ncols==len(columnnames)):
    #    pdb.set_trace()

    with open(fn,'w') as f:
        f.write( r'\begin{table*}'+'\n' )
        f.write( r'\tiny'+'\n' )
        f.write( r'\begin{tabular}{'+'l'*(ncols+1)+r'}'+'\n' )
        line = r' '
        for i in range(ncols):
            line+= r' & '+columnnames[i]
        line+=r'\\'+'\n'
        f.write(line)
        f.write(r'\hline'+'\n')
        f.write(r'\hline'+'\n')
        for j in range(nrows):
            line = rownames[j]
            for i in range(ncols):
                color='black'
                if j==np.argmax(arr[:,i]): # the maximum in this column, i.e. the target most influenced by this feature
                    color='blue'
                if i==np.argmax(arr[j,:]): # the maximum in this row, i.e. most important feature for this target
                    color='red'
                if color=='black':
                    line += r' & '+str(np.round(arr[j,i],2))
                else:
                    line += r' & { \bf \textcolor{'+color+r'}{'+str(np.round(arr[j,i],4)) + r'}}'
            line+= r'\\'+'\n'
            f.write(line)
        f.write(r'\end{tabular}'+'\n')
        f.write(r'\end{table*}'+'\n')

class fntModel:
    ### a simple wrapper for models that have been transformed/predict transformed quantities:
    def __init__(self, skModel, xtr, ytr):
        self.model = skModel
        self.xtr = xtr
        self.ytr = ytr
    def predict(self, x):
        return self.ytr.inverseTransform( self.model.predict( self.xtr.transform(x)))
    def predictFill(self,x):
        if np.shape(x)[1]<len(self.xtr.minxps):
            tofill =  len(self.xtr.minxps) - np.shape(x)[1]
            fillers = np.random.random(size=(np.shape(x)[0], tofill))
            thisx = np.vstack([x.T, fillers.T]).T
            transformed_thisx = self.xtr.transform(thisx)
            #transformed_thisx = np.vstack([transformed_thisx.T, fillers.T]).T
            transformed_thisx[:,-tofill:] = fillers[:,:]
            pred = self.model.predict( transformed_thisx ).reshape(np.shape(x)[0],1)
            yinv = self.ytr.inverseTransform(pred)
            return yinv
        else:
            return self.predict(x)

            


def estimateFeatureImportances(analyze=True, pick=True):
    ### Start with a model we like after our experience w/ searchTreeParams
    ### Train and score the model
    ### Predict Y, but shuffle values of one feature at a time in one mass bin at a time.
    ### Record score reduction at each shuffle.
    ### Plot score reduction as fn of mass for each feature.
    from sklearn.metrics import r2_score

    X_train_orig, X_validate, X_test_orig, Ys_train_orig, Ys_validate, Ys_test_orig, labels = readData(trainFrac=0.85, validateFrac=0, naccr=8, fn='broad05_to_lasso_fix.txt') # no need to feed in arr, since we're just reading the data once.
    nsamples = np.shape(X_train_orig)[0]
    nfeatures = np.shape(X_train_orig)[1]

    # 10x cross-validation
    nmasses = 1 # debug
    ntargets = np.shape(Ys_train_orig)[1] # len(labels)
    #ntargets=2 # debug
    ncv =  1 # for debuggerino
    neval = 10
    scores = np.zeros((nfeatures, nmasses, ntargets))
    r2scores = np.zeros(ntargets)
    feature_importances = np.zeros((nfeatures, ntargets))
    derivs_mass = np.zeros((nfeatures,ntargets, nmasses))
    predictions_mass = np.zeros((nfeatures, ntargets, nmasses, neval))
    ordinates_mass = np.zeros((nfeatures, ntargets, nmasses, neval))

    #feature_names = ["Mh0", "raccRvir", "rstarRed", "rgasRed", "fg0mult", "muColScaling", "muFgScaling", "muNorm", "ZIGMfac", "zmix", "eta", "Qf", "alphaMRI", "epsquench", "accCeiling", "conRF", "kZ", "xiREC"]
    feature_names = [r'$M_{h,0}$', r'$r_\mathrm{acc}/R_v$', r'$r_{\mathrm{acc},0}/r_{*,0}$',  r'$r_{\mathrm{acc},0}/r_{g,0}$',  r'$\chi_{f_{g,0}}$', r'$\alpha_\Sigma$', r'$\alpha_{f_g}$', r'$\mu_0$', r'$\chi_{Z_\mathrm{IGM}}$', r'$\xi_\mathrm{acc}$', r'$\eta$', r'$Q_f$', r'$\alpha_\mathrm{MRI}$', r'$\epsilon_\mathrm{quench}$', r'$\epsilon_\mathrm{quench}$', r'$\alpha_\mathrm{con}$', r'$k_Z$', r'$\xi$']
    texlabels = labels[:]
    for i in range(nfeatures-len(feature_names)):
        #feature_names.append("AccHist"+str(i))
        feature_names.append(r'$\dot{M}_'+str(i)+r'$')


    for k in range(ntargets):
        for cvi in range(ncv):
            ### Select a random subset of the X_train to be the validation set for this iteration
            validationIndices = np.random.uniform(0,nsamples-1, size=int(nsamples*0.1))
            vi = [int(v) for v in validationIndices]
            validationSample = np.array([ss in vi for ss in range(nsamples)])
            X_validate = X_train_orig.copy()[validationSample, :]
            X_train = X_train_orig.copy()[np.logical_not(validationSample), :]
            Ys_validate = Ys_train_orig.copy()[validationSample,:]
            Ys_train = Ys_train_orig.copy()[np.logical_not(validationSample),:]

            # Regularize the input and output vectors to give the learning algorithm an easier time
            Xtra = xtransform(X_train)
            X_train = Xtra.transform(X_train)
            X_validate = Xtra.transform(X_validate)
            X_test = Xtra.transform(X_test_orig.copy()) 

            Ytra = xtransform(Ys_train)
            Ys_train = Ytra.transform(Ys_train)
            Ys_validate = Ytra.transform(Ys_validate)
            Ys_test = Ytra.transform(Ys_test_orig.copy())


            # k=8 corresponds to z=0 sfr. Try really hard to get this right!
            errors_train_this, errors_validate_this, errors_test_this, labels_this, feature_importances_this, theModel = learnRF(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, n_estimators=100, k=k, max_depth=1000, max_features='auto', min_per_leaf=3 )
            if pick:
                pickle.dump( fntModel(theModel,Xtra,Ytra) , open('rfnt_'+str(k)+'_'+str(cvi)+'.pickle','w')) ### save the model

            if analyze: 
                feature_importances[:,k] += feature_importances_this[:]/float(ncv)
                massmin, massmax = np.min(X_train[:,0]), np.max(X_train[:,0])
                massminOrig, massmaxOrig = np.min(X_train_orig[:,0]), np.max(X_train_orig[:,0])
                massbins = np.linspace(massminOrig, massmaxOrig, num=nmasses+1)
                bincenters = (massbins[1:]+massbins[:-1])/2.0
                acc = r2_score(Ys_validate[:,k], theModel.predict(X_validate))
                r2scores[k]+=acc/float(ncv)
                print "r2 score, errors: ",cvi, acc, errors_train_this, errors_validate_this, errors_test_this

                for nm in range(nmasses):
                    binmin = massmin + (massmax-massmin)*float(nm)/float(nmasses)
                    binmax = massmin + (massmax-massmin)*float(nm+1)/float(nmasses)
                    binselec = np.logical_and( X_validate[:,0]>binmin, X_validate[:,0]<binmax )

                    #Xvm = X_validate[binselec,:]
                    #Yvm = Ys_validate[binselec,:]
                    for nfi in range(nfeatures):

                        X_construct_0 = np.mean(Xtra.inverseTransform(X_validate[binselec,:]), axis=0)
                        X_construct_min = np.min(Xtra.inverseTransform(X_validate[binselec,:]), axis=0)
                        X_construct_max = np.max(Xtra.inverseTransform(X_validate[binselec,:]), axis=0)
                        X_construct_1 = X_construct_0.copy()
                        X_construct_1[nfi] = X_construct_1[nfi]*1.1
                        y1 = Ytra.inverseTransform(theModel.predict( Xtra.transform( X_construct_1.reshape(1,nfeatures)) ))[0,k]
                        y0 = Ytra.inverseTransform(theModel.predict( Xtra.transform( X_construct_0.reshape(1,nfeatures)) ))[0,k]
                        val = (y1-y0)/(X_construct_1[nfi] - X_construct_0[nfi])
                        derivs_mass[nfi,k,nm] = val
                        
                        for ne in range(neval):
                            xnfi = X_construct_min[nfi] + (X_construct_max[nfi]-X_construct_min[nfi])*float(ne)/float(neval-1)
                            X_construct_1[nfi] = xnfi
                            predictions_mass[nfi,k,nm,ne] = Ytra.inverseTransform(theModel.predict( Xtra.transform( X_construct_1.reshape(1,nfeatures))))[0,k]
                            ordinates_mass[nfi,k,nm,ne] = xnfi


                        X_v = X_validate.copy()
                        X_vb = X_v[binselec,:]
                        X_vnb = X_v[np.logical_not(binselec),:]
                        #pdb.set_trace()
                        np.random.shuffle(X_vb[:,nfi])
                        X_v[binselec,nfi] = X_vb[:,nfi]
                        shuffle_acc = r2_score(Ys_validate[:,k], theModel.predict(X_v))
                        scores[nfi,nm,k] += ((acc-shuffle_acc)/acc)/float(ncv)



        if analyze:
            try:
                fig,ax = plt.subplots()
                #colors = [(1.0-0.9*float(i)/float(nfeatures), 0.7, 0.1+0.9*float(i)/float(nfeatures)) for i in range(nfeatures)]
                meanscores = np.mean(scores[:,:,k], axis=1)
                sortedscores = sorted(meanscores)[::-1]
                importantfeatures = []
                for impfi in range(len(meanscores)):
                    if meanscores[impfi] in sortedscores[:3]:
                        importantfeatures.append(impfi)
                #for nfi in range(nfeatures):
                    #if np.mean(scores[nfi,:,k])>0.004:
                for nfi in np.array(importantfeatures).flatten():
                    ax.plot(bincenters, scores[nfi,:,k], label=feature_names[nfi])
                ax.text( np.max(bincenters)*0.8, np.max(scores[:,:,k])*0.8, r'$R^2 = $'+str(r2scores[k])[:5] )
                ax.set_yscale('log')
                ax.set_xlabel(r'$\log_{10} M_{h,0}$')
                ax.set_ylabel('Mean decrease accuracy')
                ax.legend(loc=0)
                plt.savefig('feature_importances_'+labels[k]+'.png')
                plt.close(fig)


                fig,ax = plt.subplots()
                #for nfi in range(nfeatures):
                #    if np.mean(scores[nfi,:,k])>0.004:
                for nfi in np.array(importantfeatures).flatten():
                    y = derivs_mass[nfi,k,:]
                    print "y: ",y
                    ax.plot(bincenters, y, label=feature_names[nfi])
                ax.set_xlabel(r'$\log_{10} M_{h,0}$')
                ax.set_ylabel(r'$\partial$'+labels[k]+r'$/\partial$ feature')
                ax.legend(loc=0)
                plt.savefig('feature_derivs_'+labels[k]+'.png')
                plt.close(fig)

                #for nfi in range(nfeatures):
                for nfi in np.array(importantfeatures).flatten():
                    fig,ax = plt.subplots()
                    for nm in range(nmasses):
                        y = predictions_mass[nfi,k,nm, :]
                        x = ordinates_mass[nfi,k,nm, :]
                        ax.plot(x, y, label='Mass Bin '+str(nm))
                    ax.set_xlabel(feature_names[nfi])
                    ax.set_ylabel(labels[k])
                    ax.legend(loc=0)
                    plt.savefig('feature_predictions_'+str(nfi)+'_'+labels[k]+'.png')
                    plt.close(fig)

            except:
                print "FAILED to make plots for k=",k

    if analyze:
        printFeatureImportances(feature_importances.T*10.0, r2scores, feature_names, texlabels)
        printTable(feature_importances.T*10.0, feature_names, texlabels, fn='feature_importances.tex')
        printTable(np.mean(scores,axis=1).T, feature_names, texlabels, fn='avg_swap_scores.tex')

    #print "Average feature importances from skl: ",feature_importances

def printFeatureImportances(arr, rsquared, feature_names, texLabels, fn='feature_importance_compact.tex'):
    assert len(np.shape(arr))==2
    nrows, ncols = np.shape(arr)
    #if not (nrows==len(rownames) and ncols==len(columnnames)):
    #    pdb.set_trace()

    with open(fn,'w') as f:
        f.write( r'\begin{table*}'+'\n' )
        f.write( r'\tiny'+'\n' )
        f.write( r'\begin{tabular}{'+'l'*(5+1+1)+r'}'+'\n' )
        #line = r' '
        #for i in range(ncols):
        #    line+= r' & '+columnnames[i]
        line = r' & 0 & 1 & 2 & 3 & 4' ## the 5 most important features and their scores
        line+=r'\\'+'\n'
        f.write(line)
        f.write(r'\hline'+'\n')
        f.write(r'\hline'+'\n')
        for j in range(nrows):
            if j<=80 or ( j>80 and j%5==0 ):
                # only do the radially resolved quantities at a few points 
                line = rownames[j]
                #for i in range(ncols):
                zipped = zip(arr[j,:], range(ncols))
                zsorted = sorted(zipped, key=lambda x: x[1])
                for i in range(5):
                    line += r' & '+columnnames[zsorted[i][1]] + r': ' + str(np.round(zsorted[i][0],4)) 
                line+= r'& R$^2 = $'+str(np.round(rsquared[j],4))+r' \\'+'\n'
                f.write(line)
        f.write(r'\end{tabular}'+'\n')
        f.write(r'\end{table*}'+'\n')

def searchTreeParams(maxiter = 1000):
    ntrees, max_depth, min_per_leaf, naccr = 200, 1240, 3, 5 # initial guesses
    hyperparams = [ ntrees, max_depth, min_per_leaf, naccr ]
    hyperparambounds = [[3, 2000], [2, 30000], [1,200], [1,1000]]

    history = np.zeros( (len(hyperparams)+3, maxiter) )
    standardErrors = np.zeros( (3,maxiter) )

    arr = np.loadtxt('broad567.txt')

    def trial(hyperparams):
        print "Beginning trial with hyperparams: ", hyperparams
        X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.85, validateFrac=0, naccr=int(hyperparams[3]), arr=copy.deepcopy(arr), fn='broad567.txt')
        ### Select a random subset of the X_train to be the validation set for this iteration
        nsamples = np.shape(X_train)[0]
        validationIndices = np.random.uniform(0,nsamples-1, size=int(nsamples*0.1))
        vi = [int(v) for v in validationIndices]
        validationSample = np.array([k in vi for k in range(nsamples)])
        X_validate = X_train[validationSample, :]
        X_train = X_train[np.logical_not(validationSample), :]
        Ys_validate = Ys_train[validationSample,:]
        Ys_train = Ys_train[np.logical_not(validationSample),:]

        # Regularize the input and output vectors to give the learning algorithm an easier time
        Xtra = xtransform(X_train)
        X_train = Xtra.transform(X_train)
        X_validate = Xtra.transform(X_validate)
        X_test = Xtra.transform(X_test)

        Ytra = xtransform(Ys_train)
        Ys_train = Ytra.transform(Ys_train)
        Ys_validate = Ytra.transform(Ys_validate)
        Ys_test = Ytra.transform(Ys_test)


        # k=8 corresponds to z=0 sfr. Try really hard to get this right!
        errors_train_this, errors_validate_this, errors_test_this, labels_this, feature_importances_this, theModel = learnRF(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, n_estimators=int(hyperparams[0]), k=8, max_depth=int(hyperparams[1]), max_features='auto', min_per_leaf=int(hyperparams[2]) )

        return errors_train_this, errors_validate_this, errors_test_this, np.std(Ys_train[:,8]), np.std(Ys_validate[:,8]), np.std(Ys_test[:,8])

    # initialize history
    history[:-3,0] = hyperparams[:]
    history[-3,0] = 100
    history[-2,0] = 100
    history[-1,0] = 100
    nacceptances=0

    for i in np.arange(1,maxiter):
        proposal = np.array(hyperparams) * np.power(10.0, np.random.normal(0,1,size=len(hyperparams))*0.05) # log-normal w/ 10% variance
        # clean up the proposal: make sure it's in bounds and is composed of integers
        for k in range(len(hyperparams)):
            if proposal[k]>hyperparambounds[k][1]:
                proposal[k] = hyperparambounds[k][1]
            if proposal[k]<hyperparambounds[k][0]:
                proposal[k] = hyperparambounds[k][0]

        tre, ve, tee, sTr, sV, sTe = trial(proposal)
        standardErrors[0,i] = sTr
        standardErrors[1,i] = sV
        standardErrors[2,i] = sTe
        accept=False
        # if this trial is better (in the sense that its validation error is lower), accept it
        if ve < history[-2,i-1]:
            accept=True
        # even if this trial is worse, accept with a probability which gets drastically lower as ve gets worse
        #if np.power(history[-2,i-1]/ve, 8) > np.random.uniform():
        #    accept=True
        if np.power(10.0, (history[-2,i-1]-ve)/0.003) > np.random.uniform():
            accept=True

        if accept:
            hyperparams = copy.deepcopy(proposal)
            history[:-3,i] = copy.deepcopy(proposal)
            history[-3,i] = tre
            history[-2,i] = ve
            history[-1,i] = tee
            nacceptances += 1
        else:
            history[:,i] = copy.deepcopy(history[:,i-1])
        print "finished iteration ",i," of ",maxiter," acceptance fraction = ",float(nacceptances)/float(i)

    # now plot the results
    fig,ax = plt.subplots(5,1, figsize=(8,12))
    for i in range(len(hyperparams)):
        ax[i].plot(history[i,:], c='k')
    color=[ 'k', 'r', 'b']
    labels=[ 'train', 'validate', 'test']
    for i in range(3):
        ax[4].plot(history[4+i,1:], c=color[i], label=labels[i])
        ax[4].plot(standardErrors[i,1:], ls='--', color=color[i])
    ax[0].set_ylabel('Number of Trees')
    ax[1].set_ylabel('Maximum Depth')
    ax[2].set_ylabel('Minimum per Leaf')
    ax[3].set_ylabel('Number of Accretion Divisions')
    ax[4].set_ylabel('Errors')
    ax[4].set_xlabel('Number of Iterations')
    plt.legend()
    plt.savefig('searchTreeParams.png')
    plt.close(fig)




def validateRF():
    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=0.85, validateFrac=0.1)

    labels = labels[:80] # only look at first 80 variables for the moment
    # nfeatures = [100,1000,10000,100000] # 3 just to test something out
    ntrees = [400]
    max_depths = [ 640 ]
    max_features = ['auto'] 
    min_per_leaf = [ 6 ]
    #max_features = ['dummy']

    ## Transform the ordinates to (perhaps) improve behavior of the fitting algorithms
    Xtra = xtransform(X_train)
    X_train = Xtra.transform(X_train)
    X_validate = Xtra.transform(X_validate)
    X_test = Xtra.transform(X_test)

    Ytra = xtransform(Ys_train)
    Ys_train = Ytra.transform(Ys_train)
    Ys_validate = Ytra.transform(Ys_validate)
    Ys_test = Ytra.transform(Ys_test)


    errors_train = np.zeros((len(labels), len(ntrees), len(max_depths), len(max_features), len(min_per_leaf)))
    errors_validate = np.zeros((len(labels), len(ntrees), len(max_depths), len(max_features), len(min_per_leaf)))
    errors_test = np.zeros((len(labels), len(ntrees), len(max_depths), len(max_features), len(min_per_leaf)))
    #residuals_test = np.zeros((len(labels), len(ntrees), len(max_depths), len(max_features)))
    model_array = np.empty((len(labels), len(ntrees), len(max_depths), len(max_features), len(min_per_leaf)), dtype=object)
    #ridge_coefs = np.zeros( (len(labels), np.shape(X_train)[1]-19) )
    #lasso_coefs = np.zeros( (len(labels), 19) )
    feature_importances = np.zeros((len(labels), np.shape(X_train)[1], len(ntrees), len(max_depths), len(max_features), len(min_per_leaf)))
    print "lasso_coefs shape: ", np.shape(feature_importances)
    minIs = np.zeros( len(labels), dtype=int )
    minJs = np.zeros( len(labels), dtype=int )
    minQs = np.zeros( len(labels), dtype=int )
    minRs = np.zeros( len(labels), dtype=int )
    minValErrors = np.zeros( len(labels) ) + 1.0e30
    for i, ntree in enumerate(ntrees):
        for j, max_depth in enumerate(max_depths):
            for q, max_feature in enumerate(max_features):
                for rr, minpl in enumerate(min_per_leaf):
                    for k,label in enumerate(labels):
                        errors_train_this, errors_validate_this, errors_test_this, labels_this, feature_importances_this, theModel = learnRF(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, n_estimators=ntree, k=k, max_depth=max_depth, max_features=max_feature, min_per_leaf=minpl )
                        errors_train[k,i,j,q,rr] = errors_train_this
                        errors_validate[k,i,j,q,rr] = errors_validate_this
                        errors_test[k,i,j,q,rr] = errors_test_this
                        feature_importances[k,:, i,j,q,rr] = feature_importances_this[:]
                        model_array[k,i,j,q,rr] = theModel
    for i in range(len(ntrees)):
        for j in range(len(max_depths)):
            for q in range(len(max_features)):
                for rr in range(len(min_per_leaf)):
                    for k in range(len(labels)):
                        if errors_validate[k,i,j,q,rr] < minValErrors[k]:
                            minValErrors[k] = errors_validate[k,i,j,q,rr]
                            minIs[k] = i
                            minJs[k] = j
                            minQs[k] = q
                            minRs[k] = rr
    models = []
    for k in range(len(labels)):
        yvar = np.std(Ys_train[:,k])
        # For each dependent variable, show number of lasso coeffs, minimum validation error, corresponding training error, mean/std of ridge coeffs.
        print '--------'
        print labels[k]
        print ntrees[minIs[k]], max_depths[minJs[k]], max_features[minQs[k]], min_per_leaf[minRs[k]], errors_test[k,minIs[k],minJs[k],minQs[k],minRs[k]]/yvar, minValErrors[k]/yvar, errors_train[k,minIs[k],minJs[k],minQs[k],minRs[k]]/yvar, yvar #,  lasso_coefs[k,:,minIs[k]]
        models.append( model_array[k, minIs[k], minJs[k], minQs[k], minRs[k]] )

    pickle.dump(models, open('validatedRFdeep.pickle', 'w'))
    return models 






#def validateLasso():
#    alphas = np.power(10.0, np.linspace(-9,-1.5,30))
#    npar = 19
#    cs=[ (1.0 - float(i)/npar/2.0, 0.5, 0.5+float(i)/npar/2.0) for i in range(npar)]
#    fig,ax = plt.subplots()
#    errors_train = np.zeros((npar, len(alphas)))
#    errors_test = np.zeros((npar, len(alphas)))
#    labels=None
#    X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels = readData() 
#    for i, a in enumerate(alphas):
#        errors_train_this, errors_test_this, labels_this = learnLasso(X_train, X_test, X_validate, Ys_validate, Ys_train, Ys_test, labels, alpha=a)
#        errors_train[:,i] = errors_train_this[:]
#        errors_test[:,i] = errors_test_this[:]
#        labels=labels_this
#    for j in range(npar):
#        ax.plot(alphas, errors_train[j,:], c=cs[j], ls='--', lw=2)
#        ax.plot(alphas, errors_test[j,:], c=cs[j], ls='-', lw=2, label=labels[j])
#    #ax.legend(loc=0)
#    ax.set_xlabel(r'$\alpha$')
#    ax.set_xscale('log')
#    ax.set_ylabel(r'$\sigma_{y-y_\mathrm{pred}}$')
#    plt.savefig('lasso_validate.png')
#    plt.close(fig)

if __name__=='__main__':
    #learn(1.0e4, 1.0e-5)
    #validateSVM()
    #appendLabels()
    #ksTestAnalysis()
    #learnLasso(alpha=1.0e-3)
    #validateLasso()
    #validateRF()
    #validateLassoRidge()

    #searchTreeParams(400)
    #searchLinearModels(800)

    #runEmcee(mpi=True, continueRun=True)
    #fractionalVariancePlot()
    #ridgeCoeffsPlot()

    #validateNPR()
    #plotResiduals()

    #estimateFeatureImportances(analyze=False, pick=True) # just generate the pickled models 
    #estimateFeatureImportances(analyze=True, pick=False) # just generate the pickled models 
    
    #nuclearSearch(Nbins = 7, Niter=10000, Ninits=50)

    #appendLikelihood()

    ### analyze the fake mcmc run
    if True:
        restart={}
        updateRestart('fakemcmc_restart.pickle', restart)
        printRestart(restart)
        tracePlots(restart, 'fakemcmc_trace', burnIn=0)
        probsPlots(restart, 'fakemcmc_allProb', burnIn=0)
        #trianglePlot(restart,'fakemcmc_triangle.png', burnin=50, nspace=10)

        # Find the maximum among all models sampled so far.
        allProbs = restart['allProbs'].flatten()
        index = np.argmax(allProbs)
        print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])
        print allProbs[index]
        indices = np.unravel_index(index, np.shape(restart['allProbs']))
        xmax = restart['chain'][indices[0],indices[1],:]
        print "Favorite coordinates: ", xmax
