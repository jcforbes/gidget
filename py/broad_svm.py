import numpy as np
import time
import os
import copy
import pickle
import matplotlib.pyplot as plt
import pdb
import emcee
from sklearn import svm, linear_model, ensemble

#import pyqt_fit.nonparam_regression as smooth
#from pyqt_fit import npr_methods

# Useful information about the setup of the linear models and their fits....
logVars = [1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1]
logPreds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1]
logRadials = [1,1,1,1,1,1]
#logVars = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#logPreds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#logRadials = [0, 0, 0, 0, 0, 0]
nRadii=20
nz=4
nvars = len(logVars)

np.random.seed(1487)
randomFactorsKey01 = np.random.random(size=(50,9)) 


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

class fntModel:
    ### a simple wrapper for models that have been transformed/predict transformed quantities:
    def __init__(self, skModel, xtr, ytr, randomFactors):
        self.model = skModel
        self.xtr = xtr
        self.ytr = ytr
        self.randomFactors = randomFactors
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

def predictFill( model, x, rfkey):
    if np.shape(x)[1]<len(model.xtr.minxps):
        tofill =  len(model.xtr.minxps) - np.shape(x)[1]
        #fillers = np.random.random(size=(np.shape(x)[0], tofill))
        fillers = model.randomFactors[rfkey,:].reshape(np.shape(x)[0], tofill)
        thisx = np.vstack([x.T, fillers.T]).T
        transformed_thisx = model.xtr.transform(thisx)
        #transformed_thisx = np.vstack([transformed_thisx.T, fillers.T]).T
        transformed_thisx[:,-tofill:] = fillers[:,:]
        pred = model.model.predict( transformed_thisx ).reshape(np.shape(x)[0],1)
        yinv = model.ytr.inverseTransform(pred)
        return yinv
    else:
        return model.predict(x)



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

def fakeEmceeResiduals(emceeparams, models):
    # First transform the emceeparams into the same format used by 'X' in the fit of the linear models
    # emceeparams is the set of 18 parameters that we fit with Lasso, minus mass, so we have..
    nmh = 10 
    MhGrid = np.power(10.0, np.linspace(10.0,13,nmh))
    lnlik = 0.0

    residuals=[]
    for k,Mh in enumerate(MhGrid):
        X1 = np.array([Mh]+list(emceeparams[:-2])).reshape((1,len(emceeparams[:-2])+1))
        X1[0,1] = X1[0,1]*np.power(Mh/1.0e12, emceeparams[-2]) # modify raccRvir
        for i in range(len(logVars)):
            # Take the log of the input variable if it makes sense to do so
            if logVars[i] == 1:
                X1[:,i] = np.log10(X1[:,i])

        Y_eval = np.array( [predictFill(models[j], X1, k)[0][j] for j in range(80)] ).reshape(1,80)
        residuals.append( globalLikelihood(Y_eval, fh=emceeparams[-1], returnlikelihood=False) )

    return residuals 

def residualsFromGidget(bn, fh=0.3):
    ''' Read in models with basename bn and construct "Y_eval," then pass that to the globalLikelihood fn. You also need to specify a value for fh, the fraction of mStellarHalo to be included in mstar. '''
    import readoutput
    output = readoutput.Experiment(bn) 
    output.read(keepStars=True, fh=0)
    if len(output.models)>0:

        # borrow some code from mcmc_broad.py on hyades.
        zinds = [readoutput.Nearest( output.models[0].var['z'].sensible(), z)[0] for z in [0,1,2,3]]
        residuals=[]
        mhs = []
        variables = ['Mh', 'mstar', 'sSFR', 'sfZ', 'stZ', 'gasToStellarRatioH2', 'gasToStellarRatioHI', 'halfMassStars', 'vPhi22', 'c82', 'Sigma1', 'specificJStars', 'metallicityGradientR90', 'maxsig', 'mdotBulgeG', 'fractionGI', 'tdep', 'tDepH2', 'broeilsHI', 'mStellarHalo'] ## variables to read from the gidget models
        for model in output.models:
            Y_eval = np.zeros((1,80))
            for j in range(len(variables)):
                for k in range(nz):
                    Y_eval[0, j*nz+k ] = model.var[variables[j]].sensible(timeIndex=zinds[k])
            for ii in range(80):
                if logPreds[ ii/4 ] == 1:
                    Y_eval[0, ii] = np.log10(Y_eval[0,ii])
            residuals.append( globalLikelihood(Y_eval, fh=fh, returnlikelihood=False) )

            mhs.append(model.var['Mh'].sensible(timeIndex=zinds[0]))
        return mhs, residuals
    else:
        return [],[]

#models22 = [ pickle.load( open( 'rfnt22_'+str(k)+'_0.pickle', 'r' ) ) for k in range(80) ]
#models1718b = [ pickle.load( open( 'rfnt1718b_'+str(k)+'_0.pickle', 'r' ) ) for k in range(80) ]

def fakeEmceePlotResiduals(restart, basefn, gidgetmodels=None, xmax=None):

    xmaxOrig = copy.deepcopy(xmax)

    # Find the maximum among all models sampled so far.
    if xmax is None:
        allProbs = restart['allProbs'].flatten()
        zipped = zip(allProbs, range(len(allProbs)))
        zsorted = sorted(zipped, key=lambda x: -x[0])
        highProbInds = [zsorted[i][1] for i in range(10)]

    labels = ['SMHM 0', 'SMHM 1', 'SMHM 2', 'SMHM 3', 'sSFR 0', 'sSFR 1', 'sSFR 2', 'sSFR 3', 'Zst', 'Zg 0', 'Zg 1', 'Zg 2', 'Zg 3', 'fH2 0', 'fH2 1', 'fH2 2', 'fH2 3', 'fHI 0', 'Sigma1 0', 'Sigma1 1', 'Sigma1 2', 'Sigma1 3', 'Rst 0', 'Rst 1', 'Rst 2', 'Rst 3', 'RHI', 'c82', 'TF', 'sig 0', 'sig 1', 'sig 2', 'sig 3']

    fig,ax = plt.subplots(nrows=7, ncols=5, figsize=(14,14))
    massLim = (10,13) 

    chiSquared = np.zeros((10, 33))
    chiSquaredFh0 = np.zeros((10, 33))

    for i in range(2):
        if xmaxOrig is None:
            indices = np.unravel_index(highProbInds[i], np.shape(restart['allProbs']))
            xmaxThis = restart['chain'][indices[0],indices[1],:]
        else:
            xmaxThis=xmaxOrig[:]

        if i%2==0:
            xmaxThis[-1] = 0.0
    
        treeResiduals = np.array( fakeEmceeResiduals(xmaxThis, models17) )
        if i%2==0:
            chiSquaredFh0[i/2, :] = np.sum(np.power(treeResiduals[:,:,:],2.0), axis=0)
        else:
            chiSquared[i/2, :] = np.sum(np.power(treeResiduals[:,:,:],2.0), axis=0)

        lw=1
        c='b'
        alpha=.5
        if i%2==0:
            lw=1
            c='r'
            alpha=0.5
        for j in range(33):
            ax.flatten()[j].plot( np.power(10.0, np.linspace(massLim[0],massLim[1],np.shape(treeResiduals)[0])), treeResiduals[:,0,j], c=c, lw=lw, alpha=alpha ) 

    if gidgetmodels is not None:
        mhs, gidgetResiduals = residualsFromGidget(gidgetmodels, fh=xmaxOrig[-1])
        npGidgetResiduals = np.array(gidgetResiduals)
        for j in range(33):
            ax.flatten()[j].scatter( mhs, npGidgetResiduals[:,0,j], c='b', s=25, lw=0 )

        mhs, gidgetResidualsFh0 = residualsFromGidget(gidgetmodels, fh=0.0)
        npGidgetResidualsFh0 = np.array(gidgetResidualsFh0)
        for j in range(33):
            ax.flatten()[j].scatter( mhs, npGidgetResidualsFh0[:,0,j], c='r', s=25, lw=0 )

    chiSquaredAll = np.zeros((6,33))
    chiSquaredAll[0,:] = np.mean(chiSquared[:,:],axis=0)/10.0
    chiSquaredAll[1,:] = np.std(chiSquared[:,:],axis=0)/10.0
    chiSquaredAll[2,:] = np.mean(chiSquaredFh0[:,:],axis=0)/10.0
    chiSquaredAll[3,:] = np.std(chiSquaredFh0[:,:],axis=0)/10.0
    chiSquaredAll[4,:] = np.sum(np.power(npGidgetResiduals[:,0,:],2.0),axis=0)/33.0
    chiSquaredAll[5,:] = np.sum(np.power(npGidgetResidualsFh0[:,0,:],2.0),axis=0)/33.0

    print "treeResiduals shape: ", np.shape(treeResiduals)
    print "npGidgetResiduals shape: ", np.shape(npGidgetResiduals)
    print "npGidgetResidualsFh0 shape: ", np.shape(npGidgetResidualsFh0)

    for j in range(33):
        ax.flatten()[j].set_xlabel(r'$M_{h,0}$')
        ax.flatten()[j].set_ylabel(labels[j])
        ax.flatten()[j].set_xscale('log')
        ax.flatten()[j].plot([10**massLim[0],10**massLim[1]], [0,0], lw=2, ls='--', c='gray')
        ax.flatten()[j].fill_between( [10**massLim[0], 10**massLim[1]], -1, 1, facecolor='orange', alpha=0.1)
        #ax.flatten()[j].plot([10**massLim[0],10**massLim[1]], [1,1], lw=1, ls='-', c='gray')
        #ax.flatten()[j].plot([10**massLim[0],10**massLim[1]], [-1,-1], lw=1, ls='-', c='gray')
        ax.flatten()[j].set_xlim(10**massLim[0],10**massLim[1])
        ax.flatten()[j].set_ylim(-5,5)
        ax.flatten()[j].set_title( r'$\chi^2_{\mathrm{t}}=$ '+str(np.round(chiSquaredAll[0,j],1)) + r'$\pm$' +str(np.round(chiSquaredAll[1,j],1))+ r' $\chi^2_{\mathrm{t},f=0}=$'+str(np.round(chiSquaredAll[2,j],1)) +r'$\pm$' +str(np.round(chiSquaredAll[3,j],1))+r'$\chi^2_{\mathrm{g}}=$ '+str(np.round(chiSquaredAll[4,j],1)) +r'$\chi^2_{\mathrm{t},f=0}=$ '+str(np.round(chiSquaredAll[5,j],1)), size=7 )

        #plt.axhline(0.0, lw=2, ls='--', c='gray')
    

    plt.tight_layout()
    plt.savefig(basefn+'.pdf')
    plt.close(fig)



    
def lnlikelihood(emceeparams, models=None):
    #models = [ pickle.load( open( 'rfnt10_'+str(k)+'_0.pickle', 'r' ) ) for k in range(80) ]
    # First transform the emceeparams into the same format used by 'X' in the fit of the linear models
    # emceeparams is the set of 18 parameters that we fit with Lasso, minus mass, so we have..
    nmh = 20
    MhGrid = np.power(10.0, np.linspace(10.0,13.0,nmh)) ## try out only going up to 10^13
    lnlik = 0.0

    for k,Mh in enumerate(MhGrid):
        # Let's see.. Mh + emcee parameters, minus the last two "artificial" parameters
        X1 = np.array([Mh]+list(emceeparams[:-2])).reshape((1,len(emceeparams[:-2])+1))
        X1[0,1] = X1[0,1]*np.power(Mh/1.0e12, emceeparams[-2]) # modify raccRvir
        for i in range(len(logVars)):
            # Take the log of the input variable if it makes sense to do so
            if logVars[i] == 1:
                X1[:,i] = np.log10(X1[:,i])
        Y_eval = np.zeros((1,80))
        # 4 for each of the following:  Mhz0 mstarz0 sSFRz0 sfZz0 stZz0 gasToStellarRatioH2z0 gasToStellarRatioHIz0 halfMassStarsz0 vPhi22z0 c82z0 Sigma1z0 specificJStarsz0 metallicityGradientR90z0 maxsigz0 mdotBulgeGz0 fractionGIz0 tdepz0 tDepH2z0 broeilsHIz0 mStellarHaloz0 
        neededModels = [0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15, 16, 20,21,22,23, 24, 28,29,30,31, 32,33,34, 36, 40,41,42,43, 52,53,54,55, 72, 76,77,78,79]
        for j in range(80):
            if j in neededModels:
                Y_eval[0,j] = predictFill(models17[j], X1, k)[0][j]
        #Y_eval = np.array( [predictFill(models10[j], X1)[0][j] for j in range(80)] ).reshape(1,80)
        #lnlik += np.sum( globalLikelihood(Y_eval, fh=emceeparams[-1], returnlikelihood=True) )
        lnlik += np.sum( globalLikelihood(Y_eval, fh=0.0, returnlikelihood=True) )


    print "Returning lnlik = ", lnlik, "for emceeparams ",emceeparams
    if not np.isfinite(lnlik):
        return -np.inf
    return lnlik

def sampleFromGaussianBall():
    ## the max posterior probability estimated from the previous run.
    xmax = [0.147451569889, 2.52941811659, 2.59034734186, 2.41120695741, -0.124283831858, -0.0523679435879, 9.35680382698, -0.974822093888, 1.89905286619, 0.511551421578, 2.0337488747, 2.02929369251, 0.0642244458824, 0.00988965146683, 0.509787819545, 0.279394476293, 1.74417214913, 0.342311450585, 0.0107366934282, 0.414472066814, 1.50430105103e+12, 1.21653967757, -0.028389697679, 0.497607252288]
    draw = []
    for i in range(len(xmax)):
        draw.append( xmax[i]*(1.0 + 0.01*np.random.normal()) )
    if not np.isfinite( lnprior(draw) ):
        print "WARNING: doing recursion in broad_svm.py:sampleFromGaussianBall()"
        return sampleFromGaussianBall() ## if this is a bad draw don't use it!
    return draw


def samplefromprior( varRedFac=1.0):
    # accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormaliza     tion, kappaMassScaling = emceeparams
    # Mhz0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, ZIGMfac, zmix, eta, Qf, alphaMRI, epsqu     ench, accCeiling, conRF, kZ, xiREC = emceeparams


    return [ \
        samplefromlognormaldensity( np.log(0.141), np.log(3.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(2.0), np.log(2.0)**2.0/ varRedFac),
        samplefromnormaldensity( 0.0, 1.0**2.0 / varRedFac ),
        samplefromnormaldensity( 0.0, 1.0**2.0 / varRedFac ),
        samplefromlognormaldensity( np.log(.1), np.log(10.0)**2.0/ varRedFac),
        samplefromnormaldensity( -.5, 1.0**2.0 / varRedFac ),
        samplefromlognormaldensity( np.log(1.0), np.log(10.0)**2.0/ varRedFac ),
        samplefrombetadensity( 1.0 * varRedFac, 1.0 * varRedFac),
        samplefromlognormaldensity( np.log(1.5), np.log(2.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(1.5), np.log(2.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(0.05), np.log(2.0)**2.0/ varRedFac),
        samplefromlognormaldensity( np.log(1.0e-3), np.log(10.0)**2.0/ varRedFac),
        samplefrombetadensity( 1.0 * varRedFac, 1.0 ),
        samplefromnormaldensity( 0.3, 0.3**2.0 / varRedFac ),
        samplefromlognormaldensity( np.log(1.0), np.log(3.0)**2.0/ varRedFac),
        samplefrombetadensity( 1.0, 2.0*varRedFac ),
        samplefromlognormaldensity( np.log(1.0e-2), np.log(10.0)**2.0/varRedFac),
        samplefromnormaldensity( 0.5, 0.5**2/varRedFac), 
        samplefromlognormaldensity( np.log(1.0e12), np.log(3.0)**2.0/varRedFac ),
        samplefromlognormaldensity( np.log(1.0), np.log(3.0)**2.0/varRedFac),  
        samplefromnormaldensity( 0.3, 0.2**2.0/varRedFac) ]

def lnprob(emceeparams, models=None):
    pr = lnprior(emceeparams)
    if np.isfinite(pr):
        return lnlikelihood(emceeparams, models=models) + pr
        #return constlikelihood(emceeparams) + pr
    return pr

def printPriorOffsets(emceeparams):
    pass
#    raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, scaleAdjust, mquench, enInjFac, alpharmh, fh = emceeparams
#    print "Prior offsets: "
#    print "raccRvir: ", np.log10(raccRvir / 0.141) / np.log10(3.0)
#    print "rstarRed: ", np.log10(rstarRed / 2.0) / np.log10(2.0)
#    print "rgasRed: ", np.log10(rgasRed / 2.0) / np.log10(2.0)
#    print "fg0mult: ", np.log10(fg0mult / 2.0) / np.log10(2.0)
#    print "muColScaling: ", (muColScaling - 0.0) / 2.0
#    print "muFgScaling: ", (muFgScaling - 0.0) / 2.0
#    print "muNorm : ", np.log10(muNorm/ 1.0) / np.log10(10.0)
#    print "muMhScaling: ", (muMhScaling - -1.0) / 1.0
#    print "ZIGMfac: ", np.log10(ZIGMfac/1.0) / np.log10(3.0)
#    print "zmix", zmix-0.5
#    print "eta: ", np.log10(eta/1.5) / np.log10(2.0)
#    print "Qf: ", np.log10(Qf/1.5) / np.log10(2.0)
#    print "alphaMRI: ", np.log10(alphaMRI/0.05) / np.log10(2.0)
#    print "epsquench: ", np.log10(epsquench/1.0e-3) / np.log10(10.0)
#    print "accCeiling: ", accCeiling-0.5
#    print "conRF: ", (conRF-0.3)/0.3
#    print "kZ: ", np.log10(kZ/1.0)/np.log10(3.0)
#    print "xiREC: ", xiREC - 0.5
#    print "epsff: ", np.log10(epsff/1.0e-2)/np.log10(2.0)
#    print "scaleAdjust: ", (scaleAdjust - 0.5)/0.5
#    print "mquench: ", np.log10(mquench/1.0e12)/np.log10(3.0)
#    print "enInjFac: ", np.log10(enInjFac/1.0)/np.log10(3.0) 
#    print "alpharMh: ", alpharmh/1.0
#    print "fh: ", fh-0.5


def lnprior(emceeparams):
    raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, scaleAdjust, mquench, enInjFac, chiZslope = emceeparams
    #raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, scaleAdjust, mquench, enInjFac, alpharmh, fh = emceeparams
    #accScaleLength, muNorm, muMassScaling, muFgScaling, muColScaling, accCeiling, eta, fixedQ, Qlim, conRF, kappaNormalizat     ion, kappaMassScaling = emceeparams
    accum = 0.0

    #varRedFac = 10000.0
    varRedFac = 1.0


    accum += lnlognormaldensity( raccRvir, np.log(0.141), np.log(3.0)**2.0 )
    accum += lnlognormaldensity( rstarRed, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
    accum += lnlognormaldensity( rgasRed, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
    accum += lnlognormaldensity( fg0mult, np.log(2.0), np.log(2.0)**2.0/ varRedFac )
    accum += lnnormaldensity( muColScaling, 0.0, 1.0**2.0 / varRedFac )
    accum += lnnormaldensity( muFgScaling, 0.0, 1.0**2.0 / varRedFac )
    accum += lnlognormaldensity( muNorm, np.log(.1), np.log(10.0)**2.0/ varRedFac )
    accum += lnnormaldensity( muMhScaling, -.5, 1.0**2.0 / varRedFac )
    accum += lnlognormaldensity( ZIGMfac, np.log(1.0), np.log(10.0)**2.0/ varRedFac )
    accum += lnbetadensity( zmix, 1.0, 1.0 ) # not accurate
    accum += lnlognormaldensity( eta, np.log(1.5), np.log(2.0)**2.0/ varRedFac )
    accum += lnlognormaldensity( Qf, np.log(1.5), np.log(2.0)**2.0/ varRedFac )
    accum += lnlognormaldensity( alphaMRI, np.log(0.05), np.log(2.0)**2.0/ varRedFac )
    accum += lnlognormaldensity( epsquench, np.log(1.0e-3), np.log(10.0)**2.0/ varRedFac )
    accum += lnbetadensity( accCeiling, 1.0, 1.0 ) # not accurate
    accum += lnnormaldensity( conRF, 0.3, 0.3**2.0 / varRedFac )
    accum += lnlognormaldensity( kZ, np.log(1.0), np.log(3.0)**2.0/ varRedFac )
    accum += lnbetadensity( xiREC, 1.0, 2.0 ) # not accurate
    accum += lnlognormaldensity( epsff, np.log(1.0e-2), np.log(10.0)**2.0 / varRedFac )
    accum += lnnormaldensity( scaleAdjust, 0.5, 0.5**2 /varRedFac)
    accum += lnlognormaldensity( mquench, np.log(1.0e12), np.log(3.0)**2.0/varRedFac)
    accum += lnlognormaldensity( enInjFac, np.log(1.0), np.log(3.0)**2.0/varRedFac)
    accum += lnnormaldensity( chiZslope, 0.3, 0.2**2.0/varRedFac)
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
        ax[i,j].set_ylim( np.min(chain[:, len(chain[0,burnIn:,dim])/2 :, dim] ), np.max(chain[:, len(chain[0,burnIn:,dim])/2 :, dim] ) )

    plt.tight_layout()
    plt.savefig(fn+'.pdf')
    plt.close(fig)

def probsPlots(restart, fn, burnIn=0):
    allProbs = restart['allProbs']
    
    ndim = np.shape(allProbs)[0]
    iters = np.shape(allProbs)[1]

    # Plot the trace of the probabilities for every walker.
    fig,ax = plt.subplots()
    for walker in range(ndim):
        ax.plot(allProbs[walker,burnIn:],alpha=.3,ls='--')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Posterior probability')
    ax.set_ylim( np.min(allProbs[:, len(allProbs[0,burnIn:])/2 :] ), np.max(allProbs[:, len(allProbs[0,burnIn:])/2 :] ) )
    plt.savefig(fn+'_probs.png')
    plt.close(fig)
    print "Saved "+fn+'_probs.pdf'

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

def histplot(ax, hist, edges, c='k', vertical=True):
    for i,y in enumerate(hist):
        if vertical:
            ax.plot([edges[i],edges[i]], [0, y], lw=2, c=c, ls='-')
            ax.plot([edges[i+1],edges[i+1]], [0, y], lw=2, c=c, ls='-')
        ax.plot([edges[i],edges[i+1]], [y,y], lw=2, c=c, ls='-')

def trianglePlot(restart,fn,burnIn=0, nspace=10):
    shp = np.shape(restart['chain'])
    prs = shp[0]*(shp[1]-burnIn)*shp[2]/nspace
    shape = np.shape(restart['chain'])
    ndim = shape[2]
    labels = [r"$\eta$",r"$\epsilon_\mathrm{ff}$",r"$f_{g,0}$",r"$\mu_0$",r"$\mu_{M_h}$",r"Q",r"$r_\mathrm{acc}/r_\mathrm{vir}$", \
            r"$f_\mathrm{cool}$", r"$M_{h,0}$", r"$\sigma$ (dex)", \
            r"$x_0$",r"$x_1$",r"$x_2$",r"$x_3$",r'Error scaling',r'c offset', r'$\mu_{h_g}$']
    labels = ['raccRvir', 'rstarRed', 'rgasRed', 'fg0mult', 'muColScaling', 'muFgScaling', 'muNorm', 'muMhScaling', 'ZIGMfac', 'zmix', 'eta', 'Qf', 'alphaMRI', 'epsquench', 'accCeiling', 'conRF', 'kZ', 'xiREC', 'epsff', 'scaleAdjust', 'mquench', 'enInjFac', 'alpharmh', 'fh'] 
    
    logs = [1,1,1,1,0,0,1,0,1,0,1,1,1,1,0,0,1,0,1,0,1,1,0,0]


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
    prior = np.array([samplefromprior() for i in range(np.product(np.shape(sampleRed))/ndim * 3)])


    pickle.dump(sampleRed, open('fakemcmc22_posterior.pickle', 'w'))

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
    trifig = corner.corner(sampleRed, labels=labels, range=extents)
    trifigPrior = corner.corner(prior, color='red', plot_datapoints=False, plot_filled_contours=False, fig=trifig, range=extents)
    trifig.savefig(fn+'.pdf')


    MaP = []
    stds = []

    #### just look at the 1D distr - much more compact when the 2D marginal distributions are pretty uninteresting!
    fig,ax = plt.subplots(6,4, figsize=(14,14))
    fig.subplots_adjust(hspace=0.2, wspace=0.01)
    for i in range(np.shape(sampleRed)[1]):
        hist, edges = np.histogram( sampleRed[:,i], bins='auto', range=extents[i], density=True )
        histplot(ax.flatten()[i], hist,edges, c='k', vertical=True)
        avg = np.mean(sampleRed[:,i])
        if logs[i]==1:
            MaP.append(10.0**avg)
        else:
            MaP.append(avg)
        std = np.std(sampleRed[:,i])
        stds.append(std)
        xv = np.linspace(extents[i][0],extents[i][1],200)
        yv = 1.0/np.sqrt(2.0*np.pi*std*std) * np.exp(-(xv-avg)*(xv-avg)/(2.0*std*std))
        ax.flatten()[i].plot( xv, yv, lw=2, ls='--', c='gray')
        hist, edges = np.histogram( prior[:,i], bins='auto', range=extents[i], density=True )
        histplot(ax.flatten()[i], hist,edges, c='r', vertical=False)
        ax.flatten()[i].set_xlabel(labels[i])
    plt.savefig(fn+'_1D.pdf')
    plt.close(fig)

    msg = "Estimate of maximum posterior: ["
    msg2 = " [ "
    for i in range(len(MaP)):
        msg+= str(MaP[i])+', '
        msg2+= str(stds[i])+', '
    msg = msg[:-2] + ']'
    msg2 = msg2[:-2] + ']'
    print msg
    print msg2

    print "Shape of sampleRed: ", np.shape(sampleRed)


def runEmcee(mpi=False, continueRun=False, seedWith=None):
    import sys
    if mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        from emcee.utils import MPIPool

    #models = [ pickle.load( open( 'rfnt09_'+str(k)+'_0.pickle', 'r' ) ) for k in range(80) ]


    if mpi:
        pool = MPIPool(comm=comm, loadbalance=True)
        if not pool.is_master():
            pool.wait()
            sys.exit()
    
    ndim, nwalkers = 24, 1000 
    # fn = 'fakemcmc17a_restart.pickle' ## a ran for a long time. "Standard" result
    fn = 'fakemcmc22_restart.pickle' ## Experimentally add a term in the likelihood to reproduce Krumholz&Burkhart data on MdotSF vs. \sigma.
    restart = {}
    nsteps = 3000 
    #p0 = [ samplefromprior(varRedFac=1.0) for w in range(nwalkers) ]
    p0 = [ sampleFromGaussianBall() for w in range(nwalkers) ]

    if seedWith is not None:
        if os.path.isfile(seedWith):
            with open(seedWith,'rb') as f:
                tmp_dict = pickle.load(f)
                seedPosition = tmp_dict['currentPosition']
                p0 = copy.deepcopy(seedPosition)
                # if we need more positions for walker seeds, try just averaging.
                if len(seedPosition)<nwalkers:
                    for i in range(nwalkers-len(seedPosition)):
                        selec = np.random.choice( range(len(p0)), size=2, replace=False )
                        p0.append( (np.array(p0[selec[0]] ) + np.array(p0[selec[1]]))/2.0 )

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
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)#, args=[models])
        #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool, args=[models])
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)#, args=[models])
        #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[models])
    
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


    if mpi:
        pool.close()


def singleRelationLikelihood(x,y,datasets, fixedSigma=-1):
    import observationalData
    lnlik = np.zeros(len(x))
    distances = np.zeros((len(x),len(datasets)))
    sigmas = np.zeros((len(x),len(datasets)))
    finalResiduals = np.zeros(len(x))
    for i,ds in enumerate(datasets):
        distance, sigma = observationalData.datasets[ds].distance(x,y, fixedSigma=fixedSigma)
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
    ### neededModels = [0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15, 16, 20,21,22,23, 24, 28,29,30,31, 32,33,34, 36, 40,41,42,43, 52,53,54,55, 72, 76,77,78,79]
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
    logc82z0  = Ys_train[:,36] 
    logSigma1z0  = Ys_train[:,40] 
    logSigma1z1  = Ys_train[:,41] 
    logSigma1z2  = Ys_train[:,42] 
    logSigma1z3  = Ys_train[:,43] 
    logmaxsigz0 = Ys_train[:,52] 
    logmaxsigz1 = Ys_train[:,53] 
    logmaxsigz2 = Ys_train[:,54] 
    logmaxsigz3 = Ys_train[:,55] 
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

    lnlik = np.zeros((len(Ys_train[:,0]),33) )

    srlInd=0
    if not returnlikelihood:
        srlInd=1
    lnlik[:,0] += singleRelationLikelihood(10.0**logMh0,10.0**logMst0,['Moster13z0'])[srlInd]
    lnlik[:,1] += singleRelationLikelihood(10.0**logMh1,10.0**logMst1,['Moster13z1'])[srlInd]
    lnlik[:,2] += singleRelationLikelihood(10.0**logMh2,10.0**logMst2,['Moster13z2'])[srlInd]
    lnlik[:,3] += singleRelationLikelihood(10.0**logMh3,10.0**logMst3,['Moster13z3'])[srlInd]

    lnlik[:,4] += singleRelationLikelihood(10.0**logMst0, 10.0**logsSFRz0, ['Speagle14Specz0', 'lilly13Specz0'])[srlInd]
    lnlik[:,5] += singleRelationLikelihood(10.0**logMst1, 10.0**logsSFRz1, ['Speagle14Specz1', 'lilly13Specz1'])[srlInd]
    lnlik[:,6] += singleRelationLikelihood(10.0**logMst2, 10.0**logsSFRz2, ['Speagle14Specz2', 'lilly13Specz2'])[srlInd]
    lnlik[:,7] += singleRelationLikelihood(10.0**logMst3, 10.0**logsSFRz3, ['Speagle14Specz3', 'lilly13Specz3'])[srlInd]
    #lnlik[:,4] += singleRelationLikelihood(10.0**logMst0,10.0**logsSFRz0,['Brinchmann04Specz0'])[srlInd]
    #lnlik[:,5] += singleRelationLikelihood(10.0**logMst1,10.0**logsSFRz1,['whitaker14Specz1'])[srlInd]
    #lnlik[:,6] += singleRelationLikelihood(10.0**logMst2,10.0**logsSFRz2,['whitaker14Specz2'])[srlInd]
    #lnlik[:,7] += singleRelationLikelihood(10.0**logMst3,10.0**logsSFRz3,['lilly13Specz3'])[srlInd]

    lnlik[:,8] += singleRelationLikelihood(10.0**logMst0,10.0**logstZz0,['gallazi05','kirby13'], fixedSigma=0.17)[srlInd]

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
    lnlik[:,22] += singleRelationLikelihood(10.0**logMst0,10.0**loghalfMassStarsz0,['vdW14LTG0'])[srlInd]
    lnlik[:,23] += singleRelationLikelihood(10.0**logMst1,10.0**loghalfMassStarsz1,['vdW14LTG1'])[srlInd]
    lnlik[:,24] += singleRelationLikelihood(10.0**logMst2,10.0**loghalfMassStarsz2,['vdW14LTG2'])[srlInd]
    lnlik[:,25] += singleRelationLikelihood(10.0**logMst3,10.0**loghalfMassStarsz3,['vdW14LTG3'])[srlInd]
    #lnlik[:,22] += singleRelationLikelihood(10.0**logMst0,10.0**loghalfMassStarsz0,['vdW14ETG0','vdW14LTG0'])[srlInd]
    #lnlik[:,23] += singleRelationLikelihood(10.0**logMst1,10.0**loghalfMassStarsz1,['vdW14ETG1','vdW14LTG1'])[srlInd]
    #lnlik[:,24] += singleRelationLikelihood(10.0**logMst2,10.0**loghalfMassStarsz2,['vdW14ETG2','vdW14LTG2'])[srlInd]
    #lnlik[:,25] += singleRelationLikelihood(10.0**logMst3,10.0**loghalfMassStarsz3,['vdW14ETG3','vdW14LTG3'])[srlInd]

    lnlik[:,26] += singleRelationLikelihood(10.0**(logMst0+loggasToStellarRatioHIz0),10.0**logbroeilsHIz0,['broeils97'])[srlInd]

    lnlik[:,27] += singleRelationLikelihood(10.0**logMst0,10.0**logc82z0,['Dutton09All'])[srlInd]

    lnlik[:,28] += singleRelationLikelihood(10.0**logMst0,10.0**logvPhi22z0,['miller11'])[srlInd]

    lnlik[:,29] += singleRelationLikelihood(10.0**(logMst0+logsSFRz0-9), 10.0**logmaxsigz0, ['krumholz17'], fixedSigma=0.1)[srlInd]
    lnlik[:,30] += singleRelationLikelihood(10.0**(logMst1+logsSFRz1-9), 10.0**logmaxsigz1, ['krumholz17'], fixedSigma=0.1)[srlInd]
    lnlik[:,31] += singleRelationLikelihood(10.0**(logMst2+logsSFRz2-9), 10.0**logmaxsigz2, ['krumholz17'], fixedSigma=0.1)[srlInd]
    lnlik[:,32] += singleRelationLikelihood(10.0**(logMst3+logsSFRz3-9), 10.0**logmaxsigz3, ['krumholz17'], fixedSigma=0.1)[srlInd]

    #lnlik[:,29] += singleRelationLikelihood(10.0**(logMst0+logsSFRz0), 

    #worstfit = np.argmin(lnlik[0,:])
    #print "Worst fit ", worstfit

    if np.any(np.isnan(lnlik)):
        pdb.set_trace()

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
    import observationalData
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

def trimDownForGlue(fn='broad1718b_to_lasso.txt', naccr=8, arr=None, kscale=10, logscale=0):
    X_train, X_validate, X_test,  Ys_train, Ys_validate, Ys_test, labels = readData(trainFrac=1.0, validateFrac=0.0, fn=fn, naccr=naccr, arr=arr, kscale=kscale, logscale=logscale, includeCrossTerms=False)
    labelsFeat = ["Mh0", "raccRvir", "rstarRed", "rgasRed", "fg0mult", "muColScaling", "muFgScaling", "muNorm", "muMhScaling", "ZIGMfac", "zmix", "eta", "Qf", "alphaMRI", "epsquench", "accCeiling", "conRF", "kZ", "xiREC", "epsff", "deltabeta", "MQ", "Einj", "Mdot0KS" ]
    labelsFeatAccr = ["Mdot"+str(i).zfill(2) for i in range(naccr)]
    labelsVar = ["Mh0", "Mh1", "Mh2", "Mh3", "Mst0", "Mst1", "Mst2", "Mst3", "sSFR0", "sSFR1", "sSFR2", "sSFR3", "Zg0", "Zg1", "Zg2", "Zg3", "Zst0", "Zst1", "Zst2", "Zst3", "fgmol0", "fgmol1", "fgmol2", "fgmol3", "fgHI0", "fgHI1", "fgHI2", "fgHI3", "rst0", "rst1", "rst2", "rst3", "vphitf0", "vphitf1", "vphitf2", "vphitf3", "c0", "c1", "c2", "c3", "SigCen0", "SigCen1", "SigCen2", "SigCen3", "jst0", "jst1", "jst2", "jst3", "Zgrad0", "Zgrad1", "Zgrad2", "Zgrad3", "sigmax0", "sigmax1", "sigmax2", "sigmax3", "mdotb0", "mdotb1", "mdotb2", "mdotb3", "rGI0", "rGI1", "rGI2", "rGI3", "tdep0", "tdep1", "tdep2", "tdep3", "tdepmol0", "tdepmol1", "tdepmol2", "tdepmol3", "rHI0", "rHI1", "rHI2", "rHI3", "Mmerge0", "Mmerge1", "Mmerge2", "Mmerge3", "vphi00", "vphi01", "vphi02", "vphi03", "vphi04", "vphi05", "vphi06", "vphi07", "vphi08", "vphi09", "vphi10", "vphi11", "vphi12", "vphi13", "vphi14", "vphi15", "vphi16", "vphi17", "vphi18", "vphi19", "Sigma00", "Sigma01", "Sigma02", "Sigma03", "Sigma04", "Sigma05", "Sigma06", "Sigma07", "Sigma08", "Sigma09", "Sigma10",  "Sigma11", "Sigma12", "Sigma13", "Sigma14", "Sigma15", "Sigma16", "Sigma17", "Sigma18", "Sigma19",   "SigmaSt00", "SigmaSt01", "SigmaSt02", "SigmaSt03", "SigmaSt04", "SigmaSt05", "SigmaSt06", "SigmaSt07", "SigmaSt08", "SigmaSt09", "SigmaSt10",  "SigmaSt11", "SigmaSt12", "SigmaSt13", "SigmaSt14", "SigmaSt15", "SigmaSt16", "SigmaSt17", "SigmaSt18", "SigmaSt19",  "SigmaSFR00", "SigmaSFR01", "SigmaSFR02", "SigmaSFR03", "SigmaSFR04", "SigmaSFR05", "SigmaSFR06", "SigmaSFR07", "SigmaSFR08", "SigmaSFR09", "SigmaSFR10", "SigmaSFR11", "SigmaSFR12", "SigmaSFR13", "SigmaSFR14", "SigmaSFR15", "SigmaSFR16", "SigmaSFR17", "SigmaSFR18", "SigmaSFR19", "Zr00", "Zr01", "Zr02", "Zr03", "Zr04", "Zr05", "Zr06", "Zr07", "Zr08", "Zr09", "Zr10", "Zr11", "Zr12", "Zr13", "Zr14", "Zr15", "Zr16", "Zr17", "Zr18", "Zr19", "Age00", "Age01", "Age02", "Age03", "Age04", "Age05", "Age06", "Age07", "Age08", "Age09", "Age10", "Age11", "Age12", "Age13", "Age14", "Age15", "Age16", "Age17", "Age18", "Age19"  ]

    labels = labelsFeat+labelsFeatAccr+labelsVar
    header = ""
    for label in labels:
        header += label + " "
    #header = header.replace('0','nZEn').replace('1','nIn').replace('2','nIIn').replace('3','nIIIn').replace('4','nIVn').replace('5','nVn').replace('6','nVIn').replace('7','nVIIn').replace('8','nVIIIn').replace('9','nIXn')
    arr = np.vstack( [X_train.T, Ys_train.T] ).T


    np.savetxt('broad1718b_to_lasso_trim.txt', arr, header=header[:-1])



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

            
def quickPlot(arr, name):
    fig,ax = plt.subplots()
    hist, edges = np.histogram( arr, bins='auto', density=True )
    histplot(ax, hist,edges, c='r', vertical=False)
    plt.savefig(name)
    plt.close(fig)


def estimateFeatureImportances(analyze=True, pick=True):
    ### Start with a model we like after our experience w/ searchTreeParams
    ### Train and score the model
    ### Predict Y, but shuffle values of one feature at a time in one mass bin at a time.
    ### Record score reduction at each shuffle.
    ### Plot score reduction as fn of mass for each feature.
    from sklearn.metrics import r2_score

    X_train_orig, X_validate, X_test_orig, Ys_train_orig, Ys_validate, Ys_test_orig, labels = readData(trainFrac=0.99, validateFrac=0, naccr=8, fn='broad22partial_to_lasso.txt') # no need to feed in arr, since we're just reading the data once.
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

    #Mh0, raccRvir, rstarRed, rgasRed, fg0mult, muColScaling, muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix, eta, Qf, alphaMRI, epsquench, accCeiling, conRF, kZ, xiREC, epsff, scaleAdjust, mquench, enInjFac, chiZslope = emceeparams

    feature_names = [r'$M_{h,0}$', r'$\alpha_r$', r'$\alpha_{r,*,0}$',  r'$\alpha_{r,g,0}$',  r'$\chi_{f_{g,0}}$', r'$\alpha_\Sigma$', r'$\alpha_{f_g}$', r'$\mu_0$', r'$\alpha_{M_h}$', r'$\chi_{Z_\mathrm{IGM}}$', r'$\xi_\mathrm{acc}$', r'$\eta$', r'$Q_f$', r'$\alpha_\mathrm{MRI}$', r'$\epsilon_\mathrm{quench}$', r'$\epsilon_\mathrm{ceil}$', r'$\alpha_\mathrm{con}$', r'$k_Z$', r'$\xi$', r'$\epsilon_\mathrm{ff}$', r'$\Delta\beta$', r'$M_Q$', r'$\chi_\mathrm{inj}$, $\chi_{dlogZ/dlogM}$']
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
            
            #### store a set of random factors to be used all throughout the emcee process
            randomFactors = np.zeros(np.shape(randomFactorsKey01))
            for jj in range(50):
                randomFactors[jj,:] = np.array([X_train[ int(randomFactorsKey01[jj,ii]*len(X_train[:,0])),-9+ii] for ii in range(9)])

            Ytra = xtransform(Ys_train)
            Ys_train = Ytra.transform(Ys_train)
            Ys_validate = Ytra.transform(Ys_validate)
            Ys_test = Ytra.transform(Ys_test_orig.copy())


            # k=8 corresponds to z=0 sfr. Try really hard to get this right!
            errors_train_this, errors_validate_this, errors_test_this, labels_this, feature_importances_this, theModel = learnRF(X_train, X_validate, X_test, Ys_train, Ys_validate, Ys_test, labels, n_estimators=100, k=k, max_depth=1000, max_features='auto', min_per_leaf=3 )
            if pick:
                pickle.dump( fntModel(theModel,Xtra,Ytra,randomFactors) , open('rfnt22_'+str(k)+'_'+str(cvi)+'.pickle','w')) ### save the model

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

    ## assemble the tex labels for the tables below
    texlabels=[]
    tl0=[r'$M_h$', r'$M_*$', r'sSFR', r'$\langle Z\rangle_\mathrm{SFR}$', r'$Z_*$', r'$M_{\mathrm{H}_2}/M_*$', r'$M_{\mathrm{HI}}/M_*$', r'$r_*$', r'$v_{\phi, 2.2}$', r'$c_{82}$', r'$\Sigma_1$', r'$j_*$', r'$\partial\log_{10}Z/\partial r$', r'$\sigma_\mathrm{max}$', r'$\dot{M}_{\mathrm{bulge, gas}}$', r'$r_{GI}/r_*$', r'$t_\mathrm{dep}$', r'$t_{\mathrm{dep},\mathrm{H}_2}$', r'$r_{\mathrm{HI}}$', r'$M_{*,\mathrm{merge}}$']
    for tl in tl0:
        texlabels.append(tl)
        texlabels.append(' ')
        texlabels.append(' ')
        texlabels.append(' ')
        #texlabels.append( tl+ r'$(z=0)$')
        #texlabels.append( tl+ r'$(z=1)$')
        #texlabels.append( tl+ r'$(z=2)$')
        #texlabels.append( tl+ r'$(z=3)$')
    tl1 = [r'$v_\phi$', r'$\Sigma$', r'$\Sigma_*$', r'$\Sigma_\mathrm{SFR}$', r'$Z$', r'Age']
    for tl in tl1:
        texlabels.append(tl)
        for ri in range(1,20):
            #texlabels.append( tl+r'$(r/r_*='+ str(np.round(float(ri)/5.0, 2)) +r')$' )
            texlabels.append( '' )
    if analyze:
        printFeatureImportances(feature_importances.T*1.0, r2scores, feature_names, texlabels, normalize=False)
        printFeatureImportances(np.mean(scores,axis=1).T, r2scores, feature_names, texlabels, fn='avg_swap_scores_compact.tex', normalize=False)
        printFeatureImportances(np.mean(derivs_mass,axis=2).T, r2scores, feature_names, texlabels, fn='avg_derivs_compact.tex')
        printTable(feature_importances.T*1.0, feature_names, texlabels, fn='feature_importances.tex')
        printTable(np.mean(scores,axis=1).T, feature_names, texlabels, fn='avg_swap_scores.tex')

    #print "Average feature importances from skl: ",feature_importances

def printFeatureImportances(arr, rsquared, columnnames, rownames, fn='feature_importance_compact.tex', normalize=True):
    assert len(np.shape(arr))==2
    nrows, ncols = np.shape(arr)
    #if not (nrows==len(rownames) and ncols==len(columnnames)):
    #    pdb.set_trace()
    colorNames = ['blue', 'olive', 'CadetBlue', 'pink', 'BlueGreen', 'ForestGreen', 'LimeGreen', 'gray', 'orange', 'yellow', 'Magenta', 'cyan', 'Maroon', 'purple', 'lightslategray', 'Peach', 'Fuchsia', 'Orchid', 'ProcessBlue', 'Tan', 'Sepia', 'red', 'RawSienna']*2

    with open(fn,'w') as f:
        f.write( r'\begin{table*}'+'\n' )
        f.write( r'\tiny'+'\n' )
        f.write( r'\begin{tabular}{ll|'+'l'*5+r'|l'+r'}'+'\n' )
        line = r' & z & Feat. 0 & Feat. 1 & Feat. 2 & Feat. 3 & Feat. 4 & R$^2$' ## the 5 most important features and their scores
        line+=r'\\'+'\n'
        f.write(line)
        f.write(r'\hline'+'\n')
        f.write(r'\hline'+'\n')
        for j in range(nrows):
            if j<=80 or ( j>80 and j%5==0 ):
                if j==80:
                    ### split the table between integrated quantities and radial quantities
                    f.write(r'\end{tabular}'+'\n')
                    f.write(r'\end{table*}'+'\n')

                    f.write( r'\begin{table*}'+'\n' )
                    f.write( r'\tiny'+'\n' )
                    f.write( r'\begin{tabular}{ll|'+'l'*5+r'|l'+r'}'+'\n' )
                    line = r' & $r/r_*$ & 0 & 1 & 2 & 3 & 4 & R$^2$' ## the 5 most important features and their scores
                    line+=r'\\'+'\n'
                    f.write(line)
                    f.write(r'\hline'+'\n')
                    f.write(r'\hline'+'\n')

                # only do the radially resolved quantities at a few points 
                line = rownames[j]

                radii = ['0.1', '1.0', '1.9', '2.8']
                line += r' & '+ radii[j%4]

                #for i in range(ncols):
                zipped = zip(arr[j,:], range(ncols))
                zsorted = sorted(zipped, key=lambda x: np.abs(x[0]))
                for i in range(5):
                    if zsorted[-i-1][0]!=zsorted[-i-1][0] or zsorted[-i-1][1]!=zsorted[-i-1][1]:
                        pdb.set_trace()
                    if normalize:
                        line += r' & \cellcolor{'+colorNames[zsorted[-i-1][1]]+r'!'+str(int(min([abs(zsorted[-i-1][0])/(0.1+np.abs(zsorted[-1][0]))*100,99])))+r'} '+columnnames[zsorted[-i-1][1]] + r': ' + str(np.round(zsorted[-i-1][0],2)) 
                    else:
                        line += r' & \cellcolor{'+colorNames[zsorted[-i-1][1]]+r'!'+str(int(min([abs(zsorted[-i-1][0])*100,99])))+r'} '+columnnames[zsorted[-i-1][1]] + r': ' + str(np.round(zsorted[-i-1][0],2)) 
                line+= r'& ' +str(np.round(rsquared[j],2))+r' \\'+'\n'
                f.write(line)
                if j%4==3:
                    f.write(r'\hline'+ '\n')
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
def letter(i, caps=False):
    if caps:
        return chr(ord("A")+i)
    else:
        return chr(ord("a")+i)

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

    runEmcee(mpi=True, continueRun=False, seedWith='fakemcmc22_restart.pickle' )
    #fractionalVariancePlot()
    #ridgeCoeffsPlot()

    #validateNPR()
    #plotResiduals()

    #estimateFeatureImportances(analyze=True, pick=True) # just generate the pickled models 
    #estimateFeatureImportances(analyze=True, pick=False) # do the analysis but don't save the models
    
    #nuclearSearch(Nbins = 7, Niter=10000, Ninits=50)

    #appendLikelihood()
    
    #fakeEmceePlotResiduals(None, 'fakemcmc13_residuals', gidgetmodels='rf79', xmax=[ 1.41977729e-01,   2.03502540e+00,   2.02360753e+00,   2.01127027e+00, 2.56584739e-02,  -1.89348632e-02,   1.07203799e+00,  -9.80149073e-01, 1.00975501e+00,   4.99324933e-01,   1.50177845e+00,   1.48359195e+00, 4.94927969e-02,  9.79999389e-04,  9.98857856e-01,   2.85701876e-01, 9.59552777e-01,   9.02414277e-05,   9.91998674e-03, 0.5] )


    #fakeEmceePlotResiduals(None, 'fakemcmc13_residuals81', gidgetmodels='rf81', xmax=[ 0.3512968,   3.91202004,  3.41454232,  3.22847406, -0.52552445,  2.60129093, 4.40972677, -1.47069249,  1.25916476,  0.63954234,  1.82215254,  1.46210786, 0.05488134,  0.04106005,  0.24240757,  0.47092227,  2.25026276,  0.36920633, 0.01555641, -1.0, 2.0e12, 0.8 ] )
    #fakeEmceePlotResiduals(None, 'fakemcmc15_residuals86', xmax=  [  1.39192613e-01,   2.06852503e+00,   1.47353675e+00,   1.97365750e+00, 5.37224908e-01,   1.17595949e+00,   1.08963518e+00,  -5.08071210e-01, 1.01502597e+00,   5.17380183e-02,   5.36199804e+00,   1.25992157e+00, 6.23709052e-02,   1.72535928e-04,   9.15648977e-01,   1.06622955e-01, 5.88668979e+00,   1.79678213e-01,   6.00854464e-03,   0.00000000e+00, 9.18293025e+11, .3536  ], gidgetmodels='rf86' )



    #fakeEmceePlotResiduals(None, 'fakemcmc16_residuals88', gidgetmodels='rf88', xmax=[  1.16493377e-01,   6.85441565e-01,   2.53454237e+00,   4.66205970e+00, 7.21672391e-01,  -9.65895952e-01,   3.31207773e+00,  -8.40086055e-01, 5.10900160e+00,   4.37085835e-01,   2.12261668e+00,   1.59251012e+00, 8.54462905e-02,   3.48195437e-03,   4.07820701e-01,   1.84455124e-01, 6.85178180e-01,   2.17577982e-01,   8.94197775e-03,   2.00963528e-01, 1.18649645e+12,   5.11828699e-01] )

    #fakeEmceePlotResiduals(None, 'fakemcmc17_residuals99', gidgetmodels='rf99', xmax=[  1.76149012e-01, 2.71385761e+00, 2.43095960e+00, 2.56258037e+00, -5.22815474e-01, -7.30050429e-01, 5.96476522e+00, -1.95328561e+00, 1.13920495e+00, 6.26490351e-01, 1.57867191e+00, 1.77166604e+00, 6.13678555e-02, 6.36513703e-03, 5.42225312e-01, 1.54471933e-01, 1.69720803e+00, 4.65370474e-01, 1.28680215e-02, 4.79132819e-01, 1.64818086e+12, 9.93276414e-01, 5.41798936e-01])


    #fakeEmceePlotResiduals( None, 'fakemcmc17d_residuals112', gidgetmodels='rf112', xmax = [ 1.44656471e-01, 1.94694047e+00, 2.46300115e+00, 2.52780687e+00, 4.14257376e-02, -9.24764117e-02, 6.03022138e+00, -1.71412505e+00, 2.75050207e+00, 5.03191815e-01, 2.64146528e+00, 2.93923418e+00, 4.93805277e-02, 1.33649809e-02, 6.30847392e-01, 7.70793830e-02, 6.40127324e-01, 2.75077929e-01, 1.67309700e-02, 4.98795340e-01, 1.93406046e+12, 7.82157806e-01, -1.59310290e-01, 5.93635084e-01] )

    #fakeEmceePlotResiduals( None, 'fakemcmc17d_residuals113', gidgetmodels='rf113', xmax =[1.21990575e-01, 2.25750440e+00, 1.42283492e+00, 2.38779672e+00, -2.20443849e+00, 1.07448541e+00, 2.26876611e+00, -1.11338069e+00, 1.91939912e+00, 4.52858032e-01, 1.84182587e+00, 2.26660250e+00, 7.71741216e-02, 1.21665083e-02, 4.95477220e-01, 2.57736306e-01, 1.06866306e+00, 5.03555469e-01, 8.98255341e-03, 3.51823670e-01, 6.57947163e+11, 9.73425083e-01, -1.69906882e-01, 3.10019810e-01] )

    #### a few runners-up
    posteriors = np.array([[  1.41770654e-01,   2.45050659e+00,   2.61690983e+00,   2.63513509e+00,   3.22385410e-01,   5.25370163e-01,   4.11019928e+00,  -5.56706799e-01,   1.62770773e+00,   4.38315726e-01,   2.43862201e+00,   2.64597948e+00,   5.00168386e-02,   8.37982476e-03,   4.57828445e-01,   2.16816496e-01,   1.24032620e+00,   3.66206876e-01,   1.16197763e-02,   5.13078077e-01,   2.45445719e+12,   1.46320093e+00,  -1.88561103e-02,   3.74693962e-01],
    [  1.43471049e-01,   2.11586724e+00,   2.85565807e+00,   1.83463928e+00,   4.97074915e-02,   1.20947662e+00,   2.96755583e+00,  -7.51240372e-01,   2.73899696e+00,   6.10208444e-01,   2.38562845e+00,   2.03582376e+00,   6.21513309e-02,   4.34744165e-03,   5.12119376e-01,   1.04256577e-01,   2.12898835e+00,   4.94798946e-01,   1.02598907e-02,   3.47843175e-01,   1.36955173e+12,   1.18497000e+00,  -1.67926958e-01,   4.96547952e-01],
    [  2.04486864e-01,   2.88120002e+00,   3.25719889e+00,   2.66263057e+00,   2.27151192e-01,   5.25564583e-01,   3.42627948e+00,  -1.03353919e+00,   2.00688887e+00,   6.06262384e-01,   2.37107780e+00,   2.22078909e+00,   6.92269496e-02,   1.32278564e-02,   6.33401611e-01,   3.32810752e-01,   1.58332202e+00,   2.94782593e-01,   1.12429678e-02,   2.47279896e-01,   2.29634718e+12,   1.51504655e+00,  -5.87451991e-02,   5.61659119e-01],
    [  1.75968697e-01,   2.63077357e+00,   2.71140824e+00,   1.98690602e+00,   3.02698953e-01,   7.64297819e-01,   6.24096647e+00,  -1.15052330e+00,   2.38752909e+00,   6.95205531e-01,   2.13169580e+00,   2.27084757e+00,   7.66656260e-02,   9.82512104e-03,   6.55021191e-01,   2.21648159e-01,   1.68354039e+00,   3.57178264e-01,   8.93852873e-03,   3.99091464e-01,   1.98464589e+12,   1.45645035e+00,  -8.54418687e-02,   5.85790120e-01],
    [  1.50001597e-01,   6.63368287e-01,   3.46505862e+00,   5.29186501e-01,  -1.02397187e-01,  -1.33669354e-01,   6.30395899e+01,  -1.52932101e+00,   1.06459709e+00,   7.95051593e-01,   3.78107677e+00,   2.54264572e+00,   3.39551038e-02,   1.18020216e-03,   8.32316414e-01,   2.55190390e-01,   2.30751383e+00,   1.91278415e-01,   2.68669899e-02,   4.71418276e-01,   2.28803783e+11,   5.72407369e-01,  -1.32584911e-01,   8.54109282e-01],
    [  1.45906622e-01,   1.75985998e+00,   2.42484834e+00,   2.18460823e+00,   1.71358687e-01,  -1.26060221e+00,   3.30462760e+00,  -9.17273308e-010,   1.27074217e+00,   5.94169090e-01,   2.43907682e+00,   1.59752000e+00,   4.92771728e-02,   1.36675462e-02,   4.83738660e-01,   8.73159481e-02,   1.38621386e+00,   3.28886393e-01,   1.21186695e-02,   3.55451628e-01,   1.29547055e+12,   1.21522829e+00,  -8.24133283e-02,   4.62554955e-01],
     [  1.60878849e-01,   2.62023173e+00,   2.37961236e+00,   1.90923724e+00,  -2.17309742e+00,  -2.21250594e-01,   1.55357515e+00,  -1.26347273e+00,   2.47687330e+00,   4.64731828e-01,   1.59314023e+00,   2.39710092e+00,   7.55143918e-02,   1.38945365e-02,   5.36957565e-01,   3.50458248e-01,   1.41354094e+00,   3.66901507e-01,   8.88629028e-03,   1.66394318e-01,   1.16081408e+12,   1.04401821e+00,  -2.28896855e-01,   5.86472929e-01],
     [  1.02767464e-01,   2.82837515e+00,   2.39027386e+00,   2.28921240e+00,  -9.66654109e-01,  -2.51765156e+00,   1.78508338e+00,  -1.10373036e+00,   1.26654765e+00,   7.95578969e-01,   1.95478174e+00,   2.63760284e+00,   5.25470711e-02,   2.29941424e-03,   3.80650455e-01,   2.82853531e-01,   9.08476705e-01,   4.28674852e-01,   8.97324258e-03,   8.24456568e-01,   6.90571984e+11,   8.86126290e-01,  -1.07351512e-01,   8.39403757e-01],
     [  1.60220601e-01,   2.81426435e+00,   1.78223572e+00,   1.95563642e+00,  -5.48371939e-01,  -2.61916907e-02,   4.00486058e+00,  -1.07005557e+00,   1.57033180e+00,   4.69179408e-01,   2.26715714e+00,   2.50807400e+00,   6.75841381e-02,   6.02409794e-03,   4.79563831e-01,   3.22839595e-01,   2.12140409e+00,   3.97247225e-01,   1.01925116e-02,   5.68737281e-01,   9.74726404e+11,   1.50490111e+00,  -5.43194737e-02,   3.71735915e-01]])
    for k in range(9):
        pass
        #try:
        #    fakeEmceePlotResiduals( None, 'fakemcmc17d_residualsrf113'+letter(k,caps=True), gidgetmodels='rf113'+letter(k,caps=True), xmax=posteriors[k,:] )
        #except:
        #    print "Something bad happened for k=",k," i.e. rf113"+letter(k,caps=True), " - probably all models failed in that run"


    #fakeEmceePlotResiduals( None, 'fakemcmc17d_residuals114', gidgetmodels='rf114', xmax=[0.147451569889, 2.52941811659, 2.59034734186, 2.41120695741, -0.124283831858, -0.0523679435879, 9.35680382698, -0.974822093888, 1.89905286619, 0.511551421578, 2.0337488747, 2.02929369251, 0.0642244458824, 0.00988965146683, 0.509787819545, 0.279394476293, 1.74417214913, 0.342311450585, 0.0107366934282, 0.414472066814, 1.50430105103e+12, 1.21653967757, -0.028389697679, 0.497607252288] )

    #trimDownForGlue()

    ### analyze the fake mcmc run

    if False:
        restart={}
        updateRestart('fakemcmc1718b_restart.pickle', restart)
        printRestart(restart)
        tracePlots(restart, 'fakemcmc1718b_trace', burnIn=0)
        probsPlots(restart, 'fakemcmc1718b_allProb', burnIn=0)
        trianglePlot(restart,'fakemcmc1718b_triangle', burnIn=25, nspace=5)

        # Find the maximum among all models sampled so far.
        allProbs = restart['allProbs'].flatten() ## allprobs is presumably nwalkers*niterations
        priors = np.zeros((np.shape(restart['allProbs'])[0],np.shape(restart['allProbs'])[1]))
        for i in range(np.shape(restart['allProbs'])[0]):
            for j in range(np.shape(restart['allProbs'])[1]):
                priors[i,j] = lnprior( restart['chain'][i,j,:] )

        stlikelihoods = restart['allProbs'] - priors
        likelihoods = allProbs.flatten() - priors.flatten()
        print "likelihoods stats: ", np.max(likelihoods), np.percentile(likelihoods,[5,50,90,95,97.5,99])

        print "probs stats: ", np.max(allProbs), np.percentile(allProbs,[5,50,90,95,97.5,99])

        zipped = zip(restart['allProbs'][:,25], range(len(restart['allProbs'][:,25])))
        zsorted = sorted(zipped, key=lambda x: -x[0])
        highProbInds = [zsorted[i][1] for i in range(10)]

        zipped2 = zip(stlikelihoods[:,25], range(len(stlikelihoods[:,25])))
        zsorted2 = sorted(zipped2, key=lambda x: -x[0])
        highProbInds2 = [zsorted2[i][1] for i in range(10)]

        for i in range(2):
            #index = np.argmax(allProbs)
            #index2 = np.argmax(likelihoods)
            index = highProbInds[i]
            print 'posterior:', restart['allProbs'][index,25]
            print 'likelihood:', stlikelihoods[index,25]
            ##print 'likelihood: ', likelihoods[index2]
            #indices = np.unravel_index(index, np.shape(restart['allProbs']))
            ##indices2 = np.unravel_index(index2, np.shape(restart['allProbs']))
            #xmax = restart['chain'][indices[0],indices[1],:]
            # (nwalkers x niterations x ndim)    
            xmax = restart['chain'][index, 25, :]
            #xmax2 = restart['chain'][indices2[0],indices2[1],:]
            print "Favorite coordinates Posterior (", i, ")", xmax
            printPriorOffsets(xmax)
            #print "Favorite coordinates Likelihood: ", xmax2
            #printPriorOffsets(xmax2)

        #print " current shape of chain: (nwalkers x niterations x ndim) ",np.shape(restart['chain'])



        #fakeEmceePlotResiduals(restart, 'fakemcmc13_residuals', gidgetmodels='rf78')
