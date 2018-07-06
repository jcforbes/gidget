import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import halo
import pdb
import copy
from scipy.interpolate import interp1d
import pickle

def standardizedEpsSkewNormal(x, epsilon):
    if not hasattr(x, '__iter__'):
        if x<0:
            return 1/np.sqrt(2*np.pi) * np.exp(-x*x/(2*(1+epsilon)*(1+epsilon)))
        else:
            return 1/np.sqrt(2*np.pi) * np.exp(-x*x/(2*(1-epsilon)*(1-epsilon)))
    else:
        xThis = np.array(x)
        epsThis = np.array(epsilon)
        ltz = xThis < 0
        ret = np.zeros(np.shape(xThis))
        ret[ltz] = 1/np.sqrt(2*np.pi) * np.exp(-xThis[ltz]*xThis[ltz]/(2*(1+epsThis[ltz])*(1+epsThis[ltz])))
        nltz = np.logical_not(ltz)
        ret[nltz] = 1/np.sqrt(2*np.pi) * np.exp(-xThis[nltz]*xThis[nltz]/(2*(1-epsThis[nltz])*(1-epsThis[nltz])))
        return ret
def epsSkewNormalPDF(x, epsilon, theta, sigma):
    return standardizedEpsSkewNormal( (x-theta)/sigma, epsilon )/sigma
from scipy.stats import norm
def standardizedEpsSkewQuantile(u, epsilon):
    ret = 0
    if 0<u and u<(1.0+epsilon)/2.0:
        ret = (1.0+epsilon) * norm.ppf(u/(1.0+epsilon))
    elif (1.0+epsilon)/2.0 <= u and u<1:
        ret = (1.0-epsilon) * norm.ppf((u-epsilon)/(1.0-epsilon))
    else:
        raise ValueError
    if ret!=ret:
        pdb.set_trace()
    #print "stEpsNorm quantile: ", ret
    return ret
from scipy.optimize import brentq
def computeEpsSkewParams(q16, q50, q84):
    def to_zero(epsilon):
        ret =  (standardizedEpsSkewQuantile( 0.84, epsilon ) - standardizedEpsSkewQuantile(0.5, epsilon)) / (standardizedEpsSkewQuantile( 0.5, epsilon ) - standardizedEpsSkewQuantile(0.16, epsilon)) - (q84-q50)/(q50-q16)
        #print "eps_tozero:", ret
        return ret
    #print "Attempting brentq with ",q16,q50,q84,';', to_zero(-.99), to_zero(0), to_zero(.99)
    lower_bound = -.9
    middle_value = 0
    upper_bound = 0.9
    if( to_zero(lower_bound)*to_zero(upper_bound)>0 ):
        vlb = to_zero(lower_bound)
        vub = to_zero(upper_bound)
        vlz = to_zero(middle_value)
        if(abs(vlb)<abs(vub) and abs(vlb)<abs(vlz)):
            #print "WARNING: returning lower_bound" 
            x0 = lower_bound
        elif(abs(vub)<abs(vlb) and abs(vub)<abs(vlz)):
            #print "WARNING: returning upper_bound" 
            x0 = upper_bound
        else:
            #print "WARNING: returning epsilon=", middle_value," without minimizing"
            x0 = middle_value
    else:
        x0 = brentq(to_zero, lower_bound, upper_bound)
#    counter = 0
#    countermax = 100
#    while(to_zero(lower_bound)*to_zero(upper_bound)>0):
#        counter+=1
#        lower_bound = np.random.random()*2-1
#        upper_bound = np.random.random()*(1-lower_bound)-lower_bound 
#        if counter>countermax:
#            print "Failed to find a good bound!"
#            trial_epsilons = np.linspace(-.999,.999,50)
#            objectiveFn = [to_zero(eps) for eps in trial_epsilons]
#            print trial_epsilons
#            print objectiveFn
#            pdb.set_trace()
#            break
#    print counter
#    if counter>countermax:
#        x0 = 0
#        print "WARNING: setting epsilon to zero"
#    else:
#        x0 = brentq(to_zero, lower_bound, upper_bound)
    #print "x0: ",x0
    sigma = (q50-q16)/(standardizedEpsSkewQuantile( 0.5, x0) - standardizedEpsSkewQuantile(0.16, x0))
    #print "sigma: ", sigma
    theta = q50 - sigma*standardizedEpsSkewQuantile( 0.5, x0)
    #print "theta: ", theta
    return x0,theta,sigma


def nToC82(nn):
    # Integral of x e^(-x^(1/n)) = -n Gamma(2n, x^(1/n))
    # Note! that the incomplete gamma function in the above is NOT normalized, whereas the scipy version is.
    # ALSO note that the scipy incomplete gamma function is integrated from 0, whereas the mathematica one is integrated from infinity

    #   -n Gamma(2n, x^(1/n))
    #  = -n Gamma(2n, x^(1/n)) + n Gamma(2n, 0)
    #  = n * ( Gamma(2n, 0) - Gamma(2n, x^(1/n)) )
    #  = n * ( scipy.gamma(2n) - scipy.gamma(2n)*(1-scipy.gammainc(2n,x^(1/n))) )
    #  = n * scipy.gamma(2n) ( 1 - (1-scipy.gammainc(2n,x^(1/n))) )
    #  = n * scipy.gamma(2n)*( scipy.gammainc(2n,x^(1/n))) 
    from scipy.special import gammainc
    def cumulativeMass(x):
        return  gammainc( 2.0*nn, np.power(x,1.0/nn))
    def toZero(x):
        return cumulativeMass(x) - 0.2
    from scipy.optimize import brentq
    r20 = brentq(toZero, 0.0, 100.0**nn)
    def toZero(x):
        return cumulativeMass(x) - 0.8
    r80 = brentq(toZero, 0.0, 100.0**nn)
    return r80/r20
# r.i.p.
#def nToC82(nn):
#    ''' Convert a sersic index to a concentration index (the radius containing 80% of the mass to the radius containing 20%)'''
#    def integrand(x):
#        return x*np.exp(-np.power(x,1.0/nn))
#    from scipy.integrate import quad
#    result = quad( integrand, 0, 100)
#    totalmass = result[0]
#    from scipy.optimize import brentq
#    def toZero(x):
#        result = quad(integrand, 0, x)
#        return result[0]/totalmass - 0.8
#    r80 = brentq(toZero, .01, 10)
#    def toZero(x):
#        result = quad( integrand, 0, x) 
#        return result[0]/totalmass - 0.2
#    r20 = brentq(toZero, .01, 10)
#    return r80/r20



# This is the simplest thing that will currently do well in my setup. 
# 
class DataSet:
    def __init__(self, xvar, yvar, xval, yval, yLower=None, yUpper=None, xerr=None, yerr=None, zmin=-1, zmax=0.5, label=None, logx=True, logy=True, alpha=0.05, fixedSigma=-1, extrapolate=True):
        if label is None:
            pdb.set_trace()
        self.xvar = xvar
        self.yvar = yvar
        self.xval = xval
        self.yval = yval
        self.xerr = xerr
        if yerr is not None:
            if yUpper is None:
                yUpper = yval+yerr
            if yLower is None:
                yLower = yval-yerr
        # If nothing filled in yUpper and yLower, assume 50% errors
        if yUpper is None:
            yUpper = np.array(yval)*1.5
        if yLower is None:
            yLower = np.array(yval)/1.5
        self.yLower=yLower
        self.yUpper=yUpper

        self.alpha=alpha

        self.zmin = zmin
        self.zmax = zmax
        self.label = label
        self.logx = logx
        self.logy = logy

        self.fixedSigma = fixedSigma
        self.extrapolate = extrapolate
    def plotLikelihood(self, axIn=None):
        x0=None
        x1=None
        ex0=None
        ex1=None
        ey0=None
        ey1=None
        if self.logx:
            dynRange = np.max(self.xval)/np.min(self.xval)
            if x0 is None:
                x0 = np.min(self.xval)/dynRange**1.0
            if x1 is None:
                x1 = np.max(self.xval)*dynRange**1.0
            ex0=np.log10(x0)
            ex1=np.log10(x1)
        else:
            linRange = np.max(self.xval)-np.min(self.xval)
            if x0 is None:
                x0 = np.min(self.xval) - linRange*1.0
            if x1 is None:
                x1 = np.max(self.xval) + linRange*1.0
            ex0=x0
            ex1=x1

        if self.logx:
            # x0 and x1 are in linear space - take the log to produce even spacing in log space
            xs = np.linspace(np.log10(x0), np.log10(x1), 100)
        else:
            xs = np.linspace(x0, x1, 100)
        y0=None
        y1=None
        if self.logy:
            dynRange = np.max(self.yval)/np.min(self.yval)
            if y0 is None:
                y0 = np.min(self.yval)/dynRange**5.0
            if y1 is None:
                y1 = np.max(self.yval)*dynRange**5.0
            ey0=np.log10(y0)
            ey1=np.log10(y1)
        else:
            linRange = np.max(self.yval)-np.min(self.yval)
            if y0 is None:
                y0 = np.min(self.yval) - linRange*5.0
            if y1 is None:
                y1 = np.max(self.yval) + linRange*5.0
            ey0=y0
            ey1=y1

        if self.logy:
            # y0 and y1 are in linear space - take the log to produce even spacing in log space
            ys = np.linspace(np.log10(y0), np.log10(y1), 84)
        else:
            ys = np.linspace(y0, y1, 84)
        arr = np.zeros((len(xs),len(ys)))
        for i,x in enumerate(xs):
            xThis = x
            if self.logx:
                xThis = 10.0**x
            epsilon, theta, sigma = self.returnCachedSkewParams([xThis])

            for j,yThis in enumerate(ys):
                likThis = epsSkewNormalPDF( yThis, epsilon,theta,sigma )
                arr[i,j] = likThis

        if axIn is None:
            fig,ax=plt.subplots()
        else:
            ax=axIn
        ax.imshow( arr.T, origin='lower', interpolation='Nearest', norm=matplotlib.colors.LogNorm(), extent=(ex0,ex1,ey0,ey1))
        #ax.imshow( arr.T, origin='lower', interpolation='Nearest', extent=(ex0,ex1,ey0,ey1))
        if self.logx and self.logy:
            ax.plot( np.log10(self.xval), np.log10(self.yval), c='k', lw=2, ls='-' )
            ax.plot( np.log10(self.xval), np.log10(self.yUpper), c='k', lw=1, ls='-' )
            ax.plot( np.log10(self.xval), np.log10(self.yLower), c='k', lw=1, ls='-' )
            ax.set_xlabel(r'$\log_{10}$'+self.xvar)
            ax.set_ylabel(r'$\log_{10}$'+self.yvar)
        if self.logx and not self.logy:
            ax.plot( np.log10(self.xval), self.yval, c='k', lw=2, ls='-' )
            ax.plot( np.log10(self.xval), self.yUpper, c='k', lw=1, ls='-' )
            ax.plot( np.log10(self.xval), self.yLower, c='k', lw=1, ls='-' )
            ax.set_xlabel(r'$\log_{10}$'+self.xvar)
            ax.set_ylabel(self.yvar)
        if not self.logx and self.logy:
            ax.plot( self.xval, np.log10(self.yval), c='k', lw=2, ls='-' )
            ax.plot( self.xval, np.log10(self.yUpper), c='k', lw=1, ls='-' )
            ax.plot( self.xval, np.log10(self.yLower), c='k', lw=1, ls='-' )
            ax.set_xlabel(self.xvar)
            ax.set_ylabel(r'$\log_{10}$'+self.yvar)
        if not self.logx and not self.logy:
            ax.plot( self.xval, self.yval, c='k', lw=2, ls='-' )
            ax.plot( self.xval, self.yUpper, c='k', lw=1, ls='-' )
            ax.plot( self.xval, self.yLower, c='k', lw=1, ls='-' )
            ax.set_xlabel(self.xvar)
            ax.set_ylabel(self.yvar)


        plt.savefig('likelihood_'+self.label+'.pdf')
        if axIn is None:
            plt.close(fig)
     
    def plot(self, z, axIn=None, color='k', lw=2, scatter=False, alpha=None, showLabel=True):
        ''' Plot the requested variables at redshift z. If z is outside the z
            range (specified by zmin, zmax in __init__), plot a dotted line instead'''

        if alpha is None:
            thisAlpha = self.alpha
        else:
            thisAlpha = alpha
        if axIn is None:
            fig,ax = plt.subplots()
        else:
            ax = axIn
        label=None
        #### here's the issue lol
#        if z<0.5 and self.zmin<z and z<self.zmax:
#            label=self.label
#        elif z>0.5 and z<1.5 and self.zmin<z and z<self.zmax:
#            label=self.label
        if z>self.zmin and z<self.zmax:
            if showLabel:
                label=self.label
            if not scatter:
                ax.plot(self.xval, self.yval, c=color, lw=lw, ls='-',label=label)
            else:
                ax.scatter(self.xval, self.yval, c=color, lw=lw, label=label, s=10)
            ax.fill_between(self.xval, self.yLower, self.yUpper, facecolor=color, alpha=thisAlpha)
        else:
            ax.plot(self.xval, self.yval, c='gray', lw=lw, ls=':',label=label)
        if axIn is None:
            plt.savefig(self.label+'_'+self.xvar+'_'+self.yvar+'.png')
            plt.close(fig)

    def interior(self, x,y, sigma=1.0):
        ''' Return true if the given point is within this dataset's uncertainties'''
        xv = copy.deepcopy( self.xval )
        yv = copy.deepcopy( self.yval )
        yu = copy.deepcopy( self.yUpper )
        yl = copy.deepcopy( self.yLower )
        xEval = x
        yEval = y
        if len(xv)!=len(yv) or len(xv)!=len(yu) or len(xv)!=len(yl):
            pdb.set_trace()
        if self.logx:
            xv = np.log10(xv)
            xEval = np.log10(xEval)
        if self.logy:
            yv = np.log10(yv)
            yu = np.log10(yu)
            yl = np.log10(yl)
            yEval = np.log10(yEval)
        try:
            f = interp1d(xv,yv,kind='linear',bounds_error=False)
        except:
            pdb.set_trace()
        fu = interp1d(xv,yv+sigma*(yu-yv),kind='linear',bounds_error=False)
        fl = interp1d(xv,yv-sigma*(yv-yl),kind='linear',bounds_error=False)
        #yf = f(xEval) # don't need this
        yfu = fu(xEval)
        yfl = fl(xEval)
        #pdb.set_trace()
        inx = np.logical_and( np.min(xv) < xEval, xEval<np.max(xv) )
        ret = np.logical_and( np.logical_and( yfl<yEval, yEval < yfu ), inx)
        return ret
    def cacheQuantiles(self, x0=None, x1=None):
        cachename = 'likelihood_arrays_'+self.label+'.pickle'
        # if we already did this, just load up the data.
        #if os.path.isfile( cachename ):
        #    print "READING IN CACHED LIKELIHOOD VALUES FROM ",cachename
        #    arrays = pickle.load(open(cachename,'rb'))
        #    self.cacheq16 = arrays[0]
        #    self.cacheq50 = arrays[1]
        #    self.cacheq84 = arrays[2]
        #    self.cacheEpsilons = arrays[3]
        #    self.cacheThetas = arrays[4]
        #    self.cacheSigmas = arrays[5]
        #    return

        if self.logx:
            dynRange = np.max(self.xval)/np.min(self.xval)
            if x0 is None:
                x0 = np.min(self.xval)/dynRange**5.0
            if x1 is None:
                x1 = np.max(self.xval)*dynRange**5.0
        else:
            linRange = np.max(self.xval)-np.min(self.xval)
            if x0 is None:
                x0 = np.min(self.xval) - linRange*5.0
            if x1 is None:
                x1 = np.max(self.xval) + linRange*5.0

        if self.logx:
            # x0 and x1 are in linear space - take the log to produce even spacing in log space
            xcache = np.linspace(np.log10(x0), np.log10(x1), 100)
            q16, q50, q84, logy = self.returnQuantiles(np.power(10.0, xcache))
        else:
            xcache = np.linspace(x0, x1, 100)
            q16, q50, q84, logy = self.returnQuantiles(xcache)
        self.cacheq16 = interp1d(xcache, q16, kind='linear', bounds_error=True)
        self.cacheq50 = interp1d(xcache, q50, kind='linear', bounds_error=True)
        self.cacheq84 = interp1d(xcache, q84, kind='linear', bounds_error=True)

        epsilons=[]
        thetas=[]
        sigmas=[]
        for i in range(len(xcache)):
            epsilon, theta, sigma = computeEpsSkewParams( q16[i], q50[i], q84[i] )
            epsilons.append(epsilon)
            thetas.append(theta)
            sigmas.append(sigma)
        #print "Storing cached versions of eps skew norm parameters for ",self.label
        self.cacheEpsilons = interp1d( xcache, epsilons, kind='linear', bounds_error=False)
        self.cacheThetas = interp1d( xcache, thetas, kind='linear', bounds_error=False)
        self.cacheSigmas = interp1d( xcache, sigmas, kind='linear', bounds_error=False)

        #print "Saving cached likelihood values to ",cachename
        #arrays = [np.array(self.cacheq16), np.array(self.cacheq50), np.array(self.cacheq84), np.array(self.cacheEpsilons), np.array(self.cacheThetas), np.array(self.cacheSigmas)]
        #pickle.dump( arrays, open(cachename, 'wb'), -1) 


    def returnCachedSkewParams(self,x):
        if self.logx:
            xThis = np.log10(x)
        else:
            xThis = x

        eps = self.cacheEpsilons(xThis)
        th = self.cacheThetas(xThis)
        sig = self.cacheSigmas(xThis)

        if np.any(np.isnan(eps)) or np.any(np.isnan(th)) or np.any(np.isnan(sig)):
            invalid = np.logical_or( np.logical_or( np.isnan(eps), np.isnan(th) ), np.isnan(sig))
            valid = np.logical_not(invalid)
            if np.any(valid):
                maxvalid = np.max(sig[valid] )
                minvalid = np.min(eps[valid] )
                medvalid = np.median(th[valid] )
            else:
                # reasonable default choices I hope?
                maxvalid = 100
                medvalid = 0
                minvalid = 0
                #print "WARNING: no valid choices in returnCachedSkewParams for ", self.label, "for x=",x
            eps[invalid] = minvalid
            th[invalid] = medvalid
            sig[invalid] = maxvalid

            #print "WARNING: replacing ", np.sum(np.ones(len(eps))[invalid]), " of ",len(eps), "quantiles with extreme values owing to requested x-values lying far outside range of extrapolation. ", self.label 
        return eps,th,sig


    def returnCachedQuantiles(self,x):
        try:
            if self.logx:
                xThis = np.log10(x)
            else:
                xThis = x
            q16 = self.cacheq16(xThis)
            q50 = self.cacheq50(xThis)
            q84 = self.cacheq84(xThis)
            return q16, q50, q84, self.logx
        except:
            return self.returnQuantiles(x) 
    def returnQuantiles(self, x):
        ''' Designed to sort-of replace the workflow specified below with distance. 
            Return the 16th, 50th and 84th percentiles at the given x.'''
        sigma = 1.0
        xv = copy.deepcopy( self.xval )
        yv = copy.deepcopy( self.yval )
        yu = copy.deepcopy( self.yUpper )
        yl = copy.deepcopy( self.yLower )
        xEval = x
        if len(xv)!=len(yv) or len(xv)!=len(yu) or len(xv)!=len(yl):
            pdb.set_trace()
        if self.logx:
            xv = np.log10(xv)
            xEval = np.log10(xEval)
        if self.logy:
            yv = np.log10(yv)
            yu = np.log10(yu)
            yl = np.log10(yl)
        inx = np.logical_and( np.min(xv) < xEval, xEval<np.max(xv) )
        outx = np.logical_not(inx)
        belowx = xEval < np.min(xv)
        abovex = xEval > np.max(xv)
        inflate = np.ones(len(xEval))
        inflate[abovex] = np.exp(4.0*( xEval[abovex] - np.max(xv))**2.0 )
        inflate[belowx] = np.exp(4.0*( xEval[belowx] - np.min(xv))**2.0 )
        try:
            f = interp1d(xv,yv,kind='linear',bounds_error=False, fill_value='extrapolate')
        except:
            pdb.set_trace()
        fu = interp1d(xv,yv+sigma*(yu-yv),kind='linear',bounds_error=False)
        fl = interp1d(xv,yv-sigma*(yv-yl),kind='linear',bounds_error=False)
        yf = f(xEval) 
        yfu = fu(xEval)
        yfl = fl(xEval)

	#minvalid = np.min( np.arange(len(inx))[inx] )
	#maxvalid = np.max( np.arange(len(inx))[inx] )
        minvalid = np.argmin(xv)
        maxvalid = np.argmax(xv)
        if self.fixedSigma<=0:
            #deltaYbelowleft = yf[minvalid]-yfl[minvalid]
            #deltaYbelowright = yf[maxvalid] - yfl[maxvalid]
            #deltaYaboveleft = yfu[minvalid]-yf[minvalid]
            #deltaYaboveright = yfu[maxvalid]-yf[maxvalid]

            deltaYbelowleft = yv[minvalid]-yl[minvalid]
            deltaYbelowright = yv[maxvalid] - yl[maxvalid]
            deltaYaboveleft = yu[minvalid]-yv[minvalid]
            deltaYaboveright = yu[maxvalid]-yv[maxvalid]

            yfl[belowx] = yf[belowx] - deltaYbelowleft*inflate[belowx]
            yfl[abovex] = yf[abovex] - deltaYbelowright*inflate[abovex]

            yfu[belowx] = yf[belowx] + deltaYaboveleft*inflate[belowx]
            yfu[abovex] = yf[abovex] + deltaYaboveright*inflate[abovex]
            if( len(yfl)!=len(yf) or len(yfu)!=len(yf)):
                pdb.set_trace()
            ret = yfl, yf, yfu, self.logy
            #print "return quantiles: ",ret
 
        else:
            ret = yf-self.fixedSigma*inflate, yf, yf+self.fixedSigma*inflate, self.logy
            #print "return quantiles 2: ",ret
            if( len(inflate)!=len(yf) or len(yfu)!=len(yf)):
                pdb.set_trace()
        if np.any(np.isnan(ret[0])) or np.any(np.isnan(ret[1])) or np.any(np.isnan(ret[2])):
            invalid = np.logical_or( np.logical_or( np.isnan(ret[0]), np.isnan(ret[1]) ), np.isnan(ret[2]))
            valid = np.logical_not(invalid)
            maxvalid = np.max(ret[2][valid] )
            minvalid = np.min(ret[0][valid] )
            medvalid = np.median(ret[1][valid] )
            ret[0][invalid] = minvalid
            ret[1][invalid] = medvalid
            ret[2][invalid] = maxvalid
            #print "WARNING: replacing ", np.sum(np.ones(len(ret[0]))[invalid]), " of ",len(ret[0]), "quantiles with extreme values owing to requested x-values lying far outside range of extrapolation. ", self.label 
            #pdb.set_trace()

        
        return ret


    def distance(self, x,y):
        ''' Return the distance in sigma (above or below) the relation'''
        sigma = 1.0
        xv = copy.deepcopy( self.xval )
        yv = copy.deepcopy( self.yval )
        yu = copy.deepcopy( self.yUpper )
        yl = copy.deepcopy( self.yLower )
        xEval = x
        yEval = y
        if len(xv)!=len(yv) or len(xv)!=len(yu) or len(xv)!=len(yl):
            pdb.set_trace()
        if self.logx:
            xv = np.log10(xv)
            xEval = np.log10(xEval)
        if self.logy:
            yv = np.log10(yv)
            yu = np.log10(yu)
            yl = np.log10(yl)
            yEval = np.log10(yEval)
        try:
            f = interp1d(xv,yv,kind='linear',bounds_error=False)
        except:
            pdb.set_trace()
        fu = interp1d(xv,yv+sigma*(yu-yv),kind='linear',bounds_error=False)
        fl = interp1d(xv,yv-sigma*(yv-yl),kind='linear',bounds_error=False)
        yf = f(xEval) # don't need this
        yfu = fu(xEval)
        yfl = fl(xEval)
        if self.fixedSigma<=0:
            sigmaAbove = yfu-yf
            sigmaBelow = yf-yfl
            above = (yEval - yf)/(yfu-yf) # this will be positive for numbers above the median, and in units of the sigma above the median
            below = (yEval - yf)/(yf - yfl) # this will be negative for numbers below the median, and in units of the sigma below median
        else:
            sigmaAbove = yf*0 + self.fixedSigma
            sigmaBelow = yf*0 + self.fixedSigma
            above = (yEval - yf)/self.fixedSigma
            below = (yEval - yf)/self.fixedSigma

        isAbove = above>0
        above = np.clip(above, 0, np.inf)
        below = np.clip(below, -np.inf, 0)
        ret = above+below 
        inx = np.logical_and( np.min(xv) < xEval, xEval<np.max(xv) )
        ret[ np.logical_not(inx) ] = 0 # if we're outside the given range, don't penalize us! 
        retsigma = sigmaBelow
        retsigma[isAbove] = sigmaAbove[isAbove]
        retsigma[ np.logical_not(inx) ] = np.inf # if we're outside the given range, don't penalize us! 
        return ret, retsigma
    def test(self, nx=1000, ny=1000):
        xv = self.xval 
        yv = self.yval 
        if self.logx:
            xr = np.power(10.0, np.random.uniform( np.log10(np.min(xv)), np.log10(np.max(xv)), nx ))
        else:
            xr = np.random.uniform( np.min(xv), np.max(xv), nx )
        if self.logy:
            yr = np.power(10.0, np.random.uniform(0,1, ny)*(np.log10(np.max(self.yUpper)) -  np.log10(np.min(self.yLower))) + np.log10(np.min(self.yLower)))
        else:
            yr = np.random.uniform( np.min(self.yLower), np.max(self.yUpper), ny )
        fig,ax = plt.subplots()
        inter = self.interior(xr,yr)
        ax.scatter(xr,yr, c='k')
        ax.scatter(xr[inter], yr[inter], c='r', lw=0)
        if self.logx:
            ax.set_xscale('log')
        if self.logy:
            ax.set_yscale('log')
        ax.set_xlim( np.min(xr), np.max(xr) )
        ax.set_ylim( np.min(yr), np.max(yr) )
        plt.savefig(self.label+'_'+self.xvar+'_'+self.yvar+'_test.png')
        plt.close(fig)


def defineKravtsov13():
    Rvirs = np.power(10.0, np.linspace(-1, np.log10(25.0), 100))/0.015
    #rhalfs = 0.015* np.power(Rvirs, 0.95)
    #scatter = 10.0**0.2 # .2 dex
    rhalfs = 0.015*Rvirs
    scatter = 10.0**0.25
    datasets['Kravtsov13rvir'] = DataSet('Rvir', 'halfMassStars', Rvirs, rhalfs, yLower=rhalfs/scatter, yUpper=rhalfs*scatter, zmin=-0.5, zmax=0.5, label='Kravtsov13')


def defineMoster(z0):
    def Moster(Mh,z, mparams):
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
    logMhs = np.linspace(8.5, 14, num=100)
    Mhs = np.power(10.0, logMhs)
    effC = Moster(Mhs, z0, central)
    eff = Moster(Mhs, z0, central)
    for i in range(len(unc)):
        theseParams = copy.copy(central)
        theseParams[i] = theseParams[i]+unc[i]
        eff = np.vstack([eff, Moster(Mhs, z0, theseParams)])
        theseParams = copy.copy(central)
        theseParams[i] = theseParams[i]-unc[i]
        eff = np.vstack([eff, Moster(Mhs, z0, theseParams)])
    effM = np.min(eff, axis=0)
    effP = np.max(eff, axis=0)

    mosterDS = DataSet('Mh', 'efficiency', Mhs, effC, yLower=effC/10.0**0.15, yUpper=effC*10.0**0.15, zmin=z0-0.5, zmax=z0+0.5, label='Moster13')
    mosterDS2 = DataSet('Mh', 'mstar', Mhs, effC*Mhs, yLower=effC*Mhs/10.0**0.15, yUpper=effC*Mhs*10.0**0.15, zmin=z0-0.5, zmax=z0+0.5, label='Moster13', alpha=.5)
    datasets['Moster13effz'+str(z0)] = mosterDS
    datasets['Moster13z'+str(z0)] = mosterDS2

def defineGarrisonKimmel():
    def behrooziFn( epsilon, M1, alpha, delta, gamma, Mh ):
        def f(x):
            return - np.log10( np.power(10.0, -alpha*x) + 1.0) + delta*np.power(np.log10(1.0+np.exp(x)),gamma)/(1.0+np.exp(np.power(10.0, -x)))
        return np.log10(epsilon*M1) + f(np.log10(Mh/M1)) - f(0.0)
    
    Mh = np.power(10.0, np.linspace(10.0, 15.0, 200) )
    Mst = np.power(10.0, behrooziFn( 10.0**(-1.777), 10.0**11.514, 1.412, 3.508, 0.316, Mh))
    datasets['Behroozi13'] = DataSet( 'Mh', 'mstar', Mh, Mst, yLower=Mst/10.0**0.2, yUpper=Mst*10.0**0.2, zmin=-0.5, zmax=0.5, label='Behroozi13' )

    Mh = np.power(10.0, np.linspace(8.0, 11.5, 200) )
    for i, sigma in enumerate(np.linspace( 0.2, 1.0, 2 )):
        alpha = 0.24*sigma*sigma + 0.16*sigma + 1.99
        Mst = np.power(10.0, behrooziFn( 10.0**(-1.777), 10.0**11.514, alpha, 3.508, 0.316, Mh))
        datasets['GarrisonKimmel16_'+str(i)] = DataSet( 'Mh', 'mstar', Mh, Mst, yLower=Mst/np.power(10.0,sigma), yUpper=Mst*np.power(10.0,sigma), zmin=-0.5, zmax=0.5, label='Garrison-Kimmel16 '+r'$\sigma='+str(sigma)+r'$' )
#    for i, nu in enumerate(np.linspace( 0, 1.0, 3 )):
#        #alpha = 0.24*sigma*sigma + 0.16*sigma + 1.99
#        alpha = 0.47*nu*nu - 1.48*nu + 1.81
#        sigma = 0.2 + nu*(np.log10(Mh) - 11.514)
#        Mst = np.power(10.0, behrooziFn( 10.0**(-1.777), 10.0**11.514, alpha, 3.508, 0.316, Mh))
#        datasets['GarrisonKimmel16_'+str(i)] = DataSet( 'Mh', 'mstar', Mh, Mst, yLower=Mst/np.power(10.0,sigma), yUpper=Mst*np.power(10.0, sigma), zmin=-0.5, zmax=0.5, label='Garrison-Kimmel16 '+r'$\nu ='+str(nu)+r'$' )
   


def defineBroeils():
    mhis = np.power(10.0, np.linspace(8.0, 10.5, 100))
    DHI = np.power(10.0, (np.log10(mhis) - 6.52)/1.96)/2.0
    datasets['broeils97'] =  DataSet('MHI', 'broeilsHI', mhis, DHI, yLower=DHI/10**0.13, yUpper=DHI*10**0.13, label='Broeils97', alpha=0.5) 

def defineMetalRelations(z):
    thisMst = np.power(10.0, np.linspace(10.0, 11.5))
    ZHayward = -8.69 + 9.09*np.power(1.0+z,-0.017) - 0.0864*np.power(np.log10(thisMst) - 11.07*np.power(1.0+z,0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
            
    b = 10.4 + 4.46*np.log10(1.0+z)-1.78*np.log10(1.0+z)**2.0
    ZGenzel15 = np.power(10.0,  8.74 - 0.087 * np.power(np.log10(thisMst) -b,2.0) -8.69) # Equation 12a normalized to solar
    #bb = np.log10(mst) - 0.32*np.log10(sfr) - 10
    #ZMannucci10 = np.power(10.0, 0.21 + 0.39*bb - 0.2*bb*bb - 0.077*bb*bb*bb + 0.064*bb*bb*bb*bb
    datasets['Hayward16Z'+str(z)] = DataSet( 'mstar', 'sfZ', thisMst, ZHayward/0.02, yLower=ZHayward/0.02/10.0**0.15, yUpper=ZHayward/0.02*10.0**0.15, label='Hayward16', zmin=z-0.5, zmax=z+0.5 )
    if z>0.5:
        alphaGenzel = 0.5
    else:
        alphaGenzel = 0.2
    datasets['Genzel15Z'+str(z)] = DataSet('mstar', 'sfZ', thisMst, ZGenzel15, yLower=ZGenzel15/10.0**0.15, yUpper=ZGenzel15*10.0**0.15, label='Genzel15', zmin=z-0.5, zmax=z+0.5, alpha= alphaGenzel )
    mLee = np.power(10.0, np.linspace(5.89,9.29, 100))
    ZLee = np.power(10.0, 5.65 + 0.298*np.log10(mLee) - 8.7) # Z in Zsun
    datasets['Lee06'] = DataSet( 'mstar', 'sfZ', mLee, ZLee, yLower=ZLee/10.0**0.117, yUpper=ZLee*10.0**0.117, label='Lee06', zmax=0.5, alpha=0.5 )

    mTremonti = np.power(10.0, np.array([8.57, 8.67, 8.76, 8.86, 8.96, 9.06, 9.16, 9.26, 9.36, 9.46, 9.57, 9.66, 9.76, 9.86, 9.96, 10.06, 10.16, 10.26, 10.36, 10.46, 10.56, 10.66, 10.76, 10.86, 10.95, 11.05, 11.15, 11.25]))
    Z16Tremonti = np.power(10.0, np.array([8.25, 8.28, 8.32, 8.37, 8.46, 8.56, 8.59, 8.60,8.63, 8.66, 8.69, 8.72, 8.76, 8.80, 8.83, 8.85, 8.88, 8.92, 8.94, 8.96, 8.98, 9.00, 9.01, 9.02, 9.03, 9.03, 9.04, 9.03 ]) - 8.7 )
    Z50Tremonti = np.power(10.0, np.array([8.44, 8.48, 8.57, 8.61, 8.63, 8.66, 8.68, 8.71, 8.74, 8.78, 8.82, 8.84, 8.87, 8.90, 8.94, 8.97, 8.99, 9.01, 9.03, 9.05, 9.07, 9.08, 9.09, 9.10, 9.11, 9.11, 9.12, 9.12 ]) - 8.7)
    Z84Tremonti = np.power(10.0, np.array([8.64, 8.65, 8.70, 8.73, 8.75, 8.82, 8.82, 8.86, 8.88, 8.92, 8.94, 8.96, 8.99, 9.01, 9.05, 9.06, 9.09, 9.10, 9.11, 9.12, 9.14, 9.15, 9.15, 9.16, 9.17, 9.17, 9.18, 9.18 ]) -8.7 )
    datasets['Tremonti04'] = DataSet( 'mstar', 'sfZ', mTremonti, Z50Tremonti, yLower=Z16Tremonti, yUpper=Z84Tremonti, label='Tremonti04', zmin=-0.5, zmax=0.5, alpha=0.5 )

    def kewleyfit(a,b,c,d,rms, label, alpha=0.5):
        # 8.69 is from Asplund09
        # solar abundance is 0.0134 according to Asplund09, but throughout the rest of the codebase 0.02 is assumed.
        kmass = np.power(10, np.linspace(8.5,11,100))
        x = np.log10(kmass)
        y =np.power(10.0,  a+b*x+c*x*x+d*x*x*x - 8.69 + np.log10(0.02/0.0134))
        datasets[label] = DataSet('mstar', 'sfZ', kmass, y, yLower=y/10.0**rms, yUpper=y*10.0**rms, label=label, zmin=-0.5, zmax=0.5, alpha=alpha)
    #kewleyfit( -0.694114, 1.30207, 0.00271531, -0.00364112, 0.12, 'T04', alpha=0.1 )
    #kewleyfit( 72.0142, -20.6826, 2.22124, -0.0783089, 0.13, 'Z94', alpha=0.1 )
    #kewleyfit( 27.7911, -6.94493, 0.808097, -0.0301508, 0.1, 'KK04', alpha=0.1 )
    #kewleyfit( 28.0974, -7.23631, 0.850344, -0.0318315, 0.1, 'KD02', alpha=0.1 )
    #kewleyfit( 45.5323, -12.2469, 1.32882, -0.0471074, 0.11, 'M91', alpha=0.1 )
    #kewleyfit( -8.71120, 4.15003, -0.322156, 0.00818179, 0.08, 'D02', alpha=0.1)
    kewleyfit( 32.1488, -8.51258, 0.976384, -0.0359763, 0.10, 'PP04O3', alpha=0.1)
    #kewleyfit( 23.9049, -5.62784, 0.645142, -0.0235065, 0.09, 'PP04N2', alpha=0.1)
    #kewleyfit( 91.6457, -25.9355, 2.67172, -0.0909689, 0.12, 'P01', alpha=0.1)
    #kewleyfit( 41.9463, -10.3253, 1.04371, -0.0347747, 0.13, 'P05', alpha=0.1)

def defineStructureRelations():
    mstThis = np.power(10.0, np.linspace(9.75, 11.25, 100))
    fang = np.power(10.0, 9.29 + 0.64*(np.log10(mstThis)-10.25))  # Fang et al 2013
    datasets['Fang13'] = DataSet( 'mstar', 'Sigma1', mstThis, fang, yLower=fang/10.0**0.16, yUpper=fang*10.0**0.16, label='Fang13 Red and Green', zmin=-0.5, zmax=3.0, alpha=0.5)


    mstThis = np.power(10.0, np.linspace(8.75,11.5, 30))
    gamma, M0, logn1, logn2 = (1.95, 10.0**10.16, 0.12, 0.62) # dutton09 all galaxies
    lognAll = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    gamma, M0, logn1, logn2 = (1.70, 10.0**9.90, 0.21, 0.62) # dutton09 red galaxies
    lognRed = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    gamma, M0, logn1, logn2 = (1.70, 10.0**10.41, 0.12, 0.61) # dutton09 blue galaxies
    lognBlue = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognAll]
    datasets['Dutton09All'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 All', yLower=np.array(c82s)/1.2, yUpper=np.array(c82s)*1.2, alpha=0.5)
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognRed]
    datasets['Dutton09Red'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 Red')
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognBlue]
    datasets['Dutton09Blue'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 Blue')


    mstBarro = np.power(10.0, np.linspace(10.0, 11.0, 10))
    er = 10.0**0.14
    barrohq = np.power(10.0, 0.65*(np.log10(mstBarro)-10.5) + 9.53)
    datasets['Barro15HQ'] = DataSet( 'mstar', 'Sigma1', mstBarro, barrohq, yLower=barrohq/er, yUpper=barrohq*er, label='Barro15 Red', zmin=0.5, zmax=1.0, alpha=0.5)
    barro1q = np.power(10.0, 0.65*(np.log10(mstBarro)-10.5) + 9.64)
    datasets['Barro151Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro1q, yLower=barro1q/er, yUpper=barro1q*er, label='Barro15 Red', zmin=1.0, zmax=1.5, alpha=0.5)
    barro2q = np.power(10.0, 0.64*(np.log10(mstBarro)-10.5) + 9.76)
    datasets['Barro152Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro2q, yLower=barro2q/er, yUpper=barro2q*er, label='Barro15 Red', zmin=1.5, zmax=2.2, alpha=0.5)
    barro3q = np.power(10.0, 0.67*(np.log10(mstBarro)-10.5) + 9.80)
    datasets['Barro153Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro3q, yLower=barro3q/er, yUpper=barro3q*er, label='Barro15 Red', zmin=2.2, zmax=3.05, alpha=0.5)
    mstBarro = np.power(10.0, np.linspace(9.0, 10.5, 10))
    er = 10.0**0.25
    barrohs = np.power(10.0, 0.89*(np.log10(mstBarro)-10.5) + 9.12)
    datasets['Barro15HS'] = DataSet( 'mstar', 'Sigma1', mstBarro, barrohs, yLower=barrohs/er, yUpper=barrohs*er, label='Barro15 Blue', zmin=0.5, zmax=1.0, alpha=0.5)
    barro1s = np.power(10.0, 0.88*(np.log10(mstBarro)-10.5) + 9.16)
    datasets['Barro151S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro1s, yLower=barro1s/er, yUpper=barro1s*er, label='Barro15 Blue', zmin=1.0, zmax=1.5, alpha=0.5)
    barro2s = np.power(10.0, 0.86*(np.log10(mstBarro)-10.5) + 9.25)
    datasets['Barro152S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro2s, yLower=barro2s/er, yUpper=barro2s*er, label='Barro15 Blue', zmin=1.5, zmax=2.2, alpha=0.5)
    barro3s = np.power(10.0, 0.89*(np.log10(mstBarro)-10.5) + 9.33)
    datasets['Barro153S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro3s, yLower=barro3s/er, yUpper=barro3s*er, label='Barro15 Blue', zmin=2.2, zmax=3.05, alpha=0.5)


    obsmst = np.power(10.0, np.linspace(9.0, 12.0, 100))
    gamma = 0.1
    alpha = 0.14
    beta = 0.39
    m0 = 3.98e10
    sig1, sig2 = 0.47, 0.34
    rLTG = gamma*np.power(obsmst, alpha) * np.power(1.0 + obsmst/m0, beta-alpha)
    sigr = sig2 + (sig1-sig2)/(1.0 + (obsmst/m0)**2)
    datasets['shenLTG'] = DataSet( 'mstar', 'halfMassStars', obsmst, rLTG, yLower=rLTG/np.exp(sigr), yUpper=rLTG*np.exp(sigr), label='Shen (2003) LTG', zmax=0.5)

    obsmst = np.power(10.0, np.linspace(10.2, 12.0, 100))
    a = 0.56
    b = 2.88e-6 ### from the Erratum!
    sig1, sig2 = 0.47, 0.34
    rETG = b * np.power(obsmst, a)
    sigr = sig2 + (sig1-sig2)/(1.0 + (obsmst/m0)**2)
    #datasets['shenETG'] = DataSet( 'mstar', 'halfMassStars', obsmst, rETG, yLower=rETG/np.exp(sigr), yUpper=rETG*np.exp(sigr), label='Shen (2003) ETG', zmax=0.5)


    LTGM = [ 6.402635431918009, 6.834553440702782, 7.295754026354318, 7.808199121522694, 8.27672035139092, 8.759882869692532, 9.22108345534407, 10.231332357247437, 10.655929721815518, 11.0805270863836 ]
    LTGRmed = [2.5119091496217947, 2.5702021038821394, 2.926375959883546, 3.119661255363682, 3.202152911858091, 3.3781577674909293, 3.505694210097417, 3.687897098455056, 3.81894846469353, 4.005427082664106]
    LTG16 = [ 2.3906316042185844, 2.4489448466384207, 2.6803566658438687, 2.9048907989815342, 2.938853177971116, 3.0905883228116684, 3.2354660697439295, 3.4731164979931624, 3.669977581583761, 3.7075157486838055]
    LTG84 = [ 2.653921194025813, 2.746876468778213, 3.1515187378059704, 3.337895914979086, 3.45154511241331, 3.592938367952823, 3.689296981460003, 3.8853668268304142, 4.0336986329161855, 4.1682547787069]
           
    ETGM = [7.3469985358711565, 7.888726207906295, 8.349926793557833, 8.818448023426061, 9.279648609077597, 9.77745241581259, 10.267935578330892, 10.685212298682284, 11.087847730600291 ]
    ETGmed = [ 2.8015936349280954, 2.863198631225506, 3.066947544963633, 3.256829501688989, 3.2561904246649913, 3.3351773016071604, 3.4696015743611763, 3.6353051170119595, 3.849527793087824]
    ETG16 = [ 2.6457044894315596, 2.5618129499322033, 2.786387659388853, 2.983157446261737, 3.003303588637278, 3.099621625825474, 3.265213583599052, 3.4828497425094422, 3.693648791671034]
    ETG84 = [ 3.0718014871220904, 2.946339508823659, 3.3441243799431253, 3.454289085984601, 3.470971025126885, 3.5465038429155435, 3.6601327521902753, 3.857024268020112, 4.015819692363874]
            
    datasets['baldry12LTG0'] = DataSet( 'mstar', 'halfMassStars', np.power(10.0, LTGM), np.power(10.0, LTGRmed)/1000, yLower = np.power(10.0, LTG16)/1000, yUpper = np.power(10.0, LTG84)/1000, label='Baldry 12 LTG', zmax=0.3)
    #datasets['baldry12ETG0'] = DataSet( 'mstar', 'halfMassStars', np.power(10.0, ETGM), np.power(10.0, ETGmed)/1000, yLower = np.power(10.0, ETG16)/1000, yUpper = np.power(10.0, ETG84)/1000, label='Baldry 12 ETG', zmax=0.3)



    obsmst = np.power(10.0, np.linspace(9.0, 11.5, 100))
    obsmstETG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([0.03, 0.04, 0.13, 0.42, 0.65]))
    reffETGvdW50 = np.power(10.0, np.array([0.27, 0.27, 0.38, 0.67, 0.76]))
    reffETGvdW84 = np.power(10.0, np.array([0.46, 0.46, 0.58, 0.92, 1.08]))
    datasets['vdW14ETG0']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmax=0.5)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75 ]))
    reffLTGvdW16 = np.power(10.0, np.array([0.24, 0.36, 0.42, 0.61]))
    reffLTGvdW50 = np.power(10.0, np.array([0.49, 0.61, 0.66, 0.83]))
    reffLTGvdW84 = np.power(10.0, np.array([0.70, 0.80, 0.85, 1.01]))
    datasets['vdW14LTG0']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmax=0.5, alpha=0.5)
                
    obsmstETG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.02, -0.14, 0.02, 0.26, 0.62]))
    reffETGvdW50 = np.power(10.0, np.array([0.23, 0.21, 0.23, 0.45, 0.81]))
    reffETGvdW84 = np.power(10.0, np.array([0.43, 0.44, 0.42, 0.64, 0.97]))
    datasets['vdW14ETG1']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=0.5, zmax=1.0)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.18, 0.32, 0.39, 0.51, 0.77]))
    reffLTGvdW50 = np.power(10.0, np.array([0.43, 0.56, 0.64, 0.75, 0.90]))
    reffLTGvdW84 = np.power(10.0, np.array([0.65, 0.76, 0.83, 0.90, 1.12]))
    datasets['vdW14LTG1']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=0.5, zmax=1.0, alpha=0.5)

    obsmstETG = np.power(10.0, np.array([9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.15, -0.15, 0.07, 0.41]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.18, 0.09, 0.30, 0.58]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.42, 0.36, 0.54, 0.81]))
    datasets['vdW14ETG1h']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=1.0, zmax=1.5)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.11, 0.23, 0.33, 0.47, 0.62]))
    reffLTGvdW50 = np.power(10.0, np.array([0.37, 0.48, 0.57, 0.67, 0.82]))
    reffLTGvdW84 = np.power(10.0, np.array([0.60, 0.69, 0.77, 0.83, 0.96]))
    datasets['vdW14LTG1h']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=1.0, zmax=1.5, alpha=0.5)

    obsmstETG = np.power(10.0, np.array([9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.02, -0.27, -0.04, 0.28]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.22, 0.02, 0.19, 0.45]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.48, 0.35, 0.50, 0.74]))
    datasets['vdW14ETG2']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=1.5, zmax=2.0)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.07, 0.16, 0.28, 0.35, 0.53]))
    reffLTGvdW50 = np.power(10.0, np.array([0.33, 0.42, 0.52, 0.61, 0.70]))
    reffLTGvdW84 = np.power(10.0, np.array([0.57, 0.65, 0.72, 0.80, 0.87]))
    datasets['vdW14LTG2']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=1.5, zmax=2.0, alpha=0.5)


    obsmstETG = np.power(10.0, np.array([ 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.37, -0.20, 0.16]))
    reffETGvdW50 = np.power(10.0, np.array([-0.04, 0.08, 0.36]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.36, 0.54, 0.55]))
    datasets['vdW14ETG2h']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=2.0, zmax=2.5)

    obsmstLTG = np.power(10.0, np.array([ 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([ 0.10, 0.17, 0.26, 0.40]))
    reffLTGvdW50 = np.power(10.0, np.array([ 0.35, 0.44, 0.53, 0.64]))
    reffLTGvdW84 = np.power(10.0, np.array([ 0.57, 0.64, 0.70, 0.84]))
    datasets['vdW14LTG2h']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=2.0, zmax=2.5, alpha=0.5)

    obsmstETG = np.power(10.0, np.array([  10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.22, 0.07]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.10, 0.39]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.50, 0.68]))
    datasets['vdW14ETG3']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=2.5, zmax=3.0)

    obsmstLTG = np.power(10.0, np.array([ 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([  0.16, 0.19, 0.33]))
    reffLTGvdW50 = np.power(10.0, np.array([  0.43, 0.47, 0.55]))
    reffLTGvdW84 = np.power(10.0, np.array([  0.65, 0.71, 0.76]))
    datasets['vdW14LTG3']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=2.5, zmax=3.0, alpha=0.5)




def defineAngularMomenta():
    datasets['Fall13Disks'] = DataSet('mstar', 'specificJStars', [1.0e9,1.0e11], [10**2.3, 10**3.45], label='Fall13 Disks')
    datasets['Fall13Ellipticals'] = DataSet('mstar', 'specificJStars', [1.0e10, 5.0e11], [10**2.1, 10**2.1*50**0.6], label='Fall13 Ellipticals')


    datasets['Fall13Disks1'] = DataSet('mstar', 'specificJHI', [1.0e9,1.0e11], [10**2.3, 10**3.45], label='Fall13 Disks', zmin=0.5, zmax=0.6)
    datasets['Fall13Ellipticals1'] = DataSet('mstar', 'specificJHI', [1.0e10, 5.0e11], [10**2.1, 10**2.1*50**0.6], label='Fall13 Ellipticals', zmin=0.5, zmax=0.6)

    datasets['Fall13Disks2'] = DataSet('mstar', 'specificJH2', [1.0e9,1.0e11], [10**2.3, 10**3.45], label='Fall13 Disks', zmin=0.5, zmax=0.6)
    datasets['Fall13Ellipticals2'] = DataSet('mstar', 'specificJH2', [1.0e10, 5.0e11], [10**2.1, 10**2.1*50**0.6], label='Fall13 Ellipticals', zmin=0.5, zmax=0.6)

    datasets['Fall13Disks3'] = DataSet('mstar', 'specificJAccr', [1.0e9,1.0e11], [10**2.3, 10**3.45], label='Fall13 Disks', zmin=0.5, zmax=0.6)
    datasets['Fall13Ellipticals3'] = DataSet('mstar', 'specificJAccr', [1.0e10, 5.0e11], [10**2.1, 10**2.1*50**0.6], label='Fall13 Ellipticals', zmin=0.5, zmax=0.6)


    datasets['Burkert16'] = DataSet('mstar', 'specificJStars', [10**9.8, 10**11.4],[10.0**(3.33+2.0/3.0*(9.8-11.0)), 10.0**(3.33+2.0/3.0*(11.4-11.0))],  yLower=[10.0**(3.33+2.0/3.0*(9.8-11.0)-0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)-0.17)], yUpper=[10.0**(3.33+2.0/3.0*(9.8-11.0)+0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)+0.17)], zmin=0.8, zmax=2.6, label='Burkert16' )
    datasets['Burkert16a'] = DataSet('mstar', 'specificJHI', [10**9.8, 10**11.4],[10.0**(3.33+2.0/3.0*(9.8-11.0)), 10.0**(3.33+2.0/3.0*(11.4-11.0))],  yLower=[10.0**(3.33+2.0/3.0*(9.8-11.0)-0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)-0.17)], yUpper=[10.0**(3.33+2.0/3.0*(9.8-11.0)+0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)+0.17)], zmin=0.5, zmax=.6, label='Burkert16' )
    datasets['Burkert16b'] = DataSet('mstar', 'specificJH2', [10**9.8, 10**11.4],[10.0**(3.33+2.0/3.0*(9.8-11.0)), 10.0**(3.33+2.0/3.0*(11.4-11.0))],  yLower=[10.0**(3.33+2.0/3.0*(9.8-11.0)-0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)-0.17)], yUpper=[10.0**(3.33+2.0/3.0*(9.8-11.0)+0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)+0.17)], zmin=0.5, zmax=.6, label='Burkert16' )
    datasets['Burkert16c'] = DataSet('mstar', 'specificJAccr', [10**9.8, 10**11.4],[10.0**(3.33+2.0/3.0*(9.8-11.0)), 10.0**(3.33+2.0/3.0*(11.4-11.0))],  yLower=[10.0**(3.33+2.0/3.0*(9.8-11.0)-0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)-0.17)], yUpper=[10.0**(3.33+2.0/3.0*(9.8-11.0)+0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)+0.17)], zmin=0.5, zmax=.6, label='Burkert16' )


def defineBrinchmann(specific=True):
    z=0
    brinchmannMode = np.loadtxt(os.environ['GIDGETDIR']+'/py/brinchmann04_msmode.csv', delimiter=',')
    brinchmann975 = np.loadtxt(os.environ['GIDGETDIR']+'/py/brinchmann04_ms975.csv', delimiter=',')
    fbmode = interp1d( brinchmannMode[:,0], brinchmannMode[:,1], kind='linear' )
    fb975 = interp1d( brinchmann975[:,0], brinchmann975[:,1], kind='linear' )
    brinchmannm = np.linspace( 6.51, 11.88, 1000 ) 
    delt = fb975(brinchmannm) - fbmode(brinchmannm)
    fac = 1.0
    if specific:
        fac = np.power(10.0, brinchmannm)*1.0e-9
    if specific:
        datasets['Brinchmann04Specz'+str(z)] = DataSet('mstar', 'sSFR', np.power(10.0, brinchmannm), np.power(10.0, fbmode(brinchmannm))/fac, yLower = np.power(10.0, (fbmode(brinchmannm)-delt/2.0))/fac, yUpper=np.power(10.0, fbmode(brinchmannm)+delt/2.0)/fac, zmax=0.5, label='Brinchmann 04' )
    else:
        datasets['Brinchmann04z'+str(z)] = DataSet('mstar', 'sfr', np.power(10.0, brinchmannm), np.power(10.0, fbmode(brinchmannm))/fac, yLower = np.power(10.0, (fbmode(brinchmannm)-delt/2.0))/fac, yUpper=np.power(10.0, fbmode(brinchmannm)+delt/2.0)/fac, zmax=0.5, label='Brinchmann 04' )



def defineMS(z, specific=True):


    if z<0.5:
        ### from Cano-Di'az+ (2016) "Spatially Resolved STar Formation Main Sequence..."
        ## Note that the code uses Msun/pc^2 for colst, but the data is in Msun/kpc^2
        arr = np.loadtxt(os.environ['GIDGETDIR']+'/py/CanoDiaz16_20p.csv', delimiter=',')
        datasets['CanoDiaz16_20p'] = DataSet('colst', 'colsfr', np.power(10.0,arr[:,0]-6), np.power(10.0,arr[:,1]), zmax= 0.5, label='CanoDiaz16_20p', extrapolate=False)
        arr = np.loadtxt(os.environ['GIDGETDIR']+'/py/CanoDiaz16_40p.csv', delimiter=',')
        datasets['CanoDiaz16_40p'] = DataSet('colst', 'colsfr', np.power(10.0,arr[:,0]-6), np.power(10.0,arr[:,1]), zmax= 0.5, label='CanoDiaz16_40p', extrapolate=False)
        arr = np.loadtxt(os.environ['GIDGETDIR']+'/py/CanoDiaz16_60p.csv', delimiter=',')
        datasets['CanoDiaz16_60p'] = DataSet('colst', 'colsfr', np.power(10.0,arr[:,0]-6), np.power(10.0,arr[:,1]), zmax= 0.5, label='CanoDiaz16_60p', extrapolate=False)
        arr = np.loadtxt(os.environ['GIDGETDIR']+'/py/CanoDiaz16_80p.csv', delimiter=',')
        datasets['CanoDiaz16_80p'] = DataSet('colst', 'colsfr', np.power(10.0,arr[:,0]-6), np.power(10.0,arr[:,1]), zmax= 0.5, label='CanoDiaz16_80p', extrapolate=False)




    cos = halo.Cosmology()
    mstSpeagle = np.power(10.0, np.linspace(9.7,11.1, 10)) # just linear - no need for a large number of samples
    psiSpeagle = (0.84 - 0.026*cos.age(z))*np.log10(mstSpeagle) - (6.51 - 0.11*cos.age(z))
    fac = 1.0
    if specific:
        fac = mstSpeagle*1.0e-9
    if specific:
        datasets['Speagle14Specz'+str(z)] = DataSet('mstar', 'sSFR', mstSpeagle, np.power(10.0,psiSpeagle)/fac, yLower=np.power(10.0,psiSpeagle-0.28)/fac, yUpper=np.power(10.0,psiSpeagle+0.28)/fac, label='Speagle14', zmin=z-0.5, zmax=z+0.5, alpha=0.5)
    else:
        datasets['Speagle14z'+str(z)] = DataSet('mstar', 'sfr', mstSpeagle, np.power(10.0,psiSpeagle)/fac, yLower=np.power(10.0,psiSpeagle-0.28)/fac, yUpper=np.power(10.0,psiSpeagle+0.28)/fac, label='Speagle14', zmin=z-0.5, zmax=z+0.5)

            

    thisMst = np.power(10.0, np.linspace(10,12,100))
    whitaker12 = np.power(10.0, -1.12 + 1.14*z - 0.19*z*z - (0.3+0.13*z)*(np.log10(thisMst)-10.5)) # Gyr^-1 -- Genzel+15 eq 1
    mmin=None
    if 2.0<=z and z<2.5:
        a = -19.99
        sa = 1.87
        b = 3.44
        sb = 0.36
        c = -0.13
        sc = 0.02
        mmin= 9.3
        mmax = 11.4
        zmin = 2.0
        zmax = 2.5
    if 1.5<=z and z<2.0:
        a = -24.04
        sa = 2.08
        b = 4.17
        sb = 0.40
        c = -0.16
        sc = 0.02
        mmin= 9.2
        mmax = 11.5
        zmin = 1.5
        zmax = 2.0
    if 1.0<=z and z<1.5:
        a = -26.03
        sa = 1.69 
        b = 4.62
        sb = 0.34
        c = -0.19
        sc = 0.02
        mmin= 8.8
        mmax = 11.3
        zmin = 1.0
        zmax = 1.5
    if 0.5<=z<1.0:
        a = -27.40
        sa = 1.91
        b = 5.02
        sb = 0.39
        c = -0.22
        sc = 0.02
        mmin= 8.4
        mmax = 11.2
        zmin = 0.5
        zmax = 1.0


    if z<2:
        lilly13 = 0.117 * np.power(thisMst/3.16e10, -0.1) * np.power(1.0+z,3.0)
    else:
        lilly13 = 0.5 * np.power(thisMst/3.16e10, -0.1) * np.power(1.0+z,1.667)


    if not mmin is None:
        mwhitaker = np.power(10.0, np.linspace(mmin,mmax,100))
        if specific:
            datasets['whitaker14Specz'+str(z)] = DataSet( 'mstar', 'sSFR',mwhitaker, (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0))/mwhitaker * 1.0e9, yLower = (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 - 0.34))/mwhitaker * 1.0e9,  yUpper= (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 + 0.34))/mwhitaker * 1.0e9, label='Whitaker14', zmin=zmin, zmax=zmax )
            datasets['whitaker12Specz'+str(z)] = DataSet( 'mstar', 'sSFR',thisMst, whitaker12, yLower=whitaker12/10.0**0.34, yUpper=whitaker12*10.0**0.34, label='Whitaker12', zmin=z-0.5, zmax=z+0.5)
        else:
            datasets['whitaker14z'+str(z)] = DataSet('mstar', 'sfr', mwhitaker, (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0)), yLower = (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 - 0.34)),  yUpper= (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 + 0.34)), label='Whitaker14', zmin=zmin, zmax=zmax )
            datasets['whitaker12z'+str(z)] = DataSet( 'mstar', 'sfr',thisMst, whitaker12*thisMst*1.0e-9, yLower=whitaker12/10.0**0.34 *thisMst*1.0e-9, yUpper=whitaker12*10.0**0.34*thisMst*1.0e-9, label='Whitaker12', zmin=z-0.5, zmax=z+0.5 )
    if specific:
        datasets['lilly13Specz'+str(z)] = DataSet('mstar', 'sSFR',thisMst, lilly13, yLower=lilly13/10.0**0.34, yUpper=lilly13*10.0**0.34, label='Lilly13', zmin=z-0.5, zmax=z+0.5, alpha=0.5)
    else:
        datasets['lilly13z'+str(z)] = DataSet( 'mstar', 'sfr',thisMst, lilly13*thisMst*1.0e-9, yLower=lilly13/10.0**0.34*thisMst*1.0e-9, yUpper=lilly13*10.0**0.34*thisMst*1.0e-9, label='Lilly13', zmin=z-0.5, zmax=z+0.5 )



def defineStellarZ():
    MstKirby = np.power(10.0, np.linspace(3,9,20))
    ZKirby = np.power(10.0, -1.69 + 0.30* np.log10(MstKirby/1.0e6))

    # These are actually from Gallazzi05. (Tremonti is on the paper, so I don't feel that bad)
    MstTremonti = np.power(10.0, np.array([8.91, 9.11, 9.31, 9.51, 9.72, 9.91, 10.11, 10.31, 10.51, 10.72, 10.91, 11.11, 11.31, 11.51, 11.72, 11.91]))
    ZTremontiMed = np.power(10.0, np.array([-0.60, -0.61, -0.65, -0.61, -0.52, -0.41, -0.23, -0.11, -0.01, 0.04, 0.07, 0.10, 0.12, 0.13, 0.14, 0.15]))
    ZTremonti16 = np.power(10.0, np.array([-1.11, -1.07, -1.10, -1.03, -0.97, -0.90, -0.80, -0.65, -0.41, -0.24, -0.14, -0.09, -0.06, -0.04, -0.03, -0.03]))
    ZTremonti84 = np.power(10.0, np.array([-0.00, -0.00, -0.05, -0.01, 0.05, 0.09, 0.14, 0.17, 0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.30]))

    datasets['kirby13'] = DataSet( 'mstar', 'stZ', MstKirby, ZKirby, yLower=ZKirby*10.0**-0.17, yUpper=ZKirby*10.0**0.17, label='Kirby13', alpha=0.5)
    datasets['gallazi05'] = DataSet( 'mstar', 'stZ', MstTremonti, ZTremontiMed, yLower=ZTremonti16, yUpper=ZTremonti84, label='Gallazzi05', alpha=0.5, fixedSigma=0.17)



def defineGasFractions(z):
    peeplesMode = np.loadtxt(os.environ['GIDGETDIR']+'/py/peeples11_fgmode.csv', delimiter=',')
    peeples84 = np.loadtxt(os.environ['GIDGETDIR']+'/py/peeples11_fg68.csv', delimiter=',')
    peeplesMass = peeplesMode[:,0]
    delt = peeples84[:,1] - peeplesMode[:,1]
    datasets['peeples11'] = DataSet('mstar', 'gasToStellarRatioHI',  np.power(10.0, peeplesMass), np.power(10.0, peeplesMode[:,1]), yLower = np.power(10.0, peeplesMode[:,1]-delt), yUpper = np.power(10.0, peeplesMode[:,1]+delt), label='Peeples11' , alpha=0.5)


    thisMst = np.power(10.0, np.linspace(7.0,11.5))
    thisMHI = np.power(10.0, -0.43*np.log10(thisMst) + 3.75) * thisMst
    fgThis = thisMHI/thisMst
    datasets['papastergis12'] = DataSet( 'mstar', 'gasToStellarRatioHI', thisMst, fgThis, yLower=fgThis/10.0**0.31, yUpper=fgThis*10.0**0.31, label='Papastergis12', alpha=0.5) # 0.31 dex is the /minimal/ scatter in MHI/M* as a function of various things in Zheng2009. No scatter quoted in Papastergis '12

    thisMst = np.power(10.0, np.linspace(10.0,11.5, 100))
    def molToStarGenzel( (af2, xif2, xig2, xih2), deltaMS=0.0):
        return np.power(10.0, af2 + xif2*np.log10(1.0+z) + xig2*deltaMS + xih2*np.log10(thisMst/10.0**10.77))
    def ratioToFraction(rat):
        ## f_g = Mg/(Mg+M*) = 1/(1+M*/Mg) = 1/(1+1/rat)
        return 1.0/(1.0 + 1.0/rat)
                
    #datasets['Genzel15Lillyz'+str(z)] = DataSet( thisMst, molToStarGenzel(-0.98,2.65,0.5,-0.25), yLower=molToStarGenzel(-0.98,2.65,0.5,-0.25), yUpper=molToStarGenzel(-0.98,2.65,0.5,-0.25) )
    #datasets['Genzel15FMRz'+str(z)] = DataSet( thisMst, molToStarGenzel(-1.05,2.6,0.54,-0.41), yLower=molToStarGenzel(-1.05,2.6,0.54,-0.41), yUpper=molToStarGenzel(-1.05,2.6,0.54,-0.41) )
    #datasets['Genzel15Lillyz'+str(z)] = DataSet( thisMst, molToStarGenzel(-0.98,2.65,0.5,-0.25), yLower=molToStarGenzel(-0.98,2.65,0.5,-0.25), yUpper=molToStarGenzel(-0.98,2.65,0.5,-0.25) )

    genzelParamsCO = (-1.23, 2.71, 0.51, -0.35)
    genzelParamsCOGlobal = (-1.12, 2.71, 0.53, -0.35)
    genzelParamsDust = (-0.87, 2.26, 0.51, -0.41)
    genzelParamsDustGlobal= (-0.98, 2.32, 0.36, -0.40)
    datasets['genzel15COz'+str(z)] = DataSet( 'mstar', 'gasToStellarRatioH2', thisMst, molToStarGenzel(genzelParamsCOGlobal), yLower=molToStarGenzel(genzelParamsCOGlobal,-0.34), yUpper=molToStarGenzel(genzelParamsCOGlobal,0.34), label='Genzel15 CO', zmin=z-0.5, zmax=z+0.5, alpha=0.5)
    datasets['genzel15Dustz'+str(z)] = DataSet( 'mstar', 'gasToStellarRatioH2', thisMst, molToStarGenzel(genzelParamsDustGlobal), yLower=molToStarGenzel(genzelParamsDustGlobal,-0.34), yUpper=molToStarGenzel(genzelParamsDustGlobal,0.34), label='Genzel15 Dust', zmin=z-.5, zmax=z+0.5, alpha=0.5)
    #ax.plot(thisMst, molToStarGenzel(-0.98,2.65,0.5,-0.25), c='green', label=thisLabel1)
    #ax.plot(thisMst, molToStarGenzel(-1.05,2.6,0.54,-0.41), c='orange', label=thisLabel2)
    #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
    #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)

    saintongeM = np.power(10.0, np.linspace(10.3, 11.4, 20))
    saintongeFr = np.power(10.0, -1.607 - 0.455*(np.log10(saintongeM) - 10.70))
    uncertainty = 10.0**.331 # uncertainty on b in Table 3 (undetected@upper limits) of Saintonge11.
    datasets['saintonge11'] = DataSet( 'mstar', 'gasToStellarRatioH2', saintongeM, saintongeFr, yLower=saintongeFr/uncertainty, yUpper=saintongeFr*uncertainty, label='Saintonge11' )    

def defineTullyFisher():
    mstMiller = np.power(10.0, np.linspace(9.0, 11.5, 20))
    vMiller = np.power(10.0, (np.log10(mstMiller) - 1.718)/3.869)
    datasets['miller11'] = DataSet( 'mstar', 'vPhi22', mstMiller, vMiller, yLower=vMiller/10.0**0.058, yUpper=vMiller*10.0**0.058, label='Miller11', zmax=1.25, alpha=0.5 )


    sf = np.power(10.0, np.array([-3.1, -2.3, -1.4, -0.7, -0.0, 0.54, 0.89, 1.3, 1.67, 2.1]))
    sig = np.array([9.3, 11.5, 14.6, 14.8, 17.6, 21.6, 32.9, 53.5, 72.3, 87.6])
    arr = np.loadtxt(os.environ['GIDGETDIR']+'/py/krumholz17_sigma_compilation.txt')
    sf = np.power(10.0, arr[:,0])
    sig16 = np.power(10.0, arr[:,1])
    sig50 = np.power(10.0, arr[:,2])
    sig84 = np.power(10.0, arr[:,3])

    #datasets['krumholz17'] = DataSet('sfr', 'sfsig', sf, sig, yLower=sig/2.0, yUpper=sig*2.0, label='Krumholz17', zmin=-0.5, zmax=3.5, alpha=0.5)
    #datasets['krumholz17max'] = DataSet('sfr', 'maxsig', sf, sig, yLower=sig/2.0, yUpper=sig*2.0, label='Krumholz17', zmin=-0.5, zmax=3.5, alpha=0.5)
    datasets['krumholz17'] = DataSet('sfr', 'sfsig', sf, sig50, yLower=sig16, yUpper=sig84, label='Krumholz17', zmin=-0.5, zmax=3.5, alpha=0.5)
    datasets['krumholz17max'] = DataSet('sfr', 'maxsig', sf, sig50, yLower=sig16, yUpper=sig84, label='Krumholz17', zmin=-0.5, zmax=3.5, alpha=0.5)

#
#if xvar=='gbar' and v=='gtot':
#    gbars = np.power(10.0, np.linspace(-12,-8, 100))
#    gdagger = 1.20e-10
#    gobss = gbars/(1.0 - np.exp(-np.sqrt(gbars/gdagger)))
#    thisLabel=None
#    thisLS = '--'
#    if z[ti]<0.5:
#        thisLabel='McGaugh16'
#        thisLS = '-'
#    if histFlag:
#        ax.plot( np.log10(gbars), np.log10(gobss), label=thisLabel, ls=thisLS, c='b', lw=3 )
#        ax.plot( np.log10(gbars), np.log10(gbars), ls='--', c='gray', lw=3)
#    else:
#        ax.plot( gbars, gobss, label=thisLabel, ls=thisLS, c='b', lw=3 )
#        ax.plot( gbars, gbars, ls='--', c='gray', lw=3)
#    labelled=True


# Now let's put in some data.
datasets = {}
defineTullyFisher()
defineStellarZ()
defineStructureRelations()
defineAngularMomenta()
defineBroeils()
defineBrinchmann()
defineGarrisonKimmel()
defineKravtsov13()
# Things where the relations change as a fn of reshift 
for i in [0,1,2,3,4]:
    defineMetalRelations(i)
    defineMS(i, specific=True)
    defineMoster(i)
    defineGasFractions(i)

for dsk in datasets.keys():
    ds = datasets[dsk]
    if ds.extrapolate:
        ds.cacheQuantiles()


def identifyApplicableDatasets(xvar, yvar, z):
    res = []
    labelsInOrder = []
    uniqueLabels = []
    validRedshift = []
    for ds in sorted( datasets.keys() ):
        if datasets[ds].xvar == xvar and datasets[ds].yvar ==yvar:
            res.append( ds )
            labelsInOrder.append( datasets[ds].label )
            validRedshift.append( datasets[ds].zmin < z and z<datasets[ds].zmax )
            if not datasets[ds].label in uniqueLabels:
                uniqueLabels.append( datasets[ds].label )
    actualRes = []
    colors = []
    for il, ul in enumerate(uniqueLabels):
        # Get the datasets with the same label as ul.
        # Labels are e.g. Moster10, Miller11,... i.e. specific papers
        selec = ul == np.array(labelsInOrder)
        # For each such work, add in all of the datasets that overlap the requested reshift z.
        # I guess this should only be at most 1 dataset/redshift...
        if np.any( np.array(validRedshift)[ selec ] ):
            for ele in np.array(res)[selec][np.array(validRedshift)[selec]]:
                actualRes.append(ele)
                colors.append(il)
        else:
            actualRes.append( np.array(res)[selec][0] )
            colors.append(il)

    return actualRes, colors



if __name__=='__main__':
    for ds in datasets.keys():
        #datasets[ds].test()
        if datasets[ds].extrapolate:
            datasets[ds].plotLikelihood()



