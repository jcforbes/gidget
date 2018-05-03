import numpy as np
import scipy.stats
import copy

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
def lnparetodensity(theta, a):
    if theta<1.0:
        return -np.inf
    return np.log(a) - (a+1.0)*np.log(theta)

def samplefrombetadensityTrans(u,a,b):
    assert a>0 and b>0
    return scipy.stats.beta.ppf(u, a,b)
def samplefromgammadensityTrans(u,a,b):
    assert a>0 and b>0
    return scipy.stats.gamma.ppf(u, a, scale=1.0/b)
    #return np.random.gamma(a,1.0/b) # numpy's definition of the gamma uses e^{-x/scale}/scale^a
def samplefromlognormaldensityTrans(u,mean,var):
    assert var>0
    return np.exp(samplefromnormaldensityTrans(u,mean,var))
def samplefromnormaldensityTrans(u,mean,var):
    assert var>0
    return scipy.stats.norm.ppf(u,mean,np.sqrt(var))
    #return np.random.normal(mean,np.sqrt(var))
def samplefromloguniformdensityTrans(u,a,b):
    assert b>a
    assert a>0
    return a*(b/a)**u
    #return a*(b/a)**np.random.uniform()
def samplefromuniformdensityTrans(u,a,b):
    assert b>a
    return (b-a)*u+a
    #return np.random.uniform(a,b)
def samplefromparetodensityTrans(u,a):
    assert a>0.0
    return scipy.stats.pareto.ppf(u, a )
    #return np.random.paretoTrans(a)


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
def samplefromparetodensity(a):
    assert a>0.0
    return np.random.pareto(a)+1.0

class simpleDistribution:
    ''' A little wrapper for a 1D distribution. Give it the name and parameters of the distribution
        you want to use, and then later you can sample from this distribution or get the ln of the density at a given point.'''
    def __init__(self, token, parameters,name):
        self.token=token # specify the family of analytic distribution to use
        self.name=name # store a human-recognizable name for this variable.
        self.parameters=copy.deepcopy(parameters)
        if token=='beta':
            assert len(parameters)==2
            self.samp = samplefrombetadensity
            self.samp_trans = samplefrombetadensityTrans
            self.dens = lnbetadensity
        elif token=='gamma':
            assert len(parameters)==2
            self.samp = samplefromgammadensity
            self.samp_trans = samplefromgammadensityTrans
            self.dens = lngammadensity
        elif token=='lognormal':
            assert len(parameters)==2
            self.samp = samplefromlognormaldensity
            self.samp_trans = samplefromlognormaldensityTrans
            self.dens = lnlognormaldensity
        elif token=='normal':
            assert len(parameters)==2
            self.samp = samplefromnormaldensity
            self.samp_trans = samplefromnormaldensityTrans
            self.dens = lnnormaldensity
        elif token=='loguniform':
            assert len(parameters)==2
            self.samp = samplefromloguniformdensity
            self.samp_trans = samplefromloguniformdensityTrans
            self.dens = lnloguniformdensity
        elif token=='uniform':
            assert len(parameters)==2
            self.samp = samplefromuniformdensity
            self.samp_trans = samplefromuniformdensityTrans
            self.dens = lnuniformdensity
        elif token=='pareto':
            assert len(parameters)==1
            self.samp = samplefromparetodensity
            self.samp_trans = samplefromparetodensityTrans
            self.dens = lnparetodensity
        else:
            print "Didn't recognize the requested token", token
            raise ValueError
    def sample(self):
        return self.samp(*self.parameters)
    def sample_transform(self,u):
        return self.samp_trans(u,*self.parameters)
    def lndensity(self, x):
        return self.dens(x,*self.parameters)


class jointDistribution:
    ''' So far we've had a "prior" defined in two ways below: one function that draws samples from it, and one that 
        evaluates it at a given point. This is a little painful because any changes to the prior require updating
        both in the same way, so ideally this should be taken care of by the same object. This object.'''
    def __init__(self, list_of_distributions):
        self.list_of_dist = copy.deepcopy(list_of_distributions)
    def sample(self):
        return [distr.sample() for distr in self.list_of_dist]
    def sample_transform(self, u):
        return [self.list_of_dist[i].sample_transform(u[i])  for i in range(len(self.list_of_dist))]
    def lndensity(self, vec):
        assert len(vec) == len(self.list_of_dist)
        result = 0.0
        for i,v in enumerate(vec):
            result += self.list_of_dist[i].lndensity(v)
        if not np.isfinite(result):
            return -np.inf
        return result

