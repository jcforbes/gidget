import numpy as np
import scipy.integrate as scint
from scipy.interpolate import interp1d
import scipy.optimize
import matplotlib.pyplot as plt
import pdb
import halo

cmperkpc = 3.18e21
cmperpc = 3.18e18
sperMyr = 3.16e13
gpermsun = 1.989e33
G = 6.67e-8 * gpermsun*sperMyr**2/cmperpc**3

def dKzdz(z, r, Mh, redshift, crf):
    # z here refers to vertical height above the disk in pc
    cos = halo.Cosmology()
    h = halo.halo(Mh, redshift, cos, crf) # Assumes Mh in solar masses
    def Kz(z0):
        sphRad = np.sqrt(r*r+z0*z0/1.0e6)
        g = h.g(sphRad) # cm/s^2, assumes sphRad is in kpc
        geom = (z0/1.0e3)/sphRad
        assert geom<1
        return -g* geom * sperMyr**2/cmperpc
    hh = z*1.0e-2

    if z>0:
        res = (Kz(z+hh)-Kz(z-hh))/(2*hh) # returned in the same units as the rest of this module, 1/Myr**2
    else:
        res=0
    if res>0:
        print ("WARNING: dKzdz>0 ",res,z,hh)
    return res

def ddz(y, z0, *args):
    assert len(y) % 2 == 0
    N = len(y)/2
    sigsq = args[:N]
    r = args[N]
    Mh = args[N+1]
    redshift = args[N+2]
    crf = args[N+3]
    try:
        assert len(sigsq)==N
    except:
        print ("Failed sigsq check! ")
        print ("y=",y)
        print ("N=",N)
        print ("sigsq=",sigsq)
    deriv = np.zeros(2*N)
    for i in range(N):
        #if y[i]<0:
        #    print "Warning negative density i=",i,"z0=",z0
        #if y[N+i]>0:
        #    print "Warning setting density derivative to 0. i=",i,"z0=",z0
        deriv[i] = y[N+i]
        assert len(y[:N] == N)
        totDens=0
        for j in range(N):
            if y[j]>0:
                totDens += y[j]
        deriv[N+i] = 0
        if y[i]>0:
            densTerm = -4*np.pi*G*totDens
            dmTerm = dKzdz(z0,r,Mh,redshift,crf)
            #dmTerm = 0.0
            #print "Poisson terms: ", dmTerm, densTerm
            deriv[N+i] += y[i]/sigsq[i] * (densTerm + dmTerm) + 1.0/y[i] * y[N+i]*y[N+i]
        if y[i]<1.0e-7:
            deriv[i] = 0.0
            deriv[N+i] = 0.0
        if deriv[i]>0:
            deriv[i]=0.0
            deriv[N+i] = 0.0

    #print "Evaluating ddz. sigsq = ",sigsq
    #print " and deriv = ",deriv
    return deriv

def non_increasing(L):
    return np.all(x<=y for x,y in zip(L,L[1:]))

def singleIteration(rho0s, sigsq, r,Mh,redshift,crf):
    ''' Run the scipy integrator and see what values we get for the surface density.'''
    # Initial conditions: midplane densities and zero d rho/dz
    y0 = np.array(list(rho0s)+[0]*len(rho0s))
    # Set up grid on which we want to compute rho(z)
    zs = np.linspace( 0, 5000.0, 10000 ) # integrate the solution to much greater heights
    #   than we actually need. This is so that the column density we get in the end is accurate.
    # Feed in the args ddz needs: the square velocity dispersions and the halo parameters
    args = tuple(list(sigsq) + [r,Mh,redshift,crf])
    # Integrate using these IC's 
    y = scint.odeint(ddz, y0, zs, args=args, mxstep=50000, rtol=1.0e-13)
    cols =[]
    # Now compute the resulting column densities
    for i in range(len(rho0s)):
        if non_increasing(y[:,i]):
            pos = y[:,i]>0
            xx = zs[pos]
            yy = y[:,i][pos]
            assert len(xx)==len(yy)
            try:
                rhoInterp = interp1d( xx,yy, kind='linear')
                coli, err =  scint.quad( rhoInterp, 0, 0.99*np.max(xx)) 
                cols.append(coli)
            except:
                cols.append(0.0)
        else:
            cols.append(-100) # the solution is garbage and shouldn't be interpolated
    # Return the column densities. The calling function will then adjust rho0s
    #  until the column densities agree with the real column densities.
    return np.array(cols)

def run(colsData, sigsq, r,Mh,redshift,crf):
    def toZero(rhosAtMidplane):
        colsSoln = singleIteration(rhosAtMidplane, sigsq, r,Mh,redshift,crf)
        return (colsSoln - colsData)/colsData

    hgGuess = sigsq/(np.pi*G*colsData)
    rhoGuess = colsData / hgGuess
    sol = scipy.optimize.root(toZero, rhoGuess, tol=1.0e-11)

    rho0s = sol.x
    y0 = np.array( list(rho0s)+[0]*len(rho0s) )
    zs = np.linspace(0, 1500.0, 100)
    args = tuple(list(sigsq) + [r,Mh,redshift,crf])
    y = scint.odeint(ddz, y0, zs, args=args, mxstep=50000, rtol=1.0e-13)
    return y,zs

def weights(y,zs):
    shp = np.shape(y)
    print ("In weights, shape of y is ", shp)
    cols = np.zeros( (shp[1]/2, 4) )
    zwindows = [ [0.15,0.25], [0.25,0.5], [0.5,1.0], [1.0,1.5] ]
    for i in range(shp[1]/2):
        rhoInterp = interp1d( zs, y[:,i], kind='cubic' )
        for j in range(len(zwindows)):
            coli, err = scint.quad( rhoInterp, zwindows[j][0]*1000.0, zwindows[j][1]*1000.0)
            if coli>0:
                cols[i,j] = coli
    print ("In weights, shape of cols is ",np.shape(cols) )
    # We want each window to be normalized to 1. I.e. each column in the matrix will tell us
    #   what fraction of the mass in a given z-window is in each component
    for i in range(4):
        cols[:,i] = cols[:,i]/np.sum(cols[:,i])
    return cols
        

def plot(y,zs, axIn=None, ls='-'):
    if axIn is None:
        fig,ax = plt.subplots()
    else:
        ax = axIn
    shp = np.shape(y)
    colors=['k','r','b','g','orange','purple','lightblue','pink','grey']*5
    for i in range(shp[1]):
        ax.plot( zs, y[:,i], color=colors[i], lw=i+1, ls=ls )
    ax.set_xlabel('z (pc)')
    ax.set_ylabel(r'$\rho (M_\odot/\mathrm{pc}^3)$')
    ax.set_yscale('log')
    ax.set_ylim(1.0e-5,10.0)
    plt.savefig('rho_vs_z.png')
    if axIn is None:
        plt.close(fig)

def anEstimate(sigsq, colData):
    hGuessArr = []
    zArr = np.linspace(0,1500, 1000)
    y = np.zeros((1000,len(sigsq)*2))
    for i in range(len(sigsq)):
        fs = np.clip( np.sqrt(sigsq[i]/sigsq), 0, 1)
        colEff = np.sum(colData*fs)
        hGuessThis = sigsq[i]/(2.0*np.sqrt(2.0)*np.pi*G*colEff)
        hGuessArr.append(hGuessThis)
        rho0Guess = colData[i]/hGuessThis
        y[:,i] = rho0Guess * np.power(1.0/np.cosh(zArr/hGuessThis),2.0)
    return y,zArr


def test0():
    z = 300
    r= 8
    Mh = 10e12
    redshift=0
    crf=3.0

    print (dKzdz(z, r, Mh, redshift, crf))

def test():
    sigsq = np.power(np.array([8.0, 5.0, 10.0, 20.0, 50.0]), 2.0)
    colData = np.array([10.0, 30.0, 20.0, 20.0, 20.0])
    r, Mh, redshift, crf = (10.0, 1.0e12, 0, 3)
    y,zs = run(colData,  sigsq, r,Mh, redshift, crf)
    yAn,zAn = anEstimate(sigsq,colData)
    fig,ax = plt.subplots()
    plot(y,zs,axIn=ax)
    plot(yAn,zAn,axIn=ax, ls='--')
    plt.close(fig)
    print ("WeightsNumerical: ",weights(y,zs))
    print ("WeightsAnalytic: ",weights(yAn,zAn))

def test2():
    y,zs = run(np.array([.002, .02]), np.array([8, 20])*1.0e5 )
    plot(y,zs)
    print ("Weights: ",weights(y,zs))

def test3():
    rho0s = np.array([10.0,100.0])/500.0
    sigsq = np.power( np.array([8,30]), 2.0)

    y0 = np.array(list(rho0s)+[0]*len(rho0s))
    zs = np.linspace( 0, 1500, 100 )
    args = tuple(sigsq)
    y = scint.odeint(ddz, y0, zs, args=args, mxstep=50000, rtol=1.0e-12)


    print (y )


if __name__=='__main__':
    test()

