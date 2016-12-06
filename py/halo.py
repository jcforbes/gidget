import numpy as np
import scipy.special
import pdb


gpermsun = 1.9891e33
speryear = 3.15569e7
cmperkpc = 3.08567758e21
kB = 1.3806488e-16
mp = 1.67262178e-24 # proton mass in grams
me = 9.10938291e-28 # electron mass in grams
cmperMpc = 1000.0*cmperkpc
G = 6.67384e-8
esu = 4.80320451e-10 # electron charge in esu, where esu^2 ~ g cm^3/s^2
c = 2.99792458e10



class halo:
    def __init__(self,mass,z, cos, crf):
        zz = z
        if zz<0:
            zz=0 # Floor z at 0.
        self.mass=mass # Msun
        self.crf = crf
        self.radius = (mass * gpermsun * 3.0/(4.0*np.pi*200.0*cos.rhocrit(zz)))**(1.0/3.0) / cmperkpc # kpc
        self.vcirc = np.sqrt( G*mass*gpermsun / (self.radius*cmperkpc )) * 1.0e-5 # km/s
        self.T = mp * (self.vcirc*1.0e5)**2.0 / kB # K
        m = np.log10(mass*cos.h/1.0e12)
        nu = 10.0 ** ( -0.11 + 0.146*m + 0.0138*m*m + 0.00123*m*m*m) * (0.033 + 0.79*(1.0+zz) + 0.176*np.exp(-1.356*zz))
        a = 0.520 + (0.905 - 0.520)*np.exp(-0.617*zz**1.21)
        b = -0.101 + 0.026*zz
        self.c200 = 10.0 ** (a+b*m + crf)
        if self.c200!=self.c200:
            pdb.set_trace()
        self.alpha = 0.0095*nu*nu + 0.155
    def mInterior(self,r):
        # r in kpc
        c200 = self.c200
        rScale = self.radius*cmperkpc / c200
        alpha = self.alpha
        self.rhoScale = self.mass/( 2.0**(2.0-3.0/alpha)*np.exp(2.0/alpha)*np.pi*alpha**(+3./alpha-1.0)*rScale*rScale*rScale*scipy.special.gamma(3.0/alpha)*( scipy.special.gammainc(3.0/alpha, 2.0*(c200**alpha)/alpha))/gpermsun)
        val = self.rhoScale* 2.0**(2.0-3.0/alpha)*np.exp(2.0/alpha)*np.pi*alpha**(+3./alpha-1.0)*rScale*rScale*rScale* scipy.special.gamma(3.0/alpha)*( scipy.special.gammainc(3.0/alpha, 2.0*((r*cmperkpc/rScale)**alpha)/alpha))
        if val!=val:
            pdb.set_trace()
        return val # grams
      #   val = rhoScaleEinasto(Mh,z) * exp( -2.0/alpha * (pow(x[n]*Radius/rScale, alpha) - 1.0) );
      #   rho[n] = val;:
    def g(self, r):
        return G*self.mInterior(r)/(r*cmperkpc)**2.0 # cm/s^2
    def rho(self,r):
        ''' Return the density in g/cm**3 of the dark matter halo at the given spherical radius. Expects r in kpc.'''
        mi = self.mInterior(r)
        if r>5000:
            print "Did you accidentally request halo.rho(r) with r in cm instead of kpc?"
            assert False
        rScale = self.radius / self.c200 # kpc
        return self.rhoScale * np.exp( -2.0/self.alpha * (np.power(r/rScale, self.alpha) - 1.0) );



class Cosmology:
    def __init__(self,H0=(70.0*1.0e5/cmperMpc),OmM=0.3,OmL=0.7):
        self.H0=H0
        self.OmM = OmM
        self.OmL = OmL
        self.OmK = 1.0-OmM-OmL
        self.h = H0 / (100.0*1.0e5/cmperMpc)
    def EE(self,z):
        return np.sqrt(self.OmM*(1.0+z)**3.0 + self.OmK*(1.0+z)**2.0+self.OmL)
    def Hubble(self,z):
        return self.EE(z)*self.H0
    def rhocrit(self,z):
        H = self.Hubble(z)
        return 3.0*H*H/(8.0*np.pi*G)
    def tL(self,z):
        ''' Returns the lookback time at redshift z in years, assuming the given cosmology'''
        def toIntegrate(zz):
            return 1.0/((1.0+zz)*self.EE(zz))
        ret = scipy.integrate.quad(toIntegrate,0,z)[0]/self.H0/speryear
        return ret
    def age(self,z):
        ''' Returns the age of the universe at that redshift in Gyr. '''
        return (self.tL(10000.0) - self.tL(z))/1.0e9
    def dtdz(self,z):
        return 1.0/((1.0+z)*self.EE(z))/self.H0/speryear # units of years
    


def testhalo():
    import matplotlib.pyplot as plt
    
    cos = Cosmology()
    h = halo(1.0e12, 0.01, cos, 0)
    rs = np.power(10.0, np.linspace(-2,3,1000))
    ms = [h.mInterior(r) for r in rs]
    rhos = [h.rho(r) for r in rs]

    rho3 = rhos[-1]*np.power(rs[-1]/rs[0], 3.0)
    rho1 = rhos[0]*np.power(rs[-1]/rs[0], -1.0)

    fig,ax = plt.subplots()
    ax.plot(rs,rhos)
    ax.plot([rs[0],rs[-1]], [rho3,rhos[-1]], ls='--', c='k')
    ax.plot([rs[0],rs[-1]], [rhos[0],rho1], ls='--', c='k')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'r (kpc)')
    ax.set_ylabel(r'$\rho\ (g/cm^3)$')
    plt.savefig('test_halo_density.png')
    plt.close(fig)

    Gcgs = 6.67e-8
    cmperkpc = 3.1e21
    fig,ax = plt.subplots()
    ax.plot( rs, np.sqrt( np.array(ms) * Gcgs/ (rs*cmperkpc) )*1.0e-5 )
    ax.set_xscale('log')
    ax.set_xlabel('r (kpc)')
    ax.set_ylabel(r'Circular Velocity (km/s)')
    plt.savefig('test_halo_vphi.png')
    plt.close(fig)


if __name__=='__main__':
    testhalo()



