import numpy as np
import scipy.special


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
        self.mass=mass # Msun
        self.crf = crf
        self.radius = (mass * gpermsun * 3.0/(4.0*np.pi*200.0*cos.rhocrit(z)))**(1.0/3.0) / cmperkpc # kpc
        self.vcirc = np.sqrt( G*mass*gpermsun / (self.radius*cmperkpc )) * 1.0e-5 # km/s
        self.T = mp * (self.vcirc*1.0e5)**2.0 / kB # K
        m = np.log10(mass*cos.h/1.0e12)
        nu = 10.0 ** ( -0.11 + 0.146*m + 0.0138*m*m + 0.00123*m*m*m) * (0.033 + 0.79*(1.0+z) + 0.176*np.exp(-1.356*z))
        a = 0.520 + (0.905 - 0.520)*np.exp(-0.617*z**1.21)
        b = -0.101 + 0.026*z
        self.c200 = 10.0 ** (a+b*m + crf)
        self.alpha = 0.0095*nu*nu + 0.155
    def mInterior(self,r):
        # r in kpc
        c200 = self.c200
        rScale = self.radius*cmperkpc / c200
        alpha = self.alpha
        self.rhoScale = self.mass/( 2.0**(2.0-3.0/alpha)*np.exp(2.0/alpha)*np.pi*alpha**(+3./alpha-1.0)*rScale*rScale*rScale*scipy.special.gamma(3.0/alpha)*( scipy.special.gammainc(3.0/alpha, 2.0*(c200**alpha)/alpha))/gpermsun)
        val = self.rhoScale* 2.0**(2.0-3.0/alpha)*np.exp(2.0/alpha)*np.pi*alpha**(+3./alpha-1.0)*rScale*rScale*rScale* scipy.special.gamma(3.0/alpha)*( scipy.special.gammainc(3.0/alpha, 2.0*((r*cmperkpc/rScale)**alpha)/alpha))
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


