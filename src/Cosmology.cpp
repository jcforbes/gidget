#include "Cosmology.h"
#include "Dimensions.h"
#include <math.h>
#include <iostream>
#include <gsl/gsl_sf_gamma.h>

Cosmology Cosmology::operator=(const Cosmology& c)
{
    Cosmology toReturn(c.OmegaM(),c.OmegaL(),c.Hubble(), c.sigma8(), c.zs(), c.GetRho().size());
    return toReturn;
}


// The following two functions are inverses of each 
// other. lbt is the lookback time for a matter+lambda 
// - dominated universe.
double Cosmology::zLBT(double lbTime)
{
  double time = lbt(1.0e30)-lbTime;
  double ScaleFactor = 
    pow((OmM/OmL)*sinh(1.5*H0*time*sqrt(OmL))
	*sinh(1.5*H0*time*sqrt(OmL))
	,1./3.);
  return 1./ScaleFactor-1.;
}

// adaptation of Carroll & Ostlie Eq. 29.129
// look-back time in seconds
double Cosmology::lbt(double z)
{
  double ScaleFactor = 1./(1.+z);
  double ret= (2./(3.*H0*sqrt(OmL)) * 
	  log(sqrt(OmL/OmM)
	      +sqrt(1.0+OmL/OmM)))
    - (2./(3.*H0*sqrt(OmL)) * 
       log(sqrt(OmL*ScaleFactor*ScaleFactor*
		ScaleFactor/OmM)
	   +sqrt(1.0+OmL*ScaleFactor*ScaleFactor*
		 ScaleFactor/OmM)));
//  std::cerr << "lbt: " << ret/speryear * 1.0e-9 << ",  z: "<<z<< std::endl;
  return ret;
}
double Cosmology::Tsim(double z) 
{
  return lbt(zstart)-lbt(z);
}

double dz(double dtSim, double z, Cosmology& cos, Dimensions& dim)
{
  double dtPhys = dim.t_cgs(dtSim);
  return fabs( z - cos.zLBT(cos.lbt(z)-dtPhys)  );
}


double Cosmology::EE(double z)
{
  return sqrt(OmM*(1+z)*(1+z)*(1+z) 
	      + OmK*(1+z)*(1+z) 
	      + OmL);   
}

// the Hubble constant in 1/s
double Cosmology::Hubble(double z)
{
    return EE(z)*H0;
}

// The critical density in g/cc
double Cosmology::rhoCrit(double z)
{
    double H = Hubble(z);
    return 3.0*H*H/(8.0*M_PI*G);
}

// Give us Mh in solar masses and get back r200 in cm.
double Cosmology::r200(double Mh, double z)
{
    return pow(Mh*3.0/(4.0*M_PI*200.0*rhoCrit(z))  ,1./3.) * pow(MSol, 1./3.);
}

double Cosmology::rhoScaleEinasto(double Mh, double z)
{
    double c200 = cEinasto(Mh, z);
    double alpha = alphaEinasto(Mh, z);
    double rScale =  r200(Mh,z)/c200 ;
    return  Mh/( pow(2.0,2.0-3.0/alpha)*exp(2.0/alpha)*M_PI*pow(alpha,+3./alpha-1.0)*rScale*rScale*rScale*
        (gsl_sf_gamma(3.0/alpha) - gsl_sf_gamma_inc(3.0/alpha, 2.0*pow(c200,alpha)/alpha))/MSol);
//    return 200.0*rhoCrit(z) / (3.0*pow(2.0,-3.0/alpha)*exp(2.0/alpha)*pow(alpha,-1.0+3.0/alpha) /(c200*c200*c200) * 
//            (gsl_sf_gamma(3.0/alpha) - gsl_sf_gamma_inc(3.0/alpha, 2.0/alpha * pow(c200,alpha))));
}

// Mh is the halo mass in solar masses
double Cosmology::nu(double Mh, double z)
{
    double m = log10(Mh * h() / 1.0e12);
    return pow(10.0, -0.11 + 0.146*m + 0.0138*m*m + 0.00123*m*m*m) * (0.033 + 0.79*(1.0+z) + 0.176*exp(-1.356*z));

}

// Mh in solar masses
double Cosmology::alphaEinasto(double Mh, double z)
{
    double nuMz = nu(Mh,z);
    return 0.0095*nuMz*nuMz + 0.155;
}

// Mh in solar masses
double Cosmology::cEinasto(double Mh, double z)
{
    double m = log10(Mh * h() / 1.0e12);
    double b = -0.130 + 0.029*z;
    double a = 0.459 + (0.977 - 0.459)*exp(-0.490 * pow(z, 1.303));
    return pow(10, a+b*m);
}

// r in cm
// Mh in solar masses
double Cosmology::rhoEinasto(unsigned int n)
{
    return rho[n];
//    double c200 = cEinasto(Mh, z);
//    double rScale = r200(Mh,z) / c200;
//    double alpha = alphaEinasto(Mh, z);
//    return rhoScaleEinasto(Mh,z) * exp( -2.0/alpha * (pow(r/rScale, alpha) - 1.0) );

}

// solar masses contained within radius r given in cm.
double Cosmology::MrEinasto(unsigned int n)
{
    return mEinasto[n];
//    double c200 = cEinasto(Mh, z);
//    double rScale = r200(Mh,z) / c200;
//    double alpha = alphaEinasto(Mh, z);
//    return rhoScaleEinasto(Mh, z)* pow(2.0,2.0-3.0/alpha)*exp(2.0/alpha)*M_PI*pow(alpha,+3./alpha-1.0)*rScale*rScale*rScale*
//        (gsl_sf_gamma(3.0/alpha) - gsl_sf_gamma_inc(3.0/alpha, 2.0*pow(r/rScale,alpha)/alpha))/MSol;
}

void Cosmology::UpdateProfile(double Mh, double z, std::vector<double> & x, double Radius)
{
    double c200 = cEinasto(Mh, z);
    double rScale = r200(Mh,z) / c200;
    double alpha = alphaEinasto(Mh, z);
    double rhoScale = rhoScaleEinasto(Mh, z);
    for(unsigned int n=1; n<rho.size(); ++n) {
        double val = rhoScale* pow(2.0,2.0-3.0/alpha)*exp(2.0/alpha)*M_PI*pow(alpha,+3./alpha-1.0)*rScale*rScale*rScale*
            (gsl_sf_gamma(3.0/alpha) - gsl_sf_gamma_inc(3.0/alpha, 2.0*pow(x[n]*Radius/rScale,alpha)/alpha))/MSol;
	mEinasto[n] = val;
	val = rhoScaleEinasto(Mh,z) * exp( -2.0/alpha * (pow(x[n]*Radius/rScale, alpha) - 1.0) );
	rho[n] = val;
    }
}



