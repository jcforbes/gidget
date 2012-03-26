#include "Cosmology.h"
#include "Dimensions.h"
#include <math.h>
#include <iostream>

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
