#include "Dimensions.h"

class Cosmology {
 public:  
  // Define a cosmology given Omega_m, Omega_l (at z=0) the Hubble Constant and a starting redshift.
  Cosmology(double omm,double oml,double h0, double sigma8, double zs)
    :OmM(omm),OmL(oml),H0(h0),zstart(zs),OmK(1.-OmM-OmL),s8(sigma8) {};  
  Cosmology operator=(const Cosmology&);
  double lbt(double z); // compute lookback time given a redshift
  double zLBT(double lbTime); // compute redshift given a lookback time (in seconds)
  double Tsim(double z); // simulation time which has elapsed since zstart (seconds)
  double EE(double z); // Cosmological energy density as a function of redshift 
  double ZStart() const { return zstart;} 
  double OmegaM() const {return OmM;} 
  double OmegaL() const {return OmL;} 
  double OmegaK() const {return OmK;} 
  double Hubble() const {return H0;} 
  double zs() const {return zstart;} 
  double sigma8() const { return s8;} 
  double h() const { return H0 / (100 * 1.0e5 /(cmperkpc*1.0e3)); } 
 private:
  const double OmM, OmL, OmK, H0, zstart, s8;
};
