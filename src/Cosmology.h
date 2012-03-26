#include "Dimensions.h"

class Cosmology {
 public:  
  // Define a cosmology given Omega_m, Omega_l (at z=0) the Hubble Constant and a starting redshift.
  Cosmology(double omm,double oml,double h0, double zs)
    :OmM(omm),OmL(oml),H0(h0),zstart(zs),OmK(1.-OmM-OmL){};  
  double lbt(double z); // compute lookback time given a redshift
  double zLBT(double lbTime); // compute redshift given a lookback time (in seconds)
  double Tsim(double z); // simulation time which has elapsed since zstart (seconds)
  double EE(double z); // Cosmological energy density as a function of redshift 
  double ZStart() { return zstart;}
  double OmegaM() {return OmM;}
  double h() { return H0 / (100 * 1.0e5 /(cmperkpc*1.0e3)); }
 private:
  const double OmM, OmL, OmK, H0, zstart;
};
