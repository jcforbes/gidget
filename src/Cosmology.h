#include "Dimensions.h"
#include <vector> 

class Cosmology {
 public:  
  // Define a cosmology given Omega_m, Omega_l (at z=0) the Hubble Constant and a starting redshift.
  Cosmology(double omm,double oml,double h0, double sigma8, double zs, unsigned int nx, double crf)
    :OmM(omm),OmL(oml),H0(h0),zstart(zs),OmK(1.-OmM-OmL),s8(sigma8),
     rho(std::vector<double>(nx+1,0.)), mEinasto(std::vector<double>(nx+1,0.)),concentrationRandomFactor(crf){};  
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
  double crf() const {return concentrationRandomFactor;}
  double zs() const {return zstart;} 
  std::vector<double> GetRho() const {return rho;}
  double sigma8() const { return s8;} 
  double h() const { return H0 / (100 * 1.0e5 /(cmperkpc*1.0e3)); } 
  double Hubble(double z);
  // The critical density in g/cc
  double rhoCrit(double z);
  // Give us Mh in solar masses and get back r200 in cm.
  double r200(double Mh, double z);
  double rhoScaleEinasto(double Mh, double z);
  // Mh is the halo mass in solar masses
  double nu(double Mh, double z);
  // Mh in solar masses
  double alphaEinasto(double Mh, double z);
  // Mh in solar masses
  double cEinasto(double Mh, double z);
  // r in cm
  // Mh in solar masses
  double rhoEinasto(unsigned int n);
  double MrEinasto(unsigned int n);
  void UpdateProfile(double Mh, double z, std::vector<double> & x, double Radius);

 private:
  const double OmM, OmL, OmK, H0, zstart, s8;
  const double concentrationRandomFactor; // A factor which allows a constant multiplicative offset in the concentration of a galaxy.
  std::vector<double> rho;
  std::vector<double> mEinasto;
};
