#ifndef STPOP_H
#define STPOP_H

#include <vector>
#include <string>
class DiskContents;
class Cosmology;

class StellarPop {
 public:
  StellarPop();

  // Create a stellar Population which will form between 
  // lookback times youngest and oldest (in seconds)
  StellarPop(unsigned int nx,double youngest, double oldest);

  // Is this stellar population forming at this redshift?
  bool IsForming(Cosmology&,double reshift);

  friend class DiskContents;

  // Add the contents of sp2 to the calling StellarPop in such a way that
  // Mass, Kinetic Energy, and Mass in Metals are conserved
  void MergeStellarPops(const StellarPop& sp2, DiskContents&);

  // Over a time period dt and given a dimensionless velocity yy inwards, migrate 
  // the stars in such a way that mass, energy, and mass in metals are conserved.
  void MigrateStellarPop(double dt,std::vector<double>& yy, DiskContents&);

  // Set the contents of the current stellar population equal to 
  // some fraction f of the mass in the population sp2.
  // This is used as a numerical trick to speed up the calculation when a new
  // stellar population is formed and therefore Sigma_*/(dSigma_*/dt) is 
  // potentially very small
  void extract(StellarPop& sp2, double f);

 private:
  std::vector<double> spcol; // column density as a function of position.
  std::vector<double> spsig; // stellar velocity dispersion as a function of position.
  double ageAtz0; // i.e. lookback time at creation of these stars, in seconds
  std::vector<double> spZ; // metallicity of the stars as a function of position.
  std::vector<double> spZV; // metallicity variance
  std::vector<double> dQdS; // The partial derivative of Q wrt this population's S_*
  std::vector<double> dQds; // The partial derivative of Q wrt this population's s_*
  std::vector<double> dQdSerr; // the error in dQ/dS_*
  std::vector<double> dQdserr; // the error in dQ/ds_*
  double youngest,oldest; // stored in seconds, boundaries on the ages of stars in this population
};

#endif /* STPOP_H */
