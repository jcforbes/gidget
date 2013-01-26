#ifndef STPOP_H
#define STPOP_H

#include <vector>
#include <string>
class DiskContents;
class Cosmology;
class FixedMesh;

double ComputeVariance(double,double,double,double,double,double,double);

class StellarPop {
 public:
//  StellarPop();
//  StellarPop(const StellarPop&);
//  StellarPop & operator=(const StellarPop&);

  // Create a stellar Population which will form between 
  // lookback times youngest and oldest (in seconds)
  StellarPop(FixedMesh & mesh);

  // Is this stellar population forming at this redshift?
//  bool IsForming() const { return isForming; }; // It's probably a bad idea to do things this way.

  friend class DiskContents;

  // Add the contents of sp2 to the calling StellarPop in such a way that
  // Mass, Kinetic Energy, and Mass in Metals are conserved
  void MergeStellarPops(const StellarPop& sp2, DiskContents&);

  // Over a time period dt and given a dimensionless velocity yy inwards, migrate 
  // the stars in such a way that mass, energy, and mass in metals are conserved.
  void MigrateStellarPop(double dt,double ** tauvecStar, DiskContents&, std::vector<double>& MdotiPlusHalf);

  // Set the contents of the current stellar population equal to 
  // some fraction f of the mass in the population sp2.
  // This is used as a numerical trick to speed up the calculation when a new
  // stellar population is formed and therefore Sigma_*/(dSigma_*/dt) is 
  // potentially very small
  void extract(StellarPop& sp2, double f);

  void ComputeRecycling(DiskContents& , double z, double RfREC);

  void ComputeSpatialDerivs();

  std::vector<double> GetSpCol() const { return spcol; };
  std::vector<double> GetSpSigR() const { return spsigR; };
  std::vector<double> GetSpSigZ() const { return spsigZ; };
  std::vector<double> GetdSigRdr() const { return dSigRdr; };
  std::vector<double> GetdSigZdr() const { return dSigZdr; };
  std::vector<double> GetdColdr() const {return dColdr; };
  std::vector<double> GetSpZ() const { return spZ; };
  std::vector<double> GetSpZV() const { return spZV; };
  std::vector<double> GetdQdS() const { return dQdS; };
  std::vector<double> GetdQdsR() const { return dQdsR; };
  std::vector<double> GetdQdsZ() const { return dQdsZ; };
  std::vector<double> GetdQdSerr() const { return dQdSerr;};
  std::vector<double> GetdQdserr() const { return dQdserr;};
  FixedMesh & GetMesh() { return mesh; };
//  double GetYoungest() const { return youngest; };
//  double GetOldest() const { return oldest; };
  double GetAgeAtz0() const { return ageAtz0; };

 private:
  std::vector<double> spcol; // column density as a function of position.
  std::vector<double> spsigR; // R-direction stellar velocity dispersion as a function of position.
  std::vector<double> spsigZ; // phi- (and z-) direction stellar velocity dispersion.
  std::vector<double> dcoldtREC;
  double ageAtz0; // i.e. lookback time at creation of these stars, in seconds
  double startingAge, endingAge; // lookback time in seconds when stars started (stopped) being added to this pop
  std::vector<double> spZ; // metallicity of the stars as a function of position.
  std::vector<double> spZV; // metallicity variance
  std::vector<double> dQdS; // The partial derivative of Q wrt this population's S_*
  std::vector<double> dQdsR; // The partial derivative of Q wrt this population's s_*
  std::vector<double> dQdsZ;
  std::vector<double> dSigRdr;
  std::vector<double> dSigZdr;
  std::vector<double> dColdr;
  std::vector<double> dQdSerr; // the error in dQ/dS_*
  std::vector<double> dQdserr; // the error in dQ/ds_*
  //double youngest,oldest; // stored in seconds, boundaries on the ages of stars in this population
  FixedMesh & mesh;
};

#endif /* STPOP_H */
