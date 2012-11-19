#include <gsl/gsl_spline.h>
#include <string>
#include <vector>
#include <gsl/gsl_randist.h>

class Cosmology;
class Debug;


class AccretionHistory {
 public:
  AccretionHistory(double Mh,Debug & db):Mh0(Mh),allocated(false),linear(false),dbg(db){};
  ~AccretionHistory();  
  // Compute the normalized external accretion rate mdot(z)/mdot(zstart)
  // This is implemented by doing an interpolation of a given 
  // tabulated accretion history
  double AccOfZ(double z);

  double GetMh0() { return Mh0; }  

  // Normalized to Mh0.
  double MhOfZ(double z);

  // Read a tabulated accretion history
  void ReadTabulated(std::string filename, 
		     unsigned int zc, // column with redshifts
		     unsigned int acC, // column with accretion histories
		     unsigned int rs, // number of rows to skip
		     unsigned int nc, // number of columns in the file
		     double zstart); // starting redshift

  double GenerateNeistein08(double zstart, Cosmology & cos,
			    std::string fn, bool writeOut, unsigned long int seed,
			    double invMassRatioLimit, double zquench, int * attempts,
                            double domega, double zrelax);
  
  // Generate a realistic accretion history given a halo mass at z=0
  double AttemptToGenerateNeistein08(double zstart, 
			    Cosmology& cos, std::string fn, 
			    bool writeOut,gsl_rng * r,
			    double invMassRatioLimit, double zquench,
                            double domega, double zrelax);

  // Generate an accretion history given a halo mass Mh0 at zstart.
  // If MhAtz0 is true, the given halo mass Mh0 is assumed to be for z=0, not z=zstart
  // fn is the filename to which to write this accretion history if writeOut is true
  // Cosmology contains the assumed Cosmological parameters & functions to calculate
  // lookback time for a given redshift.
  double GenerateBoucheEtAl2009( double zstart, 
				Cosmology& cos, std::string fn, 
				bool writeOut,bool MhAtz0,double zquench);

  // Generate an accretion history which is just constant with redshift
  double GenerateConstantAccretionHistory(double rate,double zstart,Cosmology& cos,
				std::string fn,bool writeOut);

  double GenerateOscillatingAccretionHistory(
                double amp, double period, // parameters for the oscillation
		double phase, // angle by which to move the cosine to the right.
		double zst,  // starting redshift
		bool inRedshift, // the period is in redshift (otherwise the period is in time)
		Cosmology& cos, std::string fn, bool writeOut);



  // Given a vector of redshifts and the accretion rates at those
  // redshifts, fill in the GSL interpolation objects which will
  // be used to return the accretion rate at arbitrary redshifts.
  void InitializeGSLObjs(std::vector<double> z, 
			 std::vector<double> tabulatedAcc,
			 std::vector<double> haloMass);

  double epsin(double z, double Mh, Cosmology & cos,double zquench);
  void SetEfficiencyParams(double norm, double alpha_z, double alpha_Mh, double ceiling);

 private:
  // vector of tabulated dimensionless accretion - same 
  // number of elements as redshift:
  double *acc;
  
  // vector of redshifts at which acc rate has been tabulated
  double * redshift;

  double * hMass;

  // the earliest time the accretion histories know about = time the simulation will start
  double zstart;

  // objects used in the interpolation between the tabulated values
  gsl_interp_accel *accel, *accelMh;
  gsl_spline *spline, *splineMh;
 
  // a switch to determine whether we're using a cubic spline or a linear
  // interpolation. If linear is true, we're using the linear interpolation.
  bool linear;
 
  // true if the interpolator is ready to roll, i.e. if we
  // can call AccOfZ successfully
  bool allocated;

  // halo mass at redshift zero (units of solar masses)
  double Mh0;

  Debug & dbg;

  double normalization, alpha_z, alpha_Mh, ceiling;

};


void testAccretionHistory();
