#include <gsl/gsl_spline.h>
#include <string>
#include <vector>
class Cosmology;

class AccretionHistory {
 public:
  AccretionHistory():allocated(false),linear(false){};
  ~AccretionHistory();  
  // Compute the normalized external accretion rate mdot(z)/mdot(zstart)
  // This is implemented by doing an interpolation of a given 
  // tabulated accretion history
  double AccOfZ(double z);
  // Read a tabulated accretion history
  void ReadTabulated(std::string filename, 
		     unsigned int zc, // column with redshifts
		     unsigned int acC, // column with accretion histories
		     unsigned int rs, // number of rows to skip
		     unsigned int nc, // number of columns in the file
		     double zstart); // starting redshift
  // Generate a realistic accretion history given a halo mass at z=0
  double GenerateNeistein08(double Mh0, double zstart, 
			    Cosmology& cos, std::string fn, 
			    bool writeOut,unsigned long int seed);

  // Generate an accretion history given a halo mass at zstart.
  double GenerateBoucheEtAl2009(double Mh0, double zstart, 
				Cosmology& cos, std::string fn, 
				bool writeOut);
  // Generate an accretion history which is just constant with redshift
  double GenerateConstantAccretionHistory(double rate,double zstart,Cosmology& cos,std::string fn,bool writeOut);

  // Given a vector of redshifts and the accretion rates at those
  // redshifts, fill in the GSL interpolation objects which will
  // be used to return the accretion rate at arbitrary redshifts.
  void InitializeGSLObjs(std::vector<double> z, 
			 std::vector<double> tabulatedAcc);
 private:
  // vector of tabulated dimensionless accretion - same 
  // number of elements as redshift:
  double *acc;
  
  // vector of redshifts at which acc rate has been tabulated
  double * redshift;

  // the earliest time the accretion histories know about = time the simulation will start
  double zstart;

  // objects used in the interpolation between the tabulated values
  gsl_interp_accel *accel;
  gsl_spline *spline;
 
  // a switch to determine whether we're using a cubic spline or a linear
  // interpolation. If linear is true, we're using the linear interpolation.
  bool linear;
 
  // true if the interpolator is ready to roll, i.e. if we
  // can call AccOfZ successfully
  bool allocated;
};
