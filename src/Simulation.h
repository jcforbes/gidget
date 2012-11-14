class Cosmology;
class DiskContents;
class AccretionHistory;
class Debug;
class Dimensions;
#include <string.h>

// A small container to store the state variables
// and other numbers necessary to initialize a disk.
struct Initializer {
  std::vector<double> col,sig,col_st,sig_stR,sig_stZ,Z;
  unsigned int NActive,NPassive;
};

class Simulation {
 public:
  Simulation(const double tmax, const long int stepmax, 
	     const bool cosmologyOn, const unsigned int nnx,
	     const double TOL, const double zs,
	     const unsigned int na, const unsigned int np,
	     const double alphaMRI, const double sigth, const double ndecay,
	     DiskContents&,
	     AccretionHistory&,
	     Debug&, Dimensions&);
  int runToConvergence(const double fCondition,
		       const bool writeOut,
		       const std::string filename,
                       const double zrelax,
		       const unsigned int Noutputs);
  Initializer& GetInitializer() { return ini; }
 private:
  const double tmax;
  const unsigned int stepmax;
  const bool cosmologyOn;
  const unsigned int nx;
  const double TOL;
  const double zstart;
  const unsigned int NPassive,NActive;
  const double alphaMRI,sigth;
  const double ndecay;
  DiskContents& theDisk;
  AccretionHistory& accr;
  Debug& dbg;
  Dimensions& dim;
  Initializer ini;
  
};
