class Cosmology;
class DiskContents;
class Dimensions;
struct RafikovQParams;
#include <vector>
#include <string>
#include <gsl/gsl_math.h>

// Methods to compute radial derivatives in the disk.
// Each uses a min-mod slope limiter
double ddx(double,double);
double ddx(std::vector<double>&,unsigned int,std::vector<double>&);
double ddx(double *,unsigned int,std::vector<double>&);

// Print an error message and exit to the system
void errormsg(const std::string);

// The flux from cell n to cell n+1 given the velocity and column density of the material.
// Used to compute dSMigdt (i.e. the change in stellar column density = flux in - flux out.
// (note that this should almost always be negative)
double flux(unsigned int n, std::vector<double>& yy, std::vector<double>& x, std::vector<double>& col_st);

// Compute dS_*/dt owing to stellar migration, at this particular cell,
// given an inward velocity y and a column density profile col_st
double dSMigdt(unsigned int n,std::vector<double>& yy, 
	       std::vector<double>& x, std::vector<double>& col_st);

// Return the age in seconds of the oldest star in age bin i out of N.
// By age, I mean the star's age at redshift zero.
double OldIthBin(unsigned int i,Cosmology&,unsigned int N);

// Return the age in seconds of the youngest star in bin i out of N.
double YoungIthBin(unsigned int,Cosmology&,unsigned int);

// return the maximum of a and b.
double max(double a, double b);

// return I0(x)*exp(-x), where I0 is a modified Bessel function o
// of the first kind
double I0Exp(double x);

// return I1(x)*exp(-x)
double I1Exp(double x);

// Convert various numerical values into strings
std::string str(const double);
std::string str(const int);
std::string str(const unsigned int);
std::string str(const long int);

// Compute Q given a set of parameters qp and a guess at
// the most unstable dimensionless wavenumber absc.
// Upon completion, the actual abscissa will be stored in absc.
// If Q is computed approximately with the Romeo Wiegert (2011)
// prescription, absc will be 1, since that formalism does not
// require computing the particular wavenumber which is unstable.
double Q(RafikovQParams *qp, double *absc);

// Return the value Q - fixedQ. This is the function whose root 
// DiskContents::EnforceFixedQ finds. The argument sv is a factor 
// by which sigma and sigma_* are simultaneously scaled. EnforceFixedQ 
// finds the value of sv which yields Q=fixedQ at each cell.
// fixedQ is stored in DiskContents, and is copied into a RafikovQParams
// object (a pointer to which is passed to QmfQ as p)
double QmfQ(double sv, void * p);

double QmfQfst(double sv, void * p);

// Hide the process of creating the parameters necessary to compute Q,
// and simply return Q at a given cell n in the disk.
double Qsimple(unsigned int n, DiskContents& disk);

// Given a dt in units of the outer orbital period,
// how much will the redshift change?
double dz(double dtSim, double z, Cosmology&, Dimensions& );

// locate the root of F, which will be stored in guess.
int findRoot(gsl_function & F, double * guess);

// find the ~global minimum of fn using fn's derivative, F.
double minFromDeriv(gsl_function & F, 
		    gsl_function & fn, 
		    double * abcissa);

// calculate Q if var (in RafikovQParams) is -1, 
// in which case we ignore sv OR calculate 
// Q(Q_i) or Q(r_i) for the purposes of calculating dQ/dA,
// where A is a stand-in for some state variable.
double varQ(double sv, void * p);

// return the maximum element in the vector
double arrmax(std::vector<double>& arr);

// Compute Q(q), the function which must be minimized
// to find the Rafikov Q.
double Qq(double q, void * p);

// Return the derivative dQ(q)/dq
double dQdq(double q, void * p);
