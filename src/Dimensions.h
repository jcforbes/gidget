#ifndef _Dimensions_h
#define _Dimensions_h

// Some useful dimensional quantities:
const double speryear = 31556926.0;
const double cmperkpc = 3.08568025e21;
const double cmperpc = 3.08568025e18;
const double kmpermpc = 3.08568025e19;
const double Z_Sol = .02;
//const double Z_IGM = .1*Z_Sol; // this is now a user-specified parameter @ runtime
const double Z_BBN = 2.0e-5;
const double MSol = 1.98892e33; // grams
const double G=6.673e-8; // cm^3 s^-2 g^-1
const double mH = 1.67262158e-24; // grams
const double kB = 1.3806503e-16; // erg/K

// Store dimensional quantities for use in converting dimensionless quantities to physical units
class Dimensions {
 public:
  Dimensions(double,double,double);
  // convert dimensionless values -> dimensional values
  double v(double); // (per vphiR) -> velocity (km/s)
  double t(double); // (outer orbits) -> (years)
  double d(double); // (x=r/R) -> kpc

  // dimensionless, calculated from cgs values
  double chi();

  // convert dimensionless values -> cgs values
  double v_cgs(double); // cm/s
  double t_cgs(double); // s
  double d_cgs(double); // cm
  double col_cgs(double);

  // Respectively
  // Outer edge of computational domain,
  // Circular velocity at R, and
  // External accretion rate at z=zstart
  const double Radius, vphiR, MdotExt0; // stored in CGS units
  // These values set the physical scales relevant for the problem
};

#endif
