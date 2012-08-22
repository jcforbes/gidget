#include "AccretionHistory.h"
#include "Cosmology.h"
#include "Dimensions.h"
#include "DiskUtils.h"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>

#include <fstream>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>

AccretionHistory::~AccretionHistory()
{
  if(allocated) {
    gsl_spline_free(spline);
    gsl_interp_accel_free(accel);
    gsl_spline_free(splineMh);
    gsl_interp_accel_free(accelMh);
    delete[] acc;
    delete[] redshift;
    delete[] hMass;
  }
}

double AccretionHistory::GenerateOscillatingAccretionHistory(double amp, double period,
							     double phase,		     
					double zst, bool inRedshift, // as opposed to in Gyr.
					Cosmology& cosm, std::string fn, bool writeOut)
{
  zstart = zst;
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());

  double z=zstart;
  double MdotExt0=0.0; 
  double mass = 1.0e12 - (amp/2.0) * cosm.Tsim(0.0) / speryear;
  unsigned int N=10000;

  std::vector<double> redshifts(0), tabulatedAcc(0), masses(0);

  double eps = 1.0e-8; // safety factor to make sure we never get Mdot<0

  for(unsigned int i=0; i<=N; ++i) {
    z=(((double) (N-i))/((double) N)) * (zstart - 0.0);
    double MdotExt=0;
    double t=cosm.Tsim(z);
    if(inRedshift) 
      MdotExt = (amp/(2.0-eps)) + 
	(amp/2.0)*cos(2.0*M_PI*(z-zstart)/period-phase);
    else           
      MdotExt = (amp/(2.0-eps)) + 
	(amp/2.0)*cos(2.0*M_PI*(t/speryear)*1.0e-6/period-phase);
    if(i==0) MdotExt0 = MdotExt;

    if(writeOut) file << z <<" "<<t<<" "<<MdotExt<<" "<<-1.0<<std::endl;
    redshifts.push_back(z); tabulatedAcc.push_back(MdotExt); masses.push_back(mass);
    mass += MdotExt * fabs(cosm.Tsim(redshifts[N-i+1]) - cosm.Tsim(z)) / speryear;
  }
  file.close();
  InitializeGSLObjs(redshifts,tabulatedAcc,masses);
  return MdotExt0;
}


double AccretionHistory::GenerateConstantAccretionHistory(double rate, double zst,
					Cosmology& cos, std::string fn,bool writeOut)
{
  zstart = zst;
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());
  unsigned int N=1000;
  double z=zstart;
  double MdotExt0=0.0;
  double mass=1.0e12 - rate * cos.Tsim(0.0) / speryear;
  std::vector<double> redshifts(0),tabulatedAcc(0),masses(0);

  for(unsigned int i=0; i<=N; ++i) {
    z=(((double) (N-i))/((double) N))*(zstart-0.0);
    double MdotExt = rate; //solar masses /year
    if(i==0) MdotExt0 = MdotExt;
    if(writeOut) file << z <<" "<<cos.Tsim(z)<<" "<<MdotExt<<" "<<-1.0<<std::endl;
    redshifts.push_back(z); tabulatedAcc.push_back(MdotExt); masses.push_back(mass);
    mass += rate* fabs(cos.Tsim(redshifts[N-i+1]) - cos.Tsim(z)) / speryear;
  }
  file.close();
  InitializeGSLObjs(redshifts,tabulatedAcc,masses);
  return MdotExt0;
}

double omega(double z,void * rhs)
{
  double * theRHS = (double *) rhs;
  return 1.260*(1.0+z+0.09/(1.0+z) + 0.24*exp(-1.16*z))-(*theRHS);
}

double zOfOmega(double om)
{
  gsl_function F;
  F.function = &omega;
  double tempOm1=om;
  double z=om;
  F.params=&tempOm1;
  findRoot(F,&z);
  return z;

}
double u(double x)
{
  return 64.087*pow((1.0+1.074*pow(x,0.3) - 1.581*pow(x,0.4)+
			0.954*pow(x,0.5) - 0.185*pow(x,0.6)),-10.0);
}

struct SParams {
  double rhs,OmegaM,sigma8;
};

double S(double M, void * p)
{
  SParams * sp = (SParams *) p;
  double c0=3.804e-4;
  double Gamma=0.169; // power spectrum shape parameter
  double sigma8= (*sp).sigma8;
  double x=c0*Gamma*pow(M,1./3.)/pow((*sp).OmegaM,1./3.);
  return u(x)*u(x)*sigma8*sigma8/(u(32.0*Gamma)*u(32.0*Gamma)) - (*sp).rhs;
}

// Return M(S) in units of solar masses
double MofS(double theS,double s8,double OmegaM)
{
  gsl_function F;
  F.function = &S;
  SParams sp;
  sp.rhs=theS;
  sp.sigma8=s8;
  sp.OmegaM=OmegaM;
  F.params=&sp;
  double M=1.0e12;
  findRoot(F,&M);
  return M;
}
double npow(double x,double ex) 
{
  if(x<0.0) {
    return -pow(-x,ex);
  }
  else
    return pow(x,ex);
}

double AccretionHistory::GenerateNeistein08(double zst, Cosmology& cos, 
				std::string fn, bool writeOut,unsigned long int seed,
				double invMassRatioLimit,bool fatal)
{
  zstart = zst;
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());
  double sigma8=0.82;
  double OmegaM=cos.OmegaM();
  double dom = 0.1;
  double zero=0.0;
  double om = omega(0.0,&zero);
  linear=true; // use linear interpolation
  SParams sp;
  sp.rhs=0.0;
  sp.sigma8=sigma8;
  sp.OmegaM=OmegaM;
  std::vector<double> zs(0),accs(0),masses(0);
  double SS= S(Mh0*cos.h(),&sp);
  gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r,seed);
  bool first=true;
  do {
    // pick a value for our Gaussian random variable.
    double x= gsl_ran_gaussian(r,1.0);
    double s= log10(SS);
    // employ the formula in ND08 to compute the delta S for the chosen value of x.
    double deltaS= exp((1.367+0.012*s+0.234*s*s)*x + (-3.682+0.76*s-0.36*s*s));
    // update our list of redshifts, given our uniformly spaced value of omega.
    // subtract a small number so that when we evaluate some quantity at z=0, we don't get an interpolation error.

    double z=zOfOmega(om)-1.0e-4; // avoids interpolation errors at the end of the simulation...
    zs.push_back(z);
    if(!first) zs.push_back(z+.00000001);
    
    // update our tabulated value of halo mass.
    double M=MofS(SS,sp.sigma8,sp.OmegaM)/cos.h();
    masses.push_back(M);
    if(!first)  masses.push_back(M); // solar masses
    std::cout<<"z,om;S,M;ds,s,x: "<<zs[zs.size()-1]<<" "<<om<<"; "<<SS<<" "<<
	masses[masses.size()-1]<<"; "<<deltaS<<" "<<s<<" "<<x<<std::endl;

    SS+=deltaS;
    om += dom;
    first=false;
  } while(zs[zs.size()-1]<zstart*1.3);
  gsl_rng_free(r);
//  for(unsigned int i=0; i!=masses.size()-1; ++i) {
  for(unsigned int i=0; i<masses.size()-1; i+=2) {
    if((masses[i]-masses[i+1])/masses[i+1] > invMassRatioLimit && zs[i+1]<zstart) {
	errormsg("Major merger at z="+str(zs[i+1])+"<zstart; throwing out this accretion history. dM="+str(masses[i]-masses[i+1])+", M="+str(masses[i+1])+", (dM/dz)/M="+str(((masses[i]-masses[i+1])/(zs[i+1]-zs[i]))/masses[i+1]),fatal);
        return -1.0;
    }

    // The accretion rate is the difference in DM halo masses between the two time steps 
    // times the cosmic baryon fraction, divided by the time between redshift steps 
    double accr = (.18*epsin(zs[i],masses[i]*1.0e-12,cos)*(masses[i]-masses[i+1])/
		   (fabs(cos.Tsim(zs[i+1]) - cos.Tsim(zs[i]))
		    /speryear)
		   );
    // if the "first" flag above is working correctly, these two accretion rates should
    // occur at very different redshifts.
    accs.push_back(accr);
    accs.push_back(accr);
  }
  zs.pop_back();
  masses.pop_back();


  reverse(zs.begin(),zs.end());
  reverse(accs.begin(),accs.end());
  reverse(masses.begin(),masses.end());
  if(writeOut) {
    for(unsigned int i=0; i!=accs.size(); ++i) {
      file << zs[i] <<" "<<accs[i]<<std::endl;
    }
    file.close();
  }
  InitializeGSLObjs(zs,accs,masses);

  return accs[accs.size()-1] / AccOfZ(zs[zs.size()-1]);
}

// Mh in units of 10^12 solar masses
double AccretionHistory::epsin(double z, double Mh,Cosmology & cos)
{
    double fOfz;
    if(z>=2.2)
        fOfz=1.0;
    else if(z<=1.0)
        fOfz=0.5;
    else {
        fOfz = 1.0 - (cos.Tsim(z)-cos.Tsim(2.2))* 0.5 / (cos.Tsim(1.0) - cos.Tsim(2.2));
    }
    double eps;
    if(Mh < 5.0) 
        eps = 0.7*fOfz;
    else
        eps = 0.0;
    return eps;
}

double AccretionHistory::GenerateBoucheEtAl2009( double zs, Cosmology& cos, 
					std::string fn, bool writeOut, bool MhAtz0)
{
  double Mh012 = Mh0*1.0e-12;
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());
  zstart = zs;
  unsigned int N=1000; 
  double Mh=Mh012; // units of 10^12 solar masses
  double z=zstart;
  double fbp18 = 1.0; // baryon fraction / 0.18
  double MdotExt0;
  std::vector<double> redshifts(0),tabulatedAcc(0),haloMass(0);

  // Loop over redshift from z=0 to z=zstart in N increments if Mh is specified at z=0
  // Otherwise, loop over redshift from zstart to z=0.
  for(unsigned int i=0; i<=N; ++i) {
    if(!MhAtz0) // if Mh is given at z=zstart, start from high redshift and go to z=0 
	z=((double) (N-i))/((double) N)*(zstart-0.0);
    else // if Mh is given at z=0, start from low redshift and go to z=zstart.
	z=((double) i)/((double) N) * (zstart - 0.0);

    // Use the Bouche formula to tell us the dark matter accretion rate
    double dMh = 34.0 * pow(Mh,1.14)*pow(1.0+z,2.4) * 1.0e-12; // in 10^12 Msol/yr
    
    // Use the analogous formula to tell us the baryonic accretion rate.
    double MdotExt = 7.0 * epsin(z,Mh,cos) * fbp18 * pow(Mh,1.1)*pow(1+z,2.2); // in solar masses /year

    if((i==0 && !MhAtz0) || (i==N && MhAtz0)) // always set MdotExt0 to be MdotExt at z=2
	MdotExt0= MdotExt;

    haloMass.push_back(Mh); // solar masses

    // Basically compute Mh(z) by taking an Euler step, since from the above we know dMh (which is actually dMh/dt)
    if(!MhAtz0) // starting from high redshift..
        Mh+= dMh* -1.0*( cos.Tsim(z) - cos.Tsim( ((double) (N-i-1))/((double) N) * (zstart-0.0))) / speryear;
    else // starting from low redshift..
        Mh+= dMh* -1.0*( cos.Tsim(z) - cos.Tsim( ((double) (i+1))/((double) N) * (zstart-0.0))) / speryear;

    redshifts.push_back(z); tabulatedAcc.push_back(MdotExt);
    if(writeOut) file << z << " "<< cos.Tsim(z) <<" "<<MdotExt<<" "<<Mh<<std::endl;
  }

  // 
  if(MhAtz0) { // then we need to reverse redshifts and tabulatedAcc
    reverse(redshifts.begin(),redshifts.end());
    reverse(tabulatedAcc.begin(), tabulatedAcc.end());
    reverse(haloMass.begin(),haloMass.end());
  }

  file.close();
  InitializeGSLObjs(redshifts,tabulatedAcc,haloMass);
  return MdotExt0;
}

// redshifts should range from zstart to 0, tabulatedAcc, haloMass, and redshifts have 
// to have the same size; redshifts must start high and end low
void AccretionHistory::InitializeGSLObjs(std::vector<double> redshifts, std::vector<double> tabulatedAcc,
					 std::vector<double> haloMass)
{
  if(redshifts.size()!=tabulatedAcc.size() || redshifts.size()!=haloMass.size())
    errormsg("InitializeGSLObjs: tabulatedAcc and redshifts must be the same size.");

  allocated=true;
  accel=gsl_interp_accel_alloc();
  accelMh=gsl_interp_accel_alloc();
  if(!linear) {
    spline=gsl_spline_alloc(gsl_interp_cspline,tabulatedAcc.size());
    splineMh=gsl_spline_alloc(gsl_interp_cspline,haloMass.size());
  }
  else {
    spline=gsl_spline_alloc(gsl_interp_linear,tabulatedAcc.size());
    splineMh = gsl_spline_alloc(gsl_interp_linear,haloMass.size());
  }
  acc = new double[tabulatedAcc.size()];
  redshift = new double[redshifts.size()];
  hMass = new double[haloMass.size()];
  unsigned int tabAcc=tabulatedAcc.size();
//  double norm = tabulatedAcc[0];
//  // Normalize accretion rate to the highest-redshift value given.
//  for(unsigned int i=0; i<tabAcc; ++i) {
//    tabulatedAcc[i]/=norm;
//  }

  // reverse order!
  for(unsigned int i=0; i<tabAcc; ++i) {
    acc[tabAcc-i-1]=tabulatedAcc[i];
    redshift[tabAcc-i-1]=redshifts[i];
    hMass[tabAcc-i-1]=haloMass[i];
  }

  gsl_spline_init(spline,redshift,acc, tabulatedAcc.size());  
  gsl_spline_init(splineMh,redshift,hMass, haloMass.size());
}


// fn=filename
// zc= redshift column
// acC= accretion column
// rs = rows skipped before data is read
// nc = total number of columns
void AccretionHistory::ReadTabulated(std::string fn,unsigned int zc, unsigned int acC, 
                                     unsigned int rs, unsigned int nc, double zst)
{
  zstart = zst;
  std::vector<double> redshifts(0);
  std::vector<double> tabulatedAcc(0);
  std::vector<double> haloMass(0);
  std::ifstream f(fn.c_str());
  if(!f.is_open()) errormsg("Error opening file containing the tabulated accretion history!");
  std::string line;
  for(unsigned int i=0; i<rs; ++i) {
    getline(f,line); // get a line and discard it
  }
  double nm;
  while(f.good() ) {
    for(unsigned int i=0; i<nc; ++i) {  // for each column..
      f >> nm; // read the number
      if(i+1 == zc) { // if we're in the redshift column, ..
        if(redshifts.size()==0 && zstart >= zc) redshifts.push_back(zstart);
        else if(redshifts.size()==0) redshifts.push_back(nm*1.0001);
        redshifts.push_back(nm); // 
      }
      if(i+1 == acC) {
        tabulatedAcc.push_back(nm);
        if(tabulatedAcc.size()==1) tabulatedAcc.push_back(nm);
      }
    }
  }
  tabulatedAcc.pop_back();
  redshifts.pop_back();
  if(redshifts[redshifts.size()-1] != 0.0) {
    tabulatedAcc.push_back(tabulatedAcc[tabulatedAcc.size()-1]);
    redshifts.push_back(0.);
  }

  f.close();
  
  InitializeGSLObjs( redshifts,tabulatedAcc, haloMass);
}

double AccretionHistory::AccOfZ(double z)
{
  if(z > zstart) errormsg("AccOfZ: The accretion rate at this z is unspecified. z and zstart are "
	+str(z)+" "+str(zstart));
  return gsl_spline_eval(spline,z,accel)/gsl_spline_eval(spline,zstart,accel);
}

double AccretionHistory::MhOfZ(double z)
{
  if(z>zstart) errormsg("MhOfZ: The halo mass at this z is unspecified.");

  return gsl_spline_eval(splineMh,z,accelMh)/gsl_spline_eval(splineMh,0.0,accelMh);
}


void testAccretionHistory()
{
  Cosmology cos(1.0-.734, .734, 2.29e-18, 2.0);


  double dummy=0.0;
//  for(unsigned int whichAccretionHistory=10; whichAccretionHistory!=1000; ++whichAccretionHistory) {
//    AccretionHistory accr(1.0e12);
//    double Mdot0 = accr.GenerateNeistein08(2.0, cos, "", false, whichAccretionHistory, .3, false);
//    if(Mdot0>0.0) {
//      double dummy = 0.0;
//
//
//    }
//  }




}

















