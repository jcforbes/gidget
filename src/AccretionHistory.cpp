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
    delete[] acc;
    delete[] redshift;
  }
}
double AccretionHistory::GenerateConstantAccretionHistory(double rate, double zst,Cosmology& cos, std::string fn,bool writeOut)
{
  zstart = zst;
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());
  unsigned int N=1000;
  double z=zstart;
  double MdotExt0=0.0;
  std::vector<double> redshifts(0),tabulatedAcc(0);

  for(unsigned int i=0; i<=N; ++i) {
    z=(((double) (N-i))/((double) N))*(zstart-0.0);
    double MdotExt = rate; //solar masses /year
    if(i==0) MdotExt0 = MdotExt;
    if(writeOut) file << z <<" "<<cos.Tsim(z)<<" "<<MdotExt<<" "<<-1.0<<std::endl;
    redshifts.push_back(z); tabulatedAcc.push_back(MdotExt);
  }
  file.close();
  InitializeGSLObjs(redshifts,tabulatedAcc);
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
  return 64.087*pow((1.0+1.074*pow(x,0.3) - 1.581*pow(x,0.4)+0.954*pow(x,0.5) - 0.185*pow(x,0.6)),-10.0);
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
double MofS_AN(double theS,double sigma8,double OmegaM)
{
  double Gamma=0.169;
  double c0=3.804e-4;
  double uu=u(32.0*Gamma)/sigma8*sqrt(theS);
  double f =pow(uu/64.087,-.1);
  //  double yy= 0.9764842426426905*npow(f-1.0,1./3.)
  //	+0.4678825120300093*npow(f-1.0,2./3.)
  //	+0.3968694013650376*(f-1.0)
  //	+0.3807497793594493*npow(f - 1.,4./3.) 
  //	+0.36862778604069785*npow(f - 1.,5./3.) 
  //	+0.3328886440017323*(f - 1.)*(f-1.) 
  //	+0.24383270141307786*npow(f - 1.,7./3.) 
  //	+0.06282816765413031*npow(f - 1.,8./3.)
  //	-0.26049061964323156*npow(f - 1.,3.0) 
  //	-0.7878462079141473*npow(f - 1.,10./3.) 
  //	-1.5867615682478804*npow(f - 1.,11./3.) 
  //	-2.7159830733318246*npow(f - 1.,4.0);
	
  double yy= 0.9764842426426905*npow(f-1.0,1./3.) + 
    0.4678825120300093*npow(f-1.0,2./3.) + 0.3968694013650376* (f-1.) + 
    0.3807497793594493*npow(f-1.0,4./3.) + 
    0.3686277860406979*npow(f-1.0,5./3.) + 0.3328886440017323*npow(f-1.,2.) + 
    0.2438327014130779*npow(f-1.0,7./3.) + 
    0.06282816765413057*npow(f-1.0,8./3.) - 0.2604906196432312*npow(f-1.0,3.0) - 
    0.7878462079141472*npow(f-1.0,10./3.) - 
    1.5867615682478804*npow(f-1.0,11./3.) - 2.715983073331824*npow(f-1.0,4.0) - 
    4.196113952078296*npow(f-1.0,13./3.) - 
    5.957675479359978*npow(f-1.0,14./3.) - 7.757215883866726*npow(f-1.0,5.0) - 
    9.052022839308824*npow(f-1.0,16./3.) - 
    8.827493940114989*npow(f-1.0,17./3.) - 5.381614672194976*npow(f-1.0,6.0) + 
    3.906669992166177*npow(f-1.0,19./3.) + 
    22.757106826262525*npow(f-1.0,20./3.) + 55.9786623810211*npow(f-1.0,7.0) + 
    109.02755321470875*npow(f-1.0,22./3.) + 
    186.73908051945466*npow(f-1.0,23./3.) + 290.64696116977905*npow(f-1.0,8.0) + 
    414.0907383548234*npow(f-1.0,25./3.) + 
    534.1383992776838*npow(f-1.0,26./3.) + 599.3436845911236*npow(f-1.0,9.0) + 
    512.7459408463739*npow(f-1.0,28./3.) + 
    110.68562129044*npow(f-1.0,29./3.) - 859.470476255553*npow(f-1.0,10.0) - 
    2754.9181995691347*npow(f-1.0,31./3.) - 
    6030.697928001349*npow(f-1.0,32./3.) - 11182.995278112941*npow(f - 1.0,11.) - 
    18597.832691146858*npow(f-1.0,34./3.) - 
    28239.118022803046*npow(f-1.0,35./3.) - 
    39087.162378418005*npow(f - 1,12.) - 
    48222.329292856855*npow(f-1.0,37./3.) - 
    49454.5445168956*npow(f-1.0,38./3.) - 31456.98301353682*npow(f - 1.0,13.0) + 
    24481.226321912305*npow(f-1.0,40./3.) + 
    146641.89257394112*npow(f-1.0,41./3.) + 
    373969.27704908454*npow(f-1.0,14.) + 754053.6571371343*npow(f-1.0,43./3.) + 
    1.334519266949246e6 *npow(f-1.0,44./3.) + 
    2.1423016923339507e6 *npow(f-1.0,15.0) + 
    3.1428330186842247e6 *npow(f-1.0,46./3.) + 
    4.168715749788875e6 *npow(f-1.,47./3.);


  double xx=pow(yy,10.0);
  double M=OmegaM * pow(xx/(c0*Gamma),3.0);
//  std::cout<<"MofS_AN: "<<uu<<" "<<f<<" "<<xx<<" "<<M<<std::endl;
  return M;
}

double AccretionHistory::GenerateNeistein08(double Mh0, double zst, Cosmology& cos, std::string fn, bool writeOut,unsigned long int seed)
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
  do {
    double x= gsl_ran_gaussian(r,1.0);
    double s= log10(SS);
    double deltaS= exp((1.367+0.012*s+0.234*s*s)*x + (-3.682+0.76*s-0.36*s*s));
    zs.push_back(zOfOmega(om));
    masses.push_back(MofS(SS,sp.sigma8,sp.OmegaM)/cos.h());
    std::cout<<"z,om;S,M;ds,s,x: "<<zs[zs.size()-1]<<" "<<om<<"; "<<SS<<" "<<masses[masses.size()-1]<<"; "<<deltaS<<" "<<s<<" "<<x<<std::endl;
    SS+=deltaS;
    om += dom;
  } while(zs[zs.size()-1]<zstart*1.3);
  gsl_rng_free(r);
  for(unsigned int i=0; i!=masses.size()-1; ++i) {
    if((masses[i]-masses[i+1])/masses[i+1] > 0.3 && zs[i+1]<zstart) 
	errormsg("Major merger at z="+str(zs[i+1])+"<zstart; throwing out this accretion history. dM="+str(masses[i]-masses[i+1])+", M="+str(masses[i+1])+", (dM/dz)/M="+str(((masses[i]-masses[i+1])/(zs[i+1]-zs[i]))/masses[i+1]));

    // The accretion rate is the difference in DM halo masses between the two time steps 
    // times the cosmic baryon fraction, divided by the time between redshift steps 
    accs.push_back(.18*(masses[i]-masses[i+1])/
		   (fabs(cos.Tsim(zs[i+1]) - cos.Tsim(zs[i]))
		    /speryear)
		   );
  }
  zs.pop_back();
  reverse(zs.begin(),zs.end());
  reverse(accs.begin(),accs.end());
  if(writeOut) {
    for(unsigned int i=0; i!=accs.size(); ++i) {
      file << zs[i] <<" "<<accs[i]<<std::endl;
    }
    file.close();
  }
  InitializeGSLObjs(zs,accs);

  return accs[accs.size()-1] / AccOfZ(zs[zs.size()-1]);
}
double AccretionHistory::GenerateBoucheEtAl2009(double Mh0, double zs, Cosmology& cos, std::string fn, bool writeOut)
{
  std::ofstream file;
  if(writeOut) file.open(fn.c_str());
  zstart = zs;
  unsigned int N=10000; 
  double Mh=Mh0;
  double z=zstart;
  double fbp18 = 1.0; // baryon fraction / 0.18
  double MdotExt0;
  std::vector<double> redshifts(0),tabulatedAcc(0),haloMass(0);

  for(unsigned int i=0; i<=N; ++i) {
    z=((double) (N-i))/((double) N)*(zstart-0.0);
    double dMh = 34.0 * pow(Mh,1.14)*pow(1.0+z,2.4) * 1.0e-12; // 10^12 Msol/yr
    double fOfz;
    if(z>=2.2)
      fOfz=1.0;
    else if(z<=1.0)
      fOfz=0.5;
    else {
      fOfz = 1.0 - (cos.Tsim(z)-cos.Tsim(2.2)) * 0.5 / (cos.Tsim(1.0) - cos.Tsim(2.2));
    }
    double epsin;
    if(Mh < 1.5)
      epsin = 0.7*fOfz;
    else
      epsin=0.0;

    double MdotExt = 7.0 * epsin * fbp18 * pow(Mh,1.1)*pow(1+z,2.2); // solar masses /year
    if(i==0) MdotExt0= MdotExt;
    haloMass.push_back(Mh);
    Mh+=dMh* -1.0*( cos.Tsim(z) - cos.Tsim( ((double) (N-i-1))/((double) N) * (zstart-0.0))) / speryear;
      
    redshifts.push_back(z); tabulatedAcc.push_back(MdotExt);
    if(writeOut) file << z << " "<< cos.Tsim(z) <<" "<<MdotExt<<" "<<Mh<<std::endl;
  }

  file.close();
  InitializeGSLObjs(redshifts,tabulatedAcc);
  return MdotExt0;
}

// redshifts should range from zstart to 0, tabulatedAcc and redshifts have 
// to have the same size and redshifts must start high and end low
void AccretionHistory::InitializeGSLObjs(std::vector<double> redshifts, std::vector<double> tabulatedAcc)
{
  if(redshifts.size()!=tabulatedAcc.size())
    errormsg("InitializeGSLObjs: tabulatedAcc and redshifts must be the same size.");

  allocated=true;
  accel=gsl_interp_accel_alloc();
  if(!linear) spline=gsl_spline_alloc(gsl_interp_cspline,tabulatedAcc.size());
  else spline=gsl_spline_alloc(gsl_interp_linear,tabulatedAcc.size());
  acc = new double[tabulatedAcc.size()];
  redshift = new double[redshifts.size()];
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
  }

  gsl_spline_init(spline,redshift,acc, tabulatedAcc.size());  
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
  
  InitializeGSLObjs( redshifts,tabulatedAcc);
}

double AccretionHistory::AccOfZ(double z)
{
  if(z > zstart) errormsg("AccOfZ: The accretion rate at this z is unspecified. z and zstart are "
	+str(z)+" "+str(zstart));
  return gsl_spline_eval(spline,z,accel)/gsl_spline_eval(spline,zstart,accel);
}
