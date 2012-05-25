#include <math.h>
#include <vector>

#include "FixedMesh.h"
#include "DiskUtils.h"

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_spline.h>

FixedMesh::FixedMesh(double innerPowerlaw, double turnoverRadius, double sft, double xm, unsigned int nnx):
    ip(innerPowerlaw), b(turnoverRadius),nxc(nnx), soft(sft), dlnxc(exp(-log(xm)/(nnx-1.))), xminc(xm),
    xv(std::vector<double>(nnx+1,0.)),
    betav(std::vector<double>(nnx+1,0.)), betapv(std::vector<double>(nnx+1,0.)),
    uuv(std::vector<double>(nnx+1,0.)), psiv(std::vector<double>(nnx+1,0.)),
    stored(false)
{
  bsf = pow(b,sft);
  bmsf = 1.0/bsf;
  psi1=  (1.0 / (2.0*ip)) * (pow(1.0+bmsf, 2.0/sft)*pow(1.0+bsf, -2.0/sft)*
                  gsl_sf_hyperg_2F1(2./sft,2./sft,(2.+sft)/sft,-bmsf));
  
  for(unsigned int n=1; n<=nxc; ++n) {
    xv[n] = x(n);
    psiv[n] = psi(xv[n]);
    uuv[n] = uu(xv[n]);
    betav[n] = beta(xv[n]);
    betapv[n] = betap(xv[n]);
    x_gsl[n-1] = xv[n];
  }

}

FixedMesh::~FixedMesh()
{
  delete[] x_gsl;
}
double FixedMesh::x(unsigned int n)
{
  return xminc*exp(dlnxc*(((double) n) - 1.));
}

double FixedMesh::x(double n)
{
  return xminc*exp(dlnxc*(n-1.));
}

double FixedMesh::psi(double x)
{
  double xipsf = pow(x,ip*soft);
  double xip = pow(x,ip);
  return psi1 - 1.0/(2.0*ip) * xip*xip * pow(bsf + xipsf,-2.0/soft) * pow(1.0+bmsf*xipsf,2.0/soft)
          *gsl_sf_hyperg_2F1(2./soft,2./soft,(2.+soft)/soft,-bmsf*xipsf);
}

double FixedMesh::uu(double x)
{
  return pow(x,ip) / pow(pow(x,ip*soft)+bsf,1./soft);
}

double FixedMesh::beta(double x)
{
  return bsf*ip/(bsf+pow(x,ip*soft));
}

double FixedMesh::betap(double x)
{
  return -bsf*ip*ip*soft*pow(x,-1.+ip*soft) / ((bsf+pow(x,ip*soft))*(bsf+pow(x,ip*soft)));
}

unsigned int FixedMesh::necessaryN(double minsigst)
{
  if(stored) return necesN;

  double psi0 = psi(xminc);
  double u0 = uu(xminc);
  double psi1 = psi0;
  double u1 = u0;
  double x1=xminc;
  for(double mm=1.0; mm<=1000; ++mm) {
    double theMax =0.0;
    for(double n=1; n<=nxc; ++n) {
      for(double m=1.0; m<=mm; ++m) {
        theMax = max ( (psi(x(n+(m/mm))) - psi(x(n+((m-1.0)/mm))) ) / 3.0 
	        + (pow(uu(x(n+(m/mm))),2.0) - pow(uu(x(n+((m-1.0)/mm))),2.0) )/6.0 , theMax);
      }
    }
    if(theMax < .9 * minsigst*minsigst) {
      necesN = (unsigned int) mm;
      stored = true;
      return necesN;
    }
  }
  // hopefully mm=1000 is enough!
  errormsg("The given minsigst "+str(minsigst)+" requires prohibitively high resolution.");
}

