#include <math.h>
#include <vector>
#include <iostream>

#include "FixedMesh.h"
#include "DiskUtils.h"

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

double psiIntegrand(double xp, void * params)
{
  FixedMesh * m = (FixedMesh *) params;
  return m->uu(xp) * m->uu(xp)/xp;
}


FixedMesh::FixedMesh(double innerPowerlaw, double turnoverRadius, double sft, double xm, double mst, unsigned int nnx):
    ip(innerPowerlaw), b(turnoverRadius),nxc(nnx), soft(sft), dlnxc(-log(xm)/(((double) nnx)-1.)), 
    minsigst(mst), xminc(xm), xv(std::vector<double>(nnx+1,0.)),
    betav(std::vector<double>(nnx+1,0.)), betapv(std::vector<double>(nnx+1,0.)),
    xiPlusHalf(std::vector<double>(nnx+1,0.)),dxi(std::vector<double>(nnx+1,0.)),
    u1pbiPlusHalf(std::vector<double>(nnx+1,0.)),
    uuv(std::vector<double>(nnx+1,0.)),// psiv(std::vector<double>(nnx+1,0.)),
    stored(false),PsiInitialized(false)
{
  bsf = pow(b,sft);
  bmsf = 1.0/bsf;
  //  if(b!=0.0)   psi1=  (1.0 / (2.0*ip)) * (pow(1.0+bmsf, 2.0/sft)*pow(1.0+bsf, -2.0/sft)*
  //                  gsl_sf_hyperg_2F1(2./sft,2./sft,(2.+sft)/sft,-bmsf));
  //  else 
  psi1 = 0.0;
 
  x_gsl = new double[nxc];
 
  xiPlusHalf[0] = x(0.5);
  u1pbiPlusHalf[0] = uu(x(0.5)) * (1.0 + beta(x(0.5)));
  for(unsigned int n=1; n<=nxc; ++n) {
    xv[n] = x(n);
    xiPlusHalf[n] = x(((double) n) + 0.5);
    u1pbiPlusHalf[n] = uu(xiPlusHalf[n])*(1.0+beta(xiPlusHalf[n]));
    dxi[n] = xiPlusHalf[n]-xiPlusHalf[n-1];
//    psiv[n] = psi(xv[n]);
    uuv[n] = uu(xv[n]);
    betav[n] = beta(xv[n]);
    betapv[n] = betap(xv[n]);
    x_gsl[n-1] = xv[n];
  }


}

double FixedMesh::u1pbPlusHalf(unsigned int i)
{
  return u1pbiPlusHalf[i];
}

double FixedMesh::dx(unsigned int n)
{
  return dxi[n];
}

bool FixedMesh::InitializePsi()
{
//  unsigned int NN = necessaryN();
  unsigned int NN = 1000; // this is a guess...
  double dn = 1.0/(((double) NN)*((double) nxc));

  unsigned int nxcm1 = nxc-1;

  x_HR_GSL = new double[nxcm1*NN];
  psi_HR_GSL = new double[nxcm1*NN];
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc( NN*nxcm1 );

  double result;
  double error;

  gsl_function F;
  F.function = &psiIntegrand;
  F.params = &(*this);

  x_HR_GSL[nxcm1*NN-1] = 1.0;
  psi_HR_GSL[nxcm1*NN-1] = 0.0;

  for(unsigned int i=NN*nxcm1-1; i>=1; --i) {
    double n = ((double) i+NN-1) / ((double) NN);
    double currentX = x(n);
    gsl_integration_qags( &F, 1.0, currentX, 0.0, 1.0e-10, NN * nxcm1, w, &result, &error );
    x_HR_GSL[i-1] = currentX;
    psi_HR_GSL[i-1] = result;
  }
  
  gsl_integration_workspace_free(w);

  accel_psi = gsl_interp_accel_alloc();
  spline_psi = gsl_spline_alloc(gsl_interp_cspline,nxcm1*NN);

  gsl_spline_init(spline_psi,  x_HR_GSL,  psi_HR_GSL, nxcm1*NN);

  PsiInitialized = true;
  return true;

}

FixedMesh::~FixedMesh()
{
  delete[] x_gsl;

  if(PsiInitialized) {
    delete[] x_HR_GSL;
    delete[] psi_HR_GSL;

    gsl_spline_free(spline_psi);
    gsl_interp_accel_free(accel_psi);
  }
}
double FixedMesh::x(unsigned int n)
{
  return xminc*exp(dlnxc*(((double) n) - 1.));
}

double FixedMesh::xPlusHalf(unsigned int n)
{
  return xiPlusHalf[n];
}

double FixedMesh::x(double n)
{
  double val = xminc*exp(dlnxc*(n-1.));
//  if(val<xminc) 
//    return xminc;
//  if(val>1.0) 
//    return 1.0;
  return val;
}

double FixedMesh::psi(double x)
{
  if(!PsiInitialized) 
    InitializePsi();
  return gsl_spline_eval(spline_psi, x, accel_psi);
  

  // never get here
  
  double xipsf = pow(x,ip*soft);
  double xip = pow(x,ip);
  if(b!=0.0)
      return psi1 - 1.0/(2.0*ip) * xip*xip * pow(bsf + xipsf,-2.0/soft) * pow(1.0+bmsf*xipsf,2.0/soft)
          *gsl_sf_hyperg_2F1(2./soft,2./soft,(2.+soft)/soft,-bmsf*xipsf);
  else
     return log(x);
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

unsigned int FixedMesh::necessaryN()
{
  if(stored) return necesN;

  double psi0 = psi(xminc);
  double u0 = uu(xminc);
  double psi1 = psi0;
  double u1 = u0;
  double x1=xminc;
  for(double mm=1.0; mm<=100000; ++mm) {
    double theMax =0.0;
    for(double n=1; n<=nxc-1.0; ++n) {
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





