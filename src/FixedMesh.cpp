#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "FixedMesh.h"
#include "DiskUtils.h"

#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_bessel.h>

double psiIntegrand(double xp, void * params)
{
  FixedMesh * m = (FixedMesh *) params;
  return m->uu(xp) * m->uu(xp)/xp;
}

double sign(double a) {
    if(a>0.0)
        return 1;
    return -1;
}


double HypergeometricPFQ123232152(double y)
{
    //// A numerical approximation to HypergeometricPFQ[{1/2,3/2,3/2},{1,5/2},y^2] (in Mathematica parlance).
}

double HypergeometricPFQ115252223(double y)
{
    //// A numerical approximation to HypergeometricPFQ[{1,1,5/2,5/2},{2,2,3},y^2] (in Mathematica parlance).
}


double summand(double xnm, double xnp, double xj, int term)
{
    //// We assume here that the column density distribution in the sub-region of interest is described by the following:
    //// \Sigma = a x^2 + b x + c + d/x.
    //// The argument term is 0 if we are looking to retrive the terms proportional to a, 1 for b, etc. 
   gsl_mode_t mode = GSL_PREC_DOUBLE; // since we're only evaluating this once, we can probably afford the high precision.
   if(xnp<=xnm) {
       std::cerr << "Bad input in FixedMesh::summand()." << std::endl;
       return 0;
   }
   if(xj > xnp) {
       if( term==0 ) {
        return 2.0/(3.0*M_PI) * (4.0*xj*xj+1.0*xnm*xnm) * gsl_sf_ellint_Ecomp(xnm/xj, mode)
            + 2.0/(3.0*M_PI) * (4.0*xj*xj-1.0*xnp*xnp) * gsl_sf_ellint_Ecomp(xnp/xj, mode)
            - 2.0/(3.0*M_PI) * (4.0*xj*xj -1.0*xnm*xnm ) * gsl_sf_ellint_Kcomp(xnm/xj, mode)
            + 2.0/(3.0*M_PI) * (4.0*xj*xj - 1.0*xnp*xnp) * gsl_sf_ellint_Kcomp(xnp/xj, mode);
       }
       if( term==1 ) {
        return -xnm*xnm*xnm/(3.0*xj*xj) * HypergeometricPFQ123232152(xnm/xj) + xnp*xnp*xnp/(3.0*xj*xj) * HypergeometricPFQ123232152(xnp/xj);
       }
       if( term==2 ) {
        return (2.0/M_PI) * ( gsl_sf_ellint_Ecomp(xnm/xj ,mode) 
                             -gsl_sf_ellint_Ecomp(xnp/xj, mode)
                             -gsl_sf_ellint_Kcomp(xnm/xj, mode) // check this!
                             +gsl_sf_ellint_Kcomp(xnp/xj, mode) );
       }
       if( term==3 ) {
        return -2.0*xnm/(M_PI*xj*xj) * gsl_sf_ellint_Kcomp(xnm/xj, mode)
            + 2.0*xnp/(M_PI*xj*xj) * gsl_sf_ellint_Kcomp(xnp/xj, mode);
       }
   }
   else if(xj < xnm) {
       if( term==0 ) {
        return 2.0*(4.0*xj*xnm+1.0*xnm*xnm*xnm/xj)*gsl_sf_ellint_Ecomp(xj/xnm, mode) / (3.0*M_PI)
            +4.0/(3.0*M_PI) * (2.0*xj*xnp - xnp*xnp*xnp/xj) * gsl_sf_ellint_Ecomp(xj/xnp, mode)
            +2.0/(3.0*M_PI) * (2.0*xj*xnm - xnm*xnm*xnm/xj) * gsl_sf_ellint_Kcomp(xj/xnm, mode)
            +2.0/(3.0*M_PI) * (2.0*xj*xnp - xnp*xnp*xnp/xj) * gsl_sf_ellint_Kcomp(xj/xnp, mode);
       }
       if( term==1) {
        return -9.0*xj*xj*xj/(32.0*xnm*xnm) * HypergeometricPFQ115252223(xj/xnm)
            +9.0*xj*xj*xj/(32.0*xnp*xnp) * HypergeometricPFQ115252223(xj/xnp)
            +0.5*xj*log(xnm) - 0.5*xj*log(xnp);
       }
       if( term==2 ) {
        return (2.0/(M_PI*xj)) * ( xnm*gsl_sf_ellint_Ecomp(xj/xnm, mode)
                                  -xnp*gsl_sf_ellint_Ecomp(xj/xnp, mode) 
                                  -xnm*gsl_sf_ellint_Kcomp(xj/xnm, mode)
                                  +xnp*gsl_sf_ellint_Kcomp(xj/xnp, mode));
       }
       if( term==3 ) {
        return -2.0/(M_PI*xj) * gsl_sf_ellint_Kcomp(xj/xnm, mode)
            + 2.0/(M_PI*xj) * gsl_sf_ellint_Kcomp(xj/xnp, mode);
       }
   }
   else if(xnm<=xj && xj<=xnp) {
       if( term==0 ) {
        return 2.0/(3.0*M_PI) * (4.0*xj*xj+2.0*xnm*xnm)*gsl_sf_ellint_Ecomp(xnm/xj, mode)
            + 2.0/(3.0*M_PI) * (-4.0*xj*xnp - xnp*xnp*xnp/xj) * gsl_sf_ellint_Ecomp(xj/xnp, mode)
            + 2.0/(3.0*M_PI) * (-4.0*xj*xj + 1.0*xnm*xnm) * gsl_sf_ellint_Kcomp(xnm/xj, mode)
            + 2.0/(3.0*M_PI) * (2.0*xj*xnp + 1.0*xnp*xnp*xnp/xj) *gsl_sf_ellint_Kcomp(xj/xnp, mode);
       }
       if( term==1 ) {
        return 3.0*xj/4.0
            - xnm*xnm*xnm/(3.0*xj*xj) * HypergeometricPFQ123232152(xnm/xj) 
            + 9.0*xj*xj*xj/(32.0*xnp*xnp) * HypergeometricPFQ115252223(xj/xnp)
            - 0.25*xj*log(16*xnp*xnp/xj/xj);

       }
       if( term==2 ) {
        return (2.0/(M_PI*xj)) * ( xj*gsl_sf_ellint_Ecomp(xnm/xj, mode)
                                 -xnp*gsl_sf_ellint_Ecomp(xj/xnp, mode)
                                 - xj*gsl_sf_ellint_Kcomp(xnm/xj, mode)
                                 +xnp*gsl_sf_ellint_Kcomp(xj/xnp, mode) );
       }
       if( term==3 ) {
        return -2.0*xnm/(M_PI*xj*xj) * gsl_sf_ellint_Kcomp(xnm/xj, mode)
            + 2.0/(M_PI*xj) * gsl_sf_ellint_Kcomp(xj/xnp, mode);
       }
   }
   else {
       std::cerr << "Strange circumstance in summand()." << std::endl;
       return 0;
   }
}


void test_summand()
{
    std::ofstream f;
    f.open("test_summand.txt");
    std::cout << "Hello from test_summand 0" <<std::endl;
    unsigned int ntest = 1373;
    for(unsigned int i=0; i<ntest; ++i) {
        double x0 = ((double) i+1)/((double) ntest);
        f <<x0<<" "<< summand( 0.01, 0.011, x0, 2) <<" "<< summand(0.1,0.11, x0, 2)<<" "<< summand(.9,.95,x0, 2) << std::endl;
    }
    f.close();


    FixedMesh m(0, .1, 1, .01, .001, 100 );
    std::ofstream f2;
    f2.open("test2_summand.txt");
    for(unsigned int n=1; n<=m.nx(); ++n) {
        double u1=0.0;
        double u2=0.0;
        double u3=0.0;
        double u4=0.0;
        double u5=0.0;
        for(unsigned int nn=1; nn<=m.nx(); ++nn) {
            double weight = summand( m.xPlusHalf(nn-1), m.xPlusHalf(nn), m.x(n), 2 );
            u1=u1 + weight * exp(-m.x(nn)/0.01);
            u2=u2 + weight * exp(-m.x(nn)/0.3);
            u3=u3 + weight * exp(-m.x(nn)/2.0);
            u4=u4 + weight * 1.0;
            u5=u5 + weight * 1/m.x(nn);
        }
        u1=u1 * 2.0*M_PI*.1*m.x(n);
        u2=u2 * 2.0*M_PI*.1*m.x(n);
        u3=u3 * 2.0*M_PI*.1*m.x(n);
        u4=u4 * 2.0*M_PI*.1*m.x(n);
        u5=u5 * 2.0*M_PI*.1*m.x(n);
        f2 << m.x(n) << " " << u1 << " " << u2 << " " << u3 << " " << u4 << " " << u5 << std::endl;
    }
    f2.close();
}



FixedMesh::FixedMesh(double bet0, double turnoverRadius, double nRC, double xm, double mst, unsigned int nnx):
    beta0(bet0), b(turnoverRadius),nxc(nnx), nRotCurve(nRC), dlnxc(-log(xm)/(((double) nnx)-1.)), 
    minsigst(mst), xminc(xm), xv(std::vector<double>(nnx+1,0.)),
    betav(std::vector<double>(nnx+1,0.)), betapv(std::vector<double>(nnx+1,0.)),
    xiPlusHalf(std::vector<double>(nnx+1,0.)),dxi(std::vector<double>(nnx+1,0.)),
    u1pbiPlusHalf(std::vector<double>(nnx+1,0.)),
    uuv(std::vector<double>(nnx+1,0.)),// psiv(std::vector<double>(nnx+1,0.)),
    areas(std::vector<double>(nnx+1,0.)),
    sumTab(std::vector<std::vector<double> >()),
    stored(false),PsiInitialized(false)
{
  psi1 = 0.0;
  bsf=0;
  bmsf=0;

  x_gsl = new double[nxc];
 
  xiPlusHalf[0] = x(0.5);
  u1pbiPlusHalf[0] = uu(x(0.5)) * (1.0 + beta(x(0.5)));
  sumTab.push_back(std::vector<double>(nnx+1,0.));
  for(unsigned int n=1; n<=nxc; ++n) {
    xv[n] = x(((double) n));
    xiPlusHalf[n] = x(((double) n) + 0.5);
    u1pbiPlusHalf[n] = uu(xiPlusHalf[n])*(1.0+beta(xiPlusHalf[n]));
    dxi[n] = xiPlusHalf[n]-xiPlusHalf[n-1];
    areas[n] = M_PI*(xiPlusHalf[n] + xiPlusHalf[n-1])*(xiPlusHalf[n]-xiPlusHalf[n-1]);
//    psiv[n] = psi(xv[n]);
    uuv[n] = uu(xv[n]);
    betav[n] = beta(xv[n]);
    betapv[n] = betap(xv[n]);
    x_gsl[n-1] = xv[n];
    sumTab.push_back(std::vector<double>(nnx+1,0.));
  }


}
void FixedMesh::storeSummand()
{
    test_summand();
    for(unsigned int n=1; n<=nxc; ++n) {
        for(unsigned int nn=1; nn<=nxc; ++nn) {
            sumTab[n][nn] =  summand(xPlusHalf(nn-1), xPlusHalf(nn), x(n), 2);
        }
    }

}

double FixedMesh::summandTabulated(unsigned int n, unsigned int nn)
{
    return sumTab[n][nn];
}

double FixedMesh::u1pbPlusHalf(unsigned int i)
{
  return u1pbiPlusHalf[i];
}
double FixedMesh::area(unsigned int n)
{
    return areas[n];
}
double FixedMesh::dx(unsigned int n)
{
  return dxi[n];
}

bool FixedMesh::InitializePsi()
{
    errormsg("FixedMesh::InitializePsi");
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
//  delete[] x_gsl;
//
//  if(PsiInitialized) {
//    delete[] x_HR_GSL;
//    delete[] psi_HR_GSL;
//
//    gsl_spline_free(spline_psi);
//    gsl_interp_accel_free(accel_psi);
//  }
}
double FixedMesh::x(unsigned int n)
{
  if(n>=1 && n<=nxc) {
    return xv[n];
  }
  return x((double) n); //xminc*exp(dlnxc*(((double) n) - 1.));
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
    errormsg("FixedMesh::psi deprecated");
  if(!PsiInitialized) 
    InitializePsi();
  return gsl_spline_eval(spline_psi, x, accel_psi);
  

  // never get here
  
  double xipsf = 0;//pow(x,ip*soft);
  double xip = 0; //pow(x,ip);
//  if(b!=0.0)
//      return psi1 - 1.0/(2.0*ip) * xip*xip * pow(bsf + xipsf,-2.0/soft) * pow(1.0+bmsf*xipsf,2.0/soft)
//          *gsl_sf_hyperg_2F1(2./soft,2./soft,(2.+soft)/soft,-bmsf*xipsf);
//  else
     return log(x);
}

double FixedMesh::uu(double x)
{
    return pow(1.0+pow(x/b,-fabs(nRotCurve*beta0)),-sign(beta0)/nRotCurve);
//    return 1.0/pow(pow(b/x,nRotCurve*beta0)+1.0,1.0/nRotCurve);
//  return pow(x,ip) / pow(pow(x,ip*soft)+bsf,1./soft);
}

double FixedMesh::beta(double x)
{
  return beta0 - beta0/(1.0+pow(b/x,beta0*nRotCurve));
    
//  return bsf*ip/(bsf+pow(x,ip*soft));
}

double FixedMesh::betap(double x)
{
    return - beta0*beta0*nRotCurve*pow(b/x, beta0*nRotCurve) / (x*pow(1.0+pow(b/x,beta0*nRotCurve),2.0));
//  return -bsf*ip*ip*soft*pow(x,-1.+ip*soft) / ((bsf+pow(x,ip*soft))*(bsf+pow(x,ip*soft)));
}

unsigned int FixedMesh::necessaryN()
{
  errormsg("FixedMesh::necessaryN deprecated");
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
  return 0;
}





