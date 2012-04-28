#include "DiskUtils.h"
#include "RafikovQParams.h"
#include "Cosmology.h"
#include "DiskContents.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <vector>
#include <iostream>
#include <sstream>

double ddx(double left, double right)
{
  if(left*right <=0.)
    return 0.0;
  if(fabs(left)>fabs(right))
    return right;
  return left;
}

double ddx(double* arr, unsigned int n, std::vector<double>& x)
{
/*   if(n<nx && n>1) */
/*     return (arr[n+1]-arr[n-1])/(x[n+1]-x[n-1]); */
  unsigned int nx = x.size()-1;
  if(n<nx && n>1) {
    // slope limiter!
    double right = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double left = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    return ddx(left,right);
  }
  if(n==nx) {
    double right = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    double left = (arr[n-1]-arr[n-2])/(x[n-1]-x[n-2]);
    double center = (arr[n]-arr[n-2])/(x[n]-x[n-2]);
//    return ddx(center,right); 
    return ddx(left,right);
    return right;
  }
  if(n==1) {
    double left = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double right = (arr[n+2]-arr[n+1])/(x[n+2]-x[n+1]);
    double center = (arr[n+2]-arr[n])/(x[n+2]-x[n]);
    return left; // bbbbb
    //    return right; // bbbb
    return ddx(left,center); // bbb
    //    return ddx(left,right); bb
    return center;
    return ddx(left,center);
    return left;
  }

  return 0.; // never get here

}

double ddx(std::vector<double>& arr, unsigned int n, std::vector<double>& x)
{
/*   if(n<nx && n>1) */
/*     return (arr[n+1]-arr[n-1])/(x[n+1]-x[n-1]); */
  unsigned int nx = x.size()-1;
  if(n<nx && n>1) {
    // slope limiter!
    double right = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double left = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    return ddx(left,right);
    if(left*right <=0.) 
      return 0.; 
    if(fabs(left)>fabs(right))
      return right;
    return left;
  }
  if(n==nx) {
    double right = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    double left = (arr[n-1]-arr[n-2])/(x[n-1]-x[n-2]);
    double center = (arr[n]-arr[n-2])/(x[n]-x[n-2]);
    return ddx(left,right);
    return right;
  }
  if(n==1) {
    double left = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double right = (arr[n+2]-arr[n+1])/(x[n+2]-x[n+1]);
    double center = (arr[n+2]-arr[n])/(x[n+2]-x[n]);
    return left; // bbbbb
    return right; // bbbb
    return ddx(left,center);
    //    return ddx(left,right);
    return center;
    return ddx(left,center);
    return left;
  }

  return 0.; // never get here

}

void errormsg(const std::string msg)
{
  std::cerr << "Fatal problem encountered: "<<std::endl;
  std::cerr << msg << std::endl;
  std::cerr << "Your move." << std::endl;
  exit(1);
}

double flux(unsigned int n,std::vector<double>& yy, std::vector<double>& x, std::vector<double>& col_st)
{
  double fluxn;
  if(n>0 && n<x.size()-1) {
    double ym=yy[n+1];
//    if(fabs(ym) > fabs(yy[n])) ym=yy[n]; //minmod
    if(yy[n] * yy[n+1] <= 0.0) ym=0.0;
    double cst = col_st[n+1];
    if(ym > 0.0) cst = col_st[n];

//    fluxn= 2.0*M_PI*sqrt(x[n]*x[n+1]) * ym * cst;
    fluxn = 2.0*M_PI* x[n+1] * ym * cst;
  }
  else if(n==x.size()-1) {
    fluxn=0.0;
  }
  else { // n==0
    double ym = yy[n+1];
    if(ym > 0.0) {
      // something bad has happened
    }
    double cst = col_st[n+1];
//    fluxn = 2.0*M_PI*sqrt(x[n+1]*x[n+1] * (x[1]/x[2])) * ym *cst;
    fluxn = 2.0*M_PI * x[1]*ym*cst;
//    if(fluxn > 0.0) {
//      std::cout << "The bulge is leaking stars! flux, ym, cst, n;  "<<fluxn<<", "<<ym<<", "<<cst<<", "<<n<<std::endl;
//      fluxn=0.0;
//    }
  }

  return fluxn;
}
 

double dSMigdt(unsigned int n,std::vector<double>& yy, std::vector<double>& x, std::vector<double>& col_st)
{
//   double fluxn;
//   double fluxnm1;
// //  return -2*M_PI*(col_st[n]*ddx(yy,n,x) + ddx(col_st,n,x)*yy[n] +col_st[n]*yy[n]/x[n] );

//   if(n!=x.size()-1) {
//     double ym = yy[n+1];
//     if(fabs(ym) > fabs(yy[n])) ym=yy[n]; //minmod
//     if(yy[n] * yy[n+1] <= 0.0) ym=0.0;
//     double cst = col_st[n+1];
//     if(ym > 0.0) cst = col_st[n];
  
//     fluxn = 2.0*M_PI*sqrt(x[n]*x[n+1]) * ym * cst ;
//   }
//   else fluxn=0.0;

// //  if(n==1 && fluxn>1.0) {
// //    fluxn=0.0; // no stars can come out of the bulge.
// //    std::cerr << "Warning: positive flux of stars out of the bulge. flux= "+str(flux[n]) << std::endl;
// //  }


//   if(n!=1) {
//     double ymA = yy[n];
//     if(fabs(ymA) > fabs(yy[n-1])) ymA=yy[n-1]; //minmod
//     if(yy[n-1] * yy[n] <= 0.0) ymA=0.0;
//     double cstA = col_st[n];
//     if(ymA > 0.0) cstA = col_st[n-1];
  
//     fluxnm1 = 2.0*M_PI*sqrt(x[n-1]*x[n]) * ymA * cstA ;
//   }
//   else {
//     double ymA=yy[n];
//     double cstA=col_st[n];
//     fluxnm1 = 2.0*M_PI*sqrt(x[n]*x[n] * x[4]/x[5]) *ymA * cstA;


//     if(fluxnm1 > 0.0) {
//       std::cout << "Warning: positive flux of stars out of the bulge. flux= "+str(fluxnm1) << std::endl;	
//       fluxnm1=0.0; // no stars can come out of the bulge.
//     }
//   }
//  return (fluxnm1 - fluxn) / (x[n]*x[n]*log(x[5]/x[4]));

  double val = - (flux(n,yy,x,col_st) - flux(n-1,yy,x,col_st))/ (x[n]*x[n]*log(x[2]/x[1]));
  if(n==1) {
	std::cout.precision(15);
	//	double a1 = val;
	//double a2 = flux(n,yy,x,col_st);
	//double a3 = flux(n-1,yy,x,col_st);
	//	double a4 = flux(n-2,yy,x,col_st);
	//double a5 = flux(n,yy,x,col_st) - flux(n-1,yy,x,col_st);
	//double a6 = (flux(n,yy,x,col_st) - flux(n-1,yy,x,col_st)) / (x[n]*x[n]*log(x[2]/x[1]));
//	std::cout <<"val, flux(n), flux(n-1); y[n+1], y[n], col_st[n+1], col_st[n]:  "<< val << ", "<<flux(n,yy,x,col_st)<<", "<<flux(n-1,yy,x,col_st) << ";  "<< yy[n+1] << ", "<<yy[n]<<", "<<col_st[n+1]<<", "<<col_st[n]<< std::endl;
	//	std::cout<<"val, f(n), f(n-1), f(n-2), f(n)-f(n-1), dS: " <<a1<<" "<<a2<<" "<<a3<<" "<<a4<<" "<<a5<<" "<<a6<<std::endl;
  }
  return val;
}

double OldIthBin(unsigned int i,Cosmology& cos,unsigned int NAgeBins)
{
  return cos.lbt(cos.ZStart()) * ((double) NAgeBins - i + 1)/((double) NAgeBins);
}
double YoungIthBin(unsigned int i,Cosmology& cos,unsigned int NAgeBins)
{
  if(NAgeBins==1) return 0.;
  return cos.lbt(cos.ZStart()) * ((double) NAgeBins - i)/((double) NAgeBins);
}

double max(double a, double b)
{
  if(a > b)
    return a;
  return b;
}

std::string str(int val)
{
  std::stringstream s;
  s << val;
  return s.str();
}
std::string str(double val)
{
  std::stringstream s;
  s<<val;
  return s.str();
}
std::string str(unsigned int val)
{
  std::stringstream s;
  s << val;
  return s.str();
}
std::string str(long int val)
{
  std::stringstream s;
  s<<val;
  return s.str();
}
double Qsimple(unsigned int n,DiskContents& disk)
{
  RafikovQParams rqp;
  disk.ComputeRafikovQParams(&rqp,n);
  double absc=1.;
  return Q(&rqp,&absc);
}

// Same as QmfQ (below), except instead of varying sig and sig_st simultaneously by the same factor (sv),
// we just vary sig. The former is useful if we fix the ratio, phi0, in the initial conditions, and the 
// latter is useful if instead we just set the absolute value of sig_st and let the gas relax around the stars.
double QmfQfst(double sv, void *p)
{
  struct RafikovQParams * qp = (RafikovQParams *) p;
  (*qp).var=-1;
  (*qp).Qg *= sv;
  if((*qp).fixedQ < 0.0) errormsg("The fixedQ passed to QmfQfst (in DiskUtils) was not initialized.");
  double val = Q(qp,&((*qp).mostRecentq))-(*qp).fixedQ;
  (*qp).Qg /=sv;

  return val;
}

double QmfQ(double sv, void *p) 
{
  // here sv is a factor by which we will multiply all velocity dispersions (and hence Qsi - the ri will be constant since both the gas and stellar velocity dispersions will be multiplied by sv)
  struct RafikovQParams * qp = (RafikovQParams *) p;
  (*qp).var=-1;
  (*qp).Qg *= sv;
  for(unsigned int i=0; i!=(*qp).Qsi.size(); ++i) {
    (*qp).Qsi[i] *= sv;
  }
  if((*qp).fixedQ < 0.0) errormsg("The fixedQ passed to QmfQ (in DiskUtils) was not initialized.");
  double val= Q(qp,&((*qp).mostRecentq))-(*qp).fixedQ;
  //  if(sv<100) std::cout << "dbg Qm1: "<< (*qp).Qg <<" " << (*qp).Qsi[0] << " "<<(*qp).ri[0]<<" " << val<<std::endl;
  (*qp).Qg /= sv;
  for(unsigned int i=0; i!=(*qp).Qsi.size(); ++i) {
    (*qp).Qsi[i] /= sv;
  }  
  return val;
}

double Q(RafikovQParams *qp, double *absc)
{
  if((*qp).analyticQ ) { // Romeo and Wiegert 2011
    double QsinvRi3 = 0.;
    double QsinvRi = 0.;
    for(unsigned int i=0; i!= (*qp).ri.size(); ++i) {
      QsinvRi3 +=(*qp).ri[i] * (*qp).ri[i] * (*qp).ri[i] / ((*qp).Qsi[i]);
      QsinvRi += (*qp).ri[i] / (*qp).Qsi[i];
    }
    double rs= sqrt( QsinvRi3 / QsinvRi) ;
    double Qst = rs / QsinvRi;
    double W = 2./ (rs + 1./rs);
    
    if(Qst > (*qp).Qg) {
      return (*qp).thick/ (W/Qst + 1./((*qp).Qg) );
    }
    else {
      return (*qp).thick/(1./Qst + W/((*qp).Qg) );
    }
  }
  else {
    //  if((*qp).Qg<100)  std::cout << "dbg Q " << (*qp).Qg << " "<< (*qp).Qsi[0] <<" "<< (*qp).ri[0]<< std::endl;
    // Find Q given these parameters qp by minimizing Q(q).
    // absc contains an initial guess, and stores the final abscissa.
    gsl_function F,fn;
    F.function = &dQdq;
    fn.function= &Qq;
    F.params = &(*qp);
    fn.params = &(*qp);
    if (*absc <=0 ) { *absc=1;}
    return (*qp).thick*minFromDeriv(F,fn,&(*absc));
  }
  
  return -1; // never get here
}
double varQ(double sv, void * p)
{
  struct RafikovQParams * qp = (RafikovQParams *) p;
  
  double val;
  double absc= qp->mostRecentq; // guess
  if(qp->var == -1) { // do not change any values, ignore sv, update mostRecentq
    return Q(qp,&(qp->mostRecentq));
  }
  else if(qp->var==0) { // changing Q_g
    double tempQg = qp->Qg; // save the real value of Qg
    qp->Qg=sv; // now vary Qg.
    
    val = Q(qp,&absc);

    qp->Qg=tempQg; // return Qg to its original, real value.
  }
  else if(qp->var <= qp->Qsi.size() ) { // changing Qsi
    int ind = qp->var -1;
    double tempQsi = qp->Qsi[ind];
    qp->Qsi[ind]=sv;

    val = Q(qp,&absc);
    
    qp->Qsi[ind]=tempQsi;
  }
  else if(qp->var <= 2 * (qp->ri.size())) { // changing ri
    int ind = qp->var -1 - qp->ri.size();
    double tempri = qp->ri[ind];
    qp->ri[ind]=sv;

    val = Q(qp,&absc);

    qp->ri[ind]=tempri;
  }
  else {
    errormsg("Q: variable out of range");
  }
  return val;
}

double I0Exp(const double x)
{
  // return I_0(x)*exp(-x)
  gsl_sf_result res;
  // given the large value of x, I0 will be large
  // and exp will be small, so to avoid the loss in precision
  // and possible overflow, use a version which takes out the 
  // exponential dependence, 
  // i.e. gsl_sf_bessel_I0_scaled is in fact I0(x)exp(-x)
  int status = gsl_sf_bessel_I0_scaled_e(x,&res);
  // the error in this estimate may be accessed with
  // res.err
  return res.val;
}
double I1Exp(const double x) 
{
  gsl_sf_result res;
  int status = gsl_sf_bessel_I1_scaled_e(x,&res);
  return res.val;
}

int findRoot(gsl_function & F, double * guess)
{
  int status;
  int iter=0,max_iter=200;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r=0,low,high,flow,fhigh;
  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);
  
  bool large=false;

  low=.9*(*guess); high=1.1*(*guess);
  flow=GSL_FN_EVAL(&F,low); fhigh=GSL_FN_EVAL(&F,high);
  int niter=0;
  while(flow*fhigh > 0) {// no zero crossings
    //    if(niter<100) std::cout <<"dbg findroot: f("<<low<<")="<< flow <<", f("<<high<<")=" << fhigh << std::endl;
    // expand range
    if((fabs(flow)<fabs(fhigh) && !large) || (fabs(flow)>fabs(fhigh) && large)) {
//      low -= (high-low)*1.2;
      low *= 0.8;
      flow=GSL_FN_EVAL(&F,low);
    }
    else {
//      high += (high-low)*1.2;
      high *= 1.2;
      fhigh=GSL_FN_EVAL(&F,high);
    }
//    large = (++niter > 100);
    ++niter;
  }
  if(fabs(flow) > 1.0e30 && fabs(fhigh) > 1.0e30) {
    errormsg("Failed to bracket a root: low,high,flow,fhigh,niter: "+str(low)+" "+str(high)+" "+str(flow)+" "+str(fhigh)+" "+str(niter));
  }

  gsl_root_fsolver_set(s,&F,low,high);
  do {
    ++iter;
    //    status = gsl_root_fsolver_iterate(s);
    gsl_root_fsolver_iterate(s);
    r=gsl_root_fsolver_root(s);
    low =gsl_root_fsolver_x_lower(s);
    high=gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(low,high,0.0,1.0e-6);
  } while (status==GSL_CONTINUE && iter < max_iter);
  if(status==GSL_SUCCESS) {
    *guess=r;
  }
  else {
    errormsg("findRoot failed to converge");
  }
  return status;
}

// find the ~global minimum of fn using fn's derivative, F.
double minFromDeriv(gsl_function & F, gsl_function & fn, double * abcissa)
{
  int status;
  int iter=0, max_iter=100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  T=gsl_root_fsolver_brent;
  s=gsl_root_fsolver_alloc(T);
  double r=0,low,high;
  std::vector<double> lo(0), hi(0),extrema(0);
  unsigned int N=25;
  double delta = (*abcissa)*(10-.00001)/((double) N);
  for(unsigned int i=0; i<N; ++i) {
    low = .00001*(*abcissa) + delta*i;
    high= low+delta;
    //    std::cout << "minFromDeriv: values of Qq: " << GSL_FN_EVAL(&fn,low) << " " << GSL_FN_EVAL(&fn,high) << " " << low << " " << high << std::endl;
    if(GSL_FN_EVAL( &F, low) * GSL_FN_EVAL(&F,high) < 0) { // there's a zero crossing here!
      lo.push_back(low); hi.push_back(high);
    }
  }
  
  if(lo.size()==0) {
    //    errormsg("No root crossings found in minFromDeriv!  abcissa,F(low),F(high)  "+str(*abcissa)+" "+str(GSL_FN_EVAL(&F,low))+" "+str(GSL_FN_EVAL(&F,high)));
    return -1;
  }
  //  std::cout << "minFromDeriv: lo.size, hi.size, " << lo.size() <<" " <<hi.size() << std::endl;
  for(unsigned int i=0; i<lo.size(); ++i) { // for each zero,crossing...
    gsl_root_fsolver_set(s,&F,lo[i],hi[i]);
    do {
      ++iter;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      //      lo[i] = gsl_root_fsolver_x_lower(s);
      //      hi[i] = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_residual(GSL_FN_EVAL(&F,r),1.e-10);
    } while (status==GSL_CONTINUE && iter < max_iter);
    if(status==GSL_SUCCESS) {
      extrema.push_back(r);
    }
  }
  gsl_root_fsolver_free(s);

  double globalmin = 1.e10,v;
  for(unsigned int i=0; i<extrema.size(); ++i) {
    v = GSL_FN_EVAL(&fn,r);
    if(v<globalmin) {
      globalmin = v;
      *abcissa = r;
    }
  }
  return globalmin;
}


double dQdq(double q, void * p)
{
  struct RafikovQParams * qp = (RafikovQParams *) p;
  double sum=0.;
  double sum2=0.;
  double I0E,I1E;
  for(unsigned int i=0; i!=qp->Qsi.size(); ++i) {
    I0E = I0Exp(q*q*(qp->ri[i])*(qp->ri[i]));
    I1E = I1Exp(q*q*(qp->ri[i])*(qp->ri[i]));
    sum += (1./(qp->Qsi[i]))* (1.-I0E)/(q*(qp->ri[i]));
    sum2+= (2*I0E*(qp->ri[i]) - 2.* I1E*(qp->ri[i]))/(qp->Qsi[i]);
  }
  return -( (1-q*q)/((1+q*q)*(1+q*q)*(qp->Qg)) -  sum/q + sum2 ) /
    ( 2.*( q/((qp->Qg)*(1+q*q)) + sum) * ( q/((qp->Qg)*(1+q*q)) + sum));
}

double Qq(double q, void * p)
{
  if(q<=0.)
    return 1.e30; // to prevent minimization routines from crossing q=0

  struct RafikovQParams * qp = (RafikovQParams *) p;
  double sum=0.;
  for(unsigned int i=0; i!=qp->Qsi.size(); ++i) {
    sum += (1./qp->Qsi[i])* (1.-I0Exp(q*q*(qp->ri[i])*(qp->ri[i])))/(q*(qp->ri[i]));
  }
  return 1./( 2./(qp->Qg) * q/(1+q*q) + 2.*sum );
}


double arrmax(std::vector<double>& arr )
{
  double themax=0.;
  for(unsigned int n=1; n<arr.size(); ++n) {
    if(themax < fabs(arr[n]) )
      themax=fabs(arr[n]);
  }
  return themax;
}
