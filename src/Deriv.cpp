#include <math.h>
#include <vector>
#include "Deriv.h"
#include "DiskUtils.h"

double maximum(const double a, const double b)
{
  if(a>b)
    return a;
  else 
    return b;
}


// Compute the numerical derivative of a given function.
// Arguments: 
//  f    - the 1d function to be differentiated
//  h0   - the initial step size
//  etol - fractional error tolerance
//  x0   - position at which to take the derivative
//  p    - an object containing parameters needed to evaluate f
double derivDriver(double (*f)(double,void*),const double x0, double etol,void* p, double* err)
{
  double h0=.1*x0;
  double error=1.0e30;
  double red=1.0;
  double result=0.0;
  const double min = 1.0e-17;
  do {
    result=deriv(f,x0,h0*red,p,&error);
    red*=0.5;
  } while(fabs(error) > etol && red>min);
  if(red<=min) errormsg("Unable to get good derivative");
  *err=error;
  return result;
}

// An implementation of Ridders method to extrapolate the approximation to
// f'(x0) to h->0. The extrapolation uses Neville's algorithm
double deriv(double (*f)(double,void*),const double x0,const double h0,void *p,double *error)
{
  const double c_red=2.2; // factor by which to reduce h
  const unsigned int m=10; // maximum number of iterations
  const double safe=2.0; // factor by which the error is allowed to increase before termination
  
  double best;
  double red=c_red;
  volatile double temp1,temp2;
  double h=h0;

  std::vector<std::vector<double> > tab(m);
  for(unsigned int i=0; i!=tab.size(); ++i) {
    tab[i].resize(m);
  }
  temp1=x0+h;
  temp2=x0-h;
  double dif=temp1-temp2;
  tab[0][0] = ((*f)(temp1,p) - (*f)(temp2,p))/dif;
  double err=1.e30; // a large number..
  for(unsigned int i=1;i!=m;++i) {
    h=h/sqrt(c_red);
    temp1=x0+h;
    temp2=x0-h;
    dif=temp1-temp2;
    tab[0][i] = ((*f)(temp1,p)-(*f)(temp2,p))/dif;
    for(unsigned int j=1;j<=i;++j) {
      tab[j][i]=(tab[j-1][i]*red - tab[j-1][i-1])/(red-1.);
      red*=c_red;
      
      double temp_err = maximum(fabs(tab[j][i]-tab[j-1][i]),fabs(tab[j][i]-tab[j-1][i-1])); 
      if(temp_err <= err){
	err=temp_err;
	best=tab[j][i];
      }
    } // end inner loop over tableau, i.e. loop over terms needed to extrapolate to higher orders
    if(fabs(tab[i][i]-tab[i-1][i-1]) > safe*err)
      break;
  } // end outer loop over tableau, i.e. loop over orders of approximation to the derivative with successively decreasing stepsize
  *error=err/x0;
  return best;
}
