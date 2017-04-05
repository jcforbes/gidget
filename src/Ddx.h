#include <math.h>
#include <vector>
#include "Errors.h"

inline double ddxAvg(double left, double right)
{
  if(left!=left || right!=right)
      errormsg("Error computing ddx");

  return (left+right)/2.0;
}

inline double ddx(double left, double right)
{
  if(left!=left || right!=right)
      errormsg("Error computing ddx");

  if(left*right <=0.)
    return 0.0;
  if(fabs(left)>fabs(right))
    return right;
  return left;
};

template <class T>
double ddx(T & arr, unsigned int n, std::vector<double>& x, bool zeroSecondDerivIn, bool minmod)
{
  bool zeroSecondDeriv = false; //zeroSecondDerivIn;// true;
  unsigned int nx=x.size()-1;
  if(n<nx && n>1) {
    double right = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double left = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    if(minmod) return ddx(left,right);
    else return ddxAvg(left,right);
  }
  if(n==nx) {
    double right = (arr[n]-arr[n-1])/(x[n]-x[n-1]);
    double left = (arr[n-1]-arr[n-2])/(x[n-1]-x[n-2]);
    double center = (arr[n]-arr[n-2])/(x[n]-x[n-2]);
    if(minmod) return ddx(left,right); 
    else return ddxAvg(left,right);
  }
  if(n==1) {
    double left = (arr[n+1]-arr[n])/(x[n+1]-x[n]);
    double right = (arr[n+2]-arr[n+1])/(x[n+2]-x[n+1]);
    double center = (arr[n+2]-arr[n])/(x[n+2]-x[n]);
    if(zeroSecondDeriv) {
      if(minmod) return ddx(left,right); 
      else return ddxAvg(left,right);
    }
    else
      return left; 
  }
  errormsg("Unexpected input to ddx in Ddx.h");
  return 0;
};

