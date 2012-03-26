#include "ConvergenceCheck.h"
#include "Dimensions.h"

ConvergenceCheck::ConvergenceCheck(double st,double et,double lx,double ux) :
  startTime(st), endTime(et), lowerX(lx), upperX(ux),
  numerator(0.), denominator(0.)
{ }

// For single quantities which we'd like to average over time
void ConvergenceCheck::Update(double qu, double t)
{
  if(t>=startTime && t <= endTime) { 
    numerator += qu;
    denominator += 1;
  }

}

// For quantities defined at every point in the mesh wish we'd like to average over time and space
void ConvergenceCheck::Update(std::vector<double>& qu, double t,std::vector<double>& x)
{
  if(t>=startTime && t <=endTime) {
    for(unsigned int i=1; i<=x.size()-1; ++i) {
      if(x[i] >= lowerX && x[i] <= upperX) {
        numerator+=qu[i];
        denominator+=1;
      }
    }
  }
}
void ConvergenceCheck::Update(double* qu, double t,std::vector<double>& x)
{
  if(t>=startTime && t <=endTime) {
    for(unsigned int i=1; i<=x.size()-1; ++i) {
      if(x[i] >= lowerX && x[i] <= upperX) {
        numerator+=qu[i];
        denominator+=1;
      }
    }
  }
}
void ConvergenceCheck::UpdateD(double qu,double t, double dt,Dimensions& dim)
{
  if(t<= startTime  && startTime < t+dt ) {
    numerator= -qu;
    denominator = -dim.t(t); // tPhys(t)/speryear;
  }
  if(t<= endTime && endTime < t+dt) {
    numerator+=qu;
    denominator += dim.t(t); // tPhys(t)/speryear;
  }
}

double ConvergenceCheck::ratio()
{
  return numerator/denominator;
}
