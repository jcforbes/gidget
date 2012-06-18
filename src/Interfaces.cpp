#include "Interfaces.h"
#include <math.h>

Interfaces::Interfaces(std::vector<int> & keepTorqueOff,
		       std::vector<double> & xRef) :
  x(xRef),
  keepTorquesOff(keepTorqueOff),
  InterfacePositions(std::vector<int>(0)),
  InterfaceDirections(std::vector<int>(0)),
  DistanceToLeftInterface(std::vector<int>(xRef.size())),
  DistanceToRightInterface(std::vector<int>(xRef.size()))

{
  for(unsigned int n=2; n<=x.size()-1; ++n) {
    if(keepTorquesOff[n]==1 && keepTorquesOff[n-1]==0) {
      InterfacePositions.push_back(n);
      InterfaceDirections.push_back(-1);
    }
    if(keepTorquesOff[n]==0 && keepTorquesOff[n-1]==1) {
      InterfacePositions.push_back(n-1);
      InterfaceDirections.push_back(1);
    }
  }
  int lastUnstable=-1;
  for(unsigned int n=2; n<=x.size()-1; ++n) {
    if(keepTorquesOff[n]==0) {
      lastUnstable=n;
      DistanceToLeftInterface[n]=-1;
    }
    else if(lastUnstable!=-1) // there's been some unstable region
      DistanceToLeftInterface[n]=n-lastUnstable;
    else
      DistanceToLeftInterface[n]=-1;
  }
  lastUnstable=x.size()-1;
  for(unsigned int n=x.size()-2; n>=1; --n) {
    if(keepTorquesOff[n]==0) {
      lastUnstable=n;
      DistanceToRightInterface[n]=-1;
    }
    else
      DistanceToRightInterface[n]=lastUnstable-n;
  }
  return;
}


double Interfaces::weight(unsigned int n, bool cells, double dist)
{
  int dToR = DistanceToRightInterface[n];
  if(cells && dToR!=-1)
    return exp(-((double) dToR)/dist);

  if(!cells && dToR!=-1)
    return exp(-(x[dToR+n] - x[n]) / dist);

  return 0.0;
}

unsigned int Interfaces::index(unsigned int n, bool right)
{
  if(keepTorquesOff[n]==0)
    return 0;

  if(right)
    return DistanceToRightInterface[n]+n;

  return ((unsigned int) ((int) n) - DistanceToLeftInterface[n]);

}
