#include <vector>
class Dimensions;

// Each time step call a function which updates the following quantities
// function takes <struct>, quantity (or array of quantity), time, dt
// avg. from time t=st to t=et
// avg. from cell n=a to n=b
// need to track:
// numerator (i.e. sum of quantity over relevant time [& dist.])
// denominator: actual time elapsed
class ConvergenceCheck {
 public:
  // average over the time period st to et, and the x range lx to ux
  ConvergenceCheck(double st,double et,double lx,double ux);

  // Let qu contribute to the average if t is within st to et.
  void Update(double qu,double t);

  // Let qu contribute to the quantity's average value if 
  // t is within st to et AND only for the region of x 
  // from lx to ux.
  void Update(std::vector<double>& qu,double t,
	      std::vector<double>& x);

  // Same as above, except here the quantity is a double* 
  // instead of an std::vector
  void Update(double* ,double,std::vector<double>&);

  // Instead of averaging over space and time, compute an 
  // average rate of change between st and et.
  void UpdateD(double qu,double t,double dt,Dimensions&);

  // As a final step in all this averaging, compute the 
  // ratio of the accumulated numerator to the 
  // accumulated denominator
  double ratio();

 private:
  double startTime;
  double endTime;
  double lowerX;
  double upperX;
  double numerator;
  double denominator;
};
