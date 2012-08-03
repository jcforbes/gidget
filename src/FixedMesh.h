#include <vector>
#include <gsl/gsl_spline.h>

class FixedMesh {
 public:
  FixedMesh(double,double,double,double,double,unsigned int);
  ~FixedMesh();
  double * x_GSL() { return x_gsl;};
  std::vector<double> & x() {return xv;};
  std::vector<double> & beta() {return betav;};
  std::vector<double> & betap() {return betapv;};
  std::vector<double> & uu() {return uuv;};
  double x(unsigned int);
  double xPlusHalf(unsigned int);
  double dx(unsigned int);
  double x(double);
  double psi(unsigned int);
  double psi(double x);
  double beta(double x);
  double betap(double x);
  double uu(double x);
  double u1pbPlusHalf(unsigned int);
  double dlnx() { return dlnxc;};
  double xmin() { return xminc;};
  unsigned int nx() { return nxc;};
  unsigned int necessaryN();
  double n(unsigned int index,unsigned int neff);
  bool InitializePsi();
 private:
  std::vector<double> xv;
//  std::vector<double> psiv;
  std::vector<double> uuv;
  std::vector<double> betav;
  std::vector<double> betapv;
  std::vector<double> xiPlusHalf;
  std::vector<double> dxi;
  std::vector<double> u1pbiPlusHalf;
  const double ip, // power law index of vphi(R)
                b, // turnover radius.
             soft, // degree to which we soften the turnover.
            dlnxc,
            xminc;
  double * x_gsl;
  double * x_HR_GSL;
  double * psi_HR_GSL;

  gsl_spline * spline_psi;
  gsl_interp_accel * accel_psi;

  const unsigned int nxc;
  double psi1, bsf,bmsf;
  bool stored, PsiInitialized;
  unsigned int necesN;

  double minsigst;
};
