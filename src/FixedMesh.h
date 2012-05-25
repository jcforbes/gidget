#include <vector>

class FixedMesh {
 public:
  FixedMesh(double,double,double,double,unsigned int);
  ~FixedMesh();
  double * x_GSL() { return x_gsl;};
  std::vector<double> & x() {return xv;};
  std::vector<double> & beta() {return betav;};
  std::vector<double> & betap() {return betapv;};
  std::vector<double> & uu() {return uuv;};
  double x(unsigned int);
  double x(double);
  double psi(unsigned int);
  double psi(double x);
  double beta(double x);
  double betap(double x);
  double uu(double x);
  double dlnx() { return dlnxc;};
  double xmin() { return xminc;};
  unsigned int nx() { return nxc;};
  unsigned int necessaryN(double minsigst);
  double n(unsigned int index,unsigned int neff);
 private:
  std::vector<double> xv;
  std::vector<double> psiv;
  std::vector<double> uuv;
  std::vector<double> betav;
  std::vector<double> betapv;
  const double ip, // power law index of vphi(R)
                b, // turnover radius.
             soft, // degree to which we soften the turnover.
            dlnxc,
            xminc;
  double * x_gsl;
  const unsigned int nxc;
  double psi1, bsf,bmsf;
  bool stored;
  unsigned int necesN;
};
