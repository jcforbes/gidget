#include <vector>

// A very simple struct to store parameters necessary to
// evaluate the Rafikov Q.

struct RafikovQParams {
  RafikovQParams();
  int var;
  std::vector<double> Qsi; // Q_{*,i}
  std::vector<double> ri; // sigma_{*,i}/sigma
  double Qg;
  double mostRecentq; // the most recent least stable wavenumber
  double thickGas;
  double thickStars;
  double fixedQ;
  bool analyticQ;
};
