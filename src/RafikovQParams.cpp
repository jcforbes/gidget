#include "RafikovQParams.h"


// The constructor - initialize vectors to be 0 elements
// and set mostRecentq to an obviously non-physical value
RafikovQParams::RafikovQParams() :
  var(-1),Qsi(std::vector<double>(0)),ri(std::vector<double>(0)),
  Qg(-1.),mostRecentq(-1.),analyticQ(false),thickGas(1.0),fixedQ(-1.),
  thickStars(1.0)
{}
