#include <vector>

class Interfaces {
 public:
  Interfaces(std::vector<int> & keepTorquesOff, std::vector<double> & x );
  //  double DistanceToNearestInterface(bool distance);
  //  int CellsToNearestInterface(bool cells);

  // the weight to give to the nonphysical GI solution
  // n is the current cell, cells tells us whether we care 
  // about a distance in cells or a physical distance
  // dist is the characteristic decay distance.
  double weight(unsigned int n, bool cells, double dist);

  // the index of the closest cell which is unstable.
  unsigned int index(unsigned int n, bool right);
 private:
  std::vector<double> & x;
  std::vector<int> keepTorquesOff;
  std::vector<int> InterfacePositions;
  std::vector<int> InterfaceDirections;
  std::vector<int> DistanceToLeftInterface;
  std::vector<int> DistanceToRightInterface;

};
