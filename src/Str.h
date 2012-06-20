#include <sstream>
#include <string>

template <class T>
std::string str(const T val)
{
  std::stringstream s;
  s<<val;
  return s.str();
}
