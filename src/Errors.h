
#include <fstream>
#include <string>

class errormsg {
 public:
  errormsg(const std::string msg, bool fatal=true);
  static std::ofstream errorFile;
};

std::ofstream errormsg::errorFile;


