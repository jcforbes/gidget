#ifndef error_message_header
#define error_message_header

#include <fstream>
#include <string>

class errormsg {
 public:
  errormsg(const std::string msg, bool fatal=true);
  static std::ofstream errorFile;
};


#endif

