#include "Errors.h"
#include <iostream>

errormsg::errormsg(const std::string msg, bool fatal)
{
  std::cerr << msg << std::endl;
  errormsg::errorFile << msg;
  if(fatal) {
      errormsg::errorFile << std::endl;
      errormsg::errorFile.close();
      exit(1);
  }
};

std::ofstream errormsg::errorFile;



