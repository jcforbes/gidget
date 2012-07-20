#include "Errors.h"

errormsg::errormsg(const std::string msg, bool fatal)
{
  errormsg::errorFile << msg;
  if(fatal) {
      errormsg::errorFile << std::endl;
      errormsg::errorFile.close();
      exit(1);
  }
};

std::ofstream errormsg::errorFile;



