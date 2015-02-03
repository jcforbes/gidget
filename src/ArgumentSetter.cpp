#include "ArgumentSetter.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>

ArgumentSetter::ArgumentSetter(int ac, char** av,std::string fn) :
  argc(ac), argv(av), arg(2), cmtfile((fn).c_str())
{ 
  cmtfile.precision(15);
  cmtfile << "--------------" << std::endl;
  cmtfile << argv[1] << std::endl;
}
void ArgumentSetter::WriteOut()
{
  cmtfile << "--------------" << std::endl;
  for(unsigned int i=0;i!=argc;++i) {
    cmtfile << argv[i] <<" ";
  }
  cmtfile<<std::endl;
}
ArgumentSetter::~ArgumentSetter()
{
  cmtfile.close();
}
double ArgumentSetter::Set(double def,std::string name)
{
  double val;
  if(arg<argc) {
    val = atof(argv[arg]); 
  }
  else {
    val = def;
  }
  std::cout << "Setting "+name+ ": " << val << std::endl;
  cmtfile << name +": "<<val<<std::endl;  ++arg;
  return val;
}
