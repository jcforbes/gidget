#include <fstream>
#include <string>

// A small class to handle the command line arguments
class ArgumentSetter {
 public:
  // Read in argc and argv, and open a file fn which will
  // be a human-readable record of this run's parameters
  ArgumentSetter(int ac,char** av,std::string fn);  
  ~ArgumentSetter();
  /* Call this function a bunch of times, each time do the following
     - read the next element in argv, if there is one
     - return that value. If there was no provided value in argv, use the default, def
     - write a line in the comment file of the form "name: value" */
  double Set(double def,std::string name);
  void WriteOut(); // take the collected arguments and write them out.
 private:
  int argc, arg;
  char ** argv;
  std::ofstream cmtfile;
};
