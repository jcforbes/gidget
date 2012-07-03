#include "Debug.h"
#include <math.h>
#include <iostream>

int LargestPowerOfTwo(int number)
{
  if(number>0)
    return (unsigned int) floor(log((double) number) / log(2.0)+.00000001);
  return -1;
}

Debug::Debug(int number) :
  theList(std::vector<int>(LargestPowerOfTwo(number)+1))
{
  int subtracted = (int) number;
  for(int i=theList.size()-1; i>=0 ; --i) {
    int p = LargestPowerOfTwo(subtracted);
    if (  p == i ) {
      theList[i] = 1;
      subtracted-=pow(2,p);
    }
    else
      theList[i] = 0;


  }
}

bool Debug::opt(int number)
{
  if(number < theList.size() && number >= 0) {
    if(theList[number]==1) return true;
    else return false;
  }
  else return false;

}

void testDebug() 
{
  std::vector<bool>grades (0);
  
  grades.push_back(LargestPowerOfTwo(5)==2);
  grades.push_back(LargestPowerOfTwo(4)==2);
  grades.push_back(LargestPowerOfTwo(1)==0);

  Debug dbg(114);
  grades.push_back(!dbg.opt(7));
  grades.push_back(dbg.opt(6));
  grades.push_back(dbg.opt(5));
  grades.push_back(dbg.opt(4));
  grades.push_back(!dbg.opt(3));
  grades.push_back(!dbg.opt(2));
  grades.push_back(dbg.opt(1));
  grades.push_back(!dbg.opt(0));

  Debug dbg2(0);
  grades.push_back(!dbg2.opt(0));
  grades.push_back(!dbg2.opt(1));
  grades.push_back(!dbg2.opt(2));

  int total = 0;
  for(unsigned int i=0; i!=grades.size(); ++i) {
    if(grades[i]) total++;
    else std::cout << "testDebug has failed test number "<<i<<" !"<<std::endl;
  }

  std::cout << "testDebug passes "<<total<<" of "<<grades.size()<<" tests."<<std::endl;
  return;
}
