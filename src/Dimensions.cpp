#include "Dimensions.h"
#include <math.h>

Dimensions::Dimensions(double r,double v, double md) : Radius(r),vphiR(v),MdotExt0(md){}
double Dimensions::t_cgs(double tSim){ return tSim*2*M_PI*Radius/vphiR;}
double Dimensions::v_cgs(double vv){return vv*vphiR; }
double Dimensions::d_cgs(double rr){return rr*Radius;}
double Dimensions::col_cgs(double cd){return cd*MdotExt0/(vphiR*Radius);}

double Dimensions::v(double vv){return v_cgs(vv)*1.e-5;}
double Dimensions::t(double tt){return t_cgs(tt)/speryear;} // years
double Dimensions::d(double dd){return d_cgs(dd)/cmperkpc;} // kpc


double Dimensions::chi()
{
  return G*MdotExt0/(vphiR*vphiR*vphiR);
}
