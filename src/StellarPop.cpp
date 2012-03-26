#include "StellarPop.h"
#include "Cosmology.h"
#include "DiskContents.h"
#include "DiskUtils.h"
#include <math.h>
#include <vector>
#include <iostream>

StellarPop::StellarPop(unsigned int nx,double y,double o) :
  spcol(std::vector<double>(nx+1,0.)), spsig(std::vector<double>(nx+1,0.)),
  spZ(std::vector<double>(nx+1,0.)), spZV(std::vector<double>(nx+1,0.)),
  dQdS(std::vector<double>(nx+1,0.)),
  dQds(std::vector<double>(nx+1,0.)), dQdSerr(std::vector<double>(nx+1,0.)),
  dQdserr(std::vector<double>(nx+1,0.)),youngest(y),oldest(o),ageAtz0((o+y)/2.0)
{ }

StellarPop::StellarPop() :  
  youngest(-1.),oldest(-1.),ageAtz0(-1.)      
{ }

bool StellarPop::IsForming(Cosmology& cos, double redshift)
{
  double age = cos.lbt(redshift);
  return (age > youngest && age<=oldest);
}

// Merge sp2 into sp1. sp2 should be unaffected by the procedure.
void StellarPop::MergeStellarPops(const StellarPop& sp2,DiskContents& disk)
{
  for(unsigned int i=1; i<=spcol.size()-1; ++i) {
    if(sp2.spcol[i]>0.0) { // if there are actually stars to add.
      // otherwise do nothing.
      
      // sig3 ^2 = (col1*sig1^2 + col2*sig2^2)/(col1+col2):
      (*this).spsig[i]= sqrt((((*this).spcol[i])*((*this).spsig[i])*((*this).spsig[i]) 
                              + (sp2.spcol[i])*(sp2.spsig[i])*(sp2.spsig[i]))
                             /((*this).spcol[i]+sp2.spcol[i]));
      //    (*this).spZ[i] = ((*this).spZ[i] * (*this).spcol[i] + sp2.spZ[i] * sp2.spcol[i]) / ((*this).spcol[i] + sp2.spcol[i]);
      // explicitly state the moments of each metallicity distribution:
      double wtdAvg = ((*this).spZ[i]*(*this).spcol[i] + sp2.spZ[i]*sp2.spcol[i])/((*this).spcol[i] + sp2.spcol[i]);
      double avg1 = (*this).spZ[i];
      double avg2 = sp2.spZ[i];
      double var1 = (*this).spZV[i];
      double var2 = sp2.spZV[i];
      double wt1 = (*this).spcol[i];
      double wt2 = sp2.spcol[i];
      
      // Merge the two distributions:
      (*this).spZ[i] = wtdAvg;
      (*this).spZV[i] = wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
        + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);
      
      if((*this).spsig[i]!=(*this).spsig[i])
        errormsg("Error merging populations:  this spcol,spsig  sp2 spcol,spsig  "+str((*this).spcol[i])+" "+str((*this).spsig[i])+"  "+str(sp2.spcol[i])+" "+str(sp2.spsig[i]));
    }
  }
  double m1= disk.TotalWeightedByArea((*this).spcol); // mass of sp1
  double m2= disk.TotalWeightedByArea(sp2.spcol); // mass of sp2
  if(m1<0 || m2<0)
	errormsg("Error merging populations: m1 or m2 is negative: "+str(m1)+" "+str(m2));   

  double before = (*this).ageAtz0; 
  (*this).ageAtz0= (m1*((*this).ageAtz0) + m2*(sp2.ageAtz0)) / (m1+m2); // avg. age by mass.
  std::cerr.precision(15);
//  std::cerr << "Merge ageAtz0: " << before/(speryear*1.0e9) << " " << (*this).ageAtz0/(speryear*1.0e9) << " with m1, m2= "<<m1<<" "<<m2 << std::endl;

  // add the stars from sp2 to sp1.
  for(unsigned int i=1; i<=sp2.spcol.size()-1;++i) {
    (*this).spcol[i]+=sp2.spcol[i];
  }

}
void StellarPop::extract(StellarPop& sp2, double frac) 
{
  for(unsigned int n=1; n<=sp2.spcol.size()-1; ++n) {
    spsig[n] = sp2.spsig[n];
    spcol[n] = frac*sp2.spcol[n];
    spZ[n] = sp2.spZ[n];
    spZV[n] = sp2.spZV[n];
    sp2.spcol[n] -= spcol[n];;
  }
  ageAtz0 = sp2.ageAtz0;

}

void StellarPop::MigrateStellarPop(double dt, std::vector<double>& yy, DiskContents& disk)
{
  std::vector<double> dcolh((*this).spcol.size());
  std::vector<double> dsigh(dcolh.size());
  std::vector<double> dZh(dcolh.size());
  std::vector<double> incomingMass(dcolh.size());
  std::vector<double> outgoingMass(dcolh.size());
  std::vector<double> incomingZ(dcolh.size());
  std::vector<double> incomingZV(dcolh.size());
  std::vector<double> cellMass(dcolh.size());

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    dcolh[n] = -2.*M_PI*((*this).spcol[n]*ddx(yy,n,disk.GetX()) + ddx((*this).spcol,n,disk.GetX()) *yy[n] + (*this).spcol[n]*yy[n]/disk.GetX()[n]);
    dsigh[n] = -2.*M_PI*yy[n]*((1.+disk.GetBeta()[n])*disk.GetUu()[n]*disk.GetUu()[n]/(3.*(*this).spsig[n]*disk.GetX()[n]) + ddx((*this).spsig,n,disk.GetX()));
    dZh[n] = -2.*M_PI*yy[n]*ddx((*this).spZ,n,disk.GetX());

    if(n==spcol.size()-1) {
      incomingMass[n]=0.0;
      incomingZ[n]=0;
      incomingZV[n]=0;
    }
    else {
      incomingMass[n] = -4*M_PI*M_PI*disk.GetX()[n+1]*spcol[n+1]*yy[n+1]*dt;
      incomingZ[n] = spZ[n+1];
      incomingZV[n]=spZV[n+1];
    }
    outgoingMass[n] = -4*M_PI*M_PI*disk.GetX()[n]*spcol[n]*yy[n]*dt;
    cellMass[n] = 2.*M_PI*disk.GetX()[n]*disk.GetX()[n]*disk.GetDlnx()*spcol[n];
  }

  for(unsigned int n=1; n<=dZh.size()-1; ++n) {
    (*this).spcol[n] += dcolh[n]*dt;
    (*this).spsig[n] += dsigh[n]*dt;
    (*this).spZ[n] += dZh[n]*dt ;

    double wt1 = cellMass[n] - outgoingMass[n];
    double wt2 = incomingMass[n];
    double avg1 = spZ[n];
    double avg2 = incomingZ[n];
    double var1 = spZV[n];
    double var2 = incomingZV[n];
    double wtdAvg = (wt1*avg1 + wt2*avg2)/(wt1+wt2);
    spZV[n] = wt1/(wt1+wt2)* (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
      + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);

  }
}
