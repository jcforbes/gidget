#include "StellarPop.h"
#include "Cosmology.h"
#include "DiskContents.h"
#include "DiskUtils.h"
#include "FixedMesh.h"

#include <math.h>
#include <vector>
#include <iostream>

// typical constructor
StellarPop::StellarPop(unsigned int nx,double y,double o) :
  spcol(std::vector<double>(nx+1,0.)), spsigR(std::vector<double>(nx+1,0.)),
  spsigZ(std::vector<double>(nx+1,0.)),
  spZ(std::vector<double>(nx+1,0.)), spZV(std::vector<double>(nx+1,0.)),
  dQdS(std::vector<double>(nx+1,0.)),
  dQdsR(std::vector<double>(nx+1,0.)), 
  dQdsZ(std::vector<double>(nx+1,0.)),
  dQdSerr(std::vector<double>(nx+1,0.)),
  dQdserr(std::vector<double>(nx+1,0.)),
  youngest(y),oldest(o),ageAtz0((o+y)/2.0)
{ }

// default constructor
StellarPop::StellarPop() :  
  spcol(std::vector<double>(1,0.)), 
  spsigR(std::vector<double>(1,0.)),
  spsigZ(std::vector<double>(1,0.)),
  spZ(std::vector<double>(1,0.)), 
  spZV(std::vector<double>(1,0.)),
  dQdS(std::vector<double>(1,0.)),
  dQdsR(std::vector<double>(1,0.)), 
  dQdsZ(std::vector<double>(1,0.)),
  dQdSerr(std::vector<double>(1,0.)),
  dQdserr(std::vector<double>(1,0.)),
  youngest(-1.),oldest(-1.),ageAtz0(-1.)      
{ }

//// copy constructor
//StellarPop::StellarPop(const StellarPop& copy_from) :
//  spcol(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  spsig(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  spZ(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  spZV(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  dQdS(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  dQds(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  dQdSerr(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  dQdserr(std::vector<double>(copy_from.GetSpCol().size(),0.)),
//  youngest(copy_from.GetYoungest()),
//  oldest(copy_from.GetOldest()),
//  ageAtz0(copy_from.GetAgeAtz0())
//{
//  for(unsigned int n=1; n<=spcol.size(); ++n) {
//    spcol[n] = copy_from.GetSpCol()[n];
//    spsig[n] = copy_from.GetSpSig()[n];
//    spZ[n] = copy_from.GetSpZ()[n];
//    spZV[n] = copy_from.GetSpZV()[n];
//    dQdS[n] = copy_from.GetdQdS()[n];
//    dQds[n] = copy_from.GetdQds()[n];
//    dQdSerr[n] = copy_from.GetdQdSerr()[n];
//    dQdserr[n] = copy_from.GetdQdserr()[n];
//  }
//}

// copy constructor
StellarPop::StellarPop(const StellarPop& copy_from) :
  spcol(copy_from.GetSpCol()),
  spsigR(copy_from.GetSpSigR()),
  spsigZ(copy_from.GetSpSigZ()),
  spZ(copy_from.GetSpZ()),
  spZV(copy_from.GetSpZV()),
  dQdS(copy_from.GetdQdS()),
  dQdsR(copy_from.GetdQdsR()),
  dQdsZ(copy_from.GetdQdsZ()),
  dQdSerr(copy_from.GetdQdSerr()),
  dQdserr(copy_from.GetdQdserr()),
  youngest(copy_from.GetYoungest()),
  oldest(copy_from.GetOldest()),
  ageAtz0(copy_from.GetAgeAtz0())
{
  return;
}

// assignment operator
StellarPop & StellarPop::operator= (const StellarPop& copy_from) 
{
  const unsigned int nxp1= copy_from.GetSpCol().size();
  spcol.resize(nxp1);
  spsigR.resize(nxp1);
  spsigZ.resize(nxp1);
  spZ.resize(nxp1);
  spZV.resize(nxp1);
  dQdS.resize(nxp1);
  dQdsR.resize(nxp1);
  dQdsZ.resize(nxp1);
  dQdSerr.resize(nxp1);
  dQdserr.resize(nxp1);

  ageAtz0 = copy_from.GetAgeAtz0();
  youngest = copy_from.GetYoungest();
  oldest=copy_from.GetOldest();
  spcol = copy_from.GetSpCol();
  spsigR = copy_from.GetSpSigR();
  spsigZ = copy_from.GetSpSigZ();
  spZ = copy_from.GetSpZ();
  spZV = copy_from.GetSpZV();
  dQdS = copy_from.GetdQdS();
  dQdsR = copy_from.GetdQdsR();
  dQdsZ = copy_from.GetdQdsZ();
  dQdSerr = copy_from.GetdQdSerr();
  dQdserr = copy_from.GetdQdserr();
//
//  for(unsigned int n=1; n<=spcol.size(); ++n) {
//    spcol[n] = copy_from.GetSpCol()[n];
//    spsig[n] = copy_from.GetSpSig()[n];
//    spZ[n] = copy_from.GetSpZ()[n];
//    spZV[n] = copy_from.GetSpZV()[n];
//    dQdS[n] = copy_from.GetdQdS()[n];
//    dQds[n] = copy_from.GetdQds()[n];
//    dQdSerr[n] = copy_from.GetdQdSerr()[n];
//    dQdserr[n] = copy_from.GetdQdserr()[n];
//  }
  return *this;
}


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
      spsigR[i]= sqrt(((spcol[i])*((*this).spsigR[i])*((*this).spsigR[i]) 
                              + (sp2.spcol[i])*(sp2.spsigR[i])*(sp2.spsigR[i]))
                             /((*this).spcol[i]+sp2.spcol[i]));
      spsigZ[i]= sqrt((spcol[i]*spsigZ[i]*spsigZ[i] 
		       + sp2.spcol[i]*sp2.spsigZ[i]*sp2.spsigZ[i])
		      /(spcol[i]+sp2.spcol[i]));
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
      
      if((*this).spsigR[i]!=(*this).spsigR[i] || spsigZ[i]!=spsigZ[i])
        errormsg("Error merging populations:  this spcol,spsig  sp2 spcol,spsig  "+str((*this).spcol[i])+" "+str((*this).spsigR[i])+"  "+str(sp2.spcol[i])+" "+str(sp2.spsigR[i]));
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
    spsigR[n] = sp2.spsigR[n];
    spsigZ[n] = sp2.spsigZ[n];
    spcol[n] = frac*sp2.spcol[n];
    spZ[n] = sp2.spZ[n];
    spZV[n] = sp2.spZV[n];
    sp2.spcol[n] -= spcol[n];;
  }
  ageAtz0 = sp2.ageAtz0;

}

void StellarPop::MigrateStellarPop(double dt, double ** tauvecStar, DiskContents& disk)
{
  // A few convenience vectors to store data before updating the state variables.
  // These hold derivatives:
  std::vector<double> dcoldt(spcol.size());
  std::vector<double> dsigRdt(spcol.size());
  std::vector<double> dZdt(spcol.size());
  std::vector<double> MdotiPlusHalf(spcol.size());
  // These refer back to mesh variables referenced by the disk object:
  std::vector<double>& uu = disk.GetUu();
  std::vector<double>& beta = disk.GetBeta();
  std::vector<double>& x = disk.GetX();
  FixedMesh & mesh = disk.GetMesh();
  // These store information about the metal fluxes to calculate the new variance of Z, spZV.
  std::vector<double> incomingMass(spcol.size(),0.0);
  std::vector<double> outgoingMass(spcol.size(),0.0);
  std::vector<double> incomingZ(spcol.size());
  std::vector<double> incomingZV(spcol.size());
  std::vector<double> cellMass(spcol.size());

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    tauvecStar[2][n] = (tauvecStar[1][n+1]-tauvecStar[1][n-1])/(mesh.x(n+1)-mesh.x(n-1));
    MdotiPlusHalf[n] = -1.0/mesh.u1pbPlusHalf(n) * (tauvecStar[1][n+1]-tauvecStar[1][n])/(mesh.x(n+1)-x[n]);
  }
  MdotiPlusHalf[0]= -1.0/mesh.u1pbPlusHalf(0) * (tauvecStar[1][1]-tauvecStar[1][0])/(x[1]-mesh.x(0.0));
  
  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    dcoldt[n] = (MdotiPlusHalf[n]-MdotiPlusHalf[n-1])/(x[n]*mesh.dx(n));
    double MdotCentered = (-tauvecStar[2][n]*(spcol[n]/disk.activeColSt(n))
		          /(uu[n]*(1+beta[n]))); // FOR THIS COMPONENT (note the Sigma_*,i/Sigma_* term) 
    dsigRdt[n] = MdotCentered* 
      (1.0/(x[n]*spcol[n]*(spsigR[n] + spsigZ[n]))) *
      (2.0*spsigZ[n]*ddx(spsigZ,n,x)
       + 3.0* spsigR[n]*ddx(spsigR,n,x) 
       + spsigR[n]*spsigR[n]*ddx(spcol,n,x)/spcol[n] 
       + (spsigR[n]*spsigR[n] - spsigZ[n]*spsigZ[n])/x[n]);
    dZdt[n] =  MdotCentered*ddx(spZ,n,x)/(x[n]*spcol[n]);



    // Now we proceed to do what looks like a ridiculous amount of work to compute spZV.
    cellMass[n] = spcol[n]*x[n]*mesh.dx(n);
    bool fromRight = MdotiPlusHalf[n]>0.0;
    bool fromLeft = MdotiPlusHalf[n-1]<0.0;
    unsigned nx = spcol.size()-1;
    double spZp1, spZVp1, spZm1, spZVm1;
    if(n<nx) { 
      spZp1=spZ[n+1];
      spZVp1=spZ[n+1];
    }
    else { // these values shouldn't matter in theory.
      spZp1 = 0.0;
      spZVp1 = 0.0;
    }
    if(n>1) {
      spZm1=spZ[n-1];
      spZVm1=spZV[n-1];
    }
    else { // again, these values should not affect the calculation.
      spZm1= 0.0;
      spZVm1= 0.0;
    }
    if(fromRight) {
      incomingMass[n] += MdotiPlusHalf[n]*dt;
      if( !fromLeft) {
        incomingZ[n] = spZ[n+1];
        incomingZV[n] = spZV[n+1];
      }
    }
    if(!fromRight) {
      outgoingMass[n] -= MdotiPlusHalf[n]*dt; 
    }
    if(fromLeft) {
      incomingMass[n] -= MdotiPlusHalf[n-1]*dt;
      if(n>1 && !fromRight) {
        incomingZ[n] = spZ[n-1];
        incomingZV[n] = spZV[n-1];
      }
    }
    if(!fromLeft) {
      outgoingMass[n] += MdotiPlusHalf[n-1]*dt;
    }
    if(fromLeft && fromRight) {
      incomingZ[n] = (MdotiPlusHalf[n]*dt*spZp1 - MdotiPlusHalf[n-1]*dt*spZm1)/incomingMass[n];
      incomingZV[n] = ComputeVariance(MdotiPlusHalf[n]*dt,0.0,-MdotiPlusHalf[n-1]*dt,
                                      spZp1, spZm1, spZVp1, spZVm1);
    }

  }

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    spcol[n] += dcoldt[n]*dt;
    spsigR[n] += dsigRdt[n]*dt;
    spsigZ[n] += .5*dsigRdt[n]*dt;
    spZ[n] += dZdt[n]*dt;
    spZV[n] = ComputeVariance(cellMass[n],outgoingMass[n],incomingMass[n],
                              spZ[n],incomingZ[n],spZV[n],incomingZV[n]);
  }
}






double ComputeVariance(double cellMass, double outgoingMass, double incomingMass, 
                       double Z, double incomingZ, double ZV, double incomingZV)
{
    double wt1 = cellMass - outgoingMass;
    double wt2 = incomingMass;
    double avg1 = Z;
    double avg2 = incomingZ;
    double var1 = ZV;
    double var2 = incomingZV;
    double wtdAvg = (wt1*avg1 + wt2*avg2)/(wt1+wt2);
    return wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
         + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);

}

