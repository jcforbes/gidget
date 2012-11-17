#include "StellarPop.h"
#include "Cosmology.h"
#include "DiskContents.h"
#include "DiskUtils.h"
#include "FixedMesh.h"
#include "Debug.h"

#include <math.h>
#include <vector>
#include <iostream>

// typical constructor
StellarPop::StellarPop(FixedMesh & m) :
  spcol(std::vector<double>(m.nx()+1,0.)), 
  spsigR(std::vector<double>(m.nx()+1,0.)),
  spsigZ(std::vector<double>(m.nx()+1,0.)),
  dSigRdr(std::vector<double>(m.nx()+1,0.)),
  dSigZdr(std::vector<double>(m.nx()+1,0.)),
  dColdr(std::vector<double>(m.nx()+1,0.)),
  spZ(std::vector<double>(m.nx()+1,0.)), 
  spZV(std::vector<double>(m.nx()+1,0.)),
  dQdS(std::vector<double>(m.nx()+1,0.)),
  dQdsR(std::vector<double>(m.nx()+1,0.)), 
  dQdsZ(std::vector<double>(m.nx()+1,0.)),
  dQdSerr(std::vector<double>(m.nx()+1,0.)),
  dQdserr(std::vector<double>(m.nx()+1,0.)),
  ageAtz0(-1.0),
  mesh(m)
{ }


void StellarPop::ComputeSpatialDerivs()
{
  std::vector<double> & x = mesh.x();
  unsigned int nx = mesh.nx(); 
  for(unsigned int n=1; n<=nx; ++n) {
    dSigRdr[n] = ddx(spsigR,n,x,false);
    dSigZdr[n] = ddx(spsigZ,n,x,false);
    dColdr[n] = ddx(spcol,n,x,false);
  }
}



//bool StellarPop::IsForming(Cosmology& cos, double redshift)
//{
//  double age = cos.lbt(redshift);
//  return (age > youngest && age<=oldest);
//}

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
      double wtdAvg = (spZ[i]*spcol[i] + sp2.spZ[i]*sp2.spcol[i])/(spcol[i] + sp2.spcol[i]);
      double avg1 = spZ[i];
      double avg2 = sp2.spZ[i];
      double var1 = spZV[i];
      double var2 = sp2.spZV[i];
      double wt1 = spcol[i];
      double wt2 = sp2.spcol[i];
      
      // Merge the two distributions:
      (*this).spZ[i] = wtdAvg;
 //     (*this).spZV[i] = wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
 //       + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);
      spZV[i] = ComputeVariance(spcol[i],0.0, sp2.spcol[i],
				spZ[i],sp2.spZ[i],spZV[i],sp2.spZV[i]);
      
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

void StellarPop::MigrateStellarPop(double dt, double ** tauvecStar, DiskContents& disk, std::vector<double>& MdotiPlusHalf)
{
  // A few convenience vectors to store data before updating the state variables.
  // These hold derivatives:
  std::vector<double> dcoldt(spcol.size());
  std::vector<double> dsigRdt(spcol.size());
  std::vector<double> dZdt(spcol.size());
  //std::vector<double> MdotiPlusHalf(spcol.size());
  // These refer back to mesh variables referenced by the disk object:
  std::vector<double>& uu = disk.GetUu();
  std::vector<double>& beta = disk.GetBeta();
  std::vector<double>& x = disk.GetX();
  // FixedMesh & mesh = disk.GetMesh();
  // These store information about the metal fluxes to calculate the new variance of Z, spZV.
  std::vector<double> incomingMass(spcol.size(),0.0);
  std::vector<double> outgoingMass(spcol.size(),0.0);
  std::vector<double> incomingZ(spcol.size());
  std::vector<double> incomingZV(spcol.size());
  std::vector<double> cellMass(spcol.size());
  Debug& dbg = disk.GetDbg();
  double spsigZp1,spsigRp1,spcolp1;
  unsigned int nx=spcol.size()-1;

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
  //  tauvecStar[2][n] = (tauvecStar[1][n+1]-tauvecStar[1][n-1])/(mesh.x(n+1)-mesh.x(n-1));
  //  MdotiPlusHalf[n] = -1.0/mesh.u1pbPlusHalf(n) * (tauvecStar[1][n+1]-tauvecStar[1][n])/(mesh.x(n+1)-x[n]);
  }
//  MdotiPlusHalf[0]= -1.0/mesh.u1pbPlusHalf(0) * (tauvecStar[1][1]-tauvecStar[1][0])/(x[1]-mesh.x(0.0));
  
  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    double f = (spcol[n]/disk.activeColSt(n));
    dcoldt[n] = f*(MdotiPlusHalf[n]-MdotiPlusHalf[n-1])/(x[n]*mesh.dx(n));
    double MdotCentered = (-tauvecStar[2][n]*f
		          /(uu[n]*(1+beta[n]))); 
    if(dbg.opt(10)) MdotCentered=(-tauvecStar[2][n]*f/mesh.u1pbPlusHalf(n));
    if(n<nx) {
      spsigZp1=spsigZ[n+1];
      spsigRp1=spsigR[n+1];
      spcolp1 =spcol[n+1];
    }
    else {
      spsigZp1=spsigZ[nx];
      spsigRp1=spsigR[nx];
      spcolp1 =spcol[nx];
    }

    dsigRdt[n] = MdotCentered* 
      (1.0/(x[n]*spcol[n]*(spsigR[n] + spsigZ[n]))) *
      (2.0*spsigZ[n]* dSigZdr[n] //ddx(spsigZ,n,x,false)
       + 3.0* spsigR[n]* dSigRdr[n] //ddx(spsigR,n,x,false) 
       + spsigR[n]*spsigR[n]/spcol[n]* dColdr[n] //ddx(spcol,n,x,false)
       + (spsigR[n]*spsigR[n] - spsigZ[n]*spsigZ[n])/x[n]);
    dZdt[n] =  MdotCentered*ddx(spZ,n,x,true)/(x[n]*spcol[n]);



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
    if(spcol[n]<0 || spsigR[n]<0 || spsigZ[n]<0 || spZ[n]<0 || spZV[n]<0 ||
       spcol[n]!=spcol[n] || spsigR[n]!=spsigR[n] || spsigZ[n]!=spsigZ[n] ||
       spZ[n]!=spZ[n] || spZV[n]!=spZV[n])
       errormsg("Migrating the populations has produced nonsensical results!");
  }
}




double ComputeVariance(double cellMass, double outgoingMassINPUT, double incomingMass, 
                       double Z, double incomingZ, double ZV, double incomingZV)
{
    double outgoingMass = outgoingMassINPUT;
    double wt1 = cellMass - outgoingMass;
    if(wt1 < 0.0) wt1=0.0;
    double wt2 = incomingMass;
    double avg1 = Z;
    double avg2 = incomingZ;
    double var1 = ZV;
    double var2 = incomingZV;
    double wtdAvg = (wt1*avg1 + wt2*avg2)/(wt1+wt2);
    double val;
    val =  wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
         + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);

    if(val >= 0.0) return val;
    else if(val <-1.0e-10)
      errormsg("Something has gone wrong in computing the variance of metallicity. We were given the following: cellMass, outgoingMass, incomingMass, Z, incomingZ, ZV, incoming ZV:  "+str(cellMass)+" "+str(outgoingMass)+" "+str(incomingMass)+" "+str(Z)+" "+str(incomingZ)+" "+str(ZV)+" "+str(incomingZV));
    else return 0.0;
}

