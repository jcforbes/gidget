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
  spcol(std::vector<double>(nx+1,0.)), spsig(std::vector<double>(nx+1,0.)),
  spZ(std::vector<double>(nx+1,0.)), spZV(std::vector<double>(nx+1,0.)),
  dQdS(std::vector<double>(nx+1,0.)),
  dQds(std::vector<double>(nx+1,0.)), dQdSerr(std::vector<double>(nx+1,0.)),
  dQdserr(std::vector<double>(nx+1,0.)),
  youngest(y),oldest(o),ageAtz0((o+y)/2.0)
{ }

// default constructor
StellarPop::StellarPop() :  
  spcol(std::vector<double>(1,0.)), 
  spsig(std::vector<double>(1,0.)),
  spZ(std::vector<double>(1,0.)), 
  spZV(std::vector<double>(1,0.)),
  dQdS(std::vector<double>(1,0.)),
  dQds(std::vector<double>(1,0.)), 
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
  spsig(copy_from.GetSpSig()),
  spZ(copy_from.GetSpZ()),
  spZV(copy_from.GetSpZV()),
  dQdS(copy_from.GetdQdS()),
  dQds(copy_from.GetdQds()),
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
  spsig.resize(nxp1);
  spZ.resize(nxp1);
  spZV.resize(nxp1);
  dQdS.resize(nxp1);
  dQds.resize(nxp1);
  dQdSerr.resize(nxp1);
  dQdserr.resize(nxp1);

  ageAtz0 = copy_from.GetAgeAtz0();
  youngest = copy_from.GetYoungest();
  oldest=copy_from.GetOldest();
  spcol = copy_from.GetSpCol();
  spsig = copy_from.GetSpSig();
  spZ = copy_from.GetSpZ();
  spZV = copy_from.GetSpZV();
  dQdS = copy_from.GetdQdS();
  dQds = copy_from.GetdQds();
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

void StellarPop::MigrateStellarPop(double dt, double ** tauvecStar, DiskContents& disk)
{
  // A few convenience vectors to store data before updating the state variables.
  // These hold derivatives:
  std::vector<double> dcoldt(spcol.size());
  std::vector<double> dsigRdt(spcol.size());
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
                 (2.0*spsigZ[n]*ddx(spsigZ,n)
                   + 3.0* spsigR[n]*ddx(spsigR,n) 
                   + spsigR[n]*spsigR[n]*ddx(spcol,n)/spcol[n] 
                   + (spsigR[n]*spsigR[n] - spsigZ[n]*spsigZ[n])/x[n]);
    dZdt[n] =  MdotCentered*ddx(spZ,n,x)/(x[n]*spcol[n]);



    // Now we proceed to do what looks like a ridiculous amount of work to compute spZV.
    cellMass[n] = spcol[n]*x[n]*mesh.dx(n);
    bool fromRight = MdotiPlusHalf[n]>0.0;
    bool fromLeft = MdotiPlusHalf[n-1]<0.0;
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

  for(unsigned int =1; n<=spcol.size()-1; ++n) {
    spcol[n] += dcoldt[n]*dt;
    spsigR[n] += dsigRdt[n]*dt;
    spsigZ[n] += .5*dsigRdt[n]*dt;
    spZ[n] += dZdt[n]*dt;
    spZV[n] = ComputeVariance(cellMass[n],outgoingMass[n],incomingMass[n],
                              spZ[n],incomingZ[n],spZV[n],incomingZV[n]);
  }
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

  std::vector<double> flux(dcolh.size());
  for(unsigned int n=1; n<=spcol.size()-2; ++n) {
    double ym = yy[n+1];
    if(fabs(ym) > fabs(yy[n])) ym=yy[n]; //minmod
    if(yy[n] * yy[n+1] <= 0.0) ym=0.0;
    double cst = spcol[n+1];
    if(ym > 0.0) cst = spcol[n+1];

    flux[n] = 2.0*M_PI*sqrt(disk.GetX()[n]*disk.GetX()[n+1]) * ym * cst ;
    if(n==1 && flux[n]>1.0) {
	flux[n]=0.0; // no stars can come out of the bulge.
//	std::cerr << "Warning: positive flux of stars out of the bulge. flux= "+str(flux[n]) << std::endl;
    }
  }
  flux[spcol.size()-1]=0.0;

//  double sum1=0;
//  double sum2=0;
//  double sum3=0;
  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
//    dcolh[n] = -2.*M_PI*((*this).spcol[n]*ddx(yy,n,disk.GetX()) + ddx((*this).spcol,n,disk.GetX()) *yy[n] + (*this).spcol[n]*yy[n]/disk.GetX()[n]);
//    dsigh[n] = -2.*M_PI*yy[n]*((1.+disk.GetBeta()[n])*disk.GetUu()[n]*disk.GetUu()[n]/(3.*(*this).spsig[n]*disk.GetX()[n]) + ddx((*this).spsig,n,disk.GetX()));
    if(n<spcol.size()-1) {
     double sigp2 = (2./3.) * (disk.GetMesh().psi(disk.GetX()[n+1])-disk.GetMesh().psi(disk.GetX()[n])) + (1./3.) * (disk.GetUu()[n+1]*disk.GetUu()[n+1]-disk.GetUu()[n]*disk.GetUu()[n]) + spsig[n+1]*spsig[n+1];
     dsigh[n] = -2.0*M_PI/ (2.0*disk.GetX()[n]*disk.GetX()[n]*disk.GetDlnx()*spcol[n]*spsig[n]) * (disk.GetX()[n+1]*yy[n+1]*spcol[n+1]*(sigp2-spsig[n]*spsig[n]));
    }
    else dsigh[n] = 0.0;
     

    dZh[n] = -2.*M_PI*yy[n]*ddx((*this).spZ,n,disk.GetX());
    
    //    dcolh[n] = (flux[n-1]-flux[n])  / (disk.GetX()[n] * disk.GetX()[n] * disk.GetDlnx());
    dcolh[n] = dSMigdt(n,yy,disk.GetX(),(*this).spcol);


//    sum1+=-2.*M_PI*((*this).spcol[n]*ddx(yy,n,disk.GetX())) * disk.GetX()[n] * disk.GetX()[n];
//    sum2+=-2.*M_PI*ddx((*this).spcol,n,disk.GetX())*yy[n] * disk.GetX()[n] * disk.GetX()[n];
//    sum3+=-2.*M_PI*(*this).spcol[n]*yy[n]/disk.GetX()[n] * disk.GetX()[n] * disk.GetX()[n];

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

  //  if(sum1+sum2+sum3>0.0) {
  //    errormsg("Error migrating the population: mass was artificially created! sum1 sum2 sum3 sum: "+str(sum1)+" "+str(sum2)+" "+str(sum3)+" "+str(sum1+sum2+sum3));
  //  }

  for(unsigned int n=1; n<=dZh.size()-1; ++n) {
//    (*this).spcol[n] += dcolh[n]*dt;
    (*this).spsig[n] += dsigh[n]*dt;
    if(spsig[n] < disk.GetMinSigSt()*.9999999) 
	errormsg("Sigst below floor!");
    (*this).spZ[n] += dZh[n]*dt ;
    (*this).spcol[n] += dSMigdt(n,yy,disk.GetX(),(*this).spcol)*dt;
    spZV[n] = ComputeVariance(cellMass[n],outgoingMass[n],incomingMass[n],
                              spZ[n],incomingZ[n],spZV[n],incomingZV[n]);

  }
}


double ComputeVariance(double cellMass, double outgoingMass, double incomingMass, 
                       double Z, double incomingZ, double ZV, double incomingZV)
{
    double wt1 = cellMass - outgoingMass;
    double wt2 = incomingMass;
    double avg1 = spZ;
    double avg2 = incomingZ;
    double var1 = spZV;
    double var2 = incomingZV;
    double wtdAvg = (wt1*avg1 + wt2*avg2)/(wt1+wt2);
    return wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
         + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);

}

