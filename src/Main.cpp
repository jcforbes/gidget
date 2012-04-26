#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>

#include "Simulation.h"
#include "DiskContents.h"
#include "AccretionHistory.h"
#include "ArgumentSetter.h"
#include "Cosmology.h"
#include "Dimensions.h"
#include "DiskUtils.h"

/*
  This is the main method for the program GIDGET, the "Gravitational Instability-
  Dominated Galaxy Evolution Tool". The structure of the code is
	- Initialize variables from command line arguments
	- Run a simulation where the stars do nothing until the gas has converged
	- Use the end configuration of that simulation to initialize a new simulation.
	- Run this simulation and output the corresponding data.

  For more information on the internal workings of the simulation, consult the header files
  included above.
*/

int main(int argc, char **argv) { 
  void solvde(int itmax, double conv, double slowc, double scalv[], int indexv[],
	      int ne, int nb, int m, double **y, double ***c, double **s,int*);
  std::cout.precision(7);  
  if(argc<2) {
    std::vector<std::string> errmgLine(0);
    errmgLine.push_back("Need at least one argument (filename)-- Arguments are ");
    errmgLine.push_back(" \nrunName, \nnx, \neta, \nepsff, \nTmig (local orbital times), \n");
    errmgLine.push_back("analyticQ (1 or 0),\ncosmologyOn (1 or 0), \nxmin, \nNActive, \n");
    errmgLine.push_back("NPassive, \nvphiR (km/s), \nradius (kpc), \ngasTemp (K), \nQlim, \n");
    errmgLine.push_back("fg0, \ntempRatio (sig_*/sig_g), \nzstart, \ntmax, \nstepmax, \n");
    errmgLine.push_back("TOL (t_orb),\nMassLoadingFactor, \nBulgeRadius (kpc), \n");
    errmgLine.push_back("stDiskScale (kpc, or -1 for powerlaw),\nwhichAccretionHistory,\nalphaMRI");
    errmgLine.push_back(", \nthick,\nmigratePassive,\nQinit,\nkappaMetals,\nMh0\n");
    std::string msg="";
    for(unsigned int i=0; i!=errmgLine.size();++i) {
      msg+=errmgLine[i];
    }
    errormsg(msg);
  }

  std::string filename(argv[1]);
  
  // Set various parameters based on the command line arguments. If an argument is not
  // specified on the command line, set that parameter to the default value (first 
  // argument of as.Set())
  ArgumentSetter as(argc,argv,filename);
  const unsigned int nx=           as.Set(500,"nx");
  const double eta=                as.Set(1.5,"eta");
  const double epsff=              as.Set(.01,"eps_ff");
  const double tauHeat=            as.Set(2,"heating timescale (outer orbits)");
  const bool analyticQ=           (as.Set(1,"analytic Q")==1);
  const bool cosmologyOn=         (as.Set(1,"cosmological accretion history")==1);
  const double xmin=               as.Set(.01,"inner truncation radius (dimensionless)");
  const unsigned int NActive=      as.Set(1,"N Active Stellar Populations");
  const unsigned int NPassive=     as.Set(10,"N Passive Stellar Populations");
  const double vphiRatMh12 =       as.Set(220,"Circular velocity (km/s) for a 10^12 MSun halo")*1.e5;
  const double radius=             as.Set(20.,"Outer Radius (kpc)")*cmperkpc;
//  const double sigth=         sqrt(as.Set(7000,"Gas Temperature (K)")*kB/mH)/vphiR;
  const double Tgas =		   as.Set(7000,"Gas Temperature (K)");
  const double Qlim =              as.Set(2.5,"Limiting Q_*");
  const double fg0  =              as.Set(.5,"Initial gas fraction");
  const double tempRatio =         as.Set(1.,"Initial sigma_*/sigma_g");
  const double zstart =            as.Set(2.,"Initial redshift");
  const double tmax =              as.Set(1000,"Maximum Time (outer orbits)");
  const unsigned int stepmax=      as.Set(100000000,"Maximum Number of Steps");
  const double TOL =               as.Set(.0001,"TOL (outer orbits)");
  const double MassLoadingFactorAtMh12=  as.Set(1,"Mass Loading Factor at Mh=10^12 MSun");
  const double BulgeRadius      =  as.Set(0,"Velocity Curve Turnover Radius (kpc)");
  const double stScaleLength    =  as.Set(-1,"Initial Stellar Disk Scale Length (kpc)");
  const int whichAccretionHistory= as.Set(0,"Which Accretion History- 0-Bouche, 1-High, 2-Low");
  const double alphaMRI         =  as.Set(0,"alpha viscosity for the MRI");
  const double thick =             as.Set(1.5,"Thickness correction to Q");
  const bool migratePassive=      (as.Set(1,"Migrate Passive population")==1);
  const double Qinit =             as.Set(2.0,"The fixed Q");
  const double kappaMetals =       as.Set(.001,"Kappa Metals");
  const double Mh0 =  		   as.Set(1.0e12,"Halo Mass");

  // Scale the things which scale with halo mass.
  const double MassLoadingFactor = MassLoadingFactorAtMh12 * pow((Mh0/1.0e12) , -1./3.);
  const double vphiR = vphiRatMh12 * pow(Mh0/1.0e12,  1./4.);
  const double sigth = sqrt(Tgas *kB/mH)/vphiR;

  // Make an object to deal with things cosmological
  Cosmology cos(1.-.734, .734, 2.29e-18 ,zstart);

  // Make an object to deal with the accretion history
  AccretionHistory accr;
  double mdot0;

  // Based on the user's choice of accretion history, generate the appropriate
  // mapping between redshift and accretion rate.
  if(whichAccretionHistory==0)
    mdot0 = accr.GenerateBoucheEtAl2009(Mh0*1.0e-12,2.0,cos,filename+"_Bouche09.dat",true,true) * MSol/speryear;
  else if(whichAccretionHistory==2)
    mdot0 = accr.GenerateConstantAccretionHistory(2.34607,zstart,cos,filename+"_ConstAccHistory.dat",true) * MSol/speryear;
  else if(whichAccretionHistory==1)
    mdot0 = accr.GenerateConstantAccretionHistory(12.3368,zstart,cos,filename+"_ConstAccHistory2.dat",true)*MSol/speryear;
  else
    mdot0 = accr.GenerateNeistein08(Mh0,2.0,cos,filename+"_Neistein08_"+str(whichAccretionHistory)+".dat",true,whichAccretionHistory)*MSol/speryear;
  // Note that the following line does nothing but put a line in the comment file to
  // record MdotExt0 for this run.
  as.Set(mdot0/MSol*speryear,"Initial Accretion (MSol/yr)");

  // Done reading in arguments. Write out a comment file containing all of the arguments.
  as.~ArgumentSetter();
    

  // Set the dimensional quantities. 
  Dimensions dim(radius,vphiR,mdot0);

  //// Evolve a disk where the stars do not do anything and Mdot_ext=Mdot_ext,0.
  DiskContents diskIC(nx,xmin,1.0e30,eta,sigth,0.0,Qlim, // need Qlim to successfully set initial statevars
		      TOL,analyticQ,MassLoadingFactor,cos,dim,
		      thick, false,Qinit,kappaMetals);
  if(stScaleLength<0.0)  diskIC.Initialize(tempRatio,fg0,NActive,NPassive,BulgeRadius);
  else diskIC.Initialize(0.1*Z_Sol, .6, fg0, 50.0/220.0, Mh0, NActive, NPassive, BulgeRadius, stScaleLength);

  Simulation simIC(10.0,100000000,
                   false, nx,TOL,
                   zstart,NActive,NPassive,
                   alphaMRI,sigth,
                   diskIC,accr);
  int result = simIC.runToConvergence(1, false, filename+"_icgen"); // set false-> true to debug initial condition generator
  if(result!=5) // The simulation converges when the time step reaches 1*TOL.
    errormsg("Initial Condition generator failed to converge, code "+str(result));

  simIC.GetInitializer().BulgeRadius = BulgeRadius;

  // Now evolve a disk where the stars evolve as they should using the previous simulation's end condition
  DiskContents disk(nx,xmin,tauHeat,eta,sigth,epsff,Qlim,
		    TOL,analyticQ,MassLoadingFactor,cos,dim,
		    thick,migratePassive,Qinit,kappaMetals);
  disk.Initialize(simIC.GetInitializer(), stScaleLength < 0.0); // if we're using an exponential disk, don't mess with the initial conditions of the stellar disk when enforcing Q=Q_f, i.e. do not keep a fixed phi0.
  Simulation sim(tmax,stepmax,
		 cosmologyOn,nx,TOL,
		 zstart,NActive,NPassive,
		 alphaMRI,sigth,
		 disk,accr);
  sim.runToConvergence(1.0e10, true, filename); 
}
