#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "Simulation.h"
#include "DiskContents.h"
#include "AccretionHistory.h"
#include "ConvergenceCheck.h"
#include "DiskUtils.h"
#include "Cosmology.h"
#include "Debug.h"
#include "FixedMesh.h"
#include "Dimensions.h"

Simulation::Simulation(const double tm, const long int sm,
                       const bool co,   const unsigned int nnx,
		       const double tl, const double zs, 
		       const unsigned int na, const unsigned int np,
                       const double amri,const double sth, const double nd,
                       DiskContents& td,
                       AccretionHistory& ah, Debug& db, Dimensions& dm):
  tmax(tm), stepmax(sm),
  cosmologyOn(co),nx(nnx),
  TOL(tl), zstart(zs),
  NActive(na),NPassive(np),
  alphaMRI(amri),sigth(sth),
  ndecay(nd),
  theDisk(td),
  accr(ah),
  dbg(db),
  dim(dm)
{
  ini.NActive=na;
  ini.NPassive=np;
  return;
}

// Terminate when dt>TOL * fCondition
int Simulation::runToConvergence(const double fCondition, 
				 const bool writeOut, 
				 const std::string filename)
{
    // Initialize a 2 by nx vector to store the torque and its first derivative.
  double **tauvec;
  tauvec = new double *[3];
  for(unsigned int k=1; k<=2; ++k) {
    tauvec[k] = new double[nx+2];
  }

  // This is the equilibrium solution to the torque eq from Krumholz & Burkert 2010
  for(unsigned int n=1; n<=nx; ++n) {
    tauvec[1][n]=-theDisk.GetX()[n];
    tauvec[2][n]=-1.;
  }

  // Record various quantities at various time intervals at various
  // locations in the disk to check their convergence when certain quantities
  // are varied (e.g. TOL, xmin, nx)
  ConvergenceCheck DtoB4(4.,4.2,0.,0.), DtoB25(25.,25.2,0.,0.),
    ctrMdot4(4.,4.2,.4,.45),   ctrMdot25(25.,25.2,.4,.45),
    outMdot4(4.,4.2,.9,.95),   outMdot25(25.,25.2,.9,.95),
    ctrsig4(4.,4.2,.4,.45),    outsig4(4.,4.2,.9,.95),
    ctrsig25(25.,25.2,.4,.45), outsig25(25.,25.2,.9,.95),
    ctrSF4(4.,4.2,.4,.45),     ctrSF25(25.,25.2,.4,.45),
    outSF4(4.,4.2,.9,.95),     outSF25(25.,25.2,.9,.95); 

  // Initialize some variables which will change as the simulation progresses:
  // time (outer orbits), redshift, and timestep
  double t=0., z=zstart, dt=TOL*1.e-6; 

  unsigned int step=0; // number of elapsed timesteps

  // number of times data has been written to output
  int writeIndex=0;

  double AccRate=1.0; // normalized accretion rate

  // Stop the simulation if the timestep falls below this 
  // value (in general this means that something has gone 
  // terribly wrong and the simulation has no chance to 
  // run to completion)
  const double dtfloor=1.e-12; 

  // Loop over time steps until...
  while(t<tmax &&           // the simulation has run for too many orbits
	step<=stepmax &&    // the simulation has run for too many steps
	z>=0.0 &&           // the simulation has reached redshift zero
	dt>TOL*dtfloor &&   // the simulation has stalled
        dt<fCondition*TOL)  // the simulation has converged
  {

    // Every 10th of an outer orbit, write out the data from the simulation
    // Note that we'd like to use last step's timestep here to determine whether
    // to write out a file, since we would like to avoid a situation where a sudden
    // decrease in the size of dt causes us to skip a checkpoint
    // bool timeOut = (((floor(25.0*(t-dt)) < floor(25.0*t)) || step<2) && writeOut);
    //    bool timeOut = true; // write out at every time step. Use only if 
                               // things crash and burn very quickly
    double duration = theDisk.GetCos().lbt(zstart); //in seconds
//    double present = theDisk.GetCos().Tsim(z); // in seconds
//    double previous = theDisk.GetCos().Tsim(z+dz(dt,z,theDisk.GetCos(),theDisk.GetDim())); // in seconds
    double present = t * 2.0*M_PI*dim.Radius / dim.vphiR;
    double previous = (t-dt) * 2.0*M_PI*dim.Radius/dim.vphiR;

    bool timeOut = ((floor(200.0*previous/duration) < floor(200.0*present/duration) || step < 2) && writeOut) ;

    if(!cosmologyOn)
      timeOut= (((floor(25.0*(t-dt)) < floor(25.0*t)) || step<2) && writeOut);

    if(timeOut) {
      // Separate files for the active & passive stellar populations..
      theDisk.WriteOutStarsFile(filename+"_act",theDisk.active(),NActive,step);
      theDisk.WriteOutStarsFile(filename,theDisk.passive(),NPassive,step);
      
      theDisk.WriteOutStepFile(filename,accr,t,z,dt,step,tauvec);
      writeIndex++;

      std::cout << "Writing out file "<<writeIndex<<" at t = "<<present/(speryear*1.0e9)<<" Gyr, z= "<<z<<std::endl;
    }

    theDisk.ComputeColSFR();


    // Make sure that stars formed this time step have a stellar 
    // population to which they belong. If they don't, create
    // one and add it to the appropriate vector of stellar populations
    bool checkActiveStars = 
      theDisk.CheckStellarPops(dt, z, theDisk.active(),
			       NActive,true);
    bool checkPassiveStars = 
      theDisk.CheckStellarPops(dt, z, theDisk.passive(),
			       NPassive,false);
    
    // Compute the rate at which stars are moving inwards (yy)

    /*  // Uncomment this block of code and comment 'theDisk.ComputeY();' 
	// to set vr_*=vr_g artificially. Note, of course, that this will 
	// completely toss out the prescription that 
	// dQ_* /dt = max(Qlim-Q*,0)/(Tmig(2pi Omega)^-1)
	// I have not had success with this- weird things tend to 
	// happen near the center of the disk.
    for(unsigned int n=1;n<=nx;++n) {
      if(step==0)
	theDisk.GetYy()[n]=0.;
      else 
	theDisk.GetYy()[n] = tauvec[2][n] / (2.*PI*theDisk.GetX()[n]*theDisk.GetUu()[n]*theDisk.GetCol()[n]*(1.+theDisk.GetBeta()[n]));
    }
    */

    theDisk.ComputeY(ndecay);

    // Compute dQ/dA and its error, where A is a stand-in for every state variable
    theDisk.ComputePartials();

    // Set the non-dimensional accretion rate at this redshift. Used to set
    // the boundary conditions of the torque equation.
    if(cosmologyOn) AccRate = accr.AccOfZ(z);
    else AccRate = 1.;

    theDisk.ComputedSdTCos(AccRate);

    // Given the user's freedom to specify the accretion history, make sure
    // the accretion rate remains non-negative.
    if(AccRate<0.0) errormsg("Negative accretion history!");

    // Update the coefficients of the torque equation
    theDisk.UpdateCoeffs(z);

    int whichVar,whichCell; // store which variable and which cell limits the time step.

    // Solve the torque equation. The commented versions represent various choices
    // for the boundary conditions
    //    disk.ComputeTorques(tauvec,-1.*AccRate*xmin/pow(xmin,1./(1.-nx)),-1.*AccRate);
    //    disk.ComputeTorques(tauvec,0.,-1.*AccRate);
    double IBC = 2.0*M_PI*theDisk.GetX()[1]*theDisk.GetX()[1]*theDisk.GetCol()[1]*alphaMRI*sigth*theDisk.GetSig()[1]*(theDisk.GetBeta()[1]-1.0);
    if(dbg.opt(17)) {
      // IBC = -1.0*AccRate*theDisk.GetMesh().x(0.0); // b
      IBC = -1.0*.01*theDisk.GetMesh().x(0.0); // c
    }
    if(dbg.opt(18)) {
      IBC = -1.0*theDisk.GetMesh().x(0.0); // d
    }
    double OBC;
    if(dbg.opt(8))
        OBC=-1.0*theDisk.dmdtCosOuter(AccRate);
    else
        OBC=-1.0*AccRate;
    theDisk.ComputeGItorque(tauvec,IBC,OBC);
    //    disk.ComputeTorques(tauvec,-1.*AccRate*xmin,-1.*AccRate);

    // In situations where the alpha viscosity produces larger torques 
    // than the GI viscosity, use the alpha viscosity instead:
    theDisk.ComputeMRItorque(tauvec,alphaMRI,IBC,OBC,ndecay);

    // Given the solution to the torque equation, compute time 
    // derivatives of the state variables
    theDisk.ComputeDerivs(tauvec);

    // Given the derivatives, compute a time step over which none of the variables
    // change by too much. whichVar tells us which state variable is limiting the timestep
    // and whichCell tells us which cell is limiting the timestep. Both of these values
    // are printed every 5000 timesteps (see below).
    dt = theDisk.ComputeTimeStep(z,&whichVar,&whichCell); 

    // Every time step, check whether each of the convergence checks has a value
    // which needs to be updated.
    DtoB4.UpdateD(tauvec[1][1],t,dt,theDisk.GetDim()); DtoB25.UpdateD(tauvec[1][1],t,dt,theDisk.GetDim());  
    ctrMdot4.Update(tauvec[2],t,theDisk.GetX()); ctrMdot25.Update(tauvec[2],t,theDisk.GetX()); 
    outMdot4.Update(tauvec[2],t,theDisk.GetX()); outMdot25.Update(tauvec[2],t,theDisk.GetX()); 
    ctrsig4.Update(theDisk.GetSig(),t,theDisk.GetX()); outsig4.Update(theDisk.GetSig(),t,theDisk.GetX());
    ctrsig25.Update(theDisk.GetSig(),t,theDisk.GetX()); outsig25.Update(theDisk.GetSig(),t,theDisk.GetX()); 
    ctrSF4.Update(theDisk.GetColSFR(),t,theDisk.GetX()); ctrSF25.Update(theDisk.GetColSFR(),t,theDisk.GetX()); 
    outSF4.Update(theDisk.GetColSFR(),t,theDisk.GetX()); outSF25.Update(theDisk.GetColSFR(),t,theDisk.GetX()); 

    // Give the user a bit of information on the progress of the simulation
    if(step%5000==0 || step<5)
      std::cout << "Step "<<step<<", t="<<t<<" Outer Orbits, z="<<z<<", dt="<<dt
                <<", "<<whichVar<<", "<<whichCell<<std::endl;


    // And finally, update the state variables
    theDisk.UpdateStateVars(dt,z,tauvec,AccRate); 

    // update the independent variables.
    if(cosmologyOn) z-=dz(dt,z,theDisk.GetCos(),theDisk.GetDim());
    t+=dt;
    ++step;
  } // End of the loop over time steps.

  // All time steps completed! Why?
  int terminationCondition=0;
  if( dt <=TOL*dtfloor) {
    errormsg("Time step below floor");
    terminationCondition=1;
  }
  if( step >= stepmax) {
    std::cerr << "Maximum number of steps exceeded" << std::endl;
    terminationCondition=2;
  }
  if( t >= tmax ) {
    std::cout << "All outer orbits completed." << std::endl;
    terminationCondition=3;
  }
  if( z<=0.) {
    std::cout << "Reached redshift zero" << std::endl;
    terminationCondition=4;
  }
  if( dt >=fCondition*TOL) {
    std::cout << "---------------"<<std::endl;
    std::cout << "---------------"<<std::endl;
    std::cout << "Converged!"<<std::endl;
    std::cout << "---------------"<<std::endl;
    std::cout << "---------------"<<std::endl;
    terminationCondition=5;
    theDisk.store(ini);
  }

  if(writeOut) {
    std::ofstream convfile((filename+"_convergence2.dat").c_str());
    convfile << DtoB4.ratio() <<" "<<DtoB25.ratio() 
             <<" "<< ctrMdot4.ratio()  <<" " <<ctrMdot25.ratio() // 4
             <<" " <<outMdot4.ratio()  <<" "<<outMdot25.ratio()
             <<" "<<ctrsig4.ratio()    <<" "<<outsig4.ratio()    // 8
             <<" "<<ctrsig25.ratio()   <<" "<<outsig25.ratio()
             <<" "<<ctrSF4.ratio()     <<" "<<ctrSF25.ratio()
             <<" "<<outSF4.ratio()     <<" "<<outSF25.ratio()<<std::endl; //14
    convfile.close();

    theDisk.WriteOutStarsFile(filename+"_act",theDisk.active(),NActive,step);
    theDisk.WriteOutStarsFile(filename,theDisk.passive(),NPassive,step);
    theDisk.WriteOutStepFile(filename,accr,t,z,dt,step,tauvec);
  }

 
  // de-allocate memory
  for(unsigned int k=1; k<=2; ++k) {
    delete tauvec[k];
  }
  delete[] tauvec;


  // Tell the caller why the run halted.
  return terminationCondition;
}
