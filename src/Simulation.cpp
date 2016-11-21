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
#include "AccretionProfile.h"


Simulation::Simulation(const double tm, const long int sm,
        const bool co,   const unsigned int nnx,
        const double tl, const double zs, 
        const unsigned int na, const unsigned int np,
        const double amri,const double sth, 
        DiskContents& td,
        AccretionHistory& ah, Debug& db, Dimensions& dm,
        AccretionProfile& ap):
    tmax(tm), stepmax(sm),
    cosmologyOn(co),nx(nnx),
    TOL(tl), zstart(zs),
    NActive(na),NPassive(np),
    alphaMRI(amri),sigth(sth),
    theDisk(td),
    accr(ah),
    accProf(ap),
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
        const std::string filename,
        const double zrelax,
        const unsigned int Noutputs)
{

    if(zrelax < zstart)
        errormsg("zrelax < zstart. The code makes the assumption that this is not true.");

    std::vector<double> UU(nx+1,0.),LL(nx+1,0.),DD(nx+1,0.),FF(nx+1,0.), 
        UUst(nx+1,0.),LLst(nx+1,0.),DDst(nx+1,0.),FFst(nx+1,0.),
        MdotiPlusHalf(nx+1,0.), MdotiPlusHalfStar(nx+1,0.),MdotiPlusHalfMRI(nx+1,0.);
    // Initialize a 2 by nx vector to store the torque and its first derivative.
    double **tauvec;
    double **tauvecStar;
    double **tauvecMRI;
    tauvec = new double *[3];
    tauvecStar = new double *[3];
    tauvecMRI = new double *[3];
    for(unsigned int k=1; k<=2; ++k) {
        tauvec[k] = new double[nx+2];
        tauvecStar[k] = new double[nx+2];
        tauvecMRI[k] = new double[nx+2];
    }

    // This is the equilibrium solution to the torque eq from Krumholz & Burkert 2010
    for(unsigned int n=0; n<=nx+1; ++n) {
        tauvec[1][n]=0.0;
        tauvec[2][n]=0.0;
        tauvecStar[1][n] = 0.0;
        tauvecStar[2][n] = 0.0;
        tauvecMRI[1][n] = 0.0;
        tauvecMRI[2][n] = 0.0;
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

    const double dtfloor=1.e-25; 
    // Initialize some variables which will change as the simulation progresses:
    // time (outer orbits), redshift, and timestep
    double t=0., z=zrelax, dt=dtfloor*10, dtPrev = 1.0;

    int whichVar,whichCell; // store which variable and which cell limits the time step.

    unsigned int step=0; // number of elapsed timesteps

    // number of times data has been written to output
    int writeIndex=0;

    double AccRate=1.0; // normalized accretion rate

    // Stop the simulation if the timestep falls below this 
    // value (in general this means that something has gone 
    // terribly wrong and the simulation has no chance to 
    // run to completion)

    // t0 is the number of seconds between zrelax and zstart. Recording this lets us
    // put out the same snapshots for the same zstart value, regardless of the zrelax
    // we choose.
    const double t0 = theDisk.GetCos().lbt(zrelax) - theDisk.GetCos().lbt(zstart); // seconds
    const double tzs = t0*dim.vphiR/(2.0*M_PI*dim.Radius);
    const double duration = theDisk.GetCos().lbt(zstart); //in seconds
    std::vector<std::vector<int> > flagList(0);
    for(unsigned int n=0; n<=nx; ++n ) {
        flagList.push_back(std::vector<int>(3));
    }

    bool readIn = true;
    std::ifstream inputRandomFactors;
    std::string fnIn = filename+"_inputRandomFactorsY.txt";
    std::vector<double> irfY(0);
    int rfcounter=0;
    if(readIn) inputRandomFactors.open(fnIn.c_str());
    bool readFlag = true;
    while (readFlag && readIn) {
        std::cout << "Attempting to read in random factor "<<rfcounter<<std::endl;
        ++rfcounter;
        std::string line;
        readFlag = getline(inputRandomFactors, line);
        if (!readFlag)
            std::cout << "Failed to read in random factor " << rfcounter<<std::endl;
        double y = atof(line.c_str());
        irfY.push_back(y);
        std::cout << "Successfully read in y = "<<y<<" at z="<<z<<std::endl;
    }
    if(readIn) inputRandomFactors.close();




    // Loop over time steps until...
    while(t<tmax &&           // the simulation has run for too many orbits
            step<=stepmax &&    // the simulation has run for too many steps
            z>=0.0 &&           // the simulation has reached redshift zero
            (dtPrev>TOL*dtfloor || dt>TOL*dtfloor) &&   // the simulation has stalled
            dt<fCondition*TOL)  // the simulation has converged
    {

        // Every 10th of an outer orbit, write out the data from the simulation
        // Note that we'd like to use last step's timestep here to determine whether
        // to write out a file, since we would like to avoid a situation where a sudden
        // decrease in the size of dt causes us to skip a checkpoint
        // bool timeOut = (((floor(25.0*(t-dt)) < floor(25.0*t)) || step<2) && writeOut);
        //    bool timeOut = true; // write out at every time step. Use only if 
        // things crash and burn very quickly
        double present = t * 2.0*M_PI*dim.Radius / dim.vphiR - t0;
        double previous = (t-dt) * 2.0*M_PI*dim.Radius/dim.vphiR - t0;

        double NOUT = ((double) Noutputs);
        bool timeOut = (((floor(NOUT*previous/duration) < floor(NOUT*present/duration) && present>0.0) || step == 1) && writeOut ) ;

        if(!cosmologyOn)
            timeOut= (((floor(25.0*(t-dt)) < floor(25.0*t)) || step<2) && writeOut);

        if(timeOut) {
            // Separate files for the active & passive stellar populations..
            theDisk.WriteOutStarsFile(filename+"_act",theDisk.active(),NActive,step);
            theDisk.WriteOutStarsFile(filename,theDisk.passive(),NPassive,step);
            theDisk.WriteOutStepFile(filename,accr,t-tzs,z,dt,step,tauvec,tauvecStar,tauvecMRI,
                    MdotiPlusHalf,MdotiPlusHalfMRI,accProf.GetProfile(),accProf.fInner());

            writeIndex++;

            std::cout << "Writing out file "<<writeIndex<<" at t = "<<present/(speryear*1.0e9)<<" Gyr, z= "<<z<<std::endl;
        }

        // This step is kind of expensive, so only do it as frequently as we need to.
        bool updateRotationCurve = (step<5 || (step%10)==0);
        if(updateRotationCurve) {
            theDisk.GetCos().UpdateProfile(accr.MhOfZ(z)*accr.GetMh0(), z, theDisk.GetMesh().x(), dim.Radius);
            theDisk.UpdateRotationCurve(accr.MhOfZ(z)*accr.GetMh0(), z, dt);
        }
        else 
            theDisk.ZeroDuDt();

        theDisk.ComputeColSFRapprox(accr.MhOfZ(z) * accr.GetMh0(),z);
        // theDisk.ComputeColSFR(accr.MhOfZ(z) * accr.GetMh0(),z);
        theDisk.ComputeMassLoadingFactor(accr.MhOfZ(z) * accr.GetMh0(), theDisk.active()[0]->spcol);

        if(dbg.opt(6)) {
            for(unsigned int i=0; i!=theDisk.passive().size(); ++i) {
                theDisk.passive()[i]->ComputeRecycling(theDisk,z);
            }
        }


        // Make sure that stars formed this time step have a stellar 
        // population to which they belong. If they don't, create
        // one and add it to the appropriate vector of stellar populations
        if(NActive>1 && floor(((double) NActive)*previous/duration) < floor(((double) NActive)*present/duration) && present>0.0)
            theDisk.AddNewStellarPop(z,dt,theDisk.active(),true);
        if(NPassive>1 && floor(((double) NPassive)*previous/duration) < floor(((double) NPassive)*present/duration) && present>0.0)
            theDisk.AddNewStellarPop(z,dt,theDisk.passive(),false);


        // Look up the non-dimensional accretion rate at this redshift. 
        if(cosmologyOn) AccRate = accr.AccOfZ(z);
        else AccRate = 1.;
        
        double r200 = theDisk.GetCos().r200(accr.MhOfZ(z)*accr.GetMh0(), z)/dim.Radius;

        int irfInd = 0.0;
        if (irfY.size() > 0) {
            irfInd = irfY.size()-1 - ((int) floor( z / (zrelax/(((double) irfY.size())-1.0)) ));
            r200 *= pow(10.0, irfY[irfInd] );
            if( step % 1000 == 0) {
                std::cout << "At z= "<<z<<", irfInd= "<<irfInd<<" and irf= "<< irfY[irfInd] << std::endl;
            }
        }



        accProf.UpdateProfile(r200);
        // accProf.UpdateProfile(accr.MhOfZ(z));

        // Given the user's freedom to specify the accretion history, make sure
        // the accretion rate remains non-negative.
        if(AccRate<0.0) errormsg("Negative accretion history!");

        theDisk.UpdateStTorqueCoeffs(UUst,DDst,LLst,FFst);
        theDisk.ComputeGItorque(tauvecStar,0.0,0.0,UUst,DDst,LLst,FFst,MdotiPlusHalfStar);


        // Compute dQ/dA and its error, where A is a stand-in for every state variable
        theDisk.ComputePartials();


       ComputeFlagList(MdotiPlusHalf,MdotiPlusHalfMRI,flagList); 

        // Update the coefficients of the torque equation
        theDisk.UpdateCoeffs(z,UU,DD,LL,FF,tauvecStar,MdotiPlusHalfStar,tauvecMRI,MdotiPlusHalfMRI, accProf.GetProfile(), AccRate, flagList);


        // Solve the torque equation. The commented versions represent various choices
        // for the boundary conditions
        //    disk.ComputeTorques(tauvec,-1.*AccRate*xmin/pow(xmin,1./(1.-nx)),-1.*AccRate);
        //    disk.ComputeTorques(tauvec,0.,-1.*AccRate);

        double IBC=0.0;
        double OBC=0.0;
        theDisk.ComputeGItorque(tauvec,IBC,OBC,UU,DD,LL,FF,MdotiPlusHalf);
        //    disk.ComputeTorques(tauvec,-1.*AccRate*xmin,-1.*AccRate);

        // Artificially set gravitoturbulent torques to zero everywhere.
        if(dbg.opt(12)) { 
            for(unsigned int n=0; n<=nx; ++n) {
                tauvec[1][n] = 0.0;
                tauvec[2][n] = 0.0;
                MdotiPlusHalf[n]=0.0;
            }
        }


        // In situations where the alpha viscosity produces larger torques 
        // than the GI viscosity, use the alpha viscosity instead:
	if(step>1) {
	  theDisk.ComputeMRItorque(tauvecMRI,alphaMRI);
	  ComputeFluxes(tauvecMRI,MdotiPlusHalfMRI,theDisk.GetMesh());
	}

        // Given the solution to the torque equation, compute time 
        // derivatives of the state variables
        theDisk.ComputeDerivs(tauvec,MdotiPlusHalf,
                tauvecMRI,MdotiPlusHalfMRI,
                accProf.GetProfile(),AccRate,
                flagList);

        // Given the derivatives, compute a time step over which none of the variables
        // change by too much. whichVar tells us which state variable is limiting the timestep
        // and whichCell tells us which cell is limiting the timestep. Both of these values
        // are printed every 5000 timesteps (see below).
        dtPrev = dt;
//        if(step>1)
            dt = theDisk.ComputeTimeStep(z,&whichVar,&whichCell,tauvecStar,MdotiPlusHalfStar); 

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


        theDisk.UpdateStateVars(dt,dtPrev,z,tauvec,AccRate,tauvecStar,
				MdotiPlusHalf,MdotiPlusHalfStar,MdotiPlusHalfMRI,accProf.fInner(), accr.AccStOfZ(z)); 


        // update the independent variables.
        if(cosmologyOn) z-=dz(dt,z,theDisk.GetCos(),theDisk.GetDim());
        t+=dt;
        ++step;
    } // End of the loop over time steps.

    // All time steps completed! Why?
    int terminationCondition=0;
    if( dt <=TOL*dtfloor) {
        errormsg("Time step below floor in cell "+str(whichCell)+" with variable "+str(whichVar)+" at time "+str(t)+" step "+str(step)+" and dt "+str(dt));
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
        theDisk.WriteOutStepFile(filename,accr,t-tzs,z,dt,step,tauvec,tauvecStar,tauvecMRI,
                MdotiPlusHalf,MdotiPlusHalfMRI,accProf.GetProfile(),accProf.fInner());

        writeIndex++;
        std::cout << "Writing out file "<<writeIndex<<" at z= "<<z<<std::endl;
    }

    if(writeIndex != Noutputs + 2)
        errormsg("Unexpected number of outputs! "+str(writeIndex)+" vs "+str(Noutputs+2));

    // de-allocate memory
//    for(unsigned int k=1; k<=2; ++k) {
//        delete tauvec[k];
//        delete tauvecStar[k];
//        delete tauvecMRI[k];
//    }
//    delete[] tauvec;
//    delete[] tauvecStar;
//    delete[] tauvecMRI;


    // Tell the caller why the run halted.
    return terminationCondition;
}
