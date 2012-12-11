#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "Simulation.h"
#include "DiskContents.h"
#include "AccretionHistory.h"
#include "ArgumentSetter.h"
#include "Cosmology.h"
#include "Dimensions.h"
#include "DiskUtils.h"
#include "FixedMesh.h"
#include "Debug.h"


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
        errmgLine.push_back("Need at least one argument (filename)-- All possible arguments are ");
        errmgLine.push_back(" \nrunName, \nnx, \neta, \nepsff, \nTmig (local orbital times), \n");
        errmgLine.push_back("analyticQ (1 or 0),\ncosmologyOn (1 or 0), \nxmin, \nNActive, \n");
        errmgLine.push_back("NPassive, \nvphiR (km/s), \nradius (kpc), \ngasTemp (K), \nQlim, \n");
        errmgLine.push_back("fg0, \ntempRatio (sig_*/sig_g), \nzstart, \ntmax, \nstepmax, \n");
        errmgLine.push_back("TOL (t_orb),\nMassLoadingFactor, \nBulgeRadius (kpc), \n");
        errmgLine.push_back("stDiskScale (kpc, or -1 for powerlaw),\nwhichAccretionHistory,\nalphaMRI");
        errmgLine.push_back(", \nthick,\nmigratePassive,\nQinit,\nkappaMetals,\nMh0,\nminSigSt,\nnchanges");
        errmgLine.push_back("CAUTION: this message may be out of date, so check Main.cpp");
        std::string msg="";
        for(unsigned int i=0; i!=errmgLine.size();++i) {
            msg+=errmgLine[i];
        }
        //    errormsg(msg);
        std::cerr << msg << std::endl;
        return 1;
    }

    std::string filename(argv[1]);

    errormsg::errorFile.open((filename+"_stde.txt").c_str());

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
    const double vphiR =             as.Set(220,"Circular velocity (km/s)")*1.e5;
    const double radius=             as.Set(20.,"Outer Radius (kpc)")*cmperkpc;
    //  const double sigth=         sqrt(as.Set(7000,"Gas Temperature (K)")*kB/mH)/vphiR;
    const double Tgas =		   as.Set(7000,"Gas Temperature (K)");
    const double Qlim =              as.Set(2.5,"Limiting Q_*");
    const double fg0  =              as.Set(.5,"Initial gas fraction");
    const double tempRatio =         as.Set(1.,"Initial sigma_*/sigma_g");
    const double zstart =            as.Set(2.,"Initial redshift");
    const double tmax =              as.Set(1000,"Maximum Time (outer orbits)");
    const unsigned int stepmax=      as.Set(10000000,"Maximum Number of Steps");
    const double TOL =               as.Set(.0001,"TOL (outer orbits)");
    const double MassLoadingFactor=  as.Set(1,"Mass Loading Factor");
    const double BulgeRadius      =  as.Set(0,"Velocity Curve Turnover Radius (kpc)");
    const double innerPowerLaw    =  as.Set(.5,"Index of the inner power law part of the rot curve");
    const double softening        =  as.Set(2.0,"Softening of transition from flat to inner powerlaw rot curve");
    const double stScaleLength    =  as.Set(-1,"Initial Stellar Disk Scale Length (kpc)");
    const int whichAccretionHistory= as.Set(0,"Which Accretion History- 0-Bouche, 1-High, 2-Low");
    const double alphaMRI         =  as.Set(0,"alpha viscosity for the MRI");
    const double thick            =  as.Set(1.5,"Thickness correction to Q");
    const bool migratePassive=      (as.Set(1,"Migrate Passive population")==1);
    const double Qinit =             as.Set(2.0,"The fixed Q");
    const double kappaMetals =       as.Set(.001,"Kappa Metals");
    const double Mh0 =  		   as.Set(1.0e12,"Halo Mass");


    // Scale the things which scale with halo mass.
    //  const double MassLoadingFactor = MassLoadingFactorAtMh12 * pow((Mh0/1.0e12) , -1./3.);
    //  const double vphiR = vphiRatMh12 * pow(Mh0/1.0e12,  1./3.);
    const double sigth = sqrt(Tgas *kB/mH)/vphiR;

    const double minSigSt =          as.Set(1.0,"Minimum stellar velocity dispersion (km/s)")*1.e5/vphiR; 
    const double NChanges =          as.Set(6,"The number of times the lognormal accretion history should draw a new value.");
    const unsigned int Experimental= as.Set(0,"Debug parameter");
    const double accScaleLength    = as.Set(2.0,"Accretion ScaleLength (kpc)");

    const double zquench =           as.Set(-1.0,"Redshift at which accretion shuts off.");
    const double zrelax =            as.Set(zstart+1.0,"Redshift at which to start the relaxation of the disk");
    const double zetaREC =           as.Set(1.0,"Metal loading factor");
    const double RfREC =             as.Set(0.46,"Remnant fraction for inst. rec. approx.");
    const double deltaOmega =        as.Set(0.1,"Delta omega to generate Neistein10");
    const unsigned int Noutputs =    as.Set(200,"Number of outputs, excluding an initial and final output.");
    const double accNormalization=   as.Set(0.30959,"Normalization of accretion efficiency");
    const double accAlpha_z =        as.Set(0.38,"Power to which to raise (1+z) in accretion efficiency");
    const double accAlpha_Mh=        as.Set(-0.25,"Power to which to raise (M_h/10^12 Msun) in acc eff");
    const double accCeiling =        as.Set(1.0,"Maximum accretion efficiency.");
    const double fscatter =          as.Set(1.0,"Artificial factor by which to multiply scatter in ND10");
    const double invMassRatio =      as.Set(.3,"Maximum change in halo mass across 1 step in domega");


    // Make an object to deal with basic cosmological quantities.9
    // Omega_Lambda = .734, H0 = 2.29e-18 s^-1
    double OmegaLambda = 0.734;
    double OmegaMatter = 1.0-OmegaLambda;
    double H0 = 2.29e-18;
    Cosmology cos(OmegaMatter, OmegaLambda, H0 ,zrelax);

    Debug dbg(Experimental);

    // Make an object to deal with the accretion history
    AccretionHistory accr(Mh0,dbg);
    accr.SetEfficiencyParams(accNormalization, accAlpha_z, accAlpha_Mh, accCeiling);
    double mdot0;

    testAccretionHistory();

    int attempts=0;

    // Based on the user's choice of accretion history, generate the appropriate
    // mapping between redshift and accretion rate.
    if(whichAccretionHistory==0)
        mdot0 = accr.GenerateBoucheEtAl2009(zrelax,cos,filename+"_Bouche09.dat",true,true,zquench) * MSol/speryear;
    else if(whichAccretionHistory==2)
        mdot0 = accr.GenerateConstantAccretionHistory(2.34607,zstart,cos,filename+"_ConstAccHistory.dat",true) * MSol/speryear;
    else if(whichAccretionHistory==1)
        mdot0 = accr.GenerateConstantAccretionHistory(12.3368,zstart,cos,filename+"_ConstAccHistory2.dat",true)*MSol/speryear;
    else if(whichAccretionHistory<0 ) {
        mdot0=MSol/speryear * accr.GenerateLogNormal(zstart, zrelax, cos, 7.0*pow(Mh0*1.0e-12,1.1)*pow(1.+zstart,2.2),
                fscatter, NChanges, true, zquench, Mh0,-whichAccretionHistory,"_LogNormal.dat");
        //      mdot0 = accr.GenerateOscillatingAccretionHistory(10.0,-whichAccretionHistory,0.0,zstart,false,cos,filename+"_OscAccHistory.dat",true)*MSol/speryear;
    }
    else {
        mdot0 = accr.GenerateNeistein08(zstart,cos,filename+"_Neistein08_"+str(whichAccretionHistory)+".dat",
                true,whichAccretionHistory,invMassRatio,zquench,&attempts,
                deltaOmega, zrelax,fscatter)*MSol/speryear;
    }
    // Note that the following line does nothing but put a line in the comment file to
    // record MdotExt0 for this run.
    as.Set(mdot0/MSol*speryear,"Initial Accretion (MSol/yr)");
    as.Set(attempts,"Attempts to generate Neistein08");

    // Done reading in arguments. Write out a comment file containing all of the arguments.
    as.~ArgumentSetter();
    // If we're recording the convergence of the initial conditions, copy the comment file we just wrote out.
    if(dbg.opt(2) && dbg.opt(5)) {
        std::ifstream  src((filename+"_comment.txt").c_str());
        std::ofstream  dst((filename+"_icgen_comment.txt").c_str());

        dst << src.rdbuf();
    }



    // Set the dimensional quantities. 
    Dimensions dim(radius,vphiR,mdot0);
    FixedMesh mesh(innerPowerLaw,BulgeRadius/dim.d(1.0),softening,xmin,minSigSt,nx);
    double MhZs = accr.MhOfZ(zrelax)*Mh0; // this is in units of solar masses

    // don't relax the disk!
    if(!dbg.opt(5)) {
        DiskContents disk(tauHeat, eta, sigth, epsff, Qlim,
                TOL,analyticQ,MassLoadingFactor,cos,dim,mesh,dbg,
                thick,migratePassive,Qinit,kappaMetals,NActive,NPassive,
                minSigSt,accScaleLength/(radius/cmperkpc),RfREC,zetaREC);
        // double sig0 = 8.0/220.0; 
        double sig0 = sigth;
        double fcool = 1.0 ;//* sqrt(MhZs / 1.0e12);
        disk.Initialize(.1*Z_Sol,fcool,fg0,sig0,tempRatio,Mh0,MhZs,stScaleLength);

        Simulation sim(tmax,stepmax,cosmologyOn,nx,TOL,
                zstart,NActive,NPassive,alphaMRI,
                sigth,disk,accr,dbg,dim);
        int result = sim.runToConvergence(1.0e10, true, filename,zrelax,Noutputs);

    }

    // This section hasn't been used in a while and only really makes sense under certain conditions:
    if(dbg.opt(5)) {
        //// Evolve a disk where the stars do not do anything and Mdot_ext=Mdot_ext,0.
        DiskContents diskIC(1.0e30,eta,sigth,0.0,Qlim, // need Qlim to successfully set initial statevars
                TOL,analyticQ,MassLoadingFactor,cos,dim,mesh,dbg,
                thick, false,Qinit,kappaMetals,NActive,NPassive,minSigSt,accScaleLength/(radius/cmperkpc),
                RfREC,zetaREC);
        if(stScaleLength<0.0)  diskIC.Initialize(tempRatio,fg0);
        else diskIC.Initialize(0.1*Z_Sol, .6, fg0, tempRatio*50.0/220.0, Mh0, MhZs, stScaleLength);




        Simulation simIC(300.0,1000000,
                false, nx,TOL,
                zstart,NActive,NPassive,
                alphaMRI,sigth,
                diskIC,accr,dbg,dim);
        int result = simIC.runToConvergence(1, dbg.opt(2), filename+"_icgen",zstart,200); // set false-> true to debug initial condition generator
        if(result!=5) // The simulation converges when the time step reaches 1*TOL.
            errormsg("Initial Condition generator failed to converge, code "+str(result));


        // Now evolve a disk where the stars evolve as they should using the previous simulation's end condition
        DiskContents disk(tauHeat,eta,sigth,epsff,Qlim,
                TOL,analyticQ,MassLoadingFactor,cos,dim,mesh,dbg,
                thick,migratePassive,Qinit,kappaMetals,NActive,NPassive,
                minSigSt,accScaleLength/(radius/cmperkpc),
                RfREC,zetaREC);
        disk.Initialize(simIC.GetInitializer(), stScaleLength < 0.0); // if we're using an exponential disk, don't mess with the initial conditions of the stellar disk when enforcing Q=Q_f, i.e. do not keep a fixed phi0.
        Simulation sim(tmax,stepmax,
                cosmologyOn,nx,TOL,
                zstart,NActive,NPassive,
                alphaMRI,sigth,
                disk,accr,dbg,dim);
        result = sim.runToConvergence(1.0e10, true, filename,zstart,200);
    }

}
