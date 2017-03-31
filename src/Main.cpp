#include <stdio.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include "Simulation.h"
#include "DiskContents.h"
#include "AccretionHistory.h"
#include "AccretionProfile.h"
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
    std::cout.precision(15);  
    if(argc<2) {
        std::vector<std::string> errmgLine(0);
        errmgLine.push_back("Need at least one argument (filename)-- All possible arguments are ");
        errmgLine.push_back(" \nrunName, \nnx, \neta, \nepsff, \nTmig (local orbital times), \n");
        errmgLine.push_back("analyticQ (1 or 0),\ncosmologyOn (1 or 0), \nxmin, \nNActive, \n");
        errmgLine.push_back("NPassive, \nvphiR (km/s), \nradius (kpc), \ngasTemp (K), \nQlim, \n");
        errmgLine.push_back("fg0, \ntempRatio (sig_*/sig_g), \nzstart, \ntmax, \nstepmax, \n");
        errmgLine.push_back("TOL (t_orb),\nMassLoadingFactor, \nBulgeRadius (kpc), \n");
        errmgLine.push_back("stDiskScale ( -1 for powerlaw),\nwhichAccretionHistory,\nalphaMRI");
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
    ArgumentSetter as(argc,argv,filename+"_comment.txt");
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
    const double Tgas =		   as.Set(7000,"Gas Temperature (K)");
    const double sigth = sqrt(Tgas *kB/(mH))/vphiR;
    const double Qlim =              as.Set(2.5,"Limiting Q_*");
    const double fg0  =              as.Set(.5,"Initial gas fraction");
    const double tempRatio =         as.Set(1.,"Initial sigma_*/sigma_g");
    const double zstart =            as.Set(2.,"Initial redshift");
    const double tmax =              as.Set(1000,"Maximum Time (outer orbits)");
    const unsigned int stepmax=      as.Set(10000000,"Maximum Number of Steps");
    const double TOL =               as.Set(.0001,"TOL (outer orbits)");
    const double MassLoadingFactor=  as.Set(13.0,"Mass Loading Factor prefactor");
    const double MassLoadingColScaling =as.Set(-1.15,"Mass Loading Factor scaling with column density (normalized at 1 Msun/pc^2)");
    const double BulgeRadius      =  as.Set(0,"Velocity Curve Turnover Radius (kpc)");
    const double beta0            =  as.Set(.5,"Index of the inner power law part of the rot curve");
    const double nRotCurve        =  as.Set(2.0,"Sharpness of transition from flat to inner powerlaw rot curve");
    const double stScaleLength    =  as.Set(-1,"Initial Stellar Disk Scale Length (kpc)");
    const int whichAccretionHistory= as.Set(0,"Which Accretion History- 0-Bouche, 1-High, 2-Low");
    const double alphaMRI         =  as.Set(0,"alpha viscosity for the MRI");
    const double thick            =  as.Set(1.5,"Thickness correction to Q");
    const bool migratePassive=      (as.Set(1,"Migrate Passive population")==1);
    const double Qinit =             as.Set(2.0,"The fixed Q");
    const double kappaMetals =       as.Set(.001,"Kappa Metals");
    const double Mh0 =  		   as.Set(1.0e12,"Halo Mass");
    const double minSigSt =          as.Set(1.0,"Minimum stellar velocity dispersion (km/s)")*1.e5/vphiR; 
    const double NChanges =          as.Set(6,"The number of times the lognormal accretion history should draw a new value.");
    const unsigned int Experimental= as.Set(0,"Debug parameter");
    const double accScaleLength    = as.Set(.05,"Accretion ScaleLength (fraction of Vir radius)");

    const double zquench =           as.Set(-1.0,"Redshift at which accretion shuts off.");
    const double zrelax =            as.Set(zstart+1.0,"Redshift at which to start the relaxation of the disk");
    const double xiREC = as.Set(0.0,"Stellar ejecta isolation- 0 for perfect mixing, 1 for pure ejecta");
    const double RfREC =             as.Set(0.46,"Remnant fraction for inst. rec. approx.");
    const double deltaOmega =        as.Set(0.1,"Delta omega to generate Neistein10");
    const unsigned int Noutputs =    as.Set(200,"Number of outputs, excluding an initial and final output.");
    const double accNormalization=   as.Set(0.30959,"Normalization of accretion efficiency");
    const double accAlpha_z =        as.Set(0.38,"Power to which to raise (1+z) in accretion efficiency");
    const double accAlpha_Mh=        as.Set(-0.25,"Power to which to raise (M_h/10^12 Msun) in acc eff");
    const double accCeiling =        as.Set(1.0,"Maximum accretion efficiency.");
    const double fscatter =          as.Set(1.0,"Artificial factor by which to multiply scatter in NMD10");
    const double invMassRatio =      as.Set(.3,"Maximum change in halo mass across 1 step in domega");
    const double fcool =             as.Set(1.0,"Fraction of f_b M_h(zstart) that has cooled into a disk");
    const double whichAccretionProfile=as.Set(0,"Accretion profile- 0,1,2, respectively exponential, flat, or gaussian");
    const double alphaAccProf =      as.Set(0.0,"Powerlaw scaling of accretion length with (Mh/Mh0)");
    const double width =             as.Set(0.1,"The width of a gaussian accretion profile in units of the accretionScaleLength (only valid for whichAccretionProfile==2)");
    const double fH2Min =            as.Set(.03,"Minimum fH2");
    const double tDepH2SC =          as.Set(2.0,"Depletion time (Gyr)");
    const double ZIGM =              as.Set(.002,"Z of IGM in absolute units");
    const double yREC =              as.Set(.054,"yield - mass of metals produced per gas mass locked in stars");
    const double concentrationRandomFactor= as.Set(0.0, "Constant multiplicative offset from Mh-c relation (dex)");
    const double MassLoadingFgScaling=as.Set(0.16, "Scaling of the mass loading factor with the gas fraction");
    const double ksuppress =         as.Set(10.0, "The characteristic mode of a DFT to begin suppressing power in computation of vPhiDisk");
    const double kpower =            as.Set(2.0,"The power to which to raise the argument of the exponential in suppressing powe abve ksuppess");
    const double MQuench =           as.Set(1.0e12, "Halo mass at which quenching occurs");
    const double epsquench =         as.Set(0.1, "Factor by which accretion efficiency is depressed above MQuench");
    const double muQuench =          as.Set(0.0, "mass loading factor in central kpc above Mquench (overrides ordinarily-calculated mu if larger)");
    const double stScaleReduction=   as.Set(2.0, "Factor by which to reduct initial stellar scale length relative to initial exponential scale length of the accretion");
    const double gaScaleReduction=   as.Set(1.4, "Factor by which to reduct initial gas scale length relative to initial exponential scale length of the accretion");
    const double ZMix =              as.Set(0.5, "Fraction of outflowing metals mixed back in to inflowing gas");


    as.WriteOut();

    // Make an object to deal with basic cosmological quantities.9
    // Omega_Lambda = .734, H0 = 2.29e-18 s^-1
    double OmegaMatter = .258; // WMAP5
    double OmegaLambda = 1.0 - OmegaMatter;
    double H0 = 2.29e-18;
    double sigma8 = 0.796;

    // Planck!
    sigma8 = 0.826;
    OmegaLambda = 0.692;
    OmegaMatter = 1.0-OmegaLambda;
    H0 = 67.80*3.2407793e-20;


    Cosmology cos(OmegaMatter, OmegaLambda, H0 , sigma8, zrelax, nx, concentrationRandomFactor);

    Debug dbg(Experimental);

    // Make an object to deal with the accretion history
    AccretionHistory accr(Mh0,dbg,MQuench, epsquench);
    accr.SetEfficiencyParams(accNormalization, accAlpha_z, accAlpha_Mh, accCeiling);
    double mdot0;

    testAccretionHistory();

    int attempts=0;

    // Based on the user's choice of accretion history, generate the appropriate
    // mapping between redshift and accretion rate.
    if(whichAccretionHistory==0) {
        if(!dbg.opt(15))
            mdot0 = accr.GenerateBoucheEtAl2009(zrelax,cos,filename+"_Bouche09.dat",true,true,zquench) * MSol/speryear;
        else
            mdot0 = accr.GenerateAverageNMD10(
                zstart,cos,filename+"_AvgNeistein08.dat",
                true,whichAccretionHistory,invMassRatio,zquench,10000,
                deltaOmega, zrelax,fscatter)*MSol/speryear;
    }
    else if(whichAccretionHistory==2)
        mdot0 = accr.GenerateConstantAccretionHistory(2.34607,zstart,cos,filename+"_ConstAccHistory.dat",true) * MSol/speryear;
    else if(whichAccretionHistory==1)
        mdot0 = accr.GenerateConstantAccretionHistory(12.3368,zstart,cos,filename+"_ConstAccHistory2.dat",true)*MSol/speryear;
    else if(whichAccretionHistory<0 ) {
        bool constInTime = dbg.opt(14);
        bool readIn = dbg.opt(2); // Read in an external file, e.g. generated from Bolshoi outputs
        mdot0=MSol/speryear * accr.GenerateLogNormal(zstart, zrelax, cos, 
              fscatter, NChanges, true, zquench, Mh0,-whichAccretionHistory,filename+"_LogNormal.dat",constInTime,readIn,filename+"_inputRandomFactors.txt");
        //      mdot0 = accr.GenerateOscillatingAccretionHistory(10.0,-whichAccretionHistory,0.0,zstart,false,cos,filename+"_OscAccHistory.dat",true)*MSol/speryear;
    }
    else {
        mdot0 = accr.GenerateNeistein08(zstart,cos,filename+"_Neistein08_"+str(whichAccretionHistory)+".dat",
                true,whichAccretionHistory,invMassRatio,zquench,&attempts,
                deltaOmega, zrelax,fscatter)*MSol/speryear;
    }

    if(mdot0<=0)
        errormsg("Negative mdot0! The model does not know what you mean by negative accretion.");

    ArgumentSetter as2(0,argv,filename+"_aux.txt");
    // Note that the following line does nothing but put a line in the comment file to
    // record MdotExt0 for this run.
    as2.Set(mdot0/MSol*speryear,"Initial Accretion (MSol/yr)");
    as2.Set(attempts,"Attempts to generate Neistein08");
    // This is where we'll store a record of all the possiblities of dbg.opt
    as2.Set(dbg.opt(0), "Use constant mass loading factor instead of Lagos13"); // recommended
    as2.Set(dbg.opt(1), "Dont increase accr rate to account for matter accreting outside domain");
    as2.Set(dbg.opt(2), "Read in an external file for accretion rates");
    as2.Set(dbg.opt(3), "Neistein & Dekel (2008) instead of Neistein+ (2010)"); // check this
    as2.Set(dbg.opt(4), "dQ/dt ~ exp(Delta Q) - 1 instead of 0");
    as2.Set(dbg.opt(5), "Accreting metallicity changes along with ZDisk");
    as2.Set(dbg.opt(6), "Stellar recycling");
    as2.Set(dbg.opt(7), "Attenuate high-k modes when computing vPhiDisk");
    as2.Set(dbg.opt(8), "Like Bouche, but use prescription from Dekel13 WMAP5 cosmology" );
    as2.Set(dbg.opt(9), "not minmod in ddx(sig)");
    as2.Set(dbg.opt(10),"Override initialization params to enforce agreement with observations"); 
    as2.Set(dbg.opt(11), "Take non-Euler timesteps"); 
    as2.Set(dbg.opt(12), "Artificially set GI torque=0 everywhere");
    as2.Set(dbg.opt(13), "Print out rotation curve info at n==nx");
    as2.Set(dbg.opt(14), "For lognormal acc history, bursts uniform in time (otherwise uniform in z)");
    as2.Set(dbg.opt(15), "If whichAccretionHistory is 0, use the average from NMD instead of Bouche");
    as2.Set(dbg.opt(16), "Newly formed stars have full gas velocity dispersion instead of turbulent component only");
    as2.Set(dbg.opt(17), "upstream");
    as2.Set(dbg.opt(18), "overshoot");
    as2.Set(dbg.opt(19), "No longer used");


    // Done reading in arguments. Write out a comment file containing all of the arguments.
    as2.WriteOut();



    // Set the dimensional quantities. 
    Dimensions dim(radius,vphiR,mdot0);
    FixedMesh mesh(beta0,BulgeRadius/dim.d(1.0),nRotCurve,xmin,minSigSt,nx);
    mesh.storeSummand();
    double MhZs = accr.MhOfZ(zrelax)*Mh0; // this is in units of solar masses
    cos.UpdateProfile(MhZs, zrelax, mesh.x(), dim.Radius);
    double r200 = cos.r200(MhZs,zrelax); // this will be in cm

    AccretionProfile accProf(mesh, whichAccretionProfile, alphaAccProf, dbg, accScaleLength ,width);

    double ZIGMO = 0.0057/0.02 * ZIGM;
    double ZIGMFe = 0.0013/0.02 * ZIGM;
    DiskContents disk(tauHeat, eta, sigth, epsff, Qlim,
            TOL,analyticQ,MassLoadingFactor,MassLoadingColScaling,MassLoadingFgScaling,
            cos,dim,mesh,dbg,
            thick,migratePassive,Qinit,kappaMetals,NActive,NPassive,
          minSigSt,RfREC,xiREC,fH2Min,tDepH2SC,ZIGMO, ZIGMFe,yREC, ksuppress, kpower, MQuench, muQuench,
          ZMix);
    // double sig0 = 8.0/220.0; 
    double sig0 = sigth;
    double stScaleLengthA = accScaleLength*r200/cmperkpc; // accScaleLength * pow(MhZs/Mh0,alphaAccProf);
    disk.Initialize(fcool,fg0,sig0,tempRatio,Mh0,MhZs,stScaleLengthA,zrelax, stScaleReduction, gaScaleReduction);

    Simulation sim(tmax,stepmax,cosmologyOn,nx,TOL,
            zstart,NActive,NPassive,alphaMRI,
            sigth,disk,accr,dbg,dim,accProf);
    int result = sim.runToConvergence(1.0e10, true, filename,zrelax,Noutputs);


    return 0;
}
