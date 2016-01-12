#include "Errors.h"
#include "AccretionHistory.h"
#include "Cosmology.h"
#include "Dimensions.h"
#include "DiskUtils.h"
#include "Debug.h"


#include <fstream>
#include <math.h>
#include <algorithm>
#include <iostream>

AccretionHistory::~AccretionHistory()
{
    if(allocated) {
        gsl_spline_free(spline);
        gsl_interp_accel_free(accel);
        gsl_spline_free(splineMh);
        gsl_interp_accel_free(accelMh);
        //delete[] acc;
        //delete[] redshift;
        //delete[] hMass;
    }
}

double AccretionHistory::GenerateOscillatingAccretionHistory(double amp, double period,
        double phase,		     
        double zst, bool inRedshift, // as opposed to in Gyr.
        Cosmology& cosm, std::string fn, bool writeOut)
{
    zstart = zst;
    std::ofstream file;
    if(writeOut) file.open(fn.c_str());

    double z=zstart;
    double MdotExt0=0.0; 
    double mass = 1.0e12 - (amp/2.0) * cosm.Tsim(0.0) / speryear;
    unsigned int N=10000;

    std::vector<double> redshifts(0), tabulatedAcc(0), masses(0);

    double eps = 1.0e-8; // safety factor to make sure we never get Mdot<0

    for(unsigned int i=0; i<=N; ++i) {
        z=(((double) (N-i))/((double) N)) * (zstart - 0.0);
        double MdotExt=0;
        double t=cosm.Tsim(z);
        if(inRedshift) 
            MdotExt = (amp/(2.0-eps)) + 
                (amp/2.0)*cos(2.0*M_PI*(z-zstart)/period-phase);
        else           
            MdotExt = (amp/(2.0-eps)) + 
                (amp/2.0)*cos(2.0*M_PI*(t/speryear)*1.0e-6/period-phase);
        if(i==0) MdotExt0 = MdotExt;

        if(writeOut) file << z <<" "<<t<<" "<<MdotExt<<" "<<-1.0<<std::endl;
        redshifts.push_back(z); tabulatedAcc.push_back(MdotExt); masses.push_back(mass);
        mass += MdotExt * fabs(cosm.Tsim(redshifts[N-i+1]) - cosm.Tsim(z)) / speryear;
    }
    file.close();
    InitializeGSLObjs(redshifts,tabulatedAcc,masses);
    if (MdotExt0<=0)
        errormsg("Generating an oscillating accretion history has produced negative mdot.");

    return MdotExt0;
}


double AccretionHistory::GenerateConstantAccretionHistory(double rate, double zst,
        Cosmology& cos, std::string fn,bool writeOut)
{
    zstart = zst;
    std::ofstream file;
    if(writeOut) file.open(fn.c_str());
    unsigned int N=1000;
    double z=zstart;
    double MdotExt0=0.0;
    double mass=1.0e12 - rate * cos.Tsim(0.0) / speryear;
    std::vector<double> redshifts(0),tabulatedAcc(0),masses(0);

    for(unsigned int i=0; i<=N; ++i) {
        z=(((double) (N-i))/((double) N))*(zstart-0.0);
        double MdotExt = rate; //solar masses /year
        if(i==0) MdotExt0 = MdotExt;
        if(writeOut) file << z <<" "<<cos.Tsim(z)<<" "<<MdotExt<<" "<<-1.0<<std::endl;
        redshifts.push_back(z); tabulatedAcc.push_back(MdotExt); masses.push_back(mass);
        mass += rate* fabs(cos.Tsim(redshifts[N-i+1]) - cos.Tsim(z)) / speryear;
    }
    file.close();
    InitializeGSLObjs(redshifts,tabulatedAcc,masses);
    if (MdotExt0<=0)
        errormsg("Generating a constant accretion history has produced negative mdot.");
    return MdotExt0;
}

double omega(double z,void * rhs)
{
    double * theRHS = (double *) rhs;
    return 1.260*(1.0+z+0.09/(1.0+z) + 0.24*exp(-1.16*z))-(*theRHS);
}

double zOfOmega(double om)
{
    gsl_function F;
    F.function = &omega;
    double tempOm1=om;
    double z=om;
    F.params=&tempOm1;
    findRoot(F,&z);
    return z;

}
double u(double x)
{
    return 64.087*pow((1.0+1.074*pow(x,0.3) - 1.581*pow(x,0.4)+
                0.954*pow(x,0.5) - 0.185*pow(x,0.6)),-10.0);
}

struct SParams {
    double rhs,OmegaM,sigma8;
};

double S(double M, void * p)
{
    SParams * sp = (SParams *) p;
    double c0=3.804e-4;
    double Gamma=0.169; // power spectrum shape parameter
    double sigma8= (*sp).sigma8;
    double x=c0*Gamma*pow(M,1./3.)/pow((*sp).OmegaM,1./3.);
    return u(x)*u(x)*sigma8*sigma8/(u(32.0*Gamma)*u(32.0*Gamma)) - (*sp).rhs;
}

// Return M(S) in units of solar masses
double MofS(double theS,double s8,double OmegaM)
{
    gsl_function F;
    F.function = &S;
    SParams sp;
    sp.rhs=theS;
    sp.sigma8=s8;
    sp.OmegaM=OmegaM;
    F.params=&sp;
    double M=1.0e12;
    findRoot(F,&M);
    return M;
}

double npow(double x,double ex) 
{
    if(x<0.0) {
        return -pow(-x,ex);
    }
    else
        return pow(x,ex);
}

double AccretionHistory::GenerateLogNormal(double zst,double zrelax, Cosmology& cos,
                    double scatter, double Nchanges,
                    bool writeOut, double zquench,
                    double Mh0, unsigned int seed, std::string fn, bool constInTime,
                    bool readIn, std::string fnIn)
{
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r,seed);

    zstart = zst;
    double Mh012 = Mh0*1.0e-12;
    linear = true;
    std::ofstream file;
    std::ifstream inputRandomFactors;
    if(readIn) inputRandomFactors.open(fnIn.c_str());
    if(writeOut) file.open(fn.c_str());
    unsigned int N=1000; 
    if(Nchanges > N/100) N=100*Nchanges;
    double Mh=Mh012; // units of 10^12 solar masses
    // double fbp18 = 17.0/18.0; // baryon fraction / 0.18
    double MdotExt0, MdotExt;
    std::vector<double> redshifts(0),tabulatedAcc(0),haloMass(0);
    double duration = cos.Tsim(0.0);
    // double DeltaT = duration/Nchanges; // seconds
    double currentAccretionRate;
    double z,dz;
    double x=0.0;
    int rfcounter=0; // debugging 

    // Loop over redshift from z=0 to z=zstart in N increments if Mh is specified at z=0
    // Otherwise, loop over redshift from zstart to z=0.
    for(unsigned int i=0; i<=N; ++i) {
        z=((double) i)/((double) N) * (zrelax*1.01 - 0.0) -.0005*zrelax;
        dz = 1.0/((double) N) * (zrelax*1.01 - 0.0);

        // Use the Bouche formula to tell us the dark matter accretion rate
        double dMh = 34.0 * pow(Mh,1.14)*pow(1.0+z,2.4) * 1.0e-12; // in 10^12 Msol/yr

        // Use the analogous formula to tell us the baryonic accretion rate.
//        double MdotExt = 7.0 * epsin(z,Mh,cos,zquench) * fbp18 * pow(Mh,1.1)*pow(1+z,2.2); // in solar masses /year
        double present = cos.lbt(z); // lookback time (in seconds) of the current redshift.
        double next = cos.lbt(z+dz); // this will be a larger number, i.e. a number further back in time.
        bool drawNewNumber = (floor((Nchanges-1) * present/duration) < floor((Nchanges-1) * next/duration) && z<zstart && constInTime) || (floor((Nchanges-1) * z/zstart) < floor((Nchanges-1) * (z+dz)/zstart) && z<zstart && !constInTime);// || (z-zstart)*(z+dz-zstart) <= 0.0;
        if(drawNewNumber && !readIn) {
            x = gsl_ran_gaussian(r,1.0);
        }
        else if(drawNewNumber && readIn) {
            std::cout << "Attempting to read in random factor "<<rfcounter<<std::endl;
            ++rfcounter;
            std::string line;
            bool readFlag = getline(inputRandomFactors, line);
            if (!readFlag) {
                // errormsg("Failed to read in random factor!");
                std::cout<< "WARNING: failed to read in random factor. Using newly-generated random number" <<std::endl;
                x  = gsl_ran_gaussian(r,1.0);
            }
            else {
                x = atof(line.c_str());
            }
            std::cout << "Successfully read in x = "<<x<<" at z="<<z<<std::endl;
        }
        // Set dMh such that the average accretion rate == the value given above by dMh under the assumption that x is normally distributed
        // currentAccretionRate = dMh*1.0e12 * exp(-scatter*scatter*log(10.)*log(10.)/2.) * pow(10.0, x*scatter);//*epsin(z,Mh,cos,zquench);
        // In this version make no assumption about the distribution of x (because it's being generated by a simulation). x is now simply the log10 of the ratio of the simulated accretion rate to the "average" one given above.
        currentAccretionRate = dMh*1.0e12 * pow(10.0, x*scatter);//*epsin(z,Mh,cos,zquench);
        
        if(z>zstart) currentAccretionRate = dMh*1.0e12; // above zstart don't include any stochasticity
        MdotExt = currentAccretionRate * epsin(z,Mh,cos,zquench)*.17;

        if((z-zstart)*(z+dz-zstart)<=0.0)// always set MdotExt0 to be MdotExt at z=2
            MdotExt0= MdotExt;

        haloMass.push_back(Mh); // 10^12 solar masses

        if(true) { //z>zquench) {
            // Basically compute Mh(z) by taking an Euler step, since from the above we know dMh (which is actually dMh/dt)
            Mh+= currentAccretionRate*1.0e-12 * -1.0*( cos.Tsim(z) - cos.Tsim(z+dz)) / speryear;
        }

        // these small adjustments to z avoid interpolation errors without really affecting anything.
        redshifts.push_back(z); 

        if(MdotExt < 1.0e-10)
            errormsg("Very low Mdot. If you're sure it's fine edit line 233ish of AccretionHistory.cpp");
        if(MdotExt > 1.0e6)
            errormsg("Very large Mdot. If you're sure it's fine edit line 235ish of AccretionHistory.cpp");

        tabulatedAcc.push_back(MdotExt);
        if(writeOut) file << z << " "<< cos.Tsim(z) <<" "<<MdotExt<<" "<<Mh<<" "<<x<<" "<<scatter<<" "<<epsin(z,Mh,cos,zquench)<<" "<<currentAccretionRate<<std::endl;
    }

        // 
    reverse(redshifts.begin(),redshifts.end());
    reverse(tabulatedAcc.begin(), tabulatedAcc.end());
    reverse(haloMass.begin(),haloMass.end());

    if(writeOut) file.close();
    if(readIn) inputRandomFactors.close();
    InitializeGSLObjs(redshifts,tabulatedAcc,haloMass);
    if( MdotExt0<=0)
        errormsg("Generating a lognormal accretion history has produced negative mdot.");
    return MdotExt0;


}

// This is basically a wrapper function for AttemptToGenerateNeistein08 (below).
// Sometimes a given attempt to generate such an accretion history will fail because 
// the history requested will produce a major merger.
// Previously, the program would spit out an error and crash, so if you requested 1000 runs
// you would get only some somewhat unpredictable fraction of that. Now we guarantee that we'll
// keep trying until we get an acceptable accretion history, and keep track of the number of attempts we required.
double AccretionHistory::GenerateNeistein08(double zst, Cosmology& cos,
        std::string fn, bool writeOut, unsigned long int seed,
        double invMassRatioLimit, double zquench, int * nattempts,
        double domega, double zrelax,double fscatter)
{
    *nattempts = 0; // we have yet to attempt to generate an accretion history.

    // set up the random number generator.
    gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r,seed);

    double mdotext0 = -1.0;
    while(mdotext0 < 0.0) {
        mdotext0 = AttemptToGenerateNeistein08(zst,cos,fn,writeOut,r,invMassRatioLimit,zquench,domega,zrelax,fscatter);
        ++(*nattempts);
    }
    gsl_rng_free(r);
    if( mdotext0<=0 )
        errormsg("Generating a neistein08 accretion history has produced negative mdot.");

    return mdotext0;
}


double AccretionHistory::GenerateAverageNMD10(double zst, Cosmology& cos,
        std::string fn, bool writeOut, unsigned long int seed,
        double invMassRatioLimit, double zquench, int nToAvg,
        double domega, double zrelax, double fscatter)
{
    int TotalAttempts=0;
    double mdotext0;
    int NSamples = 1000; // resolution of redshift grid.
    zstart = zst;

    // Set up a grid in redshift
    std::vector<double> zs(NSamples,0);
    for(unsigned int j=0; j!=NSamples; ++j) {
        zs[j] = zrelax*1.05 - j*(zrelax*1.05 + .01)/((double) NSamples);
    }

    std::vector<double> avgMh(NSamples,0);
    std::vector<double> avgMdots(NSamples,0);
    double avgMdotExt0 = 0.0;
    for(int i=0; i!=nToAvg; ++i) {
        int nattempts=0;
        AccretionHistory accr(Mh0, dbg); // Same as for *this.
        accr.SetEfficiencyParams(normalization, alpha_z, alpha_Mh,  ceiling);
        mdotext0 = accr.GenerateNeistein08(zst,cos,fn,false,i,invMassRatioLimit,zquench,&nattempts,domega,zrelax,fscatter);
        avgMdotExt0 += mdotext0/((double) nToAvg);
        TotalAttempts += nattempts;
        // Now that we have an accretion history, accumulate its mdot's and mh's.
        for(unsigned int j=0; j!=zs.size(); ++j) {
            avgMh[j] += accr.MhOfZ(zs[j])/((double) nToAvg);
            double currentAcc = accr.AccOfZ(zs[j]) * mdotext0;
            avgMdots[j] += currentAcc/((double) nToAvg);
            //sqMdots[j] += currentAcc*currentAcc/((double) nToAvg);
        }
    }

    if(writeOut) {
        std::ofstream file(fn.c_str());
        for(unsigned int j=0; j!=zs.size()-1; ++j) {
            file << zs[j] << " " <<  // redshift
                avgMh[j] << " " <<   // <M_h(t)>
                avgMdots[j] << " " << // <Mdot_b(t)>
                //sqrt(sqLogMdots[j] - avgLogMdots*avgLogMdots) << " "<< // variance
                0.17*epsin(zs[j],avgMh[j],cos,zquench) << " " << // efficiency factor
                0.17*epsin(zs[j],avgMh[j],cos,zquench)*fabs((avgMh[j+1]-avgMh[j])*1.0e12*speryear/(cos.Tsim(zs[j+1])-cos.Tsim(zs[j]))) << " " << std::endl; // Mdot_b(t) derived from avg. acc history
        }
        file.close();
    }

    linear=true;
    InitializeGSLObjs(zs,avgMdots,avgMh);
    if( avgMdotExt0<=0 )
        errormsg("Avg NMD10 accr. history has produced negative mdot0");
    return avgMdotExt0;
}

// Attempts to generate an accretion history based on Neistein08. If it fails, it will return -1,
// otherwise it will give you the gas mass accretion rate at z=zstart.
double AccretionHistory::AttemptToGenerateNeistein08(double zst, Cosmology& cos, 
        std::string fn, bool writeOut,gsl_rng * r,
        double invMassRatioLimit,double zquench,
        double domega, double zrelax, double fscatter )
{
    zstart = zst;
    std::ofstream file;
    if(writeOut) file.open(fn.c_str());

    // The following bit of logic allows us to use either ND08a or NMD10 to
    // generate the accretion histories. If we use ND08a, (dbg.opt(3)), we need
    // to use the Millenium Run cosmology to generate the accretion histories,
    // but we don't want other functions of the cosmology to use MR cosmology
    // So, we will use hisCos to generate the histories, but simCos in the rest of the simulation.
    // If !dbg.opt(3), both will be the same (presumably WMAP5)
    Cosmology MR(0.25, 0.75, 2.36576888e-18, 0.9, zrelax, 10, 0.0); // Millenium run cosmology, use arbitrary number for nx, concentration random factor
    Cosmology simCos(cos); // 
    Cosmology hisCos(MR);
    if(!dbg.opt(3))
        hisCos = Cosmology(simCos);
//    double sigma8=0.82;
//    if(!dbg.opt(3)) sigma8=.796; // WMAP5
    double sigma8 = hisCos.sigma8();
    double OmegaM=hisCos.OmegaM();
    // double dom = 0.1;
    double dom = domega;
    double zero=0.0;
    double om = omega(0.0,&zero);
    double omTr = om;
    linear=true; // use linear interpolation
    SParams sp;
    sp.rhs=0.0;
    sp.sigma8=sigma8;
    sp.OmegaM=OmegaM;
    std::vector<double> zs(0),accs(0),masses(0),zs2(0),accs2(0),masses2(0);
    double SS= S(Mh0*hisCos.h(),&sp);
    bool first=true;
    do {
        // pick a value for our Gaussian random variable.
        double x= gsl_ran_gaussian(r,1.0);
        double s= log10(SS);
        // employ the formula in ND08a to compute the delta S for the chosen value of x.
        double deltaS;
        if(dbg.opt(3)) {
            double sigp = (1.367+0.012*s+0.234*s*s);
            double mup = (-3.682+0.76*s-0.36*s*s);
            deltaS= exp(sigp*x*fscatter + mup + sigp*sigp*(1-fscatter*fscatter)/2.0);
        }
        else { // use wmap5 cosmology, w/ fit from Neistein2010
            double a1 = -0.333;
            double a2 = -0.321;
            double a3 =  0.0807;
            double a4 =  0.622;
            double b1 =  0.132;
            double b2 =  2.404;
            double b3 =  0.585;
            double b4 = -0.436;
            double sigp = (a1*s + a2) * log10(dom) + a3*s + a4;
            double mup  = (b1*s + b2) * log10(dom) + b3*s + b4;
            // The reason for the last term below is to maintain a constant average deltaS,
            // regardless of fscatter. Note that for f=1, we get the standard lognormal distr. back,
            // while for f=0, we get a relation with no scatter but the same mean as before, since.
            // the mean of a lognormal distribution is exp(mu+sig^2/2).
            deltaS= exp(x*sigp*fscatter + mup + sigp*sigp*(1-fscatter*fscatter)/2.0);
        }
        // update our list of redshifts, given our uniformly spaced value of omega.
        // subtract a small number so that when we evaluate some quantity at z=0, we don't get an interpolation error.

        double fac = 1.0;
        // if we're using the MR to generate histories, rescale omega's as recommended
        // in Neistein, Maccio and Dekel (2010), section 5.1
        if(dbg.opt(3)) 
            fac = 0.86;
            
        double z=zOfOmega(omTr); 
        zs.push_back(z);
        if(!first) zs.push_back(z+.000001); // z's cannot have identical values or the GSL interpreter chokes.
        if(SS > 500.0) {
            return -1.0;
        } 
        // update our tabulated value of halo mass.
        double M=MofS(SS,sp.sigma8,sp.OmegaM)/hisCos.h();
        if(M < 1.0e4)
            return -1.0;
        masses.push_back(M);
        if(!first)  masses.push_back(M); // solar masses
//        std::cout<<"z,om;S,M;ds,s,x: "<<zs[zs.size()-1]<<" "<<om<<"; "<<SS<<" "<<
//            masses[masses.size()-1]<<"; "<<deltaS<<" "<<s<<" "<<x<<std::endl;

        SS+=deltaS;
        om += dom;
        omTr += dom*fac;
        first=false;
    } while(zs[zs.size()-1]<zrelax*1.1);
    //  for(unsigned int i=0; i!=masses.size()-1; ++i) {
    double md0;
    for(unsigned int i=0; i<masses.size()-1; i+=2) {
        // Check to confirm that no large change in M has occurred this step.
        if((masses[i]-masses[i+1])/masses[i+1] > invMassRatioLimit  && zs[i]<zstart) {
            return -1.0; // if too-large a change has occurred, throw out this history.
        }

        // The accretion rate is the difference in DM halo masses between the two time steps 
        // times the cosmic baryon fraction, divided by the time between redshift steps 
        // This is in units of solar masses per year
        // NOTE that this equation DOES NOT INCLUDE the efficiency factor. It will be added in
        // momentarily in the loop over intermediate redshifts.
        double accr = (.17*(masses[i]-masses[i+1])/
                (fabs(simCos.Tsim(zs[i+1]) - simCos.Tsim(zs[i]))
                 /speryear)
                );
        // if the "first" flag above is working correctly, these two accretion rates should
        // occur at very different redshifts.
        accs.push_back(accr);
        accs.push_back(accr);

    }
    zs.pop_back();
    masses.pop_back();

    // before the following lines, z goes from 0->zstart
    reverse(zs.begin(),zs.end());
    reverse(accs.begin(),accs.end());
    reverse(masses.begin(),masses.end());

    unsigned int accsSZ = accs.size();
    unsigned int massSZ = masses.size();

    zs.push_back(-.0001); zs.push_back(-.02); zs.push_back(-.1);
    accs.push_back(accs[accsSZ-1]); accs.push_back(accs[accsSZ-1]); accs.push_back(accs[accsSZ-1]);
    masses.push_back(masses[massSZ-1]); masses.push_back(masses[massSZ-1]); masses.push_back(masses[massSZ-1]);

    // The purpose of the following for loop is to add a bunch of intermediate points to our M(z), Mdot(z)
    // histories. This is because with only a few points, the interpolation may lead to annoying errors.
    // In particular, even though Mh should grow linearly between different draws from the lognormal distr.,
    // the interpolation scheme will interpolate in /redshift/, not time, so Mh will grow linearly in redshift.
    // This is problematic because between big enough interpolation steps, Mh will diverge from the sum of 
    // all accreted matter, only to reconverge at the next interpolation point.
    // At this point zs is going from high redshift to low.
    for(unsigned int i=0; i!=zs.size()-1; ++i) {
        zs2.push_back(zs[i]);
        accs2.push_back(accs[i] * epsin(zs[i],masses[i]*1.0e-12,simCos,zquench));
        masses2.push_back(masses[i]);
        // The number of points to add between two known z's.
        unsigned int NExpand = 100;
        for(unsigned int j=1; j!=NExpand; ++j ) {
            // Generate a value of z, equally spaced between the two known z's.
            double z = zs[i] - ((double) j)/((double) NExpand) * (zs[i]-zs[i+1]);
            zs2.push_back(z);
            // Now let's make Mh(z) linear in time not redshift!
            double Mhz = masses[i] + (masses[i+1]-masses[i]) * (simCos.Tsim(z) - simCos.Tsim(zs[i]))/(simCos.Tsim(zs[i+1])-simCos.Tsim(zs[i]));
            masses2.push_back(Mhz);
            // and finally multiply the accretion rate by some efficiency which changes with redshfit and halo mass.
//            accs2.push_back(fabs((Mhz-masses2[masses2.size()-2]) / 
//                        (cos.Tsim(z)-cos.Tsim(zs2[zs2.size()-2])))*
//                        0.17*epsin(z,Mhz*1.0e-12,cos,zquench)*speryear);
            accs2.push_back(accs[i]*epsin(z,Mhz*1.0e-12,simCos,zquench));
            // When we cross zstart, set md0 equal to the current accretion rate.
            if((z - zstart)*(zs2[zs2.size()-2]-zstart) <= 0.0)
                md0 = accs2[accs2.size()-1];

        }
    }


    if(writeOut) {
        for(unsigned int i=0; i!=accs.size(); ++i) {
            file << zs[i] <<" "<<accs[i]<<std::endl;
        }
        file.close();
    }


    InitializeGSLObjs(zs2,accs2,masses2);

    // return accs[accs.size()-1] / AccOfZ(zs[zs.size()-1]);
    return md0;
}

void AccretionHistory::SetEfficiencyParams(double norm, double a_z, double a_Mh, double ceil)
{
    normalization = norm;
    alpha_z = a_z;
    alpha_Mh = a_Mh;
    ceiling = ceil;

}

// IMPORTANT: Mh in units of 10^12 solar masses
double AccretionHistory::epsin(double z, double Mh,Cosmology & cos, double zquench)
{
    // try a reasonably general formula with user-controlled values.
    // Use the same functional form as CAFG.
    double val = normalization * pow(1.0+z, alpha_z) * pow(Mh, alpha_Mh);
    if(val > ceiling) val=ceiling;
    if(z<zquench) val = 0.0;
    if (val<0 || val>1)
        errormsg("Nonphysical accretion efficiency!");
    return val;
}


double AccretionHistory::GenerateBoucheEtAl2009( double zs, Cosmology& cos, 
        std::string fn, bool writeOut, bool MhAtz0, 
        double zquench)
{
    double Mh012 = Mh0*1.0e-12;
    std::ofstream file;
    if(writeOut) file.open(fn.c_str());
    zstart = zs;
    unsigned int N=1000; 
    double Mh=Mh012; // units of 10^12 solar masses
    double z=zstart;
    double fbp18 = 17.0/18.0; // baryon fraction / 0.18
    double MdotExt0;
    std::vector<double> redshifts(0),tabulatedAcc(0),haloMass(0);

    // Loop over redshift from z=0 to z=zstart in N increments if Mh is specified at z=0
    // Otherwise, loop over redshift from zstart to z=0.
    for(unsigned int i=0; i<=N; ++i) {
        if(!MhAtz0) // if Mh is given at z=zstart, start from high redshift and go to z=0 
            z=((double) (N-i))/((double) N)*(zstart*1.1-0.0)-.001;
        else // if Mh is given at z=0, start from low redshift and go to z=zstart.
            z=((double) i)/((double) N) * (zstart*1.1 - 0.0)-.001;

        double deltaz = (zstart*1.1-0.0) / ((double) N);

        // Use the Bouche formula to tell us the dark matter accretion rate
        double dMh = 34.0 * pow(Mh,1.14)*pow(1.0+z,2.4) * 1.0e-12; // in 10^12 Msol/yr

        double zero = 0.0;

        if(dbg.opt(8))
            dMh = 0.623 * pow(Mh,1.0+0.14) // so far this is dM/domega
                * fabs((omega(z,&zero)-omega(z-deltaz,&zero)) // d omega
                        * speryear / (cos.Tsim(z)-cos.Tsim(z-deltaz)));  // 1/dt

        double MdotExt = dMh*1.0e12 *0.18*fbp18*epsin(z,Mh,cos,zquench); 

        

        // Use the analogous formula to tell us the baryonic accretion rate.
//        double MdotExt = 7.0 * epsin(z,Mh,cos,zquench) * fbp18 * pow(Mh,1.1)*pow(1+z,2.2); // in solar masses /year

        if((i==0 && !MhAtz0) || (i==N && MhAtz0)) // always set MdotExt0 to be MdotExt at z=2
            MdotExt0= MdotExt;

        haloMass.push_back(Mh); // solar masses

        if(true) { //z>zquench) {
            // Basically compute Mh(z) by taking an Euler step, since from the above we know dMh (which is actually dMh/dt)
            if(!MhAtz0) // starting from high redshift..
                Mh+= dMh* -1.0*( cos.Tsim(z) - cos.Tsim( ((double) (N-i-1))/((double) N) * (zstart*1.1-0.0)-.001)) / speryear;
            else // starting from low redshift..
                Mh+= dMh* -1.0*( cos.Tsim(z) - cos.Tsim( ((double) (i+1))/((double) N) * (zstart*1.1-0.0)-.001)) / speryear;
        }

        // these small adjustments to z avoid interpolation errors without really affecting anything.
        redshifts.push_back(z); 

        tabulatedAcc.push_back(MdotExt);
        if(writeOut) file << z << " "<< cos.Tsim(z) <<" "<<MdotExt<<" "<<Mh<<std::endl;
    }

    // 
    if(MhAtz0) { // then we need to reverse redshifts and tabulatedAcc
        reverse(redshifts.begin(),redshifts.end());
        reverse(tabulatedAcc.begin(), tabulatedAcc.end());
        reverse(haloMass.begin(),haloMass.end());
    }

    file.close();
    InitializeGSLObjs(redshifts,tabulatedAcc,haloMass);
    if( MdotExt0<=0 )
        errormsg("Bouche accr. history has produced negative mdot0");
    return MdotExt0;
}

// redshifts should range from zstart to 0, tabulatedAcc, haloMass, and redshifts have 
// to have the same size; redshifts must start high and end low
void AccretionHistory::InitializeGSLObjs(std::vector<double> redshifts, std::vector<double> tabulatedAcc,
        std::vector<double> haloMass)
{
    if(redshifts.size()!=tabulatedAcc.size() || redshifts.size()!=haloMass.size())
        errormsg("InitializeGSLObjs: tabulatedAcc and redshifts must be the same size.");

    allocated=true;
    accel=gsl_interp_accel_alloc();
    accelMh=gsl_interp_accel_alloc();
    if(!linear) {
        spline=gsl_spline_alloc(gsl_interp_cspline,tabulatedAcc.size());
        splineMh=gsl_spline_alloc(gsl_interp_cspline,haloMass.size());
    }
    else {
        spline=gsl_spline_alloc(gsl_interp_linear,tabulatedAcc.size());
        splineMh = gsl_spline_alloc(gsl_interp_linear,haloMass.size());
    }
    acc = new double[tabulatedAcc.size()];
    redshift = new double[redshifts.size()];
    hMass = new double[haloMass.size()];
    unsigned int tabAcc=tabulatedAcc.size();
    //  double norm = tabulatedAcc[0];
    //  // Normalize accretion rate to the highest-redshift value given.
    //  for(unsigned int i=0; i<tabAcc; ++i) {
    //    tabulatedAcc[i]/=norm;
    //  }

    // reverse order!
    for(unsigned int i=0; i<tabAcc; ++i) {
        acc[tabAcc-i-1]=tabulatedAcc[i];
        redshift[tabAcc-i-1]=redshifts[i];
        hMass[tabAcc-i-1]=haloMass[i];
    }

    gsl_spline_init(spline,redshift,acc, tabulatedAcc.size());  
    gsl_spline_init(splineMh,redshift,hMass, haloMass.size());

}


// fn=filename
// zc= redshift column
// acC= accretion column
// rs = rows skipped before data is read
// nc = total number of columns
void AccretionHistory::ReadTabulated(std::string fn,unsigned int zc, unsigned int acC, 
        unsigned int rs, unsigned int nc, double zst)
{
    zstart = zst;
    std::vector<double> redshifts(0);
    std::vector<double> tabulatedAcc(0);
    std::vector<double> haloMass(0);
    std::ifstream f(fn.c_str());
    if(!f.is_open()) errormsg("Error opening file containing the tabulated accretion history!");
    std::string line;
    for(unsigned int i=0; i<rs; ++i) {
        getline(f,line); // get a line and discard it
    }
    double nm;
    while(f.good() ) {
        for(unsigned int i=0; i<nc; ++i) {  // for each column..
            f >> nm; // read the number
            if(i+1 == zc) { // if we're in the redshift column, ..
                if(redshifts.size()==0 && zstart >= zc) redshifts.push_back(zstart);
                else if(redshifts.size()==0) redshifts.push_back(nm*1.0001);
                redshifts.push_back(nm); // 
            }
            if(i+1 == acC) {
                tabulatedAcc.push_back(nm);
                if(tabulatedAcc.size()==1) tabulatedAcc.push_back(nm);
            }
        }
    }
    tabulatedAcc.pop_back();
    redshifts.pop_back();
    if(redshifts[redshifts.size()-1] != 0.0) {
        tabulatedAcc.push_back(tabulatedAcc[tabulatedAcc.size()-1]);
        redshifts.push_back(0.);
    }

    f.close();

    InitializeGSLObjs( redshifts,tabulatedAcc, haloMass);
}

double AccretionHistory::AccOfZ(double z)
{
    //  if(z > zstart) errormsg("AccOfZ: The accretion rate at this z is unspecified. z and zstart are "
    //	+str(z)+" "+str(zstart));
    double val =  gsl_spline_eval(spline,z,accel)/gsl_spline_eval(spline,zstart,accel);
    if(val<0.0) {
        // turned out it was, e.g. numbers on the order of 10^-200
        //    std::cout << "WARNING: setting accr rate = 0. Hopefully this negative value "<<val<<" is a small error because of the interpolation scheme + the sharp cutoff from zquench." << std::endl;
        val = 0.0;

    }
    return val;
}

double AccretionHistory::MhOfZ(double z)
{
    //if(z>zrelax) errormsg("MhOfZ: The halo mass at this z is unspecified.");

    // Return the halo mass in units of Mh(z=0). This gives callers of InitializeGSLObjs the option
    // to specify their list of halo masses in whatever units they so choose.
    return gsl_spline_eval(splineMh,z,accelMh)/gsl_spline_eval(splineMh,0.0,accelMh);
}

// This function should be implemented, but caused unforeseen problems with flagging any production model
// which called it as having failed (by having written to std::cerr; these should be resolved now.
void testAccretionHistory()
{
    //Cosmology cos(1.0-.734, .734, 2.29e-18, 2.0);


    //  for(unsigned int whichAccretionHistory=10; whichAccretionHistory!=1000; ++whichAccretionHistory) {
    //    AccretionHistory accr(1.0e12);
    //    double Mdot0 = accr.GenerateNeistein08(2.0, cos, "", false, whichAccretionHistory, .3, false);
    //    if(Mdot0>0.0) {
    //      double dummy = 0.0;
    //
    //
    //    }
    //  }


}

















