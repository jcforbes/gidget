#include "DiskContents.h"
#include "StellarPop.h"
#include "RafikovQParams.h"
#include "Cosmology.h"
#include "Dimensions.h"
#include "Deriv.h"
#include "DiskUtils.h"
#include "Simulation.h"
#include "FixedMesh.h"
#include "Interfaces.h"
#include "Debug.h"
#include "AccretionHistory.h"

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include <iostream>
#include <fstream>

double PosOnly(double input)
{
    if(input>0.0)
        return input;
    else
        return 0.0;
}
double NegOnly(double input)
{
    return -PosOnly(-input);
}

double ddxUpstream(std::vector<double>& vec, std::vector<double>& x, std::vector<int>& flags, unsigned int n)
{
    if(n>1 && n<x.size()-1) {
        return (flags[0]*vec[n-1] + flags[1]*vec[n] + flags[2]*vec[n+1]) / (flags[0]*x[n-1] + flags[1]*x[n] + flags[2]*x[n+1]);
    }

    if(n==1) {
        if(flags[0]==0) // unambiguous what to do:
            return (vec[n+1]-vec[n])/(x[n+1]-x[n]);
        else
            return 0; // shrug
    }
    if(n==x.size()-1) {
        if(flags[2]==0) // unambiguous:
            return (vec[n]-vec[n-1])/(x[n]-x[n-1]);
        else
            return 0;
    }

    std::cerr << "WARNING: unexpected result in ddxUpstream" << std::endl;
    return 0;
}

void ComputeUpstreamFlags(std::vector<int>& flags, std::vector<double>& mdot1, std::vector<double>& mdot2, unsigned int n)
{
    double mdotL = mdot1[n-1] + mdot2[n-1];
    double mdotR = mdot1[n] + mdot2[n];

    double theMdot;
    if(fabs(mdotL) > fabs(mdotR)) // the left-hand flux is dominant
        theMdot=mdotL;
    else
        theMdot=mdotR;


    if(theMdot < 0) {// going towards +r
        flags[2]=0;
        flags[1]=1;
        flags[0]=-1;
    }
    else { // going towards -r
        flags[0]=0;
        flags[1]=-1;
        flags[2]=1;
    }

}

void ComputeFlagList(std::vector<double>& mdot1, std::vector<double>& mdot2, std::vector<std::vector<int> >& flagList)
{
    for(unsigned int n=1; n<=mdot2.size()-1; ++n) {
        ComputeUpstreamFlags(flagList[n],mdot1,mdot2,n);
    }
}

// Fill an initializer object with the current state of this disk
void DiskContents::store(Initializer& in)
{
    in.col.resize(nx+1);
    in.sig.resize(nx+1);
    in.col_st.resize(nx+1);
    in.sig_stR.resize(nx+1);
    in.sig_stZ.resize(nx+1);
    for(unsigned int n=1; n<=nx; ++n) {
        in.col[n] = col[n];
        in.sig[n] = sig[n];
        in.col_st[n] = activeColSt(n);
        in.sig_stR[n] = activeSigStR(n);
        in.sig_stZ[n] = activeSigStZ(n);
    }
}

DiskContents::DiskContents(double tH, double eta,
        double sflr,double epsff,
        double ql,double tol,
        bool aq, double mlf, double mlfColScal,
        double mlfFgScal,
        Cosmology& c,Dimensions& d,
        FixedMesh& m, Debug& ddbg,
        double thk, bool migP,
        double Qinit, double km,
        unsigned int NA, unsigned int NP,
        double minSigSt,
		double rfrec, double xirec,
		double fh2min, double tdeph2sc,
        double ZIGM, double yrec,
        double ksup, double kpow,
        double mq, double muq,
        double Zmx) :
    nx(m.nx()),x(m.x()),beta(m.beta()),
    uu(m.uu()), betap(m.betap()),
    uDisk(std::vector<double>(m.nx()+1,0.)),
    uDM(std::vector<double>(m.nx()+1,0.)),
    uBulge(std::vector<double>(m.nx()+1,0.)),
    dim(d), mesh(m), dbg(ddbg),
    XMIN(m.xmin()),ZDisk(std::vector<double>(m.nx()+1,Z_IGM)),
    cos(c),tauHeat(tH),sigth(sflr),
    EPS_ff(epsff),ETA(eta),constMassLoadingFactor(mlf),
    mlfColScaling(mlfColScal),
    mlfFgScaling(mlfFgScal),
    MQuench(mq), muQuench(muq),
    //  spsActive(std::vector<StellarPop>(NA,StellarPop(m.nx(),0,c.lbt(1000)))),
    //  spsPassive(std::vector<StellarPop>(NP,StellarPop(m.nx(),0,c.lbt(1000)))),
    spsActive(std::vector<StellarPop*>(0)),
    spsPassive(std::vector<StellarPop*>(0)),
    dlnx(m.dlnx()),Qlim(ql),TOL(tol),ZBulge(ZIGM),
    yREC(yrec),RfREC(rfrec),xiREC(xirec), 
    analyticQ(aq),
    tDepH2SC(tdeph2sc),
    fH2Min(fh2min),
    ksuppress(ksup), kpower(kpow),
    thickness(thk), migratePassive(migP),
    col(std::vector<double>(m.nx()+1,0.)),
    sig(std::vector<double>(m.nx()+1,0.)),
    dsigdtTrans(std::vector<double>(m.nx()+1,0.)),
    dsigdtDdx(std::vector<double>(m.nx()+1,0.)),
    dsigdtHeat(std::vector<double>(m.nx()+1,0.)),
    dsigdtCool(std::vector<double>(m.nx()+1,0.)),
    dcoldtIncoming(std::vector<double>(m.nx()+1,0.)),
    dcoldtOutgoing(std::vector<double>(m.nx()+1,0.)),
    dQdu(std::vector<double>(m.nx()+1,0.)), 
    dudt(std::vector<double>(m.nx()+1,0.)), 
    dQdS(std::vector<double>(m.nx()+1,0.)), 
    dQds(std::vector<double>(m.nx()+1,0.)), 
    dQdSerr(std::vector<double>(m.nx()+1)),
    dQdserr(std::vector<double>(m.nx()+1,0.)),
    dcoldtPrev(std::vector<double>(m.nx()+1,0.)),
    dcoldt(std::vector<double>(m.nx()+1,0.)),
    dsigdtPrev(std::vector<double>(m.nx()+1,0.)),
    dsigdt(std::vector<double>(m.nx()+1,0.)),
    dZDiskdtPrev(std::vector<double>(m.nx()+1,0.)),
    dZDiskdtDiff(std::vector<double>(m.nx()+1,0.)),
    dZDiskdtAdv(std::vector<double>(m.nx()+1,0.)),
    dZDiskdt(std::vector<double>(m.nx()+1,0.)),
    dMZdt(std::vector<double>(m.nx()+1,0.)),
    colSFR(std::vector<double>(m.nx()+1,0.)),
    ColOutflows(std::vector<double>(m.nx()+1,0.)),
    MassLoadingFactor(std::vector<double>(m.nx()+1,0.)),
    mBubble(std::vector<double>(m.nx()+1,0.)),
    keepTorqueOff(std::vector<int>(m.nx()+1,0)),
    diffused_dcoldt(std::vector<double>(m.nx()+1,0.)),
    //yy(std::vector<double>(m.nx()+1,0.)),
    CumulativeSF(std::vector<double>(m.nx()+1,0.)),
    CumulativeTorqueErr2(std::vector<double>(m.nx()+1,0.)),
    CumulativeTorqueErr(std::vector<double>(m.nx()+1,0.)),
    d2taudx2(std::vector<double>(m.nx()+1,0.)),
    initialStellarMass(0.0),initialGasMass(0.0),
    cumulativeMassAccreted(0.0),
    cumulativeStarFormationMass(0.0),
    cumulativeGasMassThroughIB(0.0),
    cumulativeStellarMassThroughIB(0.0),
    CuStarsOut(std::vector<double>(m.nx()+1,0.)), 
    CuGasOut(std::vector<double>(m.nx()+1,0.)),
    fH2(std::vector<double>(m.nx()+1,0.5)),
    G0(std::vector<double>(m.nx()+1,1.0)),
    colvPhiDisk(std::vector<double>(m.nx()+1,0.0)),
    colstvPhiDisk(std::vector<double>(m.nx()+1,0.0)),
    dampingFactors(std::vector<double>(m.nx()+1,0.3)),
    //FF(std::vector<double>(m.nx()+1,0.)),
    // DD(std::vector<double>(m.nx()+1,0.)),
    // LL(std::vector<double>(m.nx()+1,0.)),
    // UU(std::vector<double>(m.nx()+1,0.)),
    // MdotiPlusHalf(std::vector<double>(m.nx()+1,0.)),
    fixedQ(Qinit),CumulativeTorque(0.0),
    kappaMetals(km),
    minsigst(minSigSt),
    NActive(NA),
    NPassive(NP),
    dd(exp(m.dlnx())),    // define the quantity d (a number slightly larger than 1)
    dm1(expm1(m.dlnx())), // dm1 = d-1. Use expm1 to find this to high precision.
    dmm1(-expm1(-m.dlnx())),   // dmm1 = 1 -d^-1
    dmdinv(expm1(2.*m.dlnx())/exp(m.dlnx())),  // dmdinv = d -d^-1
    sqd(exp(m.dlnx()/2.)), // sqd = square root of d
    Z_IGM(ZIGM),
    ZMix(Zmx)
{
    accel_colst = gsl_interp_accel_alloc();
    accel_sigst = gsl_interp_accel_alloc();
    spline_colst = gsl_spline_alloc(gsl_interp_cspline,nx);
    spline_sigst = gsl_spline_alloc(gsl_interp_cspline,nx);
    colst_gsl = new double[nx];
    sigst_gsl = new double[nx];

    lr=gsl_vector_alloc(nx+1); 
    diag=gsl_vector_alloc(nx+2); 
    ur=gsl_vector_alloc(nx+1);
    tau=gsl_vector_alloc(nx+2); 
    forcing=gsl_vector_alloc(nx+2);


    return;
}

DiskContents::~DiskContents()
{
    gsl_spline_free(spline_colst);
    gsl_spline_free(spline_sigst);
    gsl_interp_accel_free(accel_colst);
    gsl_interp_accel_free(accel_sigst);
    delete[] colst_gsl;
    delete[] sigst_gsl;

    gsl_vector_free(lr); gsl_vector_free(diag); gsl_vector_free(ur);
    gsl_vector_free(tau); gsl_vector_free(forcing);

    for(unsigned int i=0; i!=spsActive.size(); ++i) {
        delete spsActive[i];
    }
    spsActive.clear();
    for(unsigned int i=0; i!=spsPassive.size(); ++i) {
        delete spsPassive[i];
    }
    spsPassive.clear();
}


// Initialize the disk with an Initializer object, basically all the state variables
// from a disk, to which you've done whatever... (typically relaxed to an equilibrium configuration).
void DiskContents::Initialize(Initializer& in, bool fixedPhi0)
{
    StellarPop * initialStarsA = new StellarPop(mesh);
    StellarPop * initialStarsP = new StellarPop(mesh);


    for(unsigned int n=1; n<=nx; ++n) {
        ZDisk[n] = Z_IGM;

        col[n] = in.col[n];
        sig[n] = in.sig[n];
        initialStarsA->spcol[n] = in.col_st[n];
        initialStarsA->spsigR[n] = in.sig_stR[n];
        initialStarsA->spsigZ[n] = in.sig_stZ[n];
        initialStarsA->spZ[n]   = Z_IGM;
        initialStarsA->spZV[n]  = 0.0;

        initialStarsP->spcol[n] = initialStarsA->spcol[n];
        initialStarsP->spsigR[n] = initialStarsA->spsigR[n];
        initialStarsP->spsigZ[n] = initialStarsA->spsigZ[n];
        initialStarsP->spZ[n]   = initialStarsA->spZ[n];
        initialStarsP->spZV[n]  = initialStarsA->spZV[n];
    }

    // MBulge here is dimensionless:
    MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA->spcol[1]);
    MHalo = 0.0;
    initialStarsA->ageAtz0 = cos.lbt(cos.ZStart());
    initialStarsP->ageAtz0 = cos.lbt(cos.ZStart());

    initialStarsA->ComputeSpatialDerivs();
    initialStarsP->ComputeSpatialDerivs();

    //  spsActive.clear();
    //  spsPassive.clear();
    spsActive.push_back(initialStarsA);
    spsPassive.push_back(initialStarsP);
    EnforceFixedQ(fixedPhi0,true);

    initialStellarMass = TotalWeightedByArea(initialStarsA->spcol) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
    initialGasMass = TotalWeightedByArea(col) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}


// Simplify things. Just put in exponential disks with constant velocity dispersions.
// Set Q=Qf only when these simple initial conditions yield Q<Qf.
void DiskContents::Initialize(double fcool, double fg0,
        double sig0, double phi0, double Mh0,
        double MhZs, double stScaleLength, double zs, const double stScaleReduction, const double gaScaleReduction)
{
    StellarPop * initialStarsA = new StellarPop(mesh);
    StellarPop * initialStarsP = new StellarPop(mesh);

    double fg4 = fg0;
    double fc = fcool;
    double Z0 = Z_IGM;
    double xdstars = stScaleLength/dim.d(1.0)/stScaleReduction; // 2.0 -- make the initial size much smaller than the initial accr. radius
    double xdgas = stScaleLength/dim.d(1.0)/gaScaleReduction; // 1.4 -- make the initial size much smaller than the initial accr. radius
    // if S = f_g S0 exp(-x/xd), this is S0 such that the baryon budget is maintained, given that a fraction
    // fcool of the baryons have cooled to form a disk.
    // This is a correction, to account for the fact that for scale lengths >~ the radius of the disk,
    // most of this mass will fall off the edge of the grid. We want to put that mass back into the grid.
    double xouter = mesh.x(.5 + ((double) nx));
    double xinner = mesh.x(.5);
    xinner = 0.0; // don't change the column density depending on inner edge
    xouter = 1000.0;
    // double dummy = -exp(-xouter/xd)*xd*(xd+xouter) + exp(-xinner/xd)*xd*(xd+xinner);
    double eff;

    // override the inputs to this function for fcool and fg0 based on our knowledge of the actual halo mass in this galaxy at z=4.
    if(dbg.opt(10)) {
        double Mh = MhZs;
        double M10 = 11.590;
        double M11 = 1.195;
        double N10 = 0.0351;
        double N11 = -0.0247;
        double beta10= 1.376;
        double beta11 = -0.826;
        double gamma10 = 0.608;
        double gamma11 = 0.329;
        double zti = 4.0;
        double logM1z = M10 + M11*zti/(zti+1.0);
        double Nz = N10 + N11*zti/(zti+1.0);
        double betaz = beta10 + beta11*zti/(zti+1.0);
        double gammaz = gamma10 + gamma11*zti/(zti+1.0);
        double M1 = pow(10.0, logM1z);
        eff = 2.0*Nz / (pow(Mh/M1,-betaz) + pow(Mh/M1,gammaz));
        double mst = eff*Mh; // mstar according to the moster relation. ## at z=4 !!
        double f0 = 1.0/(1.0 + pow(mst/pow(10.0,9.15),0.4)); // from Hayward & Hopkins (2015) eq. B2
        double tau4 = 12.27/(12.27+1.60); // fractional lookback time at z=4
        double fgz4 = f0*pow(1.0 - tau4*(1.0-pow(f0,1.5)), -2.0/3.0);
        double reff4 = 5.28*pow(mst/1.0e10, 0.25)*pow(1.0+4.0,-0.6); // kpc (eq B3) at z=4
        double ZHayward = -8.69 + 9.09*pow(1.0+4.0,-0.017) - 0.0864*pow(log10(mst) - 11.07*pow(1.0+4.0,0.094),2.0);
        ZHayward = pow(10.0, ZHayward) * 0.02;
        Z0 = ZHayward;
        fc = 1.0/(1.0-fgz4) * mst/(0.17*Mh);
        fg4 = fgz4;
        xdstars = reff4/dim.d(1.0)/1.67835; // is scale length = reff? Numerical solution says reff/rd = 1.67835.
        xdgas= xdstars*2.0;
 
    }

    double S0gas = 0.17*fc*fg4*MhZs*MSol*dim.vphiR / (2.*M_PI*(-exp(-xouter/xdgas)*xdgas*(xdgas+xouter) + exp(-xinner/xdgas)*xdgas*(xdgas+xinner))*dim.MdotExt0*dim.Radius);
    double S0stars = 0.17*fc*(1-fg4)*MhZs*MSol*dim.vphiR / (2.*M_PI*(-exp(-xouter/xdstars)*xdstars*(xdstars+xouter) + exp(-xinner/xdstars)*xdstars*(xdstars+xinner))*dim.MdotExt0*dim.Radius);
    double floorval = 1.0e-15;

    for(unsigned int n=1; n<=nx; ++n) {
        ZDisk[n] = Z0;
        initialStarsA->spcol[n] = S0stars*exp(-x[n]/xdstars);
        if(initialStarsA->spcol[n] < S0stars*floorval)
            initialStarsA->spcol[n] = S0stars*floorval; // floor the initial value to avoid very small timesteps
        initialStarsA->spsigR[n] = max(sig0 * phi0, minsigst);
        initialStarsA->spsigZ[n] = max(sig0 * phi0, minsigst);
        initialStarsA->spZ[n] = Z_IGM;
        initialStarsA->spZV[n] = 0.0;

        initialStarsP->spcol[n] = initialStarsA->spcol[n];
        initialStarsP->spsigR[n] = initialStarsA->spsigR[n];
        initialStarsP->spsigZ[n] = initialStarsA->spsigZ[n];
        initialStarsP->spZ[n] = initialStarsA->spZ[n];
        initialStarsP->spZV[n] = initialStarsA->spZV[n];

        col[n] = S0gas*exp(-x[n]/xdgas);
        if(col[n] < S0gas*floorval)
            col[n] = S0gas*floorval; // floor the initial value to avoid very small timesteps
        sig[n] = max(sig0, sigth);


    }
    //std::cout << "DEBUG ICs.  "<< MhZs << " "<< eff<<" "<< xd<<" "<< " "<< stScaleLength/dim.d(1.0)<<" " << col[1]<< " "<<col[nx]<<  std::endl;
    
    ComputeMassLoadingFactor(MhZs, initialStarsA->spcol);

    // In simulation units, mass = M(dimensional) / [ Mdotext * 2pi R/vphiR ]
    // Ergo mbulge(dimensionless) = pi r0^2 \Sigma0 / [ Mdotext * 2pi R/vphiR]
    //                            = pi x0^2 R^2 S0 Mdotext/(vphiR R) / [ Mdotext *2pi R/vphiR]
    //                            = 0.5 x0^2 S0
    // MBulge = M_PI*x[1]*x[1]*(col[1]*RfREC/(MassLoadingFactor[1]+RfREC)+initialStarsA->spcol[1]); // dimensionless!
    MBulge = 0.5*x[1]*x[1]*(col[1]*RfREC/(MassLoadingFactor[1]+RfREC)+initialStarsA->spcol[1]); // dimensionless!
    MHalo = 0.0;
    initialStarsA->ageAtz0 = cos.lbt(cos.ZStart());
    initialStarsP->ageAtz0 = cos.lbt(cos.ZStart());

    initialStarsA->ComputeSpatialDerivs();
    initialStarsP->ComputeSpatialDerivs();


    spsActive.push_back(initialStarsA);
    spsPassive.push_back(initialStarsP);
    bool fixedPhi0 = true;
    bool EnforceWhenQgtrQf = false;

    // before we heat up the disk, update the rotation curve.
    UpdateRotationCurve(MhZs, zs, 1.0e-10);

    EnforceFixedQ(fixedPhi0,EnforceWhenQgtrQf);

    initialStellarMass = TotalWeightedByArea(initialStarsA->spcol) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
    initialGasMass = TotalWeightedByArea(col) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;

}


// Mh0 in solar masses
// sigst0 is in units of vphiR. stScaleLength is, as usual in units of kpc, as is BulgeRadius.
void DiskContents::Initialize(double fcool, double fg0,
        double sigst0, double Mh0,double MhZs,
        double stScaleLength)
{			     
    StellarPop * initialStarsA = new StellarPop(mesh);
    StellarPop * initialStarsP = new StellarPop(mesh);


    if(true) {
        //  InitializeGrid(BulgeRadius);

        double maxsig=0.0;
        unsigned int maxsign=1;
        bool lowQst;

        double fc = 1.0/6.0;
        double xd = stScaleLength/dim.d(1.0);
        double S0 = 0.17 * fc * (1-fg0) * MhZs*MSol/(dim.MdotExt0) * dim.vphiR/(2.0*M_PI*dim.Radius) * (1.0/(xd*xd));

        //// For each cell...
        for(unsigned int n=1; n<=nx; ++n) {
            ZDisk[n] = Z_IGM;

            // Put in the exponential stellar profile.
            initialStarsA->spcol[n] = S0 *exp(-x[n]/xd);
            initialStarsA->spsigR[n] = max(sigst0,minsigst);
            initialStarsA->spsigZ[n] = max(sigst0,minsigst);



            // if the initial conditions are such that Q_* < Q_lim, set Q_*=Q_lim by 
            // heating the stars beyond what the user requested.
            if(sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA->spsigR[n] / (initialStarsA->spcol[n]*M_PI*x[n]*dim.chi()) < Qlim) {
                lowQst = true;
                initialStarsA->spsigR[n] = max(Qlim*M_PI*x[n]*initialStarsA->spcol[n]*dim.chi()/(sqrt(2.*(beta[n]+1.))*uu[n]),minsigst);
            }

            if(initialStarsA->spsigR[n] > maxsig) { 
                maxsig=initialStarsA->spsigR[n];
                maxsign=n;
            }
            initialStarsA->spZ[n] = Z_IGM;
            initialStarsA->spZV[n] = 0.0;

            initialStarsP->spcol[n] = initialStarsA->spcol[n];
            initialStarsP->spsigR[n] = initialStarsA->spsigR[n];
            initialStarsP->spsigZ[n] = initialStarsA->spsigZ[n];
            initialStarsP->spZ[n] = initialStarsA->spZ[n];
            initialStarsP->spZV[n] = initialStarsA->spZV[n];

            sig[n] = max(pow(dim.chi() / (ETA*fg0), 1./3.)/sqrt(2.),  sigth);

            col[n] = ((thickness/fixedQ)*uu[n]*sqrt(2.*(beta[n]+1.))/(M_PI*dim.chi()*x[n]) - initialStarsA->spcol[n]/initialStarsA->spsigR[n]) *sig[n];


            if(col[n]<0 || sig[n] <0 || col[n]!=col[n] || sig[n]!=sig[n] || initialStarsA->spcol[n]<0.0 
                    || initialStarsA->spsigR[n]<0.0 || initialStarsA->spcol[n]!=initialStarsA->spcol[n] 
                    || initialStarsA->spsigR[n]!=initialStarsA->spsigR[n]) 
            {
                errormsg("Error initializing disk- nonphysical state vars: n, col, sig, spcol, spsig, Qst: "
                        +str(n)+" "+str(col[n])+" "+str(sig[n])+" "+str(initialStarsA->spcol[n])+" "
                        +str(initialStarsA->spsigR[n])+" "
                        +str(sqrt(2.*(beta[n]+1.))*uu[n]*initialStarsA->spsigR[n]
                            /(M_PI*x[n]*initialStarsA->spcol[n]*dim.chi())));
            }
        } // end loop over all cells

        if(lowQst) {
            // make sig_st monotonically increasing towards the center of the disk.
            for(unsigned int n=1; n<=maxsign; ++n ) {
                if(initialStarsA->spsigR[n] < maxsig) { 
                    initialStarsA->spsigR[n] = max(maxsig,minsigst);
                    initialStarsP->spsigR[n] = max(maxsig,minsigst);
                    initialStarsA->spsigZ[n] = max(maxsig,minsigst);
                    initialStarsP->spsigZ[n] = max(maxsig,minsigst);

                }
            }

        }

        double minQst=1.0e30;
        unsigned int minQstN=0;
        // locate the minimum value of Q_*
        for(unsigned int n=1; n<=nx; ++n) {
            if(sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA->spsigR[n]/(initialStarsA->spcol[n]*M_PI*x[n]*dim.chi()) < minQst) { 
                minQst = sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA->spsigR[n] / (initialStarsA->spcol[n]*M_PI*x[n]*dim.chi());
                minQstN = n;
            }
        }
        if(minQst < Qlim*.99999) errormsg("Minimum Qst is somehow below Qlim. "+str(Qlim)+" "+str(minQst));

        // Inwards of minQstN, set col_* such that Q_* doesn't turn back up at small radii.
        for(unsigned int n=1; n<=minQstN; ++n) {
            initialStarsA->spcol[n] = sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA->spsigR[n] / (minQst*M_PI*x[n]*dim.chi());
            initialStarsP->spcol[n] = initialStarsA->spcol[n];
        }
        for(unsigned int n=1; n<=nx; ++n) {
            initialStarsA->spsigR[n] = max(initialStarsA->spsigR[n] * Qlim / minQst,minsigst);
            initialStarsP->spsigR[n] = initialStarsA->spsigR[n];
            initialStarsA->spsigZ[n] = max(initialStarsA->spsigZ[n] * Qlim / minQst,minsigst);
            initialStarsP->spsigZ[n] = initialStarsA->spsigZ[n];
        }
    }



    MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA->spcol[1]); // dimensionless!
    MHalo = 0.0;
    initialStarsA->ageAtz0 = cos.lbt(cos.ZStart());
    initialStarsP->ageAtz0 = cos.lbt(cos.ZStart());

    initialStarsA->ComputeSpatialDerivs();
    initialStarsP->ComputeSpatialDerivs();

    spsActive.push_back(initialStarsA);
    spsPassive.push_back(initialStarsP);
    EnforceFixedQ(false,true);

    initialStellarMass = TotalWeightedByArea(initialStarsA->spcol) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
    initialGasMass = TotalWeightedByArea(col) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}

void DiskContents::Initialize(double tempRatio, double fg0)
{
    // an active and passive stellar population..
    StellarPop * initialStarsA = new StellarPop(mesh);
    StellarPop * initialStarsP = new StellarPop(mesh);

    //  InitializeGrid(BulgeRadius);

    // Metallicity
    for(unsigned int n=1; n<=nx; ++n) {
        ZDisk[n] = Z_IGM;
        sig[n] = pow(dim.chi()/(ETA*fg0),1./3.)/sqrt(2.);
        col[n] = (thickness/fixedQ)*uu[n]*sqrt(2.*(beta[n]+1.))*sig[n]*
            tempRatio/(x[n]*M_PI*dim.chi() * (tempRatio + (1.-fg0)/fg0));
        initialStarsA->spcol[n] = col[n]*(1.-fg0)/fg0;
        initialStarsA->spsigR[n] = max(tempRatio*sig[n],minsigst);
        initialStarsA->spsigZ[n] = initialStarsA->spsigR[n];
        initialStarsA->spZ[n]   = Z_IGM;
        initialStarsA->spZV[n]  = 0.0;
        initialStarsP->spcol[n] = initialStarsA->spcol[n];
        initialStarsP->spsigR[n] = initialStarsA->spsigR[n];
        initialStarsP->spsigZ[n] = initialStarsA->spsigZ[n];
        initialStarsP->spZ[n]   = initialStarsA->spZ[n];
        initialStarsP->spZV[n]  = initialStarsA->spZV[n];

        if(col[n]<0 || sig[n] <0 || col[n]!=col[n] || sig[n]!=sig[n] || initialStarsA->spcol[n]<0.0 || initialStarsA->spsigR[n]<0.0 || initialStarsA->spcol[n]!=initialStarsA->spcol[n] || initialStarsA->spsigR[n]!=initialStarsA->spsigR[n]) {
            errormsg("Error initializing disk- nonphysical state vars: n, col, sig, spcol, spsig, Qst: "+str(n)+" "+str(col[n])+" "+str(sig[n])+" "+str(initialStarsA->spcol[n])+" "+str(initialStarsA->spsigR[n])+" "+str(sqrt(2.*(beta[n]+1.))*uu[n]*initialStarsA->spsigR[n]/(M_PI*x[n]*initialStarsA->spcol[n]*dim.chi())));
        }

    }

    MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA->spcol[1]); // dimensionless!
    MHalo = 0.0;
    initialStarsA->ageAtz0 = cos.lbt(cos.ZStart());
    initialStarsP->ageAtz0 = cos.lbt(cos.ZStart());
    initialStarsA->ComputeSpatialDerivs();
    initialStarsP->ComputeSpatialDerivs();
    spsActive.push_back(initialStarsA);
    spsPassive.push_back(initialStarsP);
    //  EnforceFixedQ(true); // this is no longer reasonable given the floor, minsigst.
    bool fixedPhi0 = initialStarsA->spsigR[1] > 2.0*minsigst;
    EnforceFixedQ(fixedPhi0,true); // only allow sig and sigst to covary if initial stellar velocity dispersion is far above minsigst.
    if(!fixedPhi0)  std::cerr << "WARNING: minsigst set too high to allow intiial conditions to be set by covarying gas and stellar velocity dispersions.";

    initialStellarMass = TotalWeightedByArea(initialStarsA->spcol) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
    initialGasMass = TotalWeightedByArea(col) * 
        (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}



void DiskContents::ComputeDerivs(double ** tauvec, std::vector<double>& MdotiPlusHalf, 
                                double ** tauvecMRI, std::vector<double>& MdotiPlusHalfMRI,
                                std::vector<double>& accProf, double AccRate,
                                std::vector<std::vector<int> >& flagList)
{
    //  for(unsigned int i=0; i<=nx-1; ++i) {
    //    MdotiPlusHalf[i] = (-1.0 /  mesh.u1pbPlusHalf(i)) * (tauvec[1][i+1] - tauvec[1][i]) / (x[i+1]-mesh.x(i));
    //  }
    //  MdotiPlusHalf[nx] = -1.0 / (uu[nx] * (1.0+ beta[nx]))  * tauvec[2][nx];

    // Initialize some arrays
    for(unsigned int n=1; n<=nx; ++n) {
        dcoldtIncoming[n] = 0;
        dcoldtOutgoing[n] = 0;
    }
    // Add up outgoing metallicity
    double outflowingMetalMass = 0.0;
    double inflowingMass = 0.0;
    for(unsigned int n=1; n<=nx; ++n) {
        double mu = MassLoadingFactor[n]; //dSdtOutflows(n)/colSFR[n];
        double Zej = ZDisk[n] + xiREC * yREC*RfREC / max(mu, 1.0-RfREC);
        outflowingMetalMass += colSFR[n] *mu * Zej * mesh.area(n);
        inflowingMass += accProf[n]*AccRate * mesh.area(n);
        // double Zacc = ZMix * mu* colSFR[n]/(accProf[n]*AccRate) *Zej + Z_IGM;
    }
    double Zacc = ZMix * outflowingMetalMass/inflowingMass + Z_IGM;

    for(unsigned int n=1; n<=nx; ++n) {
        double dlnZdx; double dlnZdxL; double dlnZdxR;
        // dlnx=dx/x 
        // => dlnZ/dx = dlnZ/(x dlnx) ~ (lnZ_{n+1}-lnZ_{n-1})/(2 x[n] dlnx)
        if(n==1) {
            //      dlnZdxL=(log(ZDisk[1])-log(ZBulge))/(x[n]*dlnx);
            //      dlnZdxR=(log(ZDisk[2])-log(ZBulge))/(2.*x[n]*dlnx);
            //      dlnZdxL=0.0;
            //      dlnZdxR=0.0;
            dlnZdxL = (1.0/x[1]); // (1/Z1) (Z1-0)/x
            dlnZdxR = (log(ZDisk[2])-log(ZDisk[1]))/(x[2]-x[1]);
        }
        else if(n==nx) {
            dlnZdxR=(log(Z_IGM)-log(ZDisk[nx]))/(x[n]*dlnx);
            dlnZdxL=(log(Z_IGM)-log(ZDisk[nx-1]))/(2*x[n]*dlnx);
        }
        else {
            dlnZdxL=(log(ZDisk[n])-log(ZDisk[n-1]))/(x[n]-x[n-1]);
            dlnZdxR=(log(ZDisk[n+1])-log(ZDisk[n]))/(x[n+1]-x[n]);
        }
        dlnZdx = ddx(dlnZdxL,dlnZdxR);
        //    double taupp = ddx(tauvec[2],n,x);

        double Qg= sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);

        dcoldtPrev[n] = dcoldt[n]; 

        if(MdotiPlusHalf[n] > 0) {
            dcoldtIncoming[n] += MdotiPlusHalf[n] / (mesh.dx(n)*x[n]);
            if(n<nx)
                dcoldtOutgoing[n+1] += MdotiPlusHalf[n] / (x[n+1]*mesh.dx(n+1));
        }
        if(MdotiPlusHalfMRI[n] > 0) {
            dcoldtIncoming[n] += MdotiPlusHalfMRI[n] / (mesh.dx(n)*x[n]);
            if(n<nx)
                dcoldtOutgoing[n+1] += MdotiPlusHalfMRI[n] / (x[n+1]*mesh.dx(n+1));
        }
        if(MdotiPlusHalf[n] < 0) {
            if(n<nx)
                dcoldtIncoming[n+1] += -MdotiPlusHalf[n] / (mesh.dx(n+1)*x[n+1]);
            dcoldtOutgoing[n] += -MdotiPlusHalf[n] / (mesh.dx(n)*x[n]);
        }
        if(MdotiPlusHalfMRI[n] < 0) {
            if(n<nx)
                dcoldtIncoming[n+1] += -MdotiPlusHalfMRI[n] / (mesh.dx(n+1)*x[n+1]);
            dcoldtOutgoing[n] += -MdotiPlusHalfMRI[n] / (mesh.dx(n)*x[n]);
        }

        dcoldt[n] = (MdotiPlusHalf[n] + MdotiPlusHalfMRI[n]- MdotiPlusHalf[n-1]-MdotiPlusHalfMRI[n-1]) / (mesh.dx(n) * x[n])
            -RfREC * colSFR[n] - ColOutflows[n] + accProf[n]*AccRate;

        dsigdtPrev[n] = dsigdt[n];
        double ddxSig = ddx(sig,n,x,true,!dbg.opt(9));
        double tauvec2 = tauvec[2][n];
        std::vector<int>& flag = flagList[n];
        if(dbg.opt(17)) {
            ddxSig = ddxUpstream(sig,x,flagList[n],n);
            tauvec2 = (flag[0]*tauvec[1][n-1]+flag[1]*tauvec[1][n]+flag[2]*tauvec[1][n+1])/(flag[0]*mesh.x(n-1)+flag[1]*x[n]+flag[2]*mesh.x(n+1));
        }
    	double dsigdtMRI = (MdotiPlusHalfMRI[n] - MdotiPlusHalfMRI[n-1])*sig[n]/(3.0*x[n]*mesh.dx(n)*col[n]) - 5.0*ddxSig*tauvecMRI[2][n]/(3.0*(beta[n]+1.)*x[n]*col[n]*uu[n]) + uu[n]*(beta[n]-1.)*tauvecMRI[1][n]/(3.0*sig[n]*col[n]*x[n]*x[n]*x[n]);
	    double dsigdtGI = (MdotiPlusHalf[n] - MdotiPlusHalf[n-1])*sig[n]/(3.0*x[n]*mesh.dx(n)*col[n]) - 5.0*ddxSig*tauvec2/(3.0*(beta[n]+1.)*x[n]*col[n]*uu[n]) + uu[n]*(beta[n]-1.)*tauvec[1][n]/(3.0*sig[n]*col[n]*x[n]*x[n]*x[n]);
	
        dsigdtTrans[n] = (MdotiPlusHalf[n] + MdotiPlusHalfMRI[n] - MdotiPlusHalf[n-1]-MdotiPlusHalfMRI[n-1]) * sig[n]/(3.0*x[n]*mesh.dx(n)*col[n]);
        dsigdtDdx[n] = -5.0*ddxSig*(tauvec2+tauvecMRI[2][n]) / (3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n]);
        dsigdtHeat[n] = uu[n]*(beta[n]-1.)*(tauvec[1][n]+tauvecMRI[1][n]) / (3.0*sig[n]*col[n]*x[n]*x[n]*x[n]);

        if(sig[n] >= sigth)
            dsigdtCool[n] = -2.0*M_PI*M_PI*(ETA*pow(1.-sigth*sigth/(sig[n]*sig[n]),1.5))    
                *col[n]*dim.chi()*(1.0+activeColSt(n)/col[n] *sig[n]/activeSigStZ(n))/3.0;
        else 
            dsigdtCool[n] = 0.0;

        dsigdt[n] = (MdotiPlusHalf[n] + MdotiPlusHalfMRI[n] - MdotiPlusHalf[n-1]-MdotiPlusHalfMRI[n-1]) * sig[n] / (3.0*x[n]*mesh.dx(n)*col[n])
            - 5.0*ddxSig*(tauvec2+tauvecMRI[2][n]) / (3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n])
            +uu[n]*(beta[n]-1.)*(tauvec[1][n]+tauvecMRI[1][n]) / (3.0*sig[n]*col[n]*x[n]*x[n]*x[n]);
        if(sig[n] >= sigth) {
            dsigdt[n] -= 2.0*M_PI*M_PI*(ETA*pow(1. - sigth*sigth/(sig[n]*sig[n]),1.5))
                *col[n]*dim.chi()*(1.0+activeColSt(n)/col[n] * sig[n]/activeSigStZ(n))/3.0;
        }
        else {
            // do nothing, these terms are zero.
        }

        if(sig[n]/uu[n] > 100.0) {
            std::cout << "Large sig! n, sig, dsigdt, dsigdtCool " <<n<<" "<<sig[n]<<" "<<dsigdt[n]<<" "<<dsigdtCool[n]<< std::endl;
            errormsg("Very large velocity dispersion. Thin disk approximation definitely makes no sense here.");
        }

        //    colSFR[n] = dSSFdt(n);
        dZDiskdtAdv[n] =  -1.0/((beta[n]+1.0)*x[n]*col[n]*uu[n]) * ZDisk[n]  * dlnZdx *tauvec[2][n] ;

        double Znp1,Znm1;
        if(n<nx) 
            Znp1=ZDisk[n+1];
        else
            Znp1=Z_IGM;
        if(n>1)
            Znm1=ZDisk[n-1];
        else
            Znm1=ZBulge; // THIS SHOULD NOT MATTER, since there should be no mass exiting the bulge.
        double mu = MassLoadingFactor[n]; //dSdtOutflows(n)/colSFR[n];
        double Zej = ZDisk[n] + xiREC * yREC*RfREC / max(mu, 1.0-RfREC);
        //double Zacc = ZMix * mu* colSFR[n]/(accProf[n]*AccRate) *Zej + Z_IGM;
        //double ZAccFac = 0.5;
        dZDiskdtPrev[n] = dZDiskdt[n];
        if (dbg.opt(5)) {
            dMZdt[n] = 
                        accProf[n]*AccRate*Zacc *x[n]*x[n]*sinh(dlnx) +  // new metals from the IGM ;
                        ((yREC-ZDisk[n])*RfREC - mu*Zej)*colSFR[n]*x[n]*x[n]*sinh(dlnx) +
                        PosOnly(MdotiPlusHalf[n]+MdotiPlusHalfMRI[n])*Znp1 
                            - PosOnly(MdotiPlusHalf[n-1]+MdotiPlusHalfMRI[n-1])*ZDisk[n]
                            + NegOnly(MdotiPlusHalf[n]+MdotiPlusHalfMRI[n])*ZDisk[n]
                            - NegOnly(MdotiPlusHalf[n-1]+MdotiPlusHalfMRI[n-1])*Znm1;
            dZDiskdt[n] = -1.0/((beta[n]+1.0)*x[n]*col[n]*uu[n]) * ZDisk[n] 
                * dlnZdx *tauvec[2][n] 
                + yREC*(1.0-RfREC)*colSFR[n]/(col[n])
                + accProf[n]*AccRate * (Zacc - ZDisk[n])/col[n];
        }
        else {
            dMZdt[n] = 
                        accProf[n]*AccRate*Z_IGM *x[n]*x[n]*sinh(dlnx) +  // new metals from the IGM ;
                        ((yREC-ZDisk[n])*RfREC - mu*Zej)*colSFR[n]*x[n]*x[n]*sinh(dlnx) +
                        PosOnly(MdotiPlusHalf[n]+MdotiPlusHalfMRI[n])*Znp1 
                            - PosOnly(MdotiPlusHalf[n-1]+MdotiPlusHalfMRI[n-1])*ZDisk[n]
                            + NegOnly(MdotiPlusHalf[n]+MdotiPlusHalfMRI[n])*ZDisk[n]
                            - NegOnly(MdotiPlusHalf[n-1]+MdotiPlusHalfMRI[n-1])*Znm1;

            dZDiskdt[n] = -1.0/((beta[n]+1.0)*x[n]*col[n]*uu[n]) * ZDisk[n] 
                * dlnZdx *tauvec[2][n] 
                + yREC*(1.0-RfREC)*colSFR[n]/(col[n])
                + accProf[n]*AccRate * (Z_IGM - ZDisk[n])/col[n];

        }




        // Terms for non-instantaneous-recycling
        for(unsigned int i=0; i!=spsPassive.size(); ++i) {
            dZDiskdt[n] += yREC*spsPassive[i]->dcoldtREC[n]/col[n];
            dMZdt[n] += yREC*spsPassive[i]->dcoldtREC[n] * x[n]*x[n]*sinh(dlnx); /// THIS IS PROBABLY WRONG (the other two lines too)
            dcoldt[n] += spsPassive[i]->dcoldtREC[n];
        }




        // Check to make sure the results of this method are reasonable..
        if(dcoldt[n]!=dcoldt[n] || 
                dsigdt[n]!=dsigdt[n] ||  
                dZDiskdt[n]!=dZDiskdt[n]) {
            std::string space(" ");
            errormsg(std::string("Error computing derivatives - n,dcoldt,dsigdt,dZDiskdt,")
                    +std::string("tau[1],tau[2],  col,sig,taupp: ")+str(n)+" "+str(dcoldt[n])+" "
                    +str(dsigdt[n])+space+str(dZDiskdt[n])+space+str(tauvec[1][n])
                    +space+str(tauvec[2][n])+space+str(col[n])+space+str(sig[n])+" "
                    );
        }
    }
}

double DiskContents::ComputeTimeStep(const double redshift,int * whichVar, int * whichCell, double ** tauvecStar,std::vector<double>& MdotiPlusHalfStar)
{
    // Compute a bunch of timescales for variation in each cell, i.e. Quantity / (dQuantity/dt)
    // Find the maximum value of the inverse of all such timescales. 
    double dmax=0.;

    for(unsigned int n=1; n<=nx; ++n) {
        if(fabs(dZDiskdt[n]/ZDisk[n]) > dmax) {
            dmax = fabs(dZDiskdt[n]/ZDisk[n]);
            *whichVar=1;
            *whichCell=n;
        }
        if(fabs(dcoldt[n]/col[n]) > dmax) {
            dmax = fabs(dcoldt[n]/col[n]);
            *whichVar =2;
            *whichCell=n;
        }
        if(sig[n] > sigth) {
//            if(fabs(dsigdt[n]/sqrt(sig[n]*sig[n]-sigth*sigth)) > dmax) {
//                dmax = fabs(dsigdt[n]/sqrt(sig[n]*sig[n]-sigth*sigth));
            if(fabs(dsigdt[n]/sig[n]) > dmax) {
                dmax = fabs(dsigdt[n]/sig[n]);
                *whichVar = 3;
                *whichCell=n;
            }
        }

        unsigned int sa=spsActive.size();
        if(fabs(colSFR[n]/spsActive[sa-1]->spcol[n])>dmax) { 
            dmax=fabs(colSFR[n]/spsActive[sa-1]->spcol[n]);
            *whichCell=n;
            *whichVar=5;
        }

        for(unsigned int i=0; i!=sa; ++i) {
            double dcolst = dSMigdt(n,tauvecStar,(*this),spsActive[i]->spcol) + RfREC*colSFR[n];
            if(fabs( dcolst /spsActive[i]->spcol[n]) > dmax) {
                dmax=fabs( dcolst/spsActive[i]->spcol[n] );
                *whichVar=6;
                *whichCell=n;
            }
            if(fabs(dSigStZdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar)
                        /spsActive[i]->spsigZ[n]) > dmax) {
                dmax=fabs(dSigStZdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar)
                        /spsActive[i]->spsigZ[n]);
                *whichVar=7;
                *whichCell=n;
            }
            if(fabs(dSigStRdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar)
                        /spsActive[i]->spsigR[n]) > dmax) {
                dmax=fabs(dSigStRdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar)
                        /spsActive[i]->spsigR[n]);
                *whichVar=8;
                *whichCell=n;
            }

        } // end loop over active stellar populations




        bool restrictOnPassive = false;
        if(restrictOnPassive) {
            unsigned int sp = spsPassive.size();
            if(fabs(colSFR[n]/spsPassive[sp-1]->spcol[n])>dmax) { 
                dmax=fabs(colSFR[n]/spsPassive[sp-1]->spcol[n]);
                *whichCell=n;
                *whichVar=9;
            }

            for(unsigned int j=0; j!=sp; ++j) {

                if(fabs(dSMigdt(n,tauvecStar,(*this),spsPassive[j]->spcol)
                            /spsPassive[j]->spcol[n]) > dmax) {
                    dmax=fabs(dSMigdt(n,tauvecStar,(*this),spsPassive[j]->spcol)
                            /spsPassive[j]->spcol[n]);
                    *whichVar=10;
                    *whichCell=n;
                }
                if(fabs(dSigStZdt(n,j,spsPassive,tauvecStar,MdotiPlusHalfStar)
                            /spsPassive[j]->spsigZ[n]) > dmax) {
                    dmax=fabs(dSigStZdt(n,j,spsPassive,tauvecStar,MdotiPlusHalfStar)
                            /spsPassive[j]->spsigZ[n]);
                    *whichVar=11;
                    *whichCell=n;
                }
                if(fabs(dSigStRdt(n,j,spsPassive,tauvecStar,MdotiPlusHalfStar)
                            /spsPassive[j]->spsigR[n]) > dmax) {
                    dmax=fabs(dSigStRdt(n,j,spsPassive,tauvecStar,MdotiPlusHalfStar)
                            /spsPassive[j]->spsigR[n]);
                    *whichVar=12;
                    *whichCell=n;
                }
            } // end loop over passive stellar populations
        }


        
        if(dmax!=dmax) errormsg("Error setting timestep. n, whichVar, whichCell: "+str(n)+" "+str(*whichVar)+" "+str(*whichCell));

    } // end loop over cells
    double dt = TOL/max(dmax,10.0*TOL/x[1]); // maximum stepsize of .1 innermost orbital times


        // Code from DiffuseMetals. Despite its implicit nature, 
        // the diffusion algorithm has a limitation, namely the 
        // tridiagonal matrix can't have 30 orders of magnitude
        // in dynamic range! This puts a limit on the timestep.
        double KM;
        unsigned int n=nx;
        // Scale from Yang and Krumholz 2012 (proportional to lambda_J^2/t_orb).
        KM=(kappaMetals)*4.7e-3*sig[n]*sig[n]*sig[n]*sig[n]*sqrt(1.0+beta[n])*uu[n]/(col[n]*col[n]*dim.chi()*dim.chi()*x[n]); //KM = (kappaMetals*1.0e3)*sig[n]*sig[n]*sig[n]/(col[n]*dim.chi());
        // Scale from Yang and Krumholz 2012 (t_orb only)
        //else if(dbg.opt(8)) KM = (kappaMetals*1.0e3)*1.0e-4 * uu[n]/uu[nx]*x[nx]/x[n] * dim.Radius/(10.0*cmperkpc); //KM = (kappaMetals*1.0e3)*sig[n]*sig[n]*sig[n]/(M_PI*col[n]*dim.chi() * (1.0 +  sig[n]*activeColSt(n)/(activeSigStZ(n)*col[n]))); 
        // don't scale - kappa is constant.
        //else KM = kappaMetals;
        double sum = 4.0*M_PI*KM/(mesh.dx(n)*mesh.dx(n));
        double colnp1 = col[nx];
        if(n!=nx) colnp1=col[n+1];
        double ratio = mesh.x(n+1)*mesh.x(n+1)*colnp1/(x[n]*x[n]*col[n]);
        double eta = sum/(1.0+ratio);
        double dtZ =   (1.0+ratio)*1.0e30 / (   (kappaMetals)*4.7e-3*sig[n]*sig[n]*sig[n]*sig[n]*sqrt(1.0+beta[n])*uu[n]/(col[n]*col[n]*dim.chi()*dim.chi()*x[n]) * 4.0*M_PI/(mesh.dx(n)*mesh.dx(n))  );

//    if(dt>dtZ) {
//        dt=dtZ;
//        *whichVar=15;
//    }



    return dt;
}

void DiskContents::AddNewStellarPop(const double redshift, 
        const double dt,
        std::vector<StellarPop*>& sps, 
        bool active)
{
    unsigned int sz=sps.size();
    sps[sz-1]->endingAge = cos.lbt(redshift);
    StellarPop * currentlyForming = new StellarPop(mesh);
    currentlyForming->ageAtz0 = cos.lbt(redshift);
    currentlyForming->startingAge = cos.lbt(redshift);
    if(active) { // this is a huge kludge. Used to speed up 
        // runtimes in the case that NActive>1
        currentlyForming->extract(*(sps[sz-1]),.01);
    }
    else {
        FormNewStars(*currentlyForming, dt, redshift);
    }
    currentlyForming->ComputeSpatialDerivs();
    sps.push_back(currentlyForming);
    if(active) {
        std::cout << "Creating spsActive["<<sz<<"]"<<std::endl;
    }
    else {
        std::cout << "Creating spsPassive["<<sz<<"]"<<std::endl;
    }
}

void DiskContents::FormNewStars(StellarPop & currentlyForming, double dt, double redshift)
{
    for(unsigned int n=1; n<=nx; ++n) {
        currentlyForming.spcol[n] = RfREC * colSFR[n] * dt;
        currentlyForming.spZ[n] = ZDisk[n];
        currentlyForming.spZV[n] = 0.0;
        if(!dbg.opt(16)) {
            if(sigth*sigth+minsigst*minsigst<=sig[n]*sig[n])
                currentlyForming.spsigR[n] = sqrt(sig[n]*sig[n]-sigth*sigth); // the turbulent velocity dispersion of the gas
            else
                currentlyForming.spsigR[n] = minsigst;
        }
        else {
            currentlyForming.spsigR[n] = sig[n];
        }
        currentlyForming.spsigZ[n] = currentlyForming.spsigR[n];

        if(currentlyForming.spcol[n] < 0. || currentlyForming.spsigR[n]<0.0 || currentlyForming.spcol[n]!=currentlyForming.spcol[n] || currentlyForming.spsigR[n]!=currentlyForming.spsigR[n])
            errormsg("FormNewStars: newly formed stars are problematic: n, spcol, spsig, colSFR, dt, sigth:  "+str(n)+", "+str(currentlyForming.spcol[n])+", "+str(currentlyForming.spsigR[n])+", "+str(colSFR[n]) +", "+str(dt)+";  sig, sigth: "+str(sig[n])+", "+str(sigth));
    }
    currentlyForming.ageAtz0 = cos.lbt(redshift);

}

void DiskContents::UpdateStateVars(const double dt, const double dtPrev,
				   const double redshift,
				   double ** tauvec, double AccRate, 
				   double ** tauvecStar, 
				   std::vector<double>& MdotiPlusHalf,
				   std::vector<double>& MdotiPlusHalfStar,
				   std::vector<double>& MdotiPlusHalfMRI,
                   double fracAccInner, double accSt)
{
    // some things for debugging purposes
    std::vector<double> dColStDtMeasured(nx+1,0);
    std::vector<double> dSigStRDtMeasured(nx+1,0);
    std::vector<double> dSigStZDtMeasured(nx+1,0);
    std::vector<double> dColStDtNominal(nx+1,0);
    std::vector<double> dSigStRDtNominal(nx+1,0);
    std::vector<double> dSigStZDtNominal(nx+1,0);
    for(unsigned int n=1; n<=nx; ++n) {
        dColStDtMeasured[n] = -activeColSt(n)/dt; // we'll subtract from this number in a bit..
        dSigStRDtMeasured[n] = -activeSigStR(n)/dt;
        dSigStZDtMeasured[n] = -activeSigStZ(n)/dt;
        dColStDtNominal[n] = dSMigdt(n,tauvecStar, (*this), spsActive[0]->spcol) + RfREC*colSFR[n];
        dSigStRDtNominal[n] = dSigStRdt(n,0,spsActive,tauvecStar,MdotiPlusHalfStar);
        dSigStZDtNominal[n] = dSigStZdt(n,0,spsActive,tauvecStar,MdotiPlusHalfStar);

    }

    unsigned int szA = spsActive.size();
    unsigned int szP = spsPassive.size();
    StellarPop currentlyForming(mesh);
    FormNewStars(currentlyForming,dt,redshift);

    spsActive[szA-1]->MergeStellarPops(currentlyForming,(*this));
    spsPassive[szP-1]->MergeStellarPops(currentlyForming,(*this));

    for(unsigned int i=0; i!=spsPassive.size();++i) {
        if(migratePassive) {
            spsPassive[i]->MigrateStellarPop(dt,tauvecStar,(*this),MdotiPlusHalfStar);
            spsPassive[i]->ComputeSpatialDerivs();
        }
    }
    for(unsigned int i=0; i!=spsActive.size();++i) {
        spsActive[i]->MigrateStellarPop(dt,tauvecStar,(*this),MdotiPlusHalfStar);
        spsActive[i]->ComputeSpatialDerivs();
    }

    if(MdotiPlusHalf[0] < -1.0e-10 || MdotiPlusHalfMRI[0] < -1.0e-10 || MdotiPlusHalfStar[0] < -1.0e-10 || fracAccInner*AccRate < -1.0e-10)
        errormsg("Nonphysical mass flux at inner boundary.");

    double reduce = RfREC/(RfREC+MassLoadingFactor[1]);
    double MGasIn = dt*PosOnly(MdotiPlusHalf[0]+MdotiPlusHalfMRI[0]);
    double MStarsIn = dt*PosOnly(MdotiPlusHalfStar[0]);
    double MGasAcc = dt*fracAccInner*AccRate;
    double MstAcc = dt*accSt * (MSol/speryear) /(dim.MdotExt0 ); // accSt is in unnormalized dimensional units(!)
    MHalo += MstAcc;
    double MIn = MGasIn*reduce + MStarsIn + MGasAcc*reduce;
    //  double MIn = cumulativeMassAccreted -(MassLoadingFactor+RfREC)* cumulativeStarFormationMass - MBulge - (TotalWeightedByArea(col) - initialGasMass) - (TotalWeightedByArea());
    ZBulge = (ZBulge*MBulge + (yREC+ ZDisk[1])*MGasIn*reduce + (yREC+ Z_IGM)*MGasAcc*reduce + spsActive[0]->spZ[1]*MStarsIn)/(MBulge+MIn);
    double dummy = spsActive[0]->spZ[1]*MStarsIn;
    double dummy2 = spsActive[0]->spZ[1];
    if(ZBulge <= 0.0 || ZBulge >1.0) {
      errormsg(std::string("Nonphysical ZBulge- ZBulge,MBulge,dt,Mdot0,ZD1,dmdtCosInner,MIn:   ")+str(ZBulge)+" "+str(MBulge)+" "+str(dt)+" "+str(MdotiPlusHalf[0])+" "+str(ZDisk[1])+" "+str(fracAccInner*(AccRate)));
    }
    MBulge += MIn;
    CumulativeTorque+= tauvec[1][nx]*dt;
    for(unsigned int n=1; n<=nx; ++n) {
        if(n!=1) {
            //      for(unsigned int j=0; j!=spsActive.size(); ++j) {
            //        CuStarsOut[n] += 2.0*M_PI*
            //	  sqrt(x[n]*x[n-1] * 
            //	       spsActive[j]->spcol[n]*spsActive[j]->spcol[n-1] * 
            //	       yy[n]*yy[n-1]) * 
            //	  dt * 2.0*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR  * (1.0/MSol);
            //        if(CuStarsOut[n] != CuStarsOut[n])
            ///         errormsg("Error computing CuStarsOut. "+str(spsActive[j]->spcol[n])+" "+str(spsActive[j]->spcol[n-1])+" "+str(yy[n])+" "+str(yy[n-1]));
            //      }
            CuGasOut[n] +=  
                sqrt(max(tauvec[2][n]*tauvec[2][n-1],1.0e-20))
                / (sqrt(uu[n]*uu[n-1] * (1. + beta[n])*(1.+beta[n-1]))) 
                * dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
        }
        CuGasOut[1]+=(sqrt(max(tauvec[1][1]*tauvec[1][2],1.0e-20))-0.)
            /((XMIN*exp(dlnx/2.)*expm1(dlnx))*uu[1]*(1+beta[1])) * 
            dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
        //    for(unsigned int j=0; j!=spsActive.size(); ++j) {
        //      CuStarsOut[1] += 2*M_PI*x[1]*spsActive[j]->spcol[1]*yy[1]*
        //	dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
        //    }

        double MZ = ZDisk[n]*x[n]*x[n]*sinh(dlnx)*col[n];

        if(!dbg.opt(11))
            col[n] += dcoldt[n] * dt; //Euler step
        else
            col[n] += dt * ( (dcoldt[n] - dcoldtPrev[n])*.5*dt/dtPrev  + dcoldt[n]  );

        // The gas velocity dispersion should always be above sigth.
        // If it is not, e.g. a semi-pathological initial condition or
        // some numerical issue, set it to sigth.
        // 
        if(sig[n] < sigth) {
            sig[n]=sigth;
            keepTorqueOff[n] = 1;
        }
        else {
            if(dbg.opt(11))
                sig[n] += dt* ((dsigdt[n]-dsigdtPrev[n])*.5*dt/dtPrev + dsigdt[n]);
            else
                sig[n] += dt*dsigdt[n];
        }  
//        if(dbg.opt(11))
//            ZDisk[n] += dt* ((dZDiskdt[n]-dZDiskdtPrev[n])*.5*dt/dtPrev + dZDiskdt[n]);
//        else
//            ZDisk[n] += dZDiskdt[n] * dt; // -- Euler step
//
        MZ += dMZdt[n]*dt;
        ZDisk[n] = MZ/ (x[n]*x[n]*sinh(dlnx)*col[n]);
        CumulativeSF[n] += colSFR[n] * dt;

        double totalLost=0.0;
        for(unsigned int i=0; i!=spsPassive.size(); ++i) {
            double lostFromThisPop = spsPassive[i]->dcoldtREC[n]*dt;
            totalLost += lostFromThisPop;
            spsPassive[i]->spcol[n] -= lostFromThisPop;
        }
        spsActive[0]->spcol[n] -= totalLost;


        // Check that this method hasn't made the state variables non-physical
        if(col[n]<0. || sig[n]<0. || ZDisk[n]<0. || 
                col[n]!=col[n] || sig[n]!=sig[n] || ZDisk[n]!=ZDisk[n]) {
            std::string spc(" ");
            errormsg(std::string("Error updating statevars- dt,col,sig,ZDisk,dcoldt,dsigdt")
                    +std::string(",dZDiskdt: ")+str(dt)+spc+str(col[n])+spc+str(sig[n])+spc
                    +str(ZDisk[n])+spc+spc+str(dcoldt[n])+spc+str(dsigdt[n])+spc
                    +str(dZDiskdt[n])+spc+spc+str(spsActive[szA-1]->spcol[n])+spc
                    +str(spsActive[szA-1]->spsigR[n])+spc+str(spsActive[szA-1]->spZ[n]));
        }
    }

    // Artificially diffuse metals in the gas phase 
    // to maintain numerical stability
    //  DiffuseMetallicity(dt,.005);
    DiffuseMetals(dt);  

    // Record a few numbers to check things like mass conservation...
    cumulativeStarFormationMass += TotalWeightedByArea(currentlyForming.spcol)  
        *(2.0*M_PI*dim.Radius * dim.MdotExt0 / dim.vphiR) * (1.0/MSol);
    cumulativeGasMassThroughIB += dim.MdotExt0 
        * MdotiPlusHalf[0] * dt * (2*M_PI*dim.Radius/dim.vphiR) * (1.0/MSol);
    //  for(unsigned int j=0; j!=spsActive.size(); ++j) {
    //    cumulativeStellarMassThroughIB += 2.*M_PI*(x[1]*dim.Radius)
    //     *(spsActive[j]->spcol[1]*dim.MdotExt0/(dim.vphiR*dim.Radius))
    //     *(yy[1]*dim.vphiR) * dt * (2*M_PI*dim.Radius/dim.vphiR) * (1.0/MSol);
    //  }
    cumulativeMassAccreted += AccRate*dim.MdotExt0 * 
        dt * (2*M_PI*dim.Radius/dim.vphiR)* (1.0/MSol);

    for(unsigned int n=1; n<=nx; ++n) {
        dColStDtMeasured[n] += activeColSt(n)/dt; // we'll subtract from this number in a bit..
        dSigStRDtMeasured[n] += activeSigStR(n)/dt;
        dSigStZDtMeasured[n] += activeSigStZ(n)/dt;
    }

    for (unsigned int n=1; n<=10; ++n) {
        //    std::cout <<"col, sigR, sigZ: " << dColStDtMeasured[n]/dColStDtNominal[n]-1. << " " << dSigStRDtMeasured[n]/dSigStRDtNominal[n]-1.<<" "<<dSigStZDtMeasured[n]/dSigStZDtNominal[n]-1.<<std::endl;
    }
}

double DiskContents::TotalWeightedByArea(const std::vector<double>& perArea)
{
    double sum=0.0;
    for(unsigned int i=1; i!=x.size(); ++i) {
        sum += perArea[i]*x[i]*dlnx*x[i];
    } 
    return sum;
}

void DiskContents::ComputeRafikovQParams(RafikovQParams* p, unsigned int n)
{
    (*p).var=-1;
    (*p).analyticQ = analyticQ;
    (*p).thickGas = 1.5; // thickness;
    (*p).thickStars = 0.8 + 0.7 * spsActive[0]->spsigZ[n] / spsActive[0]->spsigR[n];
    //  (*p).mostRecentq = 1.;
    (*p).Qg = sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);
    (*p).Qsi.clear();
    (*p).ri.clear();
    for(unsigned int i=0; i!=spsActive.size(); ++i) {
        (*p).Qsi.push_back( sqrt(2.*(beta[n]+1.))*uu[n]*spsActive[i]->spsigR[n]
                /(M_PI*dim.chi()*x[n]*spsActive[i]->spcol[n]));
        (*p).ri.push_back( spsActive[i]->spsigR[n] / sig[n]);
    }
    (*p).fixedQ=fixedQ;
}

// if fixedPhi0 is true, Q=Q_f is enforced by adjusting sig and sig_st simultaneously
// otherwise we assume sig_st is fixed, and we adjust sig only.
void DiskContents::EnforceFixedQ(bool fixedPhi0, bool EnforceWhenQgrtQf)
{
    RafikovQParams rqp;
    gsl_function F;
    // the function whose root we want, f=Q(stateVars)-fixedQ.
    // QmfQ is declared in DiskUtils.h
    if(fixedPhi0)
        F.function = &QmfQ; 
    else
        F.function = &QmfQfst;
    double factor =1.;
    rqp.mostRecentq=1.;
    double absc=1.0;
    for(unsigned int n=1; n<=nx; ++n) {
        ComputeRafikovQParams(&rqp,n);

        // Only enforce Q=Qf when either
        // 1) EnforceWhenQgrtQf is true, in which case always enforce Q=Qf
        // or, 2) if we shouldn't enforce it always, only enforce it when Q<Qf
        // Basically, always set Q=Qf if Q<Qf, but not always if Q>Qf .
        if(EnforceWhenQgrtQf || (!EnforceWhenQgrtQf && Q(&rqp,&absc) < fixedQ)) {

            F.params = &rqp;
            findRoot(F,&factor);

            sig[n] *= factor;

            if(fixedPhi0) {
                for(unsigned int i=0; i!=spsActive.size(); ++i) {
                    spsActive[i]->spsigR[n] *= factor;
                    spsActive[i]->spsigZ[n] *= factor;
                }
                for(unsigned int i=0; i!=spsPassive.size(); ++i) {
                    spsPassive[i]->spsigR[n] *=factor;
                    spsPassive[i]->spsigZ[n] *=factor;
                }
            }
        }
    }
}

void DiskContents::ComputeMRItorque(double ** tauvecMRI, const double alpha) 
{
    for(unsigned int n=1; n<=nx; ++n) {
        tauvecMRI[1][n] = -2.0*M_PI*x[n]*x[n]*col[n]*sigth*sigth*alpha;
    }
    tauvecMRI[1][0]=-2.0*M_PI*mesh.x((unsigned int) 0)*mesh.x((unsigned int) 0)*col[1]*sigth*sigth*alpha;
    tauvecMRI[1][nx+1] = tauvecMRI[1][nx]; // Mdot = 0
    // tauvecMRI[1][nx+1] = tauvecMRI[1][nx] + mesh.u1pbPlusHalf(nx)/mesh.u1pbPlusHalf(nx-1) * (mesh.x(nx+1)-x[nx])/(x[nx]-x[nx-1]) * (tauvecMRI[1][nx]-tauvecMRI[1][nx-1]); // dcoldt=0
}

// Should apply to both gas and stars
void DiskContents::ComputeGItorque(double ** tauvec, const double IBC, const double OBC,
            std::vector<double>& UU, std::vector<double>& DD, 
            std::vector<double>& LL, std::vector<double>& FF, 
            std::vector<double>& MdotiPlusHalf)
{  
    // Read in the given vectors to gsl_vectors suitable for tridiagonal-solving.
    for(unsigned int n=1; n<=nx; ++n) {
        gsl_vector_set(lr,   n-1, LL[n]);
        gsl_vector_set(diag,   n, DD[n]);
        gsl_vector_set(ur,     n, UU[n]);
        gsl_vector_set(forcing,n, FF[n]);
    }

    // Set the appropriate boundary conditions
    gsl_vector_set(ur  , 0, 0.0);
    gsl_vector_set(diag, 0, 1.0);

    gsl_vector_set(lr, nx,  0.0);//-1.0/(mesh.u1pbPlusHalf(nx)*(mesh.x(nx+1)-x[nx])));
    gsl_vector_set(diag,nx+1, 1.0); //1.0/(mesh.u1pbPlusHalf(nx)*(mesh.x(nx+1)-x[nx])));
    gsl_vector_set(forcing,0,IBC);
    gsl_vector_set(forcing,nx+1,OBC);

    // solve the tridiagonal matrix equation.
    int status = gsl_linalg_solve_tridiag(diag,ur,lr,forcing,tau);
    if(status!=GSL_SUCCESS) 
        errormsg("Failed to solve the torque equation!");

    // Now read the solution back out.
    for(unsigned int n=0; n<=nx+1; ++n) {
        tauvec[1][n] = gsl_vector_get(tau,n);
    }

    // Zero out positive (non-physical) torques.
    for(unsigned int n=0; n<=nx+1; ++n) {
        if(tauvec[1][n]>0.0) {
            tauvec[1][n]=0.0;
        }
    }

    if(dbg.opt(18)) {
        std::vector<double> avgtorque(nx+1,0);
        int halfWindow=1;
        for(int n=1; n<=nx; ++n) {
            double accumulated=0.0;
            int counter = 0;
            for(int nn=max(1,n-halfWindow);nn<=min(((int) nx),n+halfWindow);++nn) {
                accumulated += tauvec[1][nn];
                ++counter;
            }
            accumulated += 10.0*tauvec[1][n];
            counter += 10;
            avgtorque[n] = accumulated/((double) counter);


        }
        for(unsigned int n=1; n<=nx; ++n) {
            tauvec[1][n] = avgtorque[n];
        }
//        std::vector<int> overshoot(nx+1,0);
//        for(unsigned int n=1; n<=nx; ++n) {
//            if(tauvec[1][n]>=0.0 && (tauvec[1][n-1]<0.0 || tauvec[1][n+1]<0.0))
//                overshoot[n]=1;
//        }
//        for(unsigned int n=1; n<=nx; ++n) {
//            if(overshoot[n]==1)
//                tauvec[1][n]=(tauvec[1][n-1]+tauvec[1][n+1])/3.0;
//        }
    }

    // Compute first derivatives and fluxes.
    ComputeFluxes(tauvec,MdotiPlusHalf,mesh);

    // All set!.
}

void ComputeFluxes(double ** tauvec, std::vector<double>& MdotiPlusHalf, FixedMesh& mesh )
{
    for(unsigned int n=1; n<=mesh.nx(); ++n) {
        tauvec[2][n] = (tauvec[1][n+1]-tauvec[1][n-1])/(mesh.x(n+1)-mesh.x(n-1));
        MdotiPlusHalf[n] = -1.0/mesh.u1pbPlusHalf(n) * (tauvec[1][n+1]-tauvec[1][n])/(mesh.x(n+1)-mesh.x(n));
    }
    MdotiPlusHalf[0]=-1.0/mesh.u1pbPlusHalf(0) * (tauvec[1][1]-tauvec[1][0])/(mesh.x((unsigned int) 1)-mesh.x((unsigned int) 0));

    if(MdotiPlusHalf[0] < 0.0)
      errormsg("Negative MdotiPlusHalf:  "+str(mesh.u1pbPlusHalf(0))+" "+str(tauvec[1][1])+" "+str(tauvec[1][0])+" "+str(mesh.x((unsigned int) 1))+" "+str(mesh.x((unsigned int) 0)) );
}

void DiskContents::DiffuseMetals(double dt)
{
    gsl_vector *lr, *diag, *ur;
    gsl_vector *MetalMass1, *MetalMass2;

    // what boundary condition do we use?

    lr=gsl_vector_alloc(nx);
    diag=gsl_vector_alloc(nx+1);
    ur=gsl_vector_alloc(nx);
    MetalMass1=gsl_vector_alloc(nx+1); MetalMass2=gsl_vector_alloc(nx+1);

    std::vector<double> etas(0), xis(0);
    etas.push_back(0.0); xis.push_back(0.0);
    // ZFlux[n] = net flux of metal mass from bin i+1 to bin i.
    for(unsigned int n=1; n<=nx; ++n) {
        double KM;
        // Scale from Yang and Krumholz 2012 (proportional to lambda_J^2/t_orb).
        KM=(kappaMetals)*4.7e-3*sig[n]*sig[n]*sig[n]*sig[n]*sqrt(1.0+beta[n])*uu[n]/(col[n]*col[n]*dim.chi()*dim.chi()*x[n]); //KM = (kappaMetals*1.0e3)*sig[n]*sig[n]*sig[n]/(col[n]*dim.chi());
        // Scale from Yang and Krumholz 2012 (t_orb only)
        //else if(dbg.opt(8)) KM = (kappaMetals*1.0e3)*1.0e-4 * uu[n]/uu[nx]*x[nx]/x[n] * dim.Radius/(10.0*cmperkpc); //KM = (kappaMetals*1.0e3)*sig[n]*sig[n]*sig[n]/(M_PI*col[n]*dim.chi() * (1.0 +  sig[n]*activeColSt(n)/(activeSigStZ(n)*col[n]))); 
        // don't scale - kappa is constant.
        //else KM = kappaMetals;
        double klim = sig[n]*x[n]; // uu[n]*1.0; // sig[n] * x[n]; //sig[n]*sig[n]/uu[n] * x[n];
        if(KM > klim ) KM = klim;
        double sum = 4.0*M_PI*KM/(mesh.dx(n)*mesh.dx(n));
        double colnp1 = col[nx];
        if(n!=nx) colnp1=col[n+1];
        double ratio = mesh.x(n+1)*mesh.x(n+1)*colnp1/(x[n]*x[n]*col[n]);
        etas.push_back(sum/(1.0+ratio));
        xis.push_back(sum*ratio/(1.0+ratio));
    }
    xis.push_back(0.0);

    for(unsigned int i=0; i<=nx-1; ++i) {
        gsl_vector_set(lr, i, -xis[i+1]*dt);
        gsl_vector_set(diag, i, 1.0+dt*(xis[i+1]+etas[i]));
        gsl_vector_set(ur,i,-etas[i+1]*dt);

    }
    gsl_vector_set(diag, nx, 1.0+dt*etas[nx]);

    std::vector<double> ZFlux(nx+1);
    for(unsigned int n=1; n<=nx; ++n) {
        gsl_vector_set(MetalMass1,n-1, ZDisk[n]*col[n]*x[n]*x[n]*dlnx);
    }

    // Now the default is to set Z_boundary = Z_IGM - dbg.opt(1) is being used for something else
    //    if(dbg.opt(1)) { // Z_boundary = Z_IGM, i.e. metals are free to flow off the edge of the disk.
        gsl_vector_set(MetalMass1,nx,Z_IGM*col[nx]*mesh.x(nx+1)*mesh.x(nx+1)*dlnx);
	//    }
	//    else { // Zero flux condition
	//        gsl_vector_set(MetalMass1,nx,ZDisk[nx]*col[nx]*mesh.x(nx+1)*mesh.x(nx+1)*dlnx);
	//    }

    for(unsigned int n=0; n<=nx; ++n) {
        double left=0;
        double right=0;
        if(n!=0) left=-xis[n]*x[n]*x[n]*col[n]*dlnx;
        if(n!=nx) right = etas[n]*x[n+1]*x[n+1]*col[n+1]*dlnx;
        else right = etas[n]*mesh.x(n+1)*mesh.x(n+1)*col[nx]*dlnx;
        ZFlux[n] = left+right;

    }

    gsl_linalg_solve_tridiag(diag,ur,lr,MetalMass1,MetalMass2);

    for(unsigned int n=1; n<=nx; ++n) {
        dZDiskdtDiff[n] = -ZDisk[n]/dt;
        ZDisk[n] = gsl_vector_get(MetalMass2,n-1)/ (col[n]*x[n]*x[n]*dlnx);
        dZDiskdtDiff[n] += ZDisk[n]/dt;
        if(ZDisk[n]!=ZDisk[n] || ZDisk[n]<0.0 || ZDisk[n]>0.5)
            errormsg("Error diffusing the metals. Printing n, ZDisk[n], col[n]:  "+str(n)+" "+str(ZDisk[n])+" "+str(col[n]));
    }
    gsl_vector_free(lr); gsl_vector_free(diag); gsl_vector_free(ur);
    gsl_vector_free(MetalMass1); gsl_vector_free(MetalMass2);
}


void DiskContents::DiffuseMetalsUnstable(double dt, double km)
{
    gsl_vector *lr, *diag, *ur;
    gsl_vector *MetalMass1, *MetalMass2;

    lr=gsl_vector_alloc(nx-1);
    diag=gsl_vector_alloc(nx);
    ur=gsl_vector_alloc(nx-1);
    MetalMass1=gsl_vector_alloc(nx); MetalMass2=gsl_vector_alloc(nx);

    unsigned int n;
    if(ZBulge<Z_BBN || ZBulge!=ZBulge) {
        std::cerr<<"Warning: ZBulge hit the metallicity floor"<<std::endl;
        ZBulge=Z_BBN;
    }
    //  gsl_vector_set(MetalMass1,0, ZDisk[n]*MBulge);
    for(n=1; n<=nx; ++n) {
        if(ZDisk[n]<Z_BBN || ZDisk[n]!=ZDisk[n]) {
            std::cerr<<"Warning:ZDisk["<<n<<"] hit the metallicity floor"<<std::endl;
            ZDisk[n]=Z_IGM;
        }
        // note that ZDisk is indexed from 1 but MetalMass1 is indexed from 0.
        gsl_vector_set(MetalMass1,n-1,  ZDisk[n]*col[n]*x[n]*x[n]);

    }
    for(n=1; n<nx-1; ++n) {
        gsl_vector_set(lr, n-1, -dt*km/(x[n]*x[n]*dlnx*dlnx));
        gsl_vector_set(diag,n, 1+dt*km*2/(x[n]*x[n]*dlnx*dlnx));
        gsl_vector_set(ur,n,   -dt*km/(x[n]*x[n]*dlnx*dlnx));
    }
    gsl_vector_set(diag,0, 1.+dt*km/(x[1]*x[1]*dlnx*dlnx));
    gsl_vector_set(ur,0,  -dt*km/(x[1]*x[1]*dlnx*dlnx));
    gsl_vector_set(lr,nx-2, -dt*km/(x[nx]*x[nx]*dlnx*dlnx));
    gsl_vector_set(diag,nx-1, 1.+dt*km/(x[nx]*x[nx]*dlnx*dlnx));

    gsl_linalg_solve_tridiag(diag,ur,lr,MetalMass1,MetalMass2);

    for(n=1; n<=nx; ++n) {
        ZDisk[n]=gsl_vector_get(MetalMass2,n-1)/ (col[n]*x[n]*x[n]);
        if(ZDisk[n]!=ZDisk[n])
            errormsg("Nonphysical metallicity: n,ZDisk2: "+str(n)+" "+str(ZDisk[n]));
    }

    gsl_vector_free(lr); gsl_vector_free(diag); gsl_vector_free(ur);
    gsl_vector_free(MetalMass1), gsl_vector_free(MetalMass2);
}

// This is not correct because logZ is not what is conserved
void DiskContents::DiffuseMetallicity(double dt,double km)
{
    gsl_vector *lr, * diag, *ur;
    gsl_vector *diffLogZ1, *diffLogZ2;

    lr=gsl_vector_alloc(nx); 
    diag=gsl_vector_alloc(nx+1); 
    ur=gsl_vector_alloc(nx);
    diffLogZ1=gsl_vector_alloc(nx+1); diffLogZ2=gsl_vector_alloc(nx+1);

    unsigned int n;
    if(ZBulge<Z_BBN || ZBulge!=ZBulge) {
        std::cerr<<"Warning: ZBulge hit the metallicity floor"<<std::endl;
        ZBulge=Z_BBN;
    }
    gsl_vector_set(diffLogZ1,0,log10(ZBulge));

    for(n=1; n<=nx; ++n) {
        if(ZDisk[n]<Z_BBN || ZDisk[n]!=ZDisk[n]) {
            std::cerr << "Warning: ZDisk["<<n
                <<"] hit the metallicity floor"<<std::endl;
            ZDisk[n]=Z_IGM;
        }

        // note that ZDisk is indexed from1 but diffLogZ1 is indexed from 0.
        gsl_vector_set(diffLogZ1,n,log10(ZDisk[n])); 
    }

    for(n=1; n<nx; ++n) {
        gsl_vector_set(lr,n-1,  -dt*km/(x[n]*x[n]*dlnx*dlnx));
        gsl_vector_set(diag,n,  1+dt*km*2/(x[n]*x[n]*dlnx*dlnx));
        gsl_vector_set(ur,n,  -dt*km/(x[n]*x[n]*dlnx*dlnx));
    }
    gsl_vector_set(diag,0,  1.+dt*km/(x[1]*x[1]*dlnx*dlnx));
    gsl_vector_set(ur,0,  -dt*km/(x[1]*x[1]*dlnx*dlnx));
    gsl_vector_set(lr,nx-1,  -dt*km/(x[nx]*x[nx]*dlnx*dlnx));
    gsl_vector_set(diag,nx,  1.+dt*km/(x[nx]*x[nx]*dlnx*dlnx));

    gsl_linalg_solve_tridiag(diag,ur,lr,diffLogZ1,diffLogZ2);

    for(n=1; n<=nx; ++n) {
        ZDisk[n]=pow(10.0,gsl_vector_get(diffLogZ2,n));
        if(ZDisk[n]!=ZDisk[n])
            errormsg("Nonphysical metallicity: n,ZDisk2,ZBulge: "+str(n)
                    +" "+str(ZDisk[n])+" "+str(ZBulge));
    }

    ZBulge=pow(10.0,gsl_vector_get(diffLogZ2,0));

    gsl_vector_free(lr); gsl_vector_free(diag); gsl_vector_free(ur);
    gsl_vector_free(diffLogZ1); gsl_vector_free(diffLogZ2);
}


// Note that the KMT09 model sets ch as in the commented line below,
//  but in the iterative model of K13, it's an input.
double DiskContents::ComputeH2Fraction(double ch, double thisCol, double thisZ)
{

    // Krumholz & Dekel 2011
    double Z0 = thisZ/Z_Sol;
    double clumping = 5.0;
    double Sig0 = dim.col_cgs(thisCol);
    //    double ch = 3.1 * (1.0 + 3.1*pow(Z0,0.365))/4.1;
    double tauc = 320.0 * clumping * Sig0 * Z0;
    double ss = log(1.0 + 0.6 * ch + .01*ch*ch)/(0.6*tauc);
    double val = 1.0 - 0.75 * ss/(1.0+0.25*ss);
    // if(val<fH2Min ) val = fH2Min; // This "feature" was an old hack to get around KMT09's deficiencies. It's solved by going to K13.
    if(val<0) val=0.0;
    if(val<0. || val>1.0 || val!=val)
        errormsg("Nonphysical H2 Fraction :" + str(val) + 
                ", ch,tauc,ss,ZDisk,ZBulge,col= " +
                " "+str(ch)+" "+str(tauc)+" "+str(ss)+" "
                +str(thisZ)+" "+str(ZBulge)+" "+str(thisCol));
    return val;
}

double DiskContents::ComputeRhoSD(unsigned int n)
{
    //// The dark matter density at a given cell indexed by n.
    // double r = x[n]*dim.Radius; // cm
    return GetCos().rhoEinasto(n) + activeColSt(n)*dim.col_cgs(1.0) / hStars(n);

}

void DiskContents::ZeroDuDt()
{
    for(unsigned int n=1; n<=nx; ++n)
        dudt[n] = 0.0;
}

// Take a pointer to a vector (0th element assumed blank), and average it over a window from -k to k.
void windowAverage(std::vector<double>& arr, int k)
{
    int nx=arr.size()+1;
    std::vector<double> arrcopy(nx+1);
    // First copy the array
    for(int n=1; n<=nx; ++n) {
        arrcopy[n] = arr[n];
    }
    for(int n=1; n<=nx; ++n) {
        double accum=0;
        double counter=0;
        for(int nn=n-k; nn<=n+k; ++nn) {
            if(nn>=1 && nn<=nx) {
                accum += arrcopy[nn];
                counter += 1.0;
            }
        }
        double val = accum/counter;
        if(val<0 || val!=val) {
            errormsg("Nonphysical value of smoothed quantity. n, val, accum, counter: "+str(n)+" "+str(val)+" "+str(accum)+" "+str(counter));
        }
        arr[n] = val;
    }
}


void DiskContents::UpdateRotationCurve(double Mh, double z, double dt)
{
    std::vector<unsigned int> thin(nx+1,0.);
    for(unsigned int n=1; n<=nx; ++n) {
        // initialize
        dudt[n] = -uu[n]*sqrt(2.*(beta[n]+1.))/dt;
        uu[n] = 0.0;

        // If the scale height is larger than the radius, the thin disk approximation makes little sense.
        if(hStars(n) > x[n]*dim.Radius * 0.5) { 
            thin[n] = 0;
        }
        else {
            thin[n]=1;
        }

    }

    double * col_copy = new double[nx];
    double * colst_copy = new double[nx];
    //std::cout << "colst_copy before: ";
    for(unsigned int n=0; n<nx; ++n) {
        col_copy[n] = log(col[n+1]);
        colst_copy[n] = log(activeColSt(n+1));
        //std::cout << colst_copy[n] << " ";
    }
    //std::cout << std::endl;
    gsl_fft_real_radix2_transform( col_copy, 1, nx );
    gsl_fft_real_radix2_transform( colst_copy, 1, nx );
    //std::cout << "colst_copy during: ";
    //double ksuppress = 10.0;
    //double kpower = 2.0;
    for(unsigned int n=1; n<nx/2; ++n) {
        if(dbg.opt(7)) { 
            double exparg = -pow( ((double) n) /ksuppress, kpower) ;
            double theexp = exp(exparg);
            //std::cout << "Debug DFT. n, exparg, ksuppress, kpower, exp: "<<n<<" "<<exparg<<" "<<ksuppress<<" "<<kpower<<" "<<exp(exparg)<<std::endl;
            col_copy[n] = col_copy[n] * theexp;
            colst_copy[n] = colst_copy[n] * theexp;
            col_copy[nx-n] = col_copy[nx-n] * theexp;
            colst_copy[nx-n] = colst_copy[nx-n] * theexp;
        }
    }
    if( dbg.opt(7) ) {
        col_copy[nx/2] = col_copy[nx/2] * exp(-pow(((double) nx)/(2.0*ksuppress), kpower));
        colst_copy[nx/2] = colst_copy[nx/2] * exp(-pow( ((double) nx)/(2.0*ksuppress), kpower));
    }
    //std::cout << std::endl;

    gsl_fft_halfcomplex_radix2_inverse(col_copy, 1, nx);
    gsl_fft_halfcomplex_radix2_inverse(colst_copy, 1, nx);


    //std::cout << "colst_copy after: ";
    for(unsigned int n=1; n<=nx; ++n) {
        //std::cout << colst_copy[n-1] << " ";
        colvPhiDisk[n] = exp(col_copy[n-1]);
        colstvPhiDisk[n] = exp(colst_copy[n-1]);
        double haloContrib;
        double bulgeContrib;
        double thinDiskContrib = 0.0;
        double thickDiskContrib = 0.0;
        for(unsigned int nn=1; nn<=nx; ++nn) {
            //thinDiskContrib += (activeColSt(nn)*((double) thin[nn])+col[nn])*mesh.summandTabulated(n,nn);
            thinDiskContrib += (exp(colst_copy[nn-1])*((double) thin[nn])+exp(col_copy[nn-1]))*mesh.summandTabulated(n,nn);
        }
        if(thinDiskContrib < 0) {
            // std::cout << "WARNING: Negative thinDiskContrib in DiskContents::UpdateRotationCurve! z,n,val: "<<z<<" "<<n<<" "<<thinDiskContrib<<std::endl;
            thinDiskContrib=0.0;
        }
        thinDiskContrib *= 2.0*M_PI*x[n]*dim.chi();
        for(unsigned int nn=1; nn<=n; ++nn) {
                thickDiskContrib += activeColSt(nn)*((double) (1-thin[nn]))*mesh.area(nn);
        }
        thickDiskContrib *= dim.chi()/x[n];
        // v^2_bulge = G MB/r
        // u^2_bulge vphiR^2 = G MBulge(dimensionless) * Mdotext *2pi R/vphiR / (x R)
        // u^2_bulge = MBulge(dimensionless) /x  * 2 pi chi
        bulgeContrib = MBulge/x[n] * 2.0*M_PI*dim.chi();
        haloContrib = G*GetCos().MrEinasto(n)*MSol/(x[n]*dim.Radius*dim.vphiR*dim.vphiR);
//        std::cout << "n, halo contribution: "<<n<<" "<<haloContrib/(haloContrib+uu[n]) << " "<<uu[n]<<" "<<haloContrib<<" "<<GetCos().MrEinasto(x[n]*dim.Radius, Mh, z)/Mh << std::endl;
        uDM[n] = sqrt(haloContrib);
        uBulge[n] = sqrt(bulgeContrib);
        uDisk[n] = sqrt(thickDiskContrib+thinDiskContrib);
    }
    // std::cout << std::endl;

    delete[] col_copy;
    delete[] colst_copy;
//    std::vector<double> uDiskCopy(nx+1,0.0);
//    for(unsigned int n=1; n<=nx; ++n) {
//        uDiskCopy[n] = uDisk[n];
//    }
//    for(unsigned int n=1; n<=nx; ++n) {
//        double sum=0.0;
//        double nsum=0.0;
//        for(int k=n-5; k<=n+5; ++k) {
//            if(k>=1 && k<=nx) {
//                sum += uDiskCopy[n]*uDiskCopy[n];
//                nsum += 1.0;
//            }
//        }
//        uDisk[n] = sqrt(sum/nsum);
//    }
    for(unsigned int n=1; n<=nx; ++n) {
        uu[n] = sqrt(uDM[n]*uDM[n] + uBulge[n]*uBulge[n] + uDisk[n]*uDisk[n]);//sqrt(haloContrib+bulgeContrib+thickDiskContrib+thinDiskContrib);

//        if(n==nx/3 && dbg.opt(13)) {
 //           std::cout << "Debug vrot "<< haloContrib<<" "<<bulgeContrib<<" "<<thinDiskContrib<<" "<<thickDiskContrib <<std::endl;
  //      }
    }
    // Before computing beta, smooth u to prevent sudden changes in the first derivative!
    windowAverage( uu, 20 );

    for(unsigned int n=2; n<=nx-1; ++n) {
        beta[n] = log(uu[n+1]/uu[n-1])/log(x[n+1]/x[n-1]);
        if(beta[n]!=beta[n]) {
            errormsg("Nonphysical value of beta! u[n-1],u[n],u[n+1], x[n-1],x[n],x[n+1]: "+str(uu[n-1])+" "+str(uu[n])+
                    " "+str(uu[n+1])+" "+str(x[n-1])+" "+str(x[n])+" "+str(x[n+1]));;
        }
        if(beta[n]>0.99) {
            // Prevent non-physical situation wherein rotation curve is steeper than solid body
            //std::cout << "WARNING: steep ascending rotation curve in DiskContents::UpdateRotationCurve! z,n,val: "<<z<<" "<<n<<" "<<beta[n]<<std::endl;
            beta[n]=.99;
        }
        if(beta[n]<-.5) {
            // Prevent non-physical situation wherein rotation curve is much steeper than Keplerian
            //std::cout << "WARNING: steep declining rotation curve in DiskContents::UpdateRotationCurve! z,n,val: "<<z<<" "<<n<<" "<<beta[n]<<std::endl;
            beta[n]=-0.5;
        }
    }
    beta[1]=beta[2];
    beta[nx]=beta[nx-1];
    for(unsigned int n=2; n<=nx-1; ++n) {
        betap[n] = (beta[n+1]-beta[n-1])/(x[n+1]-x[n-1]);
    }
    betap[1]=0.0;
    betap[nx]=0.0;

    for(unsigned int n=1; n<=nx; ++n) {
        dudt[n] += uu[n]*sqrt(2.*(beta[n]+1.))/dt;
    }
    
}

// return the scale height of the stellar disk in cm
double DiskContents::hStars(unsigned int n)
{
    double sigz = activeSigStZ(n);
    return sigz*sigz*dim.vphiR*dim.vphiR/(M_PI*G*((col[n]+activeColSt(n))*dim.col_cgs(1.0)));
}

// return the scale height of the gas disk in cm.
double DiskContents::hGas(unsigned int n)
{
    return dim.Radius* sig[n]*sig[n]/(dim.chi()*M_PI*(col[n] + sig[n]/activeSigStZ(n) * activeColSt(n)));
}

double DiskContents::ComputeColSFRapprox(double Mh, double z)
 {
     double colConv = dim.col_cgs(1.0);
     double thirdTermConv = 32.0* 0.33 * 5.0 * 0.5 * 8.0e5*8.0e5 / (M_PI*G);
 //    double tSC = tDepH2SC*1.0e9*speryear; // seconds
     double colGMCcubed = pow(85.0*MSol / (cmperpc*cmperpc), 3.0);
     for(unsigned int n=1; n<=nx; ++n) {
         double tff = pow(M_PI,.25)/sqrt(8.0) * sig[n]*dim.vphiR / (G*pow(colGMCcubed * col[n]*dim.MdotExt0/(dim.vphiR*dim.Radius),0.25));
         double tToomre = 1.0 / (2.*M_PI
             * sqrt(M_PI)*dim.chi()*col[n]/sig[n]
             * sqrt(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigStZ(n))
             * sqrt(32.0 / (3.0*M_PI))) *2.0*M_PI*dim.Radius/dim.vphiR; // seconds
         double tdyn = min(tToomre,tff);
         double nCNM_2p_perG0p_1 = 2.30 *4.1 / (1.0+3.1*pow(ZDisk[n]/Z_Sol, 0.365)); // n/(10/cc)  --- this won't change during the iterations
         double colHI = (col[n])*colConv;
         double RH2 = 1.0;
         double rhosd = ComputeRhoSD(n);
         double thirdTerm =  thirdTermConv * rhosd / (colHI*colHI);
         double nCNM_Hydro_1_HIdom = M_PI*G*colHI*colHI/20.0 * (1.0 + 2.0*RH2 + sqrt((1.0+2.0*RH2)*(1.0+2.0*RH2) + thirdTerm)) / (11.0*( kB )*243.0);
         double n1;

         double tdepHD = 3.0*tdyn/(2.0*EPS_ff) + 22.0*speryear*1.0e9/(5.0*nCNM_Hydro_1_HIdom*ZDisk[n]/Z_Sol);
         double chGuess = 7.2/nCNM_2p_perG0p_1;
         double tdep2p = tdyn/(ComputeH2Fraction(chGuess, col[n], ZDisk[n]) * EPS_ff);

         double colSF;
         if(tdep2p < tdepHD) {
              colSF = col[n]*colConv/tdep2p;
         }
         else {
              colSF = col[n]*colConv/tdepHD;
         }
         colSFR[n] = colSF / (dim.MdotExt0  / ( 2.0*M_PI* dim.Radius*dim.Radius)); // convert back to code units
         G0[n] = colSF / (2.5e-3 *MSol / (cmperpc*cmperpc*1.0e6*speryear));
         fH2[n] = tdyn*colSF/(colHI*EPS_ff);
 
 
     }
     return -1;
 
 }


double DiskContents::ComputeColSFR(double Mh, double z)
{
    // According to Krumholz (2013), here's what we should be doing to get things right.
    // We have the following equations:
    //  nCNM = max(nCNM_2p, nCNM_hydro)
    //  chi = 7.2 G0p/n1
    //  fH2 = (formula above)
    //  G0p = colSFR / colSFR0
    //  colSFR = fH2*epsff*col/tDyn
    //
    // Now, 
    //  Guess a value for G0p.
    //  Solve for fH2
    //  Solve for colSFR
    //  Adjust G0p and iterate.
    

    double colConv = dim.col_cgs(1.0);
    double thirdTermConv = 32.0* 0.33 * 5.0 * 0.5 * 8.0e5*8.0e5 / (M_PI*G);
//    double tSC = tDepH2SC*1.0e9*speryear; // seconds
    double colGMCcubed = pow(85.0*MSol / (cmperpc*cmperpc), 3.0);
    for(unsigned int n=1; n<=nx; ++n) {
        double tff = pow(M_PI,.25)/sqrt(8.0) * sig[n]*dim.vphiR / (G*pow(colGMCcubed * col[n]*dim.MdotExt0/(dim.vphiR*dim.Radius),0.25));
        double tToomre = 1.0 / (2.*M_PI
            * sqrt(M_PI)*dim.chi()*col[n]/sig[n]
            * sqrt(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigStZ(n))
            * sqrt(32.0 / (3.0*M_PI))) *2.0*M_PI*dim.Radius/dim.vphiR; // seconds
        double tdyn = min(tToomre,tff);
        double nCNM_2p_perG0p_1 = 2.30 *4.1 / (1.0+3.1*pow(ZDisk[n]/Z_Sol, 0.365)); // n/(10/cc)  --- this won't change during the iterations
        double nCNM_Hydro_1;
        double n1;
        int nIter = 0; // number of guesses of G0p.
        bool continueFlag1 = true;
        bool continueFlag2;
        double colSF;
        double rhosd = ComputeRhoSD(n);
        double chGuess = 7.2/nCNM_2p_perG0p_1;

        double fH2_currentIteration=fH2[n];
//        double fH2_currentIteration = ComputeH2Fraction( chGuess, col[n], ZDisk[n]);
        double fH2_proposed;
        double fH2_nextIteration;

        double G0p_previous = G0[n];
        double G0p_currentIteration = G0[n];
        double G0p_proposed;
        double G0p_nextIteration;
        double damp = min(dampingFactors[n]*1.5,.5);
        while(continueFlag1 || continueFlag2) {
            // For this value of G0p, we need to solve for fH2. This again requires iteration. Take an initial guess.
//            double chGuess = 7.2/nCNM_2p_perG0p_1;
//            fH2Guess = ComputeH2Fraction( chGuess, col[n], ZDisk[n] );
            continueFlag2 = true;
            int nIter2 = 0;
            while(continueFlag2) {
                double colH2 = fH2_currentIteration* col[n] *colConv;
                double colHI = (col[n] - colH2)*colConv;
                double RH2 = colH2/colHI;
                double thirdTerm =  thirdTermConv * rhosd / (colHI*colHI);
                // dimensionless:  n/(10/cc)
                nCNM_Hydro_1 = M_PI*G*colHI*colHI/20.0 * (1.0 + 2.0*RH2 + sqrt((1.0+2.0*RH2)*(1.0+2.0*RH2) + thirdTerm)) / (11.0*( kB )*243.0);
                n1 = max(nCNM_2p_perG0p_1*G0p_currentIteration, nCNM_Hydro_1);
                chGuess = 7.2*G0p_currentIteration/n1;
                fH2_proposed = ComputeH2Fraction(chGuess, col[n], ZDisk[n]);
                //continueFlag2 = (fabs((newfH2Guess - fH2Guess)/newfH2Guess ) > 1.0e-3);
                fH2_nextIteration = exp( 0.5 * log(fH2_currentIteration) + 0.5* log(fH2_proposed) );
                continueFlag2 = (fabs((fH2_proposed- fH2_currentIteration)/fH2_currentIteration) > 1.0e-3);
                fH2_currentIteration = fH2_nextIteration;
                ++nIter2;
            }
            // }
            // Now we know fH2. Thus we can guess:
            colSF = fH2_currentIteration*EPS_ff*colConv*col[n]/tdyn;  // g/s/cm^2
            G0p_proposed = colSF/ (2.5e-3 * MSol / (cmperpc*cmperpc*1.0e6*speryear));
            G0p_nextIteration = exp((1.0-damp)*log(G0p_currentIteration) + (damp)*log(G0p_proposed));
            // are we oscillating?
            if( (G0p_nextIteration-G0p_currentIteration) * (G0p_currentIteration - G0p_previous) < 0.0 ) {
                damp*=.5;
            }
            else {
                // 1-damp_new -> .83*(1-damp_old)
                // damp_new = 1 - .83*(1-damp_old)
                // damp = 1.0 - 0.83*(1.0-damp);
            }
            // G0p_nextIteration = .95*G0p_currentIteration + .05*G0p_proposed;
            continueFlag1 = (fabs((G0p_proposed- G0p_currentIteration)/G0p_currentIteration) > 1.0e-2);
//            if(nIter++>1000) {
            if((n==199 && fH2_currentIteration>.1 && !continueFlag1) || nIter++>5000) {
                std::cerr << "*****" << std::endl;
                std::cerr << "nIter, nIter2, n, damp: " << nIter << " "<<nIter2<<" "<< n<<" "<<damp<<std::endl;
                // std::cerr << "colH2, colHI: " << colH2 << " " << colHI << std::endl;
                std::cerr << "n1, hydro, cnm, rhosd: " << n1 <<" "<<nCNM_Hydro_1<<" "<<nCNM_2p_perG0p_1*G0p_currentIteration <<" "<<rhosd<<std::endl;
                std::cerr << "rhon, rhohydro, rhocnm, rhosd: " << n1*mH*10.0 <<" "<<nCNM_Hydro_1*mH*10.0<<" "<<nCNM_2p_perG0p_1*G0p_currentIteration*mH*10.0 <<" "<<rhosd<<std::endl;
                std::cerr << "fH2_current, fH2_proposed, fH2_next (prop-old)/old" << std::endl;
                std::cerr << fH2_currentIteration<<" "<<fH2_proposed<<" "<<fH2_nextIteration<<" "<<(fH2_proposed-fH2_currentIteration)/fH2_currentIteration<<std::endl;
                std::cerr << "G0p_current, G0p_proposed, G0p_next (prop-old)/old" << std::endl;
                std::cerr << G0p_currentIteration<<" "<<G0p_proposed<<" "<<G0p_nextIteration<<" "<<(G0p_proposed-G0p_currentIteration)/G0p_currentIteration<<std::endl;
                std::cerr << "*****" << std::endl;
//                errormsg("Too many iterations computing fH2 -- top layer");
            }
//            G0p = G0p + (newG0p-G0p)/10.0; // only take half a step
          G0p_previous = G0p_currentIteration;
          G0p_currentIteration = G0p_nextIteration;

        }
        colSFR[n] = colSF / (dim.MdotExt0  / ( 2.0*M_PI* dim.Radius*dim.Radius)); // convert back to code units
        G0[n] = G0p_currentIteration;
        fH2[n] = fH2_currentIteration;
    }
    return -1;
}
void DiskContents::ComputeMassLoadingFactor(double Mh, std::vector<double>& colst)
{
    double theCurrentMLF = constMassLoadingFactor; // *pow(Mh/1.0e12, mlfMhScaling);
    double col0 = 1.0 * dim.MdotExt0/(dim.vphiR*dim.Radius) * cmperpc*cmperpc/MSol; // 1 Msun/pc^2 in code units
    for(unsigned int n=1; n<=nx; ++n) {
        double hg = hGas(n); // gas scale height in cm
        double fr = 1.5;
        mBubble[n] = (1.0-fH2[n]) * M_PI * (fr*fr - 1.0/3.0) * col[n] * hg*hg * dim.col_cgs(1.0) / MSol; // Msun
        double mGMC = 1.0e7; // solar masses
        double tauGMC = 30.0e6 * speryear / dim.t_cgs(1.0)  ; // dimensionless
        double rMerge = 67.0 * cmperpc *
            pow(.1*mGMC*.01,.032) * pow(0.5* col[n]*dim.col_cgs(1.0) / (hg * mH), -0.37) * pow(sig[n]*dim.vphiR/1.0e6,-0.4);
        // mBubble[n] * fH2[n]*col[n]*2pirdr/mGMC > (1-fH2)*2pirdr*col => more mass in bubbles than in ISM!
        // => mBubble[n] * fH2[n]/mGMC > (1-fH2) => more mass in bubbles than in ISM!
        if(mBubble[n] > (1-fH2[n])*mGMC/fH2[n] ) {
            mBubble[n] = (1.0-fH2[n])*mGMC/fH2[n];
        }
        double pConfined = 0.0; 
        if(rMerge < hg) {
            double x = 1.0 - rMerge/hg;
            //double cdf = .5 * (1.0 + gsl_sf_erf(x/sqrt(2.0)));
            //pConfined = 2.0*cdf - 1.0;
            pConfined = x;
        }
        if(pConfined<0.0 || pConfined>1.0)
            errormsg("probability not between 0 and 1");
        if(dbg.opt(0)) {
//            ColOutflows[n] = constMassLoadingFactor*colSFR[n];
            double colDM = ComputeRhoSD(n) * hg * dim.vphiR*dim.Radius/dim.MdotExt0; // g/cm**3 * cm / (Mdotext0 (g/s)/ (vPhiR(cm/s) r(cm))) -- A rough estimate of the dark matter column density
            ColOutflows[n] = theCurrentMLF*colSFR[n] * pow(col[n] / col0, mlfColScaling) * pow(col[n]/(col[n]+ colDM ), mlfFgScaling);
            if( Mh>MQuench && muQuench>ColOutflows[n]/colSFR[n] ) {
                //ColOutflows[n] = muQuench*colSFR[n] ;
                double weight = 1-exp(-((double) n)/20.0);
                ColOutflows[n] += 1.0*weight + (1-weight)*muQuench*colSFR[n];
            }
        }
        else {
            ColOutflows[n] = (1-1.0/fr) * mBubble[n] *fH2[n]*col[n]/(tauGMC * mGMC) * (1.0 - pConfined)
                + colSFR[n]*constMassLoadingFactor; // dimensionless
        }
        if(ColOutflows[n] < 0.0)
            errormsg("Negative outflow column...");

        if(colSFR[n]>0)
            MassLoadingFactor[n] = ColOutflows[n]/colSFR[n];
        else
            MassLoadingFactor[n] = theCurrentMLF;


//        if(MassLoadingFactor[n] > 1.0/(EPS_ff * fH2[n])) {
//            MassLoadingFactor[n] = 1.0/(EPS_ff* fH2[n]);
//            ColOutflows[n] = colSFR[n]*MassLoadingFactor[n];
//        }
        if(MassLoadingFactor[n]!=MassLoadingFactor[n] || MassLoadingFactor[n]<0)
            errormsg("Problem computing MLF");
    }
}



// The only place this function is used is in computing the coefficient for the torque equation which requires
// knowing d s_*/dt. This function is NOT used for actually updating s_*.
double DiskContents::dSigStRdt(unsigned int n, unsigned int sp,std::vector<StellarPop*>& sps,double ** tauvecStar,std::vector<double>& MdotiPlusHalfStar)
{
    std::vector<double>& col_st = sps[sp]->spcol ;
    std::vector<double>& sig_stR = sps[sp]->spsigR ;
    std::vector<double>& sig_stZ = sps[sp]->spsigZ ;

    double val = 0.0;


    //  val = (-tauvecStar[2][n]*(col_st[n]/activeColSt(n))/mesh.u1pbPlusHalf(n)) * (1.0/(x[n]*col_st[n]*(sig_stR[n] + sig_stZ[n]))) *
    //           (2.0*sig_stZ[n]*sps[sp]->dSigZdr[n]  //ddx(sig_stZ,n,x,false)
    //          + 3.0* sig_stR[n]*sps[sp]->dSigRdr[n]//ddx(sig_stR,n,x,false) 
    //	    + sig_stR[n]*sig_stR[n]/col_st[n] * sps[sp]->dColdr[n] //ddx(col_st,n,x,false) 
    //          + (sig_stR[n]*sig_stR[n] - sig_stZ[n]*sig_stZ[n])/x[n]);
    val =  1.0/(x[n]*col_st[n]*(sig_stR[n]+sig_stZ[n])) * ((beta[n]-1.)*uu[n]*tauvecStar[1][n]/(x[n]*x[n]) + (3.0*sig_stR[n]*sps[sp]->dSigRdr[n] + 2.0*sig_stZ[n]*sps[sp]->dSigZdr[n]) *(-tauvecStar[2][n]*(col_st[n]/activeColSt(n))/(uu[n]*(1.+beta[n]))) + sig_stR[n]*sig_stR[n]*(MdotiPlusHalfStar[n]-MdotiPlusHalfStar[n])/mesh.dx(n));




    if(sps.size()-1==sp) { // i.e. this population is having stars added presently.
        if(!dbg.opt(16)) {
            if(sigth*sigth+minsigst*minsigst <= sig[n]*sig[n]) {
                val += (sig[n]*sig[n] - sigth*sigth - sig_stR[n]*sig_stR[n])*RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stR[n]);
            }
            else { // in this case, the new stellar population will have velocity dispersion = minsigst
                val += (minsigst*minsigst - sig_stR[n]*sig_stR[n] ) *RfREC * colSFR[n]
                    /(2.0*col_st[n]*sig_stR[n]);
            }
        }
        else {
            if(minsigst*minsigst <= sig[n]*sig[n]) {
                val += (sig[n]*sig[n] - sig_stR[n]*sig_stR[n])*RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stR[n]);
            }
            else {
                val += (minsigst*minsigst - sig_stR[n]*sig_stR[n])*RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stR[n]);
            }
    
        }
    }
    return val;
}

// The only place this function is used is in computing the coefficient for the torque equation which requires
// knowing d s_*/dt. This function is NOT used for actually updating s_*.
double DiskContents::dSigStZdt(unsigned int n, unsigned int sp,std::vector<StellarPop*>& sps,double ** tauvecStar,std::vector<double>& MdotiPlusHalfStar)
{
    std::vector<double>& col_st = sps[sp]->spcol;
    std::vector<double>& sig_stR = sps[sp]->spsigR;
    std::vector<double>& sig_stZ = sps[sp]->spsigZ;

    double val = 0.0;


    val =  1.0/(x[n]*col_st[n]*(sig_stR[n]+sig_stZ[n])) * ((beta[n]-1.)*uu[n]*tauvecStar[1][n]/(x[n]*x[n]) + (3.0*sig_stR[n]*sps[sp]->dSigRdr[n] + 2.0*sig_stZ[n]*sps[sp]->dSigZdr[n]) *(-tauvecStar[2][n]*(col_st[n]/activeColSt(n))/(uu[n]*(1.+beta[n]))) + sig_stR[n]*sig_stR[n]*(MdotiPlusHalfStar[n]-MdotiPlusHalfStar[n])/mesh.dx(n));

    val *= .5;
    //  val =0.5* (-tauvecStar[2][n]*(col_st[n]/activeColSt(n))/mesh.u1pbPlusHalf(n)) * (1.0/(x[n]*col_st[n]*(sig_stR[n] + sig_stZ[n]))) *
    //           (2.0*sig_stZ[n]* sps[sp]->dSigZdr[n] // ddx(sig_stZ,n,x,false)
    //          + 3.0* sig_stR[n]* sps[sp]->dSigRdr[n] // ddx(sig_stR,n,x,false) 
    //          + sig_stR[n]*sig_stR[n]/col_st[n] * sps[sp]->dColdr[n] // ddx(col_st,n,x,false)
    //          + (sig_stR[n]*sig_stR[n] - sig_stZ[n]*sig_stZ[n])/x[n]);


    if(sps.size()-1==sp) { // i.e. this population is forming stars
        if(!dbg.opt(16)) {
            if(sigth*sigth+minsigst*minsigst <= sig[n]*sig[n]) {
                val += (sig[n]*sig[n] - sigth*sigth - sig_stZ[n]*sig_stZ[n])*RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stZ[n]);
            }
            else { // in this case, the new stellar population will have velocity dispersion = minsigst
                val += (minsigst*minsigst - sig_stZ[n]*sig_stZ[n] ) *RfREC * colSFR[n]
                    /(2.0*col_st[n]*sig_stZ[n]);
            }
        }
        else {
            if(minsigst*minsigst <= sig[n]*sig[n]) {
                val += (sig[n]*sig[n] - sig_stZ[n]*sig_stZ[n])*RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stZ[n]);
            }
            else {
                val += (minsigst*minsigst - sig_stZ[n]*sig_stZ[n] ) *RfREC*colSFR[n]
                    /(2.0*col_st[n]*sig_stZ[n]);
            }
        }
    }
    return val;
}

void DiskContents::UpdateStTorqueCoeffs(std::vector<double>& UUst, std::vector<double>& DDst, std::vector<double>& LLst, std::vector<double>& FFst)
{
    std::vector<double>& sigStR = spsActive[0]->spsigR;
    std::vector<double>& sigStZ = spsActive[0]->spsigZ;
    std::vector<double>& colst = spsActive[0]->spcol;

    for(unsigned int n=1; n<=nx; ++n) {


        //    UUst[n] = -1.0/(mesh.x(n+1)-mesh.x(n-1)) * 1.0/(uu[n]*(1+beta[n])) * 1.0/(x[n]*colst[n]) * 1.0/(sigStZ[n]+sigStR[n]) * (2.0*sigStZ[n] * spsActive[0]->dSigZdr[n] + 3.0*sigStR[n]*spsActive[0]->dSigRdr[n] + sigStR[n]*sigStR[n]*spsActive[0]->dColdr[n]/colst[n] + (sigStR[n]*sigStR[n]-sigStZ[n]*sigStZ[n])/x[n]) * 1/sigStR[n]
        //	- 1.0/colst[n] * 1.0/(x[n+1]-x[n]) * -1.0/mesh.u1pbPlusHalf(n) * 1.0/(x[n]*mesh.dx(n));
        //    DDst[n] = -1.0/colst[n] * -1.0/(mesh.x(n+1)-x[n]) * -1.0/mesh.u1pbPlusHalf(n)  * 1.0/(x[n]*mesh.dx(n))
        //	+1.0/colst[n] *  1.0/(x[n]-mesh.x(n-1)) * -1.0/mesh.u1pbPlusHalf(n-1)* 1.0/(x[n]*mesh.dx(n));
        //    LLst[n] = 1.0/(mesh.x(n+1.0)-mesh.x(n-1.0)) * 1.0/(uu[n]*(1+beta[n])) * 1.0/(x[n]*colst[n]) * 1.0/(sigStZ[n]+sigStR[n]) * (2.0*sigStZ[n] * spsActive[0]->dSigZdr[n] + 3.0*sigStR[n]*spsActive[0]->dSigRdr[n] + sigStR[n]*sigStR[n]*spsActive[0]->dColdr[n]/colst[n] + (sigStR[n]*sigStR[n]-sigStZ[n]*sigStZ[n])/x[n]) * 1/sigStR[n]
        //	- 1.0/colst[n] * 1.0/(x[n]-mesh.x(n-1)) * -1.0/mesh.u1pbPlusHalf(n-1) * 1.0/(x[n]*mesh.dx(n));


        //    UUst[n] = 1.0/(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n]) * ( -(3.0*sigStR[n]* spsActive[0]->dSigRdr[n] + 2.0*sigStZ[n]*spsActive[0]->dSigZdr[n])/((mesh.x(n+1)-mesh.x(n-1)) * uu[n]*(1.+beta[n]))  - sigStR[n]*sigStR[n]/(x[n]*(mesh.x(n+1)-x[n])*mesh.u1pbPlusHalf(n)) ) + 1.0/(colst[n]*x[n]*mesh.dx(n)*(mesh.x(n+1)-x[n])*mesh.u1pbPlusHalf(n));
        //    DDst[n] = 1.0/(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n]) * ( (beta[n]-1)*uu[n]/(x[n]*x[n]) + sigStR[n]*sigStR[n]/mesh.dx(n) * ( 1.0/((mesh.x(n+1)-x[n])*mesh.u1pbPlusHalf(n)) + 1.0/((x[n]-mesh.x(n-1))*mesh.u1pbPlusHalf(n-1)))) - 1.0/(colst[n]*x[n]*mesh.dx(n)) * ( 1.0/((mesh.x(n+1)-x[n])*mesh.u1pbPlusHalf(n)) + 1.0/((x[n]-mesh.x(n-1))*mesh.u1pbPlusHalf(n-1)));
        //    LLst[n] = 1.0/(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n]) * ( (3.0*sigStR[n]*spsActive[0]->dSigRdr[n] + 2.0*sigStZ[n]*spsActive[0]->dSigZdr[n] )/((mesh.x(n+1)-mesh.x(n-1))*uu[n]*(1.+beta[n])) - sigStR[n]*sigStR[n]/(mesh.dx(n)*(x[n]-mesh.x(n-1)) * mesh.u1pbPlusHalf(n-1))) + 1.0/(colst[n]*x[n]*mesh.dx(n)*(x[n]-mesh.x(n-1))*mesh.u1pbPlusHalf(n-1));


        UUst[n] = 1.0/(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n]) * -1./(uu[n]*(1.+beta[n])) * 1./(mesh.x(n+1)-mesh.x(n-1)) * (3.*sigStR[n]*spsActive[0]->dSigRdr[n] + 2.*sigStZ[n]*spsActive[0]->dSigZdr[n]) 
            + sigStR[n]/(x[n]*mesh.dx(n)*colst[n]*(sigStR[n]+sigStZ[n]))*(-1./mesh.u1pbPlusHalf(n)) * 1./(mesh.x(n+1)-x[n]) 
            - 1./(colst[n]*x[n]*mesh.dx(n))*-1./mesh.u1pbPlusHalf(n)*1./(mesh.x(n+1)-x[n]);
        DDst[n] = sigStR[n]/(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*mesh.dx(n))*(-1./mesh.u1pbPlusHalf(n) * -1./(mesh.x(n+1)-x[n]) - -1./mesh.u1pbPlusHalf(n-1)*1./(x[n]-mesh.x(n-1)))
            - uu[n]*(1.-beta[n])/(x[n]*x[n]*x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n])
            -1./(colst[n]*x[n]*mesh.dx(n))*(-1./mesh.u1pbPlusHalf(n)*-1./(mesh.x(n+1)-x[n]) - -1./mesh.u1pbPlusHalf(n-1)*1./(x[n]-mesh.x(n-1)));
        LLst[n] = 1./(x[n]*colst[n]*(sigStR[n]+sigStZ[n])*sigStR[n])*-1./(uu[n]*(1.+beta[n]))*-1./(mesh.x(n+1)-mesh.x(n-1))*(3.*sigStR[n]*spsActive[0]->dSigRdr[n]+2.*sigStZ[n]*spsActive[0]->dSigZdr[n])
            + sigStR[n]/(x[n]*mesh.dx(n)*colst[n]*(sigStR[n]+sigStZ[n])) * 1./mesh.u1pbPlusHalf(n-1)* -1./(x[n]-mesh.x(n-1))
            - 1./(colst[n]*x[n]*mesh.dx(n)) * 1./mesh.u1pbPlusHalf(n-1) * -1./(x[n]-mesh.x(n-1));


        double Qst = sqrt(2.*(beta[n]+1.))*uu[n]*sigStR[n]/(M_PI*dim.chi()*x[n]*colst[n]);
        FFst[n] = (Qlim - Qst)*uu[n] / (x[n]*tauHeat*Qst);

        // Instead of just allowing Q to slowly approach Qlim, make the process more vigorous.
        if(dbg.opt(4)) {
            FFst[n] = exp(FFst[n]);
        }

        if(Qst > Qlim) {
            FFst[n] = 0.0;
            UUst[n] = 0.0;
            LLst[n] = 0.0;
            DDst[n] = 1.0;
        }
        double dummy = 1.0;
    }

}

void DiskContents::UpdateCoeffs(double redshift, std::vector<double>& UU, std::vector<double>& DD,
            std::vector<double>& LL, std::vector<double>& FF,
            double ** tauvecStar,std::vector<double>& MdotiPlusHalfStar, 
            double ** tauvecMRI, std::vector<double>& MdotiPlusHalfMRI,
            std::vector<double> & accProf, double AccRate,
            std::vector<std::vector<int> >& flagList)
{
    double absc = 1.;
    RafikovQParams rqp;
    //  std::vector<double> LL(nx+1,0.), DD(nx+1,0.0), UU(nx+1,0.0), FF(nx+1,0.0);
    std::vector<double> Qs(nx+1,fixedQ);
    for(unsigned int n=1; n<=nx; ++n) {
        ComputeRafikovQParams(&rqp,n);
        Qs[n] = Q(&rqp,&absc);
    }



    for(unsigned int n=1; n<=nx; ++n) {


        double ddxSig = ddx(sig,n,x,true,!dbg.opt(9));
        std::vector<int> & flags = flagList[n];
        double baseDdxTerm = dQds[n] * (-5.0/(3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n]));
        double UUddx = baseDdxTerm * ddxSig * 1.0/(mesh.x(n+1)-mesh.x(n-1));
        double LLddx = baseDdxTerm * ddxSig *-1.0/(mesh.x(n+1)-mesh.x(n-1));
        double DDddx = 0.0;
        if(dbg.opt(17)) {
            ddxSig = ddxUpstream(sig, x, flags, n);
            double denom = (flags[2]*mesh.x(n+1) +flags[1]*x[n] +flags[0]*mesh.x(n-1));
            UUddx = baseDdxTerm * ddxSig * flags[2]/denom;
            LLddx = baseDdxTerm * ddxSig * flags[0]/denom;
            DDddx = baseDdxTerm * ddxSig * flags[1]/denom;
        }
        UU[n] = (1.0/(x[n]*mesh.dx(n))) *(-1.0/mesh.u1pbPlusHalf(n))*(1.0/(mesh.x(n+1)-mesh.x(n))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + UUddx;
        LL[n] = (1.0/(x[n]*mesh.dx(n))) * (-1.0/mesh.u1pbPlusHalf(n-1))*(1.0/(mesh.x(n)-mesh.x(n-1))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + LLddx ;
        DD[n] = ( 1.0/(mesh.u1pbPlusHalf(n)*(mesh.x(n+1)-x[n])) + 1.0/(mesh.u1pbPlusHalf(n-1)*(x[n]-mesh.x(n-1)))) * (1.0/(x[n]*mesh.dx(n))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + (uu[n]*(beta[n]-1.0)/(3.0*sig[n]*col[n]*x[n]*x[n]*x[n])) * dQds[n] + DDddx;
        // That was the easy stuff. Now we need to compute the forcing term.
        // We start with terms related to obvious source terms, i.e. the terms which 
        // appear in the continuity equation unrelated to transport through the disk.
        if(!dbg.opt(4)) {
            FF[n] = RfREC*dQdS[n]*colSFR[n] + dQdS[n]*ColOutflows[n] - dQdS[n]*diffused_dcoldt[n] - dQdS[n]*AccRate*accProf[n] - dQdu[n]*dudt[n];
            // Now we add the contribution from the changing velocity dispersion
            if(sigth<=sig[n]) {
                double Qg=  sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);
                FF[n] += dQds[n] * 2*M_PI*M_PI*(ETA*
                        pow(1.- sigth*sigth/(sig[n]*sig[n]),1.5)
                        )*col[n]*dim.chi()*(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigStZ(n))/(3.);
            }
            else {
                // do nothing - this term is zero, even though pow(<neg>,1.5) is nan.
            }
    
            // Add the contributions to FF from the stellar components. Here we see 
            // the difference between an active and a passive stellar population- 
            // only the active populations affect the gravitational stability of the disk.
            unsigned int sa = spsActive.size();
            for(unsigned int i=0; i!=sa; ++i) {
                if(sa-1 == i) { // i.e. population i is forming stars
                    FF[n] -= spsActive[i]->dQdS[n] * RfREC * colSFR[n];
                }
                FF[n] -= spsActive[i]->dQdsR[n] * dSigStRdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar) 
                    + spsActive[i]->dQdsZ[n] *dSigStZdt(n,i,spsActive,tauvecStar,MdotiPlusHalfStar)
                    + spsActive[i]->dQdS[n] * dSMigdt(n,tauvecStar,(*this),spsActive[i]->spcol);
            }
    
            // add the contribution from the MRI torques.
            FF[n] -= dQdS[n] * ((MdotiPlusHalfMRI[n] - MdotiPlusHalfMRI[n-1]) / (mesh.dx(n) * x[n]));
    
            FF[n] -= dQds[n] * ( (MdotiPlusHalfMRI[n] - MdotiPlusHalfMRI[n-1]) * sig[n] / (3.0*x[n]*mesh.dx(n)*col[n])
                    - 5.0*ddxSig*tauvecMRI[2][n] / (3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n]) 
                    +uu[n]*(beta[n]-1.)*tauvecMRI[1][n] / (3.0*sig[n]*col[n]*x[n]*x[n]*x[n]) );

            // and finally add in the contribution from non-instantaneously-recycled stellar winds:
            unsigned int sp = spsPassive.size();
            if(spsActive.size()==1) {
                for(unsigned int i=0; i!=sp; ++i) {
                    FF[n] += spsPassive[i]->dcoldtREC[n] * (-dQdS[n] + spsActive[0]->dQdS[n]);
                }
            }
            else {
                for(unsigned int i=0; i!=sa; ++i) {
                    FF[n] += spsActive[i]->dcoldtREC[n] * (-dQdS[n] + spsActive[0]->dQdS[n]);
                }
            }
        }
        
        
        //if(redshift <.1 && n>=20 && n<=40) {
        //    std::cout << "Components of torque eq: n, S,sR,sZ "<<n<<" "<<spsActive[i].dQdsR[n] * dSigStRdt(n,i,redshift,spsActive,tauvecStar)<<" "<<
        //}
        double QQ = Qs[n];// Q(&rqp,&absc);

        // Forget all of the contributions to forcing we just computed, and
        // set it to an exponential.
        if(dbg.opt(4)) {
          double faster = 10.0;
    	  FF[n] = expm1((fixedQ-QQ)*uu[n]*faster / (x[n]));

        }


        // When this cell is stable, set the torque equal to zero. This may imply a non-zero
        // mass flux if tau' is nonzero, in which case mass may flow into this cell, but that's
        // exactly what we want to happen.
        if( QQ>fixedQ ) {
            FF[n]=0.0;
            UU[n]=0.0;
            DD[n]=1.0;
            LL[n]=0.0;
        }



        if(FF[n]!=FF[n] || DD[n]!=DD[n] || LL[n]!=LL[n] || UU[n]!=UU[n]) {
            std::string spc(" ");
            errormsg(std::string("Error calculating torque eq. coefficients: FF,DD,LL,RR")
                    +std::string("   col,sig  dQdS,dQds,dQdSst,dQdsst ")+str(FF[n])+" "+str(DD[n])
                    +spc+str(LL[n])+spc+str(UU[n])+spc+spc+str(col[n])+spc+str(sig[n])
                    +spc+spc+str(dQdS[n])+" "+str(dQds[n])+spc+str(spsActive[0]->dQdS[n])
                    +spc+str(spsActive[0]->dQdsR[n])+spc+spc+str(dSigStRdt(n,0,spsActive,tauvecStar,MdotiPlusHalfStar))+spc
                    +str(dSMigdt(n,tauvecStar,*this,spsActive[0]->spcol)));
        }

    }
}


void DiskContents::WriteOutStarsFile(std::string filename,
        std::vector<StellarPop *>& sps,
        unsigned int NAgeBins,unsigned int step)
{
    std::ofstream starsFile;
    if(step != 0) {
        starsFile.open((filename+"_stars.dat").c_str(),
                std::ios::binary | std::ios::app);
    }
    else { // if this is the first step, don't append..
        starsFile.open((filename+"_stars.dat").c_str(),std::ios::binary);
    }

    int NABp1,sz,nnx;
    NABp1=NAgeBins+1; sz=sps.size(); nnx=nx;
    starsFile.write((char *) &NABp1, sizeof(NABp1));
    starsFile.write((char *) &sz,sizeof(sz));
    starsFile.write((char *) &nnx,sizeof(nnx));

    for(unsigned int n=1; n<=nx; ++n) {
        starsFile.write((char *) &(x[n]),sizeof(x[n]));
    }

    for(unsigned int i=0; i!=sps.size(); ++i) {
        double yrs = sps[i]->ageAtz0/speryear;
        double start = sps[i]->startingAge/speryear;
        double end = sps[i]->endingAge/speryear;
        starsFile.write((char *) &yrs,sizeof(yrs));
        starsFile.write((char *) &start,sizeof(start));
        starsFile.write((char *) &end,sizeof(end));
        //    std::cerr << "i, step, Age in Ga (resp): "<<i<<" "<<step<<" "<<yrs*1.0e-9<<std::endl;
        int nError = -1;
	int kError = -1;
        for(unsigned int n=1; n<=nx; ++n) {
            if(sps[i]->spcol[n]!=sps[i]->spcol[n]) {
                nError = n;
                kError = 0;
            }
            starsFile.write((char *) &(sps[i]->spcol[n]),sizeof(sps[i]->spcol[n]));
        }
        for(unsigned int n=1; n<=nx; ++n) {
            if(sps[i]->spsigR[n]!=sps[i]->spsigR[n]) {
                nError = n;
                kError = 1;
            }
            starsFile.write((char *) &(sps[i]->spsigR[n]),sizeof(sps[i]->spsigR[n]));
        }
        for(unsigned int n=1; n<=nx; ++n) {
            if(sps[i]->spsigZ[n]!=sps[i]->spsigZ[n]) {
                nError = n;
                kError = 2;
            }
            starsFile.write((char *) &(sps[i]->spsigZ[n]),sizeof(sps[i]->spsigZ[n]));
        }
        for(unsigned int n=1; n<=nx; ++n) {
            if(sps[i]->spZ[n]!=sps[i]->spZ[n]) {
                nError = n;
                kError = 3;
            }
            starsFile.write((char *) &(sps[i]->spZ[n]),sizeof(sps[i]->spZ[n]));
        }

        for(unsigned int n=1; n<=nx; ++n) {
            double zv = sqrt(sps[i]->spZV[n]);
            if(zv!=zv) {
                nError = n;
                kError = 4;
            }
            starsFile.write((char *) &(zv),sizeof(zv));
        }
	if(kError >= 0)
	  errormsg("Attempted to write out NaN to starsfile. n,k: "+str(nError)+" "+str(kError));
        //    for(unsigned int n=1; n<=nx; ++n) {
        //      starsFile.write((char *) &(sps[i].dQdS[n]),sizeof(sps[i].dQdS[n]));
        //    }
        //    for(unsigned int n=1; n<=nx; ++n) {
        //      starsFile.write((char *) &(sps[i].dQds[n]),sizeof(sps[i].dQds[n]));
        //    }
        //    for(unsigned int n=1; n<=nx; ++n) {
        //      starsFile.write((char *) &(sps[i].dQdSerr[n]),sizeof(sps[i].dQdSerr[n]));
        //    }
        //    for(unsigned int n=1; n<=nx; ++n) {
        //      starsFile.write((char *) &(sps[i].dQdserr[n]),sizeof(sps[i].dQdserr[n]));
        //    }
    }
    starsFile.close();
}
double DiskContents::ComputeQst(unsigned int n)
{
    return sqrt(2.*(beta[n]+1.))*uu[n]*activeSigStR(n)/(M_PI*dim.chi()*x[n]*activeColSt(n));
}
void DiskContents::WriteOutStepFile(std::string filename, AccretionHistory & accr,
        double t, double z, double dt, 
        unsigned int step,double **tauvec, double **tauvecStar, double ** tauvecMRI,
        std::vector<double>& MdotiPlusHalf, std::vector<double>& MdotiPlusHalfMRI,
        std::vector<double>& accProf, double fAccInner)
{
    std::ofstream file;
    if(step==0) {
        file.open((filename+"_radial.dat").c_str(),std::ios::binary);
    }
    else {
        file.open((filename+"_radial.dat").c_str(),
                std::ios::binary | std::ios::app);
    }

    RafikovQParams rqp;

    // Pretend the stars are all in one population..
    std::vector<double> col_st(nx+1),sig_stR(nx+1),sig_stZ(nx+1),Mts(nx+1);
    for(unsigned int n=1;n<=nx;++n) {
        //    col_st[n]=0.; sig_st[n]=0.;
        //    for(unsigned int i=0;i!=spsActive.size();++i) {
        //      col_st[n]+=spsActive[i]->spcol[n];
        //      sig_st[n]+=spsActive[i]->spcol[n]
        //                *spsActive[i]->spsig[n]*spsActive[i]->spsig[n];
        //    }
        //    sig_st[n]=sqrt(sig_st[n]/col_st[n]);
        col_st[n]=activeColSt(n);
        sig_stR[n]=activeSigStR(n);
        sig_stZ[n]=activeSigStZ(n);
    }

    int kError=-1;
    int nError=-1;

    // loop over each cell.
    // Print out a bunch of quantities, some of which we'll have to do
    // a few calculations to figure out.
    for(unsigned int n=1;n<=nx;++n){
        double dcol_stdt,dsig_stdt,currentQ,mrq,Qst,Qg,Q_WS,Q_RW,
               torqueErr,vrg,Q_R,lambdaT,Mt,verify;
        double alpha,fh2,taupp;
        mrq=1.;
        ComputeRafikovQParams(&rqp,n);
        currentQ=Q(&rqp,&mrq);
        //    verify =Qq(mrq,&rqp);
        rqp.analyticQ=false;
        double temp=1.0;
        Q_R = Q(&rqp,&temp);
        double temp2=temp;
        verify = Qq(temp,&rqp);
        rqp.analyticQ=true;
        Q_RW = Q(&rqp,&temp);
        //    Qst = sqrt(2.*(beta[n]+1.))*uu[n]*sig_st[n]/(M_PI*dim.chi()*x[n]*col_st[n]);
        Qst=ComputeQst(n);
        Qg  = sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);
        Q_WS = 1./(1./Qg + 1./Qst);
        //    torqueErr=h2[n]*ddx(tauvec[2],n,x) + h1[n]*tauvec[2][n] + 
        //      h0[n]*tauvec[1][n] - H[n];
        torqueErr=0;
        //dsig_stdt = -2.*M_PI*yy[n]*((1+beta[n])*uu[n]*uu[n]/(3.*sig_st[n]*x[n]) 
        //      + ddx(sig_st,n,x))
        //  + (sig[n]*sig[n] -  sig_st[n]*sig_st[n])*RfREC*colSFR[n]
        //      /(2.*col_st[n]*sig_st[n]);
        //dcol_stdt = -2.*M_PI*(col_st[n]*ddx(yy,n,x) + ddx(col_st,n,x)*yy[n] 
        //  + col_st[n]*yy[n]/x[n]) + RfREC*colSFR[n];
        vrg = (tauvec[2][n]+tauvecMRI[2][n]) / (2.*M_PI*x[n]*uu[n]*col[n]*(1.+beta[n]));
        fh2 = fH2[n]; //ComputeH2Fraction(n);
        //    taupp = (H[n] - h1[n]*tauvec[2][n] - h0[n]*tauvec[1][n])/h2[n];
        taupp = d2taudx2[n];
        if(mrq<=0) mrq=1.;
        // lambdaT is the dimensionless Toomre length
        lambdaT = 2.*M_PI*sig[n]*x[n]/(temp2*sqrt(2.*(beta[n]+1.))*uu[n]); 
        Mt = lambdaT*lambdaT*col[n];
        Mts.push_back(Mt);
        double yy = tauvecStar[2][n]/(2.0*M_PI*x[n]*col_st[n]*uu[n]*(1+beta[n]));
        double dcolst = dSMigdt(n,tauvecStar,(*this),spsActive[0]->spcol) + RfREC*colSFR[n];
  
        // actually this might not be the correct definition:
        //alpha = (-tauvec[2][nx])* dim.chi()/(3. * sig[n]*sig[n]*sig[n]);
        alpha = (-tauvec[1][n]) / (2.0*M_PI*x[n]*x[n]*sig[n]*sig[n]*col[n]);
        std::vector<double> wrt(0);
        wrt.push_back(x[n]);wrt.push_back(tauvec[1][n]);wrt.push_back(tauvec[2][n]);  // 1..3
        wrt.push_back(col[n]);wrt.push_back(sig[n]);wrt.push_back(col_st[n]);         // 4..6
        wrt.push_back(sig_stR[n]);wrt.push_back(dcoldt[n]);wrt.push_back(dsigdt[n]);   // 7..9
        wrt.push_back(dcolst);wrt.push_back(0.0);wrt.push_back(currentQ);    // 10..12
        wrt.push_back(dZDiskdtAdv[n]);wrt.push_back(MdotiPlusHalfMRI[n]*dim.MdotExt0*speryear/MSol);wrt.push_back(beta[n]);   // 13..15
        wrt.push_back(uu[n]);wrt.push_back(col[n]/(col[n]+col_st[n]));wrt.push_back(temp2); // 16..18
        wrt.push_back(lambdaT);wrt.push_back(Mt);wrt.push_back(dZDiskdt[n]); // 19..21
        wrt.push_back(ZDisk[n]);wrt.push_back(Qst);wrt.push_back(Qg);        // 22..24
        wrt.push_back(Q_R);wrt.push_back(Q_WS);wrt.push_back(Q_RW);          // 25..27
        wrt.push_back(sig_stZ[n]);wrt.push_back(colSFR[n]);wrt.push_back(accProf[n]*accr.AccOfZ(z)); // 28..30
        wrt.push_back(dQdS[n]);wrt.push_back(dQds[n]);wrt.push_back(dQdSerr[n]); // 31..33
        wrt.push_back(dQdserr[n]);wrt.push_back(yy);wrt.push_back(torqueErr); // 34..36
        wrt.push_back(vrg);wrt.push_back(CuStarsOut[n]);wrt.push_back((MdotiPlusHalf[n]+MdotiPlusHalfMRI[n])*dim.MdotExt0*speryear/MSol); // 37..39
        wrt.push_back(dsigdtTrans[n]);wrt.push_back(dsigdtDdx[n]);wrt.push_back(dsigdtHeat[n]);//40..42
        wrt.push_back(ddx(tauvec[2],n,x,false,true));wrt.push_back(ddx(sig,n,x,true,!dbg.opt(9)));wrt.push_back(dsigdtCool[n]); // 43..45
        wrt.push_back(dZDiskdtDiff[n]);wrt.push_back(alpha);wrt.push_back(fh2); // 46..48
        wrt.push_back(CumulativeTorqueErr[n]); wrt.push_back(CumulativeTorqueErr2[n]);// 49..50
        wrt.push_back(d2taudx2[n]); wrt.push_back(CumulativeSF[n]); // 51..52
        wrt.push_back(dcoldtIncoming[n]); wrt.push_back(dcoldtOutgoing[n]); // 53..54
        wrt.push_back(MassLoadingFactor[n]); wrt.push_back(colvPhiDisk[n]); wrt.push_back(colstvPhiDisk[n]); // 55..57
        
        if(n==1 ) {
            int ncol = wrt.size();
            int nrow = nx;
            file.write((char*) &ncol,sizeof(ncol));
            file.write((char*) &nrow,sizeof(nrow));
        }
        for(unsigned int k=0;k!=wrt.size();++k) {
            double a=wrt[k];
            if(a!=a) {
                kError = k;
                nError = n;
            }
            file.write((char *) &a,sizeof(a));
        }
    }
    file.close();

    if(kError!=-1)
        errormsg("Error writing file!  k,n: "+str(kError)+" "+str(nError));

    std::ofstream file2;
    if(step==0) {
        // overwrite if it already exists.
        file2.open((filename+"_evolution.dat").c_str(),std::ios::binary);
    }
    else {
        file2.open((filename+"_evolution.dat").c_str(),
                std::ios::app | std::ios::binary);
    }
    std::vector<double> wrt2(0);
    double totalMass = TotalWeightedByArea(col);
    double gasMass=totalMass;
    for(unsigned int aa=0; aa!=spsActive.size(); ++aa) {
        totalMass+= TotalWeightedByArea(spsActive[aa]->spcol);
    }
    wrt2.push_back((double) step);wrt2.push_back(t);wrt2.push_back(dt); // 1..3
    wrt2.push_back(MBulge);wrt2.push_back(ZBulge);wrt2.push_back(MHalo); // 4..6
    wrt2.push_back(gasMass/totalMass);wrt2.push_back(arrmax(Mts));wrt2.push_back(MdotiPlusHalf[0]+MdotiPlusHalfMRI[0]); // 7..9
    wrt2.push_back(z); wrt2.push_back(TotalWeightedByArea(colSFR)); // 10..11
    double currentStellarMass=0.0;
    double currentGasMass = TotalWeightedByArea(col)* (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) * (1.0/MSol);
    for(unsigned int i=0; i!=spsActive.size(); ++i) {
        currentStellarMass += TotalWeightedByArea(spsActive[i]->spcol) * 
            (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) * (1.0/MSol);
    }
    wrt2.push_back(currentGasMass-initialGasMass); //12
    wrt2.push_back(currentStellarMass-initialStellarMass);//13
    wrt2.push_back(cumulativeGasMassThroughIB);//14
    wrt2.push_back(cumulativeStellarMassThroughIB);//15
    wrt2.push_back(cumulativeStarFormationMass);//16
    wrt2.push_back(cumulativeMassAccreted);//17
    wrt2.push_back(CumulativeTorque);//18
    wrt2.push_back(accr.GetMh0() * accr.MhOfZ(z)); // 19
    wrt2.push_back(accr.AccOfZ(z) * dim.MdotExt0/MSol * speryear); // 20
    wrt2.push_back(fAccInner*accr.AccOfZ(z)); // 21
    wrt2.push_back(accr.AccStOfZ(z)); // 22
    if(step==1) {
        int ncol = wrt2.size();
        file2.write((char *) &ncol,sizeof(ncol));
        //////    file2 << wrt2.size()<<std::endl;
    }
    for(unsigned int k=0; k!=wrt2.size(); ++k) {
        double a=wrt2[k];
        file2.write((char *) &a,sizeof(a));
        ///////    file2 << wrt2[k] << " ";
    }
    ///////  file2<<std::endl;
    file2.close();
}


void DiskContents::ComputePartials()
{
    if(!analyticQ) {
        double dfridr(double (*func)(double,void*),double x,
                double h, void* p,double *err);
        RafikovQParams rqp;
        rqp.analyticQ = analyticQ;
        gsl_function F;
        double result,error,val,result2,error2;
        std::vector<double> auxpartials(0), auxpartialerrors(0);
        F.function = &varQ;
        for(unsigned int n=1; n<=nx; ++n) {
            auxpartials.clear(); auxpartialerrors.clear();
            ComputeRafikovQParams(&rqp,n);
            for(unsigned int i=0; i<=rqp.ri.size() * 2; ++i) {
                rqp.var = i; // pick the quantity to vary        F.params= &(rqp);
                if(i==0) val=(rqp).Qg; 
                else if(i<=spsActive.size()) val=(rqp).Qsi[i-1];
                else if(i<=spsActive.size()*2) val=(rqp).ri[i-1-spsActive.size()];

                //gsl_deriv_central(&F, val, val*1.e-11*red, &result2,&error2);

                result = derivDriver(varQ,val,1.0e-8,&rqp,&error);

                auxpartials.push_back(result);
                auxpartialerrors.push_back(error);
            }

            // Compute dQ/dS
            dQdS[n] = auxpartials[0] * (-1)* rqp.Qg/col[n];
            dQdSerr[n] = auxpartialerrors[0] * dQdS[n] / auxpartials[0];
            double sum=0.;
            double errsum=0.;
            for(unsigned int j=0; j<rqp.ri.size(); ++j) {
                sum += auxpartials[j+1+spsActive.size()]*-1.*rqp.ri[j]/sig[n];
                //// fabs(auxpartials[j+1+spsActive.size()] * (-1.)* rqp.ri[j]/sig[n]);
                errsum = sqrt(errsum*errsum 
                        + pow(auxpartialerrors[j+1+spsActive.size()]*rqp.ri[j]/sig[n],2.));
            }
            dQds[n] = auxpartials[0] * (rqp.Qg/sig[n]) + sum;
            dQdserr[n] = sqrt(pow(auxpartialerrors[0] * (rqp.Qg/sig[n]),2.)
                    +errsum*errsum);

            for(unsigned int k=0; k<rqp.ri.size(); ++k) {
                spsActive[k]->dQdS[n] = auxpartials[k+1] * 
                    (-1)*(rqp.Qsi[k])/(spsActive[k]->spcol[n]);
                spsActive[k]->dQdsR[n] = auxpartials[k+1] * 
                    rqp.Qsi[k]/(spsActive[k]->spsigR[n]) 
                    + auxpartials[k+1+spsActive.size()] * rqp.ri[k] 
                    / (spsActive[k]->spsigR[n]);

                spsActive[k]->dQdSerr[n] = fabs(auxpartialerrors[k+1]
                        *(-rqp.Qsi[k])/(spsActive[k]->spcol[n]));
                spsActive[k]->dQdserr[n] = sqrt(pow(auxpartialerrors[k+1] 
                            *(rqp.Qsi[k]/(spsActive[k]->spsigR[n])),2.) 
                        + pow(auxpartialerrors[k+1+spsActive.size()]
                            *rqp.ri[k]/(spsActive[k]->spsigR[n]),2.));
            }
        }
    }

    else {
        if(spsActive.size()>1)
            std::cerr << "WARNING: More active stellar populations than assumed: " << spsActive.size()<< std::endl;
        for(unsigned int n=1; n<=nx; ++n) {
            double col_st=activeColSt(n);
            double sigStR=activeSigStR(n);
            double sigStZ=activeSigStZ(n);

            double Qst = sqrt(2.*(beta[n]+1.))*uu[n]*sigStR
                /(M_PI*dim.chi()*x[n]*col_st);
            double Qg = sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]
                /(M_PI*dim.chi()*x[n]*col[n]);
            double rs = sigStR/sig[n];
            double W = 2./(rs + 1./rs);

            double thickGas = 1.5; // .8 + .7 * s_z/s_r
            double thickStars = 0.8 + 0.7 * sigStZ/sigStR;
            // Q_RW = 1./(W/Qst + 1./Qg) if Qst>Qg   or  1./(1./Qst + W/Qg) otherwise

            double uTimesfb = uu[n]*sqrt(2.*(beta[n]+1.));

            if(Qst*thickStars>Qg*thickGas) {
                dQdS[n] = -(2./3.)/(Qg*col[n]*pow((2./3.)/Qg + 2.0*sig[n]/(Qst/sigStR  * (sig[n]*sig[n]+sigStR*sigStR)*thickStars),2.0));

                dQds[n] = - ( -(2./3.)/(Qg*sig[n]) - 4.0*sig[n]*sig[n]*sigStR/(Qst*pow(sig[n]*sig[n]+sigStR*sigStR,2.0)*thickStars)   + 2.0*sigStR/(Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars)) /
                    pow( (2./3.)/Qg  + 2.0*sig[n]*sigStR/(Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars),2.0);

                spsActive[0]->dQdS[n] = - 2.0*sig[n]*sig[n] / (Qg*col[n]*(sig[n]*sig[n]+sigStR*sigStR)*thickStars*pow((2./3.)/Qg + 2.0 *sig[n] * sigStR/(Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars),2.0));

                spsActive[0]->dQdsR[n] = -( 1.4*sig[n]*sigStZ/(Qst*sigStR*(sig[n]*sig[n]+sigStR*sigStR)*thickStars*thickStars)  - 4.0*sig[n]*sigStR*sigStR/(Qst*pow(sig[n]*sig[n]+sigStR*sigStR,2.0)*thickStars)  ) /
                    pow(2./(3.*Qg)  + 2*sig[n]*sigStR/(Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars),2.0);

                spsActive[0]->dQdsZ[n] = 1.4 * sig[n] / (Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars*thickStars*pow(2./(3.*Qg) + 2.0*sig[n]*sigStR/(Qst*(sig[n]*sig[n]+sigStR*sigStR)*thickStars),2.0));

                dQdu[n] = 1.0/(W/(Qst*thickStars/uu[n]) + 1.0/(Qg*thickGas/uu[n]));
            }

            else {
                dQdS[n] = (-4./3.)*sigStR*sig[n]/(Qg*col[n]*(sig[n]*sig[n]+sigStR*sigStR)*pow((4./3.)*sigStR*sig[n]/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 1.0/(Qst*thickStars),2.0));

                dQds[n] = (8./3.)*sig[n]*sig[n]*sigStR/(Qg*pow(sig[n]*sig[n]+sigStR*sigStR,2.0)*pow((4./3.)*sig[n]*sigStR/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 1.0/(Qst*thickStars),2.0));

                spsActive[0]->dQdS[n] = -1.0/( Qst*col_st*thickStars*pow( (4./3.)*sigStR*sig[n]/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 1.0/(Qst*thickStars),2.0));

                spsActive[0]->dQdsR[n] = -(-(8./3.)*sig[n]*sigStR*sigStR/(Qg*pow(sig[n]*sig[n]+sigStR*sigStR,2.0)) + (4./3.)*sig[n]/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 0.7*sigStZ/(Qst*sigStR*sigStR*thickStars*thickStars) - 1.0/(Qst*sigStR*thickStars)) / pow((4./3.)*sigStR*sig[n]/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 1.0/(Qst*thickStars),2.0); 

                spsActive[0]->dQdsZ[n] = 0.7/(Qst*sigStR*thickStars*thickStars*pow((4./3.)*sigStR*sig[n]/(Qg*(sig[n]*sig[n]+sigStR*sigStR)) + 1.0/(Qst*thickStars),2.0));


                dQdu[n] = 1.0/(1.0/(Qst*thickStars/uTimesfb) + W/(Qg*thickGas/uTimesfb));



            }
            dQdSerr[n]=dQdserr[n]=spsActive[0]->dQdSerr[n]=spsActive[0]->dQdserr[n]=0.;

            if(dQdS[n]!=dQdS[n] || dQds[n]!=dQds[n] 
                    || spsActive[0]->dQdS[n]!=spsActive[0]->dQdS[n] 
                    || spsActive[0]->dQdsR[n]!=spsActive[0]->dQdsR[n]) {
                std::string spc(" ");
                errormsg(std::string("Error computing partials:  dQdS,dQds,dQdSst,dQdsst  ")
                        +std::string("Qst,Qg   W,rs  ")+str(dQdS[n])+spc+str(dQds[n])+spc
                        +str(spsActive[0]->dQdS[n])+spc+str(spsActive[0]->dQdsR[n])
                        +spc+spc+str(Qst)+spc+str(Qg)+spc+spc+str(W)+spc+str(rs) 
                        +std::string("  col, sig, uu, beta: ")+str(col[n])+spc+str(sig[n])
                        +spc+str(uu[n])+spc+str(beta[n])+spc+std::string("n=")+str(n));
            }

        }
    }

}

double DiskContents::activeColSt(unsigned int n)
{
    double val = 0.0;
    bool alert=false;
    for(unsigned int i=0; i!=spsActive.size(); ++i) {
        double tval = spsActive[i]->spcol[n];
        val += tval;
        if(tval < 0.0) alert=true;
    }
    if(alert) {
        std::string msg=("Negative column density! n= " + str(n)+"; ");
        for(unsigned int i=0; i!=spsActive.size(); ++i) {
            msg+=str(spsActive[i]->spcol[n]) + ", ";
        }
        errormsg(msg,true);
    }
    return val;
}

double DiskContents::activeSigStR(unsigned int n)
{
    double val = 0.0;
    for(unsigned int i=0; i!=spsActive.size(); ++i ) {
        val += spsActive[i]->spcol[n] * spsActive[i]->spsigR[n]*spsActive[i]->spsigR[n];
    }
    val/=activeColSt(n);
    val = sqrt(val);
    return val;
}

double DiskContents::activeSigStZ(unsigned int n)
{
    double val = 0.0;
    for(unsigned int i=0; i!=spsActive.size(); ++i ) {
        val += spsActive[i]->spcol[n] * spsActive[i]->spsigZ[n]*spsActive[i]->spsigZ[n];
    }
    val/=activeColSt(n);
    val = sqrt(val);
    return val;

}



