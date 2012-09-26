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

#include <iostream>
#include <fstream>

// Fill an initializer object with the current state of this disk
void DiskContents::store(Initializer& in)
{
  in.col.resize(nx+1);
  in.sig.resize(nx+1);
  in.col_st.resize(nx+1);
  in.sig_st.resize(nx+1);
  for(unsigned int n=1; n<=nx; ++n) {
    in.col[n] = col[n];
    in.sig[n] = sig[n];
    in.col_st[n] = activeColSt(n);
    in.sig_st[n] = activeSigSt(n);
  }
}

DiskContents::DiskContents(double tH, double eta,
                           double sflr,double epsff,
			   double ql,double tol,
                           bool aq, double mlf, 
			   Cosmology& c,Dimensions& d,
			   FixedMesh& m, Debug& ddbg,
			   double thk, bool migP,
                           double Qinit, double km,
                           unsigned int NA, unsigned int NP,
			   double minSigSt, double accSL) :
  nx(m.nx()),x(m.x()),beta(m.beta()),
  uu(m.uu()), betap(m.betap()),
  dim(d), mesh(m), dbg(ddbg),
  XMIN(m.xmin()),ZDisk(std::vector<double>(m.nx()+1,Z_IGM)),
  cos(c),tauHeat(tH),sigth(sflr),
  EPS_ff(epsff),ETA(eta),MassLoadingFactor(mlf),
//  spsActive(std::vector<StellarPop>(NA,StellarPop(m.nx(),0,c.lbt(1000)))),
//  spsPassive(std::vector<StellarPop>(NP,StellarPop(m.nx(),0,c.lbt(1000)))),
  spsActive(std::vector<StellarPop>(0)),
  spsPassive(std::vector<StellarPop>(0)),
  dlnx(m.dlnx()),Qlim(ql),TOL(tol),ZBulge(Z_IGM),
  yREC(.054),RfREC(0.46),zetaREC(1.0), 
  analyticQ(aq),
  thickness(thk), migratePassive(migP),
  col(std::vector<double>(m.nx()+1,0.)),
  sig(std::vector<double>(m.nx()+1,0.)),  
  dQdS(std::vector<double>(m.nx()+1,0.)), 
  dQds(std::vector<double>(m.nx()+1,0.)), 
  dQdSerr(std::vector<double>(m.nx()+1)),
  dQdserr(std::vector<double>(m.nx()+1,0.)),
  dcoldt(std::vector<double>(m.nx()+1,0.)),
  dsigdt(std::vector<double>(m.nx()+1,0.)),
  dZDiskdt(std::vector<double>(m.nx()+1,0.)),
  colSFR(std::vector<double>(m.nx()+1,0.)),
  dcoldtCos(std::vector<double>(m.nx()+1,0.)),
  keepTorqueOff(std::vector<int>(m.nx()+1,0)),
  diffused_dcoldt(std::vector<double>(m.nx()+1,0.)),
  yy(std::vector<double>(m.nx()+1,0.)),
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
  FF(std::vector<double>(m.nx()+1,0.)),
  DD(std::vector<double>(m.nx()+1,0.)),
  LL(std::vector<double>(m.nx()+1,0.)),
  UU(std::vector<double>(m.nx()+1,0.)),
  MdotiPlusHalf(std::vector<double>(m.nx()+1,0.)),
  fixedQ(Qinit),CumulativeTorque(0.0),
  kappaMetals(km),
  accScaleLength(accSL),
  minsigst(minSigSt),
  NActive(NA),
  NPassive(NP),
  dd(exp(m.dlnx())),    // define the quantity d (a number slightly larger than 1)
  dm1(expm1(m.dlnx())), // dm1 = d-1. Use expm1 to find this to high precision.
  dmm1(-expm1(-m.dlnx())),   // dmm1 = 1 -d^-1
  dmdinv(expm1(2.*m.dlnx())/exp(m.dlnx())),  // dmdinv = d -d^-1
  sqd(exp(m.dlnx()/2.)) // sqd = square root of d
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

  if(dbg.opt(13))
    zetaREC = .7;
  if(dbg.opt(14))
    zetaREC = .3;

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

}


// Initialize the disk with an Initializer object, basically all the state variables
// from a disk, to which you've done whatever... (typically relaxed to an equilibrium configuration).
void DiskContents::Initialize(Initializer& in, bool fixedPhi0)
{
  StellarPop initialStarsA(nx,
         YoungIthBin(0,cos,in.NActive),
         OldIthBin(0,cos,in.NActive));
  StellarPop initialStarsP(nx,
         YoungIthBin(0,cos,in.NPassive),
         OldIthBin(0,cos,in.NPassive));

//  InitializeGrid(in.BulgeRadius);

  double Z_Init = 0.1 * Z_Sol;
  for(unsigned int n=1; n<=nx; ++n) {
    ZDisk[n] = Z_Init;
    
    col[n] = in.col[n];
    sig[n] = in.sig[n];
    initialStarsA.spcol[n] = in.col_st[n];
    initialStarsA.spsig[n] = in.sig_st[n];
    initialStarsA.spZ[n]   = Z_Init;
    initialStarsA.spZV[n]  = 0.0;

    initialStarsP.spcol[n] = initialStarsA.spcol[n];
    initialStarsP.spsig[n] = initialStarsA.spsig[n];
    initialStarsP.spZ[n]   = initialStarsA.spZ[n];
    initialStarsP.spZV[n]  = initialStarsA.spZV[n];
  }

  // MBulge here is dimensionless:
  MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA.spcol[1]);
  initialStarsA.ageAtz0 = cos.lbt(cos.ZStart());
  initialStarsP.ageAtz0 = cos.lbt(cos.ZStart());

//  spsActive.clear();
//  spsPassive.clear();
  spsActive.push_back(initialStarsA);
  spsPassive.push_back(initialStarsP);
  EnforceFixedQ(fixedPhi0,true);

  initialStellarMass = TotalWeightedByArea(initialStarsA.spcol) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
  initialGasMass = TotalWeightedByArea(col) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}


// Simplify things. Just put in exponential disks with constant velocity dispersions.
// Set Q=Qf only when these simple initial tries yield Q<Qf.
void DiskContents::Initialize(double Z_Init, double fcool, double fg0,
			      double sig0, double phi0, double Mh0,
                              double MhZs, double stScaleLength)
{
    StellarPop initialStarsA(nx,YoungIthBin(0,cos,NActive),
			   OldIthBin(0,cos,NActive));
    StellarPop initialStarsP(nx,YoungIthBin(0,cos,NPassive),
			   OldIthBin(0,cos,NPassive));

    double xd = stScaleLength/dim.d(1.0);
    double S0 = 0.18 * fcool * MhZs*MSol / (dim.MdotExt0) * dim.vphiR/(2.0*M_PI*dim.Radius) * (1.0/(xd*xd));

    for(unsigned int n=1; n<=nx; ++n) {
      ZDisk[n] = Z_Init;
      initialStarsA.spcol[n] = S0*(1-fg0)*exp(-x[n]/xd);
      initialStarsA.spsig[n] = max(sig0 * phi0, minsigst);
      initialStarsA.spZ[n] = Z_Init;
      initialStarsA.spZV[n] = 0.0;

      initialStarsP.spcol[n] = initialStarsA.spcol[n];
      initialStarsP.spsig[n] = initialStarsA.spsig[n];
      initialStarsP.spZ[n] = initialStarsA.spZ[n];
      initialStarsP.spZV[n] = initialStarsA.spZV[n];

      col[n] = S0*fg0*exp(-x[n]/xd);
      sig[n] = max(sig0, sigth);

      
    }

  MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA.spcol[1]); // dimensionless!
  initialStarsA.ageAtz0 = cos.lbt(cos.ZStart());
  initialStarsP.ageAtz0 = cos.lbt(cos.ZStart());


  spsActive.push_back(initialStarsA);
  spsPassive.push_back(initialStarsP);
  EnforceFixedQ(true,false);

  initialStellarMass = TotalWeightedByArea(initialStarsA.spcol) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
  initialGasMass = TotalWeightedByArea(col) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;

}


// Z_Init is in absolute units (i.e solar metallicity would be ~0.02), Mh0 in solar masses
// sigst0 is in units of vphiR. stScaleLength is, as usual in units of kpc, as is BulgeRadius.
void DiskContents::Initialize(double Z_Init,double fcool, double fg0,
			      double sigst0, double Mh0,double MhZs,
			      double stScaleLength)
{			     
  StellarPop initialStarsA(nx,YoungIthBin(0,cos,NActive),
			   OldIthBin(0,cos,NActive));
  StellarPop initialStarsP(nx,YoungIthBin(0,cos,NPassive),
			   OldIthBin(0,cos,NPassive));


  if(true) {
  //  InitializeGrid(BulgeRadius);
  
    double maxsig=0.0;
    unsigned int maxsign=1;
    bool lowQst;
  
    double fc = 1.0/6.0;
    double xd = stScaleLength/dim.d(1.0);
    double S0 = 0.18 * fc * (1-fg0) * MhZs*MSol/(dim.MdotExt0) * dim.vphiR/(2.0*M_PI*dim.Radius) * (1.0/(xd*xd));
  
    //// For each cell...
    for(unsigned int n=1; n<=nx; ++n) {
      ZDisk[n] = Z_Init;

      // Put in the exponential stellar profile.
      initialStarsA.spcol[n] = S0 *exp(-x[n]/xd);
      initialStarsA.spsig[n] = max(sigst0,minsigst);



      // if the initial conditions are such that Q_* < Q_lim, set Q_*=Q_lim by 
      // heating the stars beyond what the user requested.
      if(sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA.spsig[n] / (initialStarsA.spcol[n]*M_PI*x[n]*dim.chi()) < Qlim) {
        lowQst = true;
        initialStarsA.spsig[n] = max(Qlim*M_PI*x[n]*initialStarsA.spcol[n]*dim.chi()/(sqrt(2.*(beta[n]+1.))*uu[n]),minsigst);
      }
 
      if(initialStarsA.spsig[n] > maxsig) { 
  	maxsig=initialStarsA.spsig[n];
  	maxsign=n;
      }
      initialStarsA.spZ[n] = Z_Init;
      initialStarsA.spZV[n] = 0.0;
  
      initialStarsP.spcol[n] = initialStarsA.spcol[n];
      initialStarsP.spsig[n] = initialStarsA.spsig[n];
      initialStarsP.spZ[n] = initialStarsA.spZ[n];
      initialStarsP.spZV[n] = initialStarsA.spZV[n];

      sig[n] = max(pow(dim.chi() / (ETA*fg0), 1./3.)/sqrt(2.),  sigth);
    
      col[n] = ((thickness/fixedQ)*uu[n]*sqrt(2.*(beta[n]+1.))/(M_PI*dim.chi()*x[n]) - initialStarsA.spcol[n]/initialStarsA.spsig[n]) *sig[n];


      if(col[n]<0 || sig[n] <0 || col[n]!=col[n] || sig[n]!=sig[n] || initialStarsA.spcol[n]<0.0 
          || initialStarsA.spsig[n]<0.0 || initialStarsA.spcol[n]!=initialStarsA.spcol[n] 
          || initialStarsA.spsig[n]!=initialStarsA.spsig[n]) 
      {
  	errormsg("Error initializing disk- nonphysical state vars: n, col, sig, spcol, spsig, Qst: "
                  +str(n)+" "+str(col[n])+" "+str(sig[n])+" "+str(initialStarsA.spcol[n])+" "
                  +str(initialStarsA.spsig[n])+" "
                  +str(sqrt(2.*(beta[n]+1.))*uu[n]*initialStarsA.spsig[n]
                             /(M_PI*x[n]*initialStarsA.spcol[n]*dim.chi())));
      }
    } // end loop over all cells
  
    if(lowQst) {
     // make sig_st monotonically increasing towards the center of the disk.
      for(unsigned int n=1; n<=maxsign; ++n ) {
        if(initialStarsA.spsig[n] < maxsig) { 
          initialStarsA.spsig[n] = max(maxsig,minsigst);
          initialStarsP.spsig[n] = max(maxsig,minsigst);
        }
      }
  
    }

    double minQst=1.0e30;
    unsigned int minQstN=0;
    // locate the minimum value of Q_*
    for(unsigned int n=1; n<=nx; ++n) {
      if(sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA.spsig[n]/(initialStarsA.spcol[n]*M_PI*x[n]*dim.chi()) < minQst) { 
  	minQst = sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA.spsig[n] / (initialStarsA.spcol[n]*M_PI*x[n]*dim.chi());
  	minQstN = n;
      }
    }
    if(minQst < Qlim*.99999) errormsg("Minimum Qst is somehow below Qlim. "+str(Qlim)+" "+str(minQst));

    // Inwards of minQstN, set col_* such that Q_* doesn't turn back up at small radii.
    for(unsigned int n=1; n<=minQstN; ++n) {
      initialStarsA.spcol[n] = sqrt(2.*(beta[n]+1.0))*uu[n]*initialStarsA.spsig[n] / (minQst*M_PI*x[n]*dim.chi());
      initialStarsP.spcol[n] = initialStarsA.spcol[n];
    }
    for(unsigned int n=1; n<=nx; ++n) {
      initialStarsA.spsig[n] = max(initialStarsA.spsig[n] * Qlim / minQst,minsigst);
      initialStarsP.spsig[n] = initialStarsA.spsig[n];
    }
  }



  MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA.spcol[1]); // dimensionless!
  initialStarsA.ageAtz0 = cos.lbt(cos.ZStart());
  initialStarsP.ageAtz0 = cos.lbt(cos.ZStart());


  spsActive.push_back(initialStarsA);
  spsPassive.push_back(initialStarsP);
  EnforceFixedQ(false,true);

  initialStellarMass = TotalWeightedByArea(initialStarsA.spcol) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
  initialGasMass = TotalWeightedByArea(col) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}

void DiskContents::Initialize(double tempRatio, double fg0)
{
  // an active and passive stellar population..
  StellarPop initialStarsA(nx,
      YoungIthBin(0,cos,NActive),
      OldIthBin(0,cos,NActive));
  StellarPop initialStarsP(nx,
      YoungIthBin(0,cos,NPassive),
      OldIthBin(0,cos,NPassive));

//  InitializeGrid(BulgeRadius);

  // Metallicity
  double Z_Init = 0.1*Z_Sol;
  for(unsigned int n=1; n<=nx; ++n) {
    ZDisk[n] = Z_Init;
    sig[n] = pow(dim.chi()/(ETA*fg0),1./3.)/sqrt(2.);
    col[n] = (thickness/fixedQ)*uu[n]*sqrt(2.*(beta[n]+1.))*sig[n]*
         tempRatio/(x[n]*M_PI*dim.chi() * (tempRatio + (1.-fg0)/fg0));
    initialStarsA.spcol[n] = col[n]*(1.-fg0)/fg0;
    initialStarsA.spsig[n] = max(tempRatio*sig[n],minsigst);
    initialStarsA.spZ[n]   = Z_Init;
    initialStarsA.spZV[n]  = 0.0;
    initialStarsP.spcol[n] = initialStarsA.spcol[n];
    initialStarsP.spsig[n] = initialStarsA.spsig[n];
    initialStarsP.spZ[n]   = initialStarsA.spZ[n];
    initialStarsP.spZV[n]  = initialStarsA.spZV[n];

    if(col[n]<0 || sig[n] <0 || col[n]!=col[n] || sig[n]!=sig[n] || initialStarsA.spcol[n]<0.0 || initialStarsA.spsig[n]<0.0 || initialStarsA.spcol[n]!=initialStarsA.spcol[n] || initialStarsA.spsig[n]!=initialStarsA.spsig[n]) {
	errormsg("Error initializing disk- nonphysical state vars: n, col, sig, spcol, spsig, Qst: "+str(n)+" "+str(col[n])+" "+str(sig[n])+" "+str(initialStarsA.spcol[n])+" "+str(initialStarsA.spsig[n])+" "+str(sqrt(2.*(beta[n]+1.))*uu[n]*initialStarsA.spsig[n]/(M_PI*x[n]*initialStarsA.spcol[n]*dim.chi())));
    }

  }

  MBulge = M_PI*x[1]*x[1]*(col[1]+initialStarsA.spcol[1]); // dimensionless!
  initialStarsA.ageAtz0 = cos.lbt(cos.ZStart());
  initialStarsP.ageAtz0 = cos.lbt(cos.ZStart());
  spsActive.push_back(initialStarsA);
  spsPassive.push_back(initialStarsP);
//  EnforceFixedQ(true); // this is no longer reasonable given the floor, minsigst.
  bool fixedPhi0 = initialStarsA.spsig[1] > 2.0*minsigst;
  EnforceFixedQ(fixedPhi0,true); // only allow sig and sigst to covary if initial stellar velocity dispersion is far above minsigst.
  if(!fixedPhi0)  std::cerr << "WARNING: minsigst set too high to allow intiial conditions to be set by covarying gas and stellar velocity dispersions.";

  initialStellarMass = TotalWeightedByArea(initialStarsA.spcol) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
  initialGasMass = TotalWeightedByArea(col) * 
    (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) / MSol;
}

void DiskContents::ComputeDerivs(double ** tauvec)
{
  //  for(unsigned int i=0; i<=nx-1; ++i) {
  //    MdotiPlusHalf[i] = (-1.0 /  mesh.u1pbPlusHalf(i)) * (tauvec[1][i+1] - tauvec[1][i]) / (x[i+1]-mesh.x(i));
  //  }
  //  MdotiPlusHalf[nx] = -1.0 / (uu[nx] * (1.0+ beta[nx]))  * tauvec[2][nx];


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
    double taupp = ddx(tauvec[2],n,x);
////    double taupp = d2taudx2[n];
//    
//    if(taupp!=taupp) {
//      taupp=0.;
//      std::cerr << "WARNING: torque equation may be ill-posed here- n,"<<
//          "tauvec[1],tauvec[2],H,h0,h1,h2: " << n << ", "<<tauvec[1][n]<<
//          ", "<<tauvec[2][n]<<", "<<H[n]<<", "<<h0[n]<<", "<<h1[n]<<", "
//          <<h2[n]<<std::endl;
//    }
    
    double Qg= sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);


    if(dbg.opt(11)) {
      dcoldt[n] = -taupp/((beta[n]+1.)*uu[n]*x[n]) 
        + (beta[n]*beta[n]+beta[n]+x[n]*betap[n])*tauvec[2][n]
             /((beta[n]+1.)*(beta[n]+1.)*uu[n]*x[n]*x[n])
        - RfREC * colSFR[n] - dSdtOutflows(n);
  
      dsigdt[n] = uu[n]*(beta[n]-1.)*tauvec[1][n]/
                                        (3.*sig[n]*col[n]*x[n]*x[n]*x[n])
                  + (sig[n]*(beta[n]+beta[n]*beta[n]+x[n]*betap[n])
                                 /(3.*(beta[n]+1.)*(beta[n]+1.)*col[n]*uu[n]*x[n]*x[n])
                          -5.*ddx(sig,n,x)/(3.*(beta[n]+1.)*col[n]*uu[n]*x[n]))
                    *tauvec[2][n]
                 - sig[n]*taupp/(3.*(beta[n]+1.)*col[n]*uu[n]*x[n]);
      if(sigth<=sig[n]) {
        dsigdt[n]-= 2.*M_PI*M_PI*(ETA*pow(1. - sigth*sigth/(sig[n]*sig[n]),1.5))
           *col[n]*dim.chi()*(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigSt(n))/(3.);
      }
      else {
        // do nothing - this term is zero.
      }
    }
    else {
      dcoldt[n] = (MdotiPlusHalf[n] - MdotiPlusHalf[n-1]) / (mesh.dx(n) * x[n])
                    -RfREC * colSFR[n] - dSdtOutflows(n) + dcoldtCos[n];
      
      dsigdt[n] = (MdotiPlusHalf[n] - MdotiPlusHalf[n-1]) * sig[n] / (3.0*x[n]*mesh.dx(n)*col[n])
                    - 5.0*ddx(sig,n,x)*tauvec[2][n] / (3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n])
                    +uu[n]*(beta[n]-1.)*tauvec[1][n] / (3.0*sig[n]*col[n]*x[n]*x[n]*x[n]);
      if(sig[n] >= sigth) {
        dsigdt[n] -= 2.0*M_PI*M_PI*(ETA*pow(1. - sigth*sigth/(sig[n]*sig[n]),1.5))
                        *col[n]*dim.chi()*(1.0+activeColSt(n)/col[n] * sig[n]/activeSigSt(n))/3.0;
      }
      else {
        // do nothing.
      }
    }
    
//    colSFR[n] = dSSFdt(n);
    dZDiskdt[n] = -1.0/((beta[n]+1.0)*x[n]*col[n]*uu[n]) * ZDisk[n] 
                 * dlnZdx *tauvec[2][n] 
                 + yREC*(1.0-RfREC)*zetaREC*colSFR[n]/(col[n])
		 + dcoldtCos[n] * (Z_IGM - ZDisk[n])/col[n];

    // Check to make sure the results of this method are reasonable..
    if(dcoldt[n]!=dcoldt[n] || 
       dsigdt[n]!=dsigdt[n] ||  
       dZDiskdt[n]!=dZDiskdt[n]) {
      std::string space(" ");
      errormsg(std::string("Error computing derivatives - n,dcoldt,dsigdt,dZDiskdt,")
        +std::string("tau[1],tau[2],  col,sig,taupp: ")+str(n)+" "+str(dcoldt[n])+" "
        +str(dsigdt[n])+space+str(dZDiskdt[n])+space+str(tauvec[1][n])
        +space+str(tauvec[2][n])+space+str(col[n])+space+str(sig[n])+" "
        +str(taupp));
    }
  }
}

double DiskContents::ComputeTimeStep(const double redshift,int * whichVar, int * whichCell)
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
      if(fabs(dsigdt[n]/sqrt(sig[n]*sig[n]-sigth*sigth)) > dmax) {
	dmax = fabs(dsigdt[n]/sqrt(sig[n]*sig[n]-sigth*sigth));
	*whichVar = 3;
	*whichCell=n;
      }
    }
//    else { 
//      if(fabs(dsigdt[n]/sqrt(sigth*sigth-sig[n]*sig[n])) > dmax) {
//	dmax = fabs(dsigdt[n]/sqrt(sigth*sigth - sig[n]*sig[n]));
//	*whichVar = 4;
//	*whichCell=n;
//      }
//    }

//    if(fabs((ZFlux[n]-ZFlux[n-1])/(x[n]*x[n]*dlnx*col[n]*ZDisk[n])) > dmax) {
//      dmax = fabs((ZFlux[n]-ZFlux[n-1])/(col[n]*ZDisk[n]));
//      *whichVar=4;
//      *whichCell=n;
//    }

    for(unsigned int i=0; i!=spsActive.size(); ++i) {
      if(spsActive[i].IsForming(cos,redshift)) {
        if(fabs(colSFR[n]/spsActive[i].spcol[n])>dmax) { 
          dmax=fabs(colSFR[n]/spsActive[i].spcol[n]);
	  *whichCell=n;
          *whichVar=5;
        }
      }
      if(fabs(dSMigdt(n,yy,x,spsActive[i].spcol)
               /spsActive[i].spcol[n]) > dmax) {
        dmax=fabs(dSMigdt(n,yy,x,spsActive[i].spcol)
          /spsActive[i].spcol[n]);
        *whichVar=6;
	*whichCell=n;
      }
      if(fabs(dSigstdt(n,i,redshift,spsActive)
                /spsActive[i].spsig[n]) > dmax) {
        dmax=fabs(dSigstdt(n,i,redshift,spsActive)
           /spsActive[i].spsig[n]);
        *whichVar=7;
	*whichCell=n;
      }
//      if(fabs(flux(n,yy,x,spsActive[i].spcol) / (2.0*M_PI*x[n]*yy[n]*spsActive[i].spcol[n])) > dmax) {
//	dmax=fabs(flux(n,yy,x,spsActive[i].spcol)/(2.0*M_PI*x[n]*yy[n]*spsActive[i].spcol[n]));
//	*whichVar=8;
//	*whichCell=n;
//      }
//      if(fabs(flux(n-1,yy,x,spsActive[i].spcol) /(2.0*M_PI*x[n]*yy[n]* spsActive[i].spcol[n])) > dmax) {
//	dmax=fabs(flux(n-1,yy,x,spsActive[i].spcol)/(2.0*M_PI*x[n]*yy[n]*spsActive[i].spcol[n]));
//	*whichVar=9;
//	*whichCell=n;
//      }
    } // end loop over stellar populations

//    double Qst = ComputeQst(n);
//    if(fabs(max(Qlim-Qst,0.0)*uu[n]/(2.0*M_PI*x[n])) > dmax) {
//	dmax=fabs(max(Qlim-Qst,0.0)*uu[n]/(2.0*M_PI*x[n]));
//	*whichVar=10;
//	*whichCell=n;
//    }

    //    if(fabs(2*PI*yy[n]/(x[n]*dlnx))>dmax) {
    //      dmax = fabs(2*PI*yy[n]/(x[n]*dlnx)) ;
    //      *which=7;
    //    }


    if(dmax!=dmax) errormsg("Error setting timestep. n, whichVar, whichCell: "+str(n)+" "+str(*whichVar)+" "+str(*whichCell));

  } // end loop over cells
  return TOL/max(dmax,10.0*TOL/x[1]); // maximum stepsize of .1 innermost orbital times
}

bool DiskContents::CheckStellarPops(const double dt, const double redshift,
				    std::vector<StellarPop>& sps, 
				    unsigned int NAB,bool active)
{
  bool inc=false;
  for(unsigned int i=0; i!=sps.size(); ++i) {
    inc=(inc || sps[i].IsForming(cos,redshift));
  }
  if(inc)
    return true;
  if(!inc) { // add in the new population!
    unsigned int sz=sps.size();
    StellarPop currentlyForming(nx,
        YoungIthBin(sz,cos,NAB),
        OldIthBin(sz,cos,NAB));
    currentlyForming.ageAtz0 = cos.lbt(redshift);
    if(active) { // this is a huge kludge. Used to speed up 
                 // runtimes in the case that NActive>1
      currentlyForming.extract(sps[sz-1],.01);
    }
    else {
      for(unsigned int n=1; n<=nx; ++n) {
        currentlyForming.spcol[n] = RfREC*colSFR[n]*dt;
        if(sigth*sigth+minsigst*minsigst <= sig[n]*sig[n]) currentlyForming.spsig[n] = sqrt(sig[n]*sig[n]-sigth*sigth);
	else currentlyForming.spsig[n] = minsigst;
        currentlyForming.spZ[n]   = ZDisk[n];
        currentlyForming.spZV[n]  = 0.0;

	if( currentlyForming.spcol[n] <0.0 || currentlyForming.spsig[n]<0.0 
	 || currentlyForming.spZ[n]<0.0    || currentlyForming.spZV[n] <0.0
	 || currentlyForming.spcol[n]!=currentlyForming.spcol[n]
	 || currentlyForming.spsig[n]!=currentlyForming.spsig[n]
	 || currentlyForming.spZ[n]!=currentlyForming.spZ[n]
	 || currentlyForming.spZV[n]!=currentlyForming.spZV[n])
	     errormsg("Error forming new stellar population: "+str(currentlyForming.spcol[n])+" "+str(colSFR[n])+" "+str(dt));
      }
    }
    sps.push_back(currentlyForming);
    if(active) {
      std::cout << "Creating spsActive["<<sz<<"]"<<std::endl;
    }
    else {
      std::cout << "Creating spsPassive["<<sz<<"]"<<std::endl;
    }
  }
  return false;
}

void DiskContents::UpdateStateVars(const double dt, const double redshift,
                                   double ** tauvec, double AccRate)
{
  double ostars1=spsActive[0].spcol[200];
  unsigned int szA = spsActive.size();
  unsigned int szP = spsPassive.size();
  StellarPop currentlyForming(nx,
        YoungIthBin(szA,cos,1),
        OldIthBin(szA,cos,1));
  for(unsigned int n=1; n<=nx; ++n) {
    // The stars being formed this time step have..
    currentlyForming.spcol[n] = RfREC* colSFR[n] * dt; // col. density of SF*dt 
    currentlyForming.spZ[n]=ZDisk[n];             // the metallicity of the gas
    currentlyForming.spZV[n]=0.0;
    if(sigth*sigth+minsigst*minsigst<=sig[n]*sig[n])
      currentlyForming.spsig[n] = sqrt(sig[n]*sig[n]-sigth*sigth); // the velocity dispersion of the gas
    else
	currentlyForming.spsig[n] = minsigst;
//      currentlyForming.spsig[n] = sigth; // what the velocity dispersion of the gas should be!

    if(currentlyForming.spcol[n] < 0. || currentlyForming.spsig[n]<0.0 || currentlyForming.spcol[n]!=currentlyForming.spcol[n] || currentlyForming.spsig[n]!=currentlyForming.spsig[n])
      errormsg("UpdateStateVars: newly formed stars are problematic: n, spcol, spsig, colSFR, dt, sigth:  "+str(n)+", "+str(currentlyForming.spcol[n])+", "+str(currentlyForming.spsig[n])+", "+str(colSFR[n]) +", "+str(dt)+";  sig, sigth: "+str(sig[n])+", "+str(sigth));
  }
  currentlyForming.ageAtz0 = cos.lbt(redshift);

  // are the currently-forming stars included in a pre-existing population?
  bool incA=false;
  for(unsigned int i=0; i!=spsActive.size();++i) {
    incA = (incA || spsActive[i].IsForming(cos,redshift));
    spsActive[i].MigrateStellarPop(dt,yy,(*this));
  }
  bool incP=false;
  for(unsigned int i=0; i!=spsPassive.size();++i) {
    incP = (incP || spsPassive[i].IsForming(cos,redshift));
    if(migratePassive) spsPassive[i].MigrateStellarPop(dt,yy,(*this));
  }
  double ostars2=spsActive[0].spcol[200];

  if(!incA || !incP) {
    errormsg(std::string("UpdateStateVars: currently forming stars not included in the extant stellar populations!"));
  }
  else {
    spsActive[szA-1].MergeStellarPops(currentlyForming,(*this));
    spsPassive[szP-1].MergeStellarPops(currentlyForming,(*this));
  }
  double ostars3=spsActive[0].spcol[200];

  //  std::cout << "Stars at n=200: " << ostars1 << " " << ostars2 << " "<< ostars3 <<std::endl;

//  double MIn = - dt*tauvec[2][1]/(uu[1]*(1+beta[1]));
  double MIn = dt*MdotiPlusHalf[0]+dt*dmdtCosInner(AccRate);
//  double MIn = cumulativeMassAccreted -(MassLoadingFactor+RfREC)* cumulativeStarFormationMass - MBulge - (TotalWeightedByArea(col) - initialGasMass) - (TotalWeightedByArea());
  ZBulge = (ZBulge*MBulge +MIn*ZDisk[1])/(MBulge + MIn);
  MBulge += MIn;
  CumulativeTorque+= tauvec[1][nx]*dt;
  for(unsigned int n=1; n<=nx; ++n) {
    if(n!=1) {
      for(unsigned int j=0; j!=spsActive.size(); ++j) {
        CuStarsOut[n] += 2.0*M_PI*
	  sqrt(x[n]*x[n-1] * 
	       spsActive[j].spcol[n]*spsActive[j].spcol[n-1] * 
	       yy[n]*yy[n-1]) * 
	  dt * 2.0*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR  * (1.0/MSol);
        if(CuStarsOut[n] != CuStarsOut[n])
          errormsg("Error computing CuStarsOut. "+str(spsActive[j].spcol[n])+" "+str(spsActive[j].spcol[n-1])+" "+str(yy[n])+" "+str(yy[n-1]));
      }
      CuGasOut[n] +=  
	sqrt(max(tauvec[2][n]*tauvec[2][n-1],1.0e-20))
	/ (sqrt(uu[n]*uu[n-1] * (1. + beta[n])*(1.+beta[n-1]))) 
	* dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
    }
    CuGasOut[1]+=(sqrt(max(tauvec[1][1]*tauvec[1][2],1.0e-20))-0.)
      /((XMIN*exp(dlnx/2.)*expm1(dlnx))*uu[1]*(1+beta[1])) * 
      dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
    for(unsigned int j=0; j!=spsActive.size(); ++j) {
      CuStarsOut[1] += 2*M_PI*x[1]*spsActive[j].spcol[1]*yy[1]*
	dt * 2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR * (1.0/MSol);
    }
    col[n] += dcoldt[n] * dt;

    // The gas velocity dispersion should always be above sigth.
    // If it is not, e.g. a semi-pathological initial condition or
    // some numerical issue, set it to sigth.
    // 
    if(sig[n] < sigth) {
      sig[n]=sigth;
      keepTorqueOff[n] = 1;
    }
    else {
      sig[n] += dsigdt[n] * dt;
    }  
    ZDisk[n] += dZDiskdt[n] * dt;
    CumulativeSF[n] += colSFR[n] * dt;

    double Q_RW = Qsimple(n,*this);
 
    if(keepTorqueOff[n]==1) {
      double dummy=0.0;
    }

    // Check that this method hasn't made the state variables non-physical
    if(col[n]<0. || sig[n]<0. || ZDisk[n]<0. || 
       col[n]!=col[n] || sig[n]!=sig[n] || ZDisk[n]!=ZDisk[n]) {
	std::string spc(" ");
      errormsg(std::string("Error updating statevars- dt,col,sig,ZDisk,dcoldt,dsigdt")
         +std::string(",dZDiskdt: ")+str(dt)+spc+str(col[n])+spc+str(sig[n])+spc
         +str(ZDisk[n])+spc+spc+str(dcoldt[n])+spc+str(dsigdt[n])+spc
         +str(dZDiskdt[n])+spc+spc+str(spsActive[szA-1].spcol[n])+spc
         +str(spsActive[szA-1].spsig[n])+spc+str(spsActive[szA-1].spZ[n]));
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
  for(unsigned int j=0; j!=spsActive.size(); ++j) {
    cumulativeStellarMassThroughIB += 2.*M_PI*(x[1]*dim.Radius)
     *(spsActive[j].spcol[1]*dim.MdotExt0/(dim.vphiR*dim.Radius))
     *(yy[1]*dim.vphiR) * dt * (2*M_PI*dim.Radius/dim.vphiR) * (1.0/MSol);
  }
  cumulativeMassAccreted += AccRate*dim.MdotExt0 * 
    dt * (2*M_PI*dim.Radius/dim.vphiR)* (1.0/MSol);

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
  (*p).thick = thickness;
  //  (*p).mostRecentq = 1.;
  (*p).Qg = sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);
  (*p).Qsi.clear();
  (*p).ri.clear();
  for(unsigned int i=0; i!=spsActive.size(); ++i) {
    (*p).Qsi.push_back( sqrt(2.*(beta[n]+1.))*uu[n]*spsActive[i].spsig[n]
           /(M_PI*dim.chi()*x[n]*spsActive[i].spcol[n]));
    (*p).ri.push_back( spsActive[i].spsig[n] / sig[n]);
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
    if(EnforceWhenQgrtQf || !EnforceWhenQgrtQf && Q(&rqp,&absc) < fixedQ) {

      F.params = &rqp;
      findRoot(F,&factor);
    
      sig[n] *= factor;

      if(fixedPhi0) {
        for(unsigned int i=0; i!=spsActive.size(); ++i) {
          spsActive[i].spsig[n] *= factor;
        }
        for(unsigned int i=0; i!=spsPassive.size(); ++i) {
          spsPassive[i].spsig[n] *=factor;
        }
      }
    }
  }
}

void DiskContents::ComputeMRItorque(double ** tauvec, const double alpha, 
                                    const double IBC, const double OBC, const double ndecay)
{
  std::vector<double> tauMRI(nx+1,0.0);
  std::vector<double> taupMRI(nx+1,0.0);
  std::vector<int> replaceWithMRI(nx+1,0);
  for(unsigned int n=1; n<=nx; ++n) {
    // Compute a torque and its derivative given an alpha.
    //    tauMRI[n] = 2.0*M_PI*x[n]*x[n]*col[n]*alpha*sigth*sig[n]*(beta[n]-1);
    //    taupMRI[n] = tauMRI[n] * (2.0/x[n] + ddx(col,n,x)/col[n] 
    //       + ddx(sig,n,x)/sig[n] + betap[n]/(beta[n]-1));
    tauMRI[n] = IBC;


    // Where GI has shut down and no longer transports mass inwards, 
    // allow another source of viscosity to drive gas inwards.
//    if(tauMRI[n] < tauvec[1][n] || taupMRI[n] < (tauvec[1][n+1]-tauvec[1][n-1])/(mesh.x(n+1)-mesh.x(n-1)) ) { // larger negative value
    if(tauMRI[n] < tauvec[1][n]) {
      tauvec[1][n]=tauMRI[n];
//      tauvec[2][n]=taupMRI[n];
        replaceWithMRI[n] = 1;
    }
  }


  TauPrimeFromTau(tauvec);
}


void DiskContents::ComputeGItorque(double ** tauvec, const double IBC, const double OBC)
{  
  for(unsigned int n=1; n<=nx; ++n) {
    gsl_vector_set(lr,   n-1, LL[n]);
    gsl_vector_set(diag,   n, DD[n]);
    gsl_vector_set(ur,     n, UU[n]);
    gsl_vector_set(forcing,n, FF[n]);
  }
  gsl_vector_set(ur  , 0, 0.0);
  gsl_vector_set(diag, 0, 1.0);

  gsl_vector_set(lr, nx,  -1.0/(mesh.u1pbPlusHalf(nx)*(mesh.x(nx+1)-x[nx])));
  gsl_vector_set(diag,nx+1, 1.0/(mesh.u1pbPlusHalf(nx)*(mesh.x(nx+1)-x[nx])));
  gsl_vector_set(forcing,0,IBC);
  gsl_vector_set(forcing,nx+1,OBC);

  int status = gsl_linalg_solve_tridiag(diag,ur,lr,forcing,tau);
  if(status!=GSL_SUCCESS) 
    errormsg("Failed to solve the torque equation!");

  for(unsigned int n=0; n<=nx+1; ++n) {
    tauvec[1][n] = gsl_vector_get(tau,n);
  }

}


void DiskContents::TauPrimeFromTau(double ** tauvec)
{
  for(unsigned int n=1; n<=nx; ++n) {
    tauvec[2][n] = (tauvec[1][n+1]-tauvec[1][n-1])/(mesh.x(n+1)-mesh.x(n-1));
    MdotiPlusHalf[n] = -1.0/mesh.u1pbPlusHalf(n) * (tauvec[1][n+1]-tauvec[1][n])/(mesh.x(n+1)-x[n]);

    if(dbg.opt(4) && MdotiPlusHalf[n] <0.0) {
      tauvec[2][n] = 0.0;
      MdotiPlusHalf[n] = 0.0;
    }
  }
  MdotiPlusHalf[0]=-1.0/mesh.u1pbPlusHalf(0) * (tauvec[1][1]-tauvec[1][0])/(x[1]-mesh.x(0.0));

  if(dbg.opt(4) && MdotiPlusHalf[0] < 0.0) {
      MdotiPlusHalf[0] = 0.0;
  }

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
      if(dbg.opt(15)) KM = (kappaMetals*1.0e3)*sig[n]*sig[n]*sig[n]/(col[n]*dim.chi());
      else KM = kappaMetals;
      double sum = 4.0*M_PI*KM/(mesh.dx(n)*mesh.dx(n));
      double colnp1 = col[nx];
      if(n!=nx) colnp1=col[n+1];
      double ratio = mesh.x(n+1)*mesh.x(n+1.0)*colnp1/(x[n]*x[n]*col[n]);
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

    if(dbg.opt(1)) { // Z_boundary = Z_IGM
      gsl_vector_set(MetalMass1,nx,Z_IGM*col[nx]*mesh.x(nx+1.0)*mesh.x(nx+1.0)*dlnx);
    }
    else { // Zero flux condition
      gsl_vector_set(MetalMass1,nx,ZDisk[nx]*col[nx]*mesh.x(nx+1.0)*mesh.x(nx+1.0)*dlnx);
    }

    for(unsigned int n=0; n<=nx; ++n) {
      double left=0;
      double right=0;
      if(n!=0) left=-xis[n]*x[n]*x[n]*col[n]*dlnx;
      if(n!=nx) right = etas[n]*x[n+1]*x[n+1]*col[n+1]*dlnx;
      else right = etas[n]*mesh.x(n+1.0)*mesh.x(n+1.0)*col[nx]*dlnx;
      ZFlux[n] = left+right;
      
    }

    gsl_linalg_solve_tridiag(diag,ur,lr,MetalMass1,MetalMass2);

    for(unsigned int n=1; n<=nx; ++n) {
      ZDisk[n] = gsl_vector_get(MetalMass2,n-1)/ (col[n]*x[n]*x[n]*dlnx);
      if(ZDisk[n]!=ZDisk[n] || ZDisk[n]<0.0 || ZDisk[n]>1.0)
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


double DiskContents::ComputeH2Fraction(unsigned int n)
{
  // McKee & Krumholz 2009
  ////  double ch = 3.1 * (1 + 3.1 * pow(ZDisk[n]/Z_Sol,.365)) / 4.1;
  ////  double tauc = 0.066 * dim.coldensity(col[n]) * (ZDisk[n]/Z_Sol);
  ////  double ss = log(1.+.6*ch+.01*ch*ch)/(0.6*tauc);
  ////  double val = 1.0 - 0.75*ss/(1.+0.25*ss);

  // Krumholz & Dekel 2011
  double Z0 = ZDisk[n]/Z_Sol;
  double Sig0 = dim.col_cgs(col[n]);
  double ch = 3.1 * (1.0 + 3.1*pow(Z0,0.365))/4.1;
  double tauc = 320.0 * 5.0 * Sig0 * Z0;
  double ss = log(1.0 + 0.6 * ch + .01*ch*ch)/(0.6*tauc);
  double val = 1.0 - 0.75 * ss/(1.0+0.25*ss);

  if(val<0.03) val = 0.03;
  if(val<0. || val>1.0 || val!=val)
    errormsg("Nonphysical H2 Fraction :" + str(val) + 
         ", n,ch,tauc,ss,ZDisk,ZBulge,col= " +str(n)+
         " "+str(ch)+" "+str(tauc)+" "+str(ss)+" "
         +str(ZDisk[n])+" "+str(ZBulge)+" "+str(col[n]));
  return val;
}
double DiskContents::ComputeColSFR()
{
  for(unsigned int n=1; n<=nx; ++n) {
    double fH2 = ComputeH2Fraction(n);
    double val = fH2 * 2.*M_PI*EPS_ff//*sqrt(
  //      uu[n]*col[n]*col[n]*col[n]*dim.chi()/(sig[n]*x[n]));
  	* sqrt(M_PI)*dim.chi()*col[n]*col[n]/sig[n]
  	* sqrt(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigSt(n))
  	* sqrt(32.0 / (3.0*M_PI));
    if(dbg.opt(19)) { // constant SF depletion time.
      double tdep = 2.0; // in Ga
      val = (EPS_ff/.01) * fH2 * col[n] / (tdep * 1.0e9 * speryear * dim.vphiR/ (2.0*M_PI*dim.Radius));
    }
    if(val < 0 || val!=val)
      errormsg("Error computing colSFR:  n, val, fH2, col, sig   "
         +str(n)+" "+str(val)+" "+str(ComputeH2Fraction(n))+" "+str(col[n])
         +" "+str(sig[n]));
    colSFR[n]=val;
  }
}
double DiskContents::dSdtOutflows(unsigned int n)
{
  return colSFR[n]*MassLoadingFactor;
}

double DiskContents::ComputedSdTCos(double AccRate)
{
    if(!dbg.opt(8)) return 0.0;

    double sum = 0.0;

    for(unsigned int n=1; n<=nx; ++n) {
      if(dbg.opt(8)) {
        double nD = ((double) n);
        double xlo = mesh.x(nD -0.5);
        double xhi = mesh.x(nD +0.5);
      
        dcoldtCos[n] =  (  (accScaleLength+xlo)*exp(-xlo/accScaleLength) 
                          -(accScaleLength+xhi)*exp(-xhi/accScaleLength))
                       *  AccRate*2.0 / (accScaleLength*(xhi-xlo)*(xhi+xlo));

	sum += dcoldtCos[n] * x[n]*x[n]*2.0*sinh(dlnx/2.0); 

      }
      else dcoldtCos[n] = 0.0;
    }

    sum += dmdtCosOuter(AccRate) + dmdtCosInner(AccRate);
 
}
double DiskContents::dmdtCosOuter(double AccRate)
{
    if(!dbg.opt(8)) return AccRate;

    double xb = mesh.x(.5 + ((double) nx));
    return AccRate *  (1.0 + xb/accScaleLength)*exp(-xb/accScaleLength);
}
double DiskContents::dmdtCosInner(double AccRate)
{
    if(!dbg.opt(8)) return 0.0;
    double xb = mesh.x(0.5);
    return AccRate * (1 - exp(-xb/accScaleLength) * (1.0+xb/accScaleLength));
}

double DiskContents::dSigstdt(unsigned int n, unsigned int sp,double redshift,std::vector<StellarPop>& sps)
{
  std::vector<double>& col_st = sps[sp].spcol;
  std::vector<double>& sig_st = sps[sp].spsig;

//  if(sig_st[n] <=0. || col_st[n]<=0)
//    return 0.;

//  double val =  - 2.*M_PI*yy[n]*((1.+beta[n])*uu[n]*uu[n]
//             /(3.*sig_st[n]*x[n]) + ddx(sig_st,n,x)  );

  double val = 0.0;
  
  if(n<nx) { 
    double sigp2 = (2./3.) * (mesh.psi(mesh.x(n+1))-mesh.psi(mesh.x(n))) + (1./3.) * (uu[n+1]*uu[n+1]-uu[n]*uu[n]) + sig_st[n+1]*sig_st[n+1];
    val = -2.0*M_PI/ (2.0*x[n]*x[n]*dlnx*col_st[n]*sig_st[n]) * (x[n+1]*yy[n+1]*col_st[n+1]*(sigp2-sig_st[n]*sig_st[n]));
  }
  if(sps[sp].IsForming(cos,redshift)) {
    if(sigth*sigth+minsigst*minsigst <= sig[n]*sig[n]) {
      val += (sig[n]*sig[n] - sigth*sigth - sig_st[n]*sig_st[n])*RfREC*colSFR[n]
             /(2.0*col_st[n]*sig_st[n]);
    }
    else { // in this case, the new stellar population will have velocity dispersion = minsigst
      val += (minsigst*minsigst - sig_st[n]*sig_st[n] ) *RfREC * colSFR[n]
	    /(2.0*col_st[n]*sig_st[n]);
    }
  }
  return val;
}

void DiskContents::UpdateCoeffs(double redshift)
{
  double absc = 1.;
  RafikovQParams rqp;
//  std::vector<double> LL(nx+1,0.), DD(nx+1,0.0), UU(nx+1,0.0), FF(nx+1,0.0);
  for(unsigned int n=1; n<=nx; ++n) {
    ComputeRafikovQParams(&rqp,n);
    UU[n] = (1.0/(x[n]*mesh.dx(n))) *(-1.0/mesh.u1pbPlusHalf(n))*(1.0/(mesh.x(n+1.0)-mesh.x(n))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + dQds[n] * (-5.0*ddx(sig,n,x)/(3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n]))*(1.0/(mesh.x(n+1.0)-mesh.x(n-1.0)));
    LL[n] = (1.0/(x[n]*mesh.dx(n))) * (-1.0/mesh.u1pbPlusHalf(n-1))*(1.0/(mesh.x(n)-mesh.x(n-1))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + dQds[n]*(-5.0*ddx(sig,n,x)/(3.0*(beta[n]+1.0)*x[n]*col[n]*uu[n]))*(-1.0/(mesh.x(n+1.0)-mesh.x(n-1.0)));
    DD[n] = ( 1.0/(mesh.u1pbPlusHalf(n)*(mesh.x(n+1.0)-x[n])) + 1.0/(mesh.u1pbPlusHalf(n-1)*(x[n]-mesh.x(n-1.0)))) * (1.0/(x[n]*mesh.dx(n))) * (dQdS[n] + dQds[n]*sig[n]/(3.0*col[n])) + (uu[n]*(beta[n]-1.0)/(3.0*sig[n]*col[n]*x[n]*x[n]*x[n])) * dQds[n];
    FF[n] = RfREC*dQdS[n]*colSFR[n] + dQdS[n]*dSdtOutflows(n) - dQdS[n]*diffused_dcoldt[n] - dQdS[n]*dcoldtCos[n];
    if(sigth<=sig[n]) {
      double Qg=  sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]/(M_PI*dim.chi()*x[n]*col[n]);
      FF[n] += dQds[n] * 2*M_PI*M_PI*(ETA*
				      pow(1.- sigth*sigth/(sig[n]*sig[n]),1.5)
				      )*col[n]*dim.chi()*(1.0 + activeColSt(n)/col[n] * sig[n]/activeSigSt(n))/(3.);
    }
    else {
      // do nothing - this term is zero, even though pow(<neg>,1.5) is nan.
    }
    
    // Add the contributions to H from the stellar components. Here we see 
    // the difference between an active and a passive stellar population- 
    // only the active populations contribute to H:
    for(unsigned int i=0; i!=spsActive.size(); ++i) {
      if(spsActive[i].IsForming(cos,redshift)) {
        FF[n] -= spsActive[i].dQdS[n] * RfREC * colSFR[n];
      }
      FF[n] -= spsActive[i].dQds[n] * dSigstdt(n,i,redshift,spsActive) 
	+ spsActive[i].dQdS[n] * dSMigdt(n,yy,x,spsActive[i].spcol);
    }

    double QQ = Q(&rqp,&absc);

    if(QQ>fixedQ) {
	FF[n]=0.0;
	UU[n]=0.0;
	DD[n]=1.0;
	LL[n]=-1.0;
    }



    if(FF[n]!=FF[n] || DD[n]!=DD[n] || LL[n]!=LL[n] || UU[n]!=UU[n]) {
      std::string spc(" ");
      errormsg(std::string("Error calculating torque eq. coefficients: FF,DD,LL,RR")
         +std::string("   col,sig  dQdS,dQds,dQdSst,dQdsst ")+str(FF[n])+" "+str(DD[n])
         +spc+str(LL[n])+spc+str(UU[n])+spc+spc+str(col[n])+spc+str(sig[n])
         +spc+spc+str(dQdS[n])+" "+str(dQds[n])+spc+str(spsActive[0].dQdS[n])
         +spc+str(spsActive[0].dQds[n])+spc+spc+str(dSigstdt(n,0,redshift,spsActive))+spc
	+str(dSMigdt(n,yy,x,spsActive[0].spcol)));
    }
  
  }
}


void DiskContents::WriteOutStarsFile(std::string filename,
				     std::vector<StellarPop>& sps,
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
    double yrs = sps[i].ageAtz0/speryear;
    starsFile.write((char *) &yrs,sizeof(yrs));
//    std::cerr << "i, step, Age in Ga (resp): "<<i<<" "<<step<<" "<<yrs*1.0e-9<<std::endl;
    for(unsigned int n=1; n<=nx; ++n) {
      starsFile.write((char *) &(sps[i].spcol[n]),sizeof(sps[i].spcol[n]));
    }
    for(unsigned int n=1; n<=nx; ++n) {
      starsFile.write((char *) &(sps[i].spsig[n]),sizeof(sps[i].spsig[n]));
    }
    for(unsigned int n=1; n<=nx; ++n) {
      starsFile.write((char *) &(sps[i].spZ[n]),sizeof(sps[i].spZ[n]));
    }

    for(unsigned int n=1; n<=nx; ++n) {
      double zv = sqrt(sps[i].spZV[n]);
      starsFile.write((char *) &(zv),sizeof(zv));
    }
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
  return sqrt(2.*(beta[n]+1.))*uu[n]*activeSigSt(n)/(M_PI*dim.chi()*x[n]*activeColSt(n));
}
void DiskContents::WriteOutStepFile(std::string filename, AccretionHistory & accr,
                                    double t, double z, double dt, 
                                    unsigned int step,double **tauvec)
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
  std::vector<double> col_st(nx+1),sig_st(nx+1),Mts(nx+1);
  for(unsigned int n=1;n<=nx;++n) {
//    col_st[n]=0.; sig_st[n]=0.;
//    for(unsigned int i=0;i!=spsActive.size();++i) {
//      col_st[n]+=spsActive[i].spcol[n];
//      sig_st[n]+=spsActive[i].spcol[n]
//                *spsActive[i].spsig[n]*spsActive[i].spsig[n];
//    }
//    sig_st[n]=sqrt(sig_st[n]/col_st[n]);
    col_st[n]=activeColSt(n);
    sig_st[n]=activeSigSt(n);
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
    dsig_stdt = -2.*M_PI*yy[n]*((1+beta[n])*uu[n]*uu[n]/(3.*sig_st[n]*x[n]) 
      + ddx(sig_st,n,x))
      + (sig[n]*sig[n] -  sig_st[n]*sig_st[n])*RfREC*colSFR[n]
          /(2.*col_st[n]*sig_st[n]);
    dcol_stdt = -2.*M_PI*(col_st[n]*ddx(yy,n,x) + ddx(col_st,n,x)*yy[n] 
      + col_st[n]*yy[n]/x[n]) + RfREC*colSFR[n];
    vrg = tauvec[2][n] / (2.*M_PI*x[n]*uu[n]*col[n]*(1.+beta[n]));
    fh2 = ComputeH2Fraction(n);
//    taupp = (H[n] - h1[n]*tauvec[2][n] - h0[n]*tauvec[1][n])/h2[n];
    taupp = d2taudx2[n];
    if(mrq<=0) mrq=1.;
    // lambdaT is the dimensionless Toomre length
    lambdaT = 2.*M_PI*sig[n]*x[n]/(temp2*sqrt(2.*(beta[n]+1.))*uu[n]); 
    Mt = lambdaT*lambdaT*col[n];
    Mts.push_back(Mt);

    // actually this might not be the correct definition:
    //alpha = (-tauvec[2][nx])* dim.chi()/(3. * sig[n]*sig[n]*sig[n]);
    alpha = (-tauvec[1][n]) / (2.0*M_PI*x[n]*x[n]*sig[n]*sig[n]*col[n]);
    std::vector<double> wrt(0);
    wrt.push_back(x[n]);wrt.push_back(tauvec[1][n]);wrt.push_back(tauvec[2][n]);  // 1..3
    wrt.push_back(col[n]);wrt.push_back(sig[n]);wrt.push_back(col_st[n]);         // 4..6
    wrt.push_back(sig_st[n]);wrt.push_back(dcoldt[n]);wrt.push_back(dsigdt[n]);   // 7..9
    wrt.push_back(dcol_stdt);wrt.push_back(dsig_stdt);wrt.push_back(currentQ);    // 10..12
    wrt.push_back(DD[n]);wrt.push_back(LL[n]);wrt.push_back(UU[n]);   // 13..15
    wrt.push_back(FF[n]);wrt.push_back(col[n]/(col[n]+col_st[n]));wrt.push_back(temp2); // 16..18
    wrt.push_back(lambdaT);wrt.push_back(Mt);wrt.push_back(dZDiskdt[n]); // 19..21
    wrt.push_back(ZDisk[n]);wrt.push_back(Qst);wrt.push_back(Qg);        // 22..24
    wrt.push_back(Q_R);wrt.push_back(Q_WS);wrt.push_back(Q_RW);          // 25..27
    wrt.push_back(verify);wrt.push_back(colSFR[n]);wrt.push_back(taupp); // 28..30
    wrt.push_back(dQdS[n]);wrt.push_back(dQds[n]);wrt.push_back(dQdSerr[n]); // 31..33
    wrt.push_back(dQdserr[n]);wrt.push_back(yy[n]);wrt.push_back(torqueErr); // 34..36
    wrt.push_back(vrg);wrt.push_back(CuStarsOut[n]);wrt.push_back(MdotiPlusHalf[n]*dim.MdotExt0*speryear/MSol); // 37..39
    wrt.push_back(flux(n-1,yy,x,col_st));wrt.push_back(0);wrt.push_back(0);//40..42
    wrt.push_back(ddx(tauvec[2],n,x));wrt.push_back(ddx(sig,n,x));wrt.push_back(0); // 43..45
    wrt.push_back(0);wrt.push_back(alpha);wrt.push_back(fh2); // 46..48
    wrt.push_back(CumulativeTorqueErr[n]); wrt.push_back(CumulativeTorqueErr2[n]);// 49..50
    wrt.push_back(d2taudx2[n]); wrt.push_back(CumulativeSF[n]); // 51..52

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
    totalMass+= TotalWeightedByArea(spsActive[aa].spcol);
  }
  wrt2.push_back((double) step);wrt2.push_back(t);wrt2.push_back(dt); // 1..3
  wrt2.push_back(MBulge);wrt2.push_back(ZBulge);wrt2.push_back(gasMass); // 4..6
  wrt2.push_back(gasMass/totalMass);wrt2.push_back(arrmax(Mts));wrt2.push_back(MdotiPlusHalf[0]); // 7..9
  wrt2.push_back(z); wrt2.push_back(TotalWeightedByArea(colSFR)); // 10..11
  double currentStellarMass=0.0;
  double currentGasMass = TotalWeightedByArea(col)* (2*M_PI*dim.Radius*dim.MdotExt0/dim.vphiR) * (1.0/MSol);
  for(unsigned int i=0; i!=spsActive.size(); ++i) {
    currentStellarMass += TotalWeightedByArea(spsActive[i].spcol) * 
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
  if(step==0) {
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


void DiskContents::ComputeY(const double ndecay)
{

  yy[nx]=0.;
  std::vector<double> col_st(nx+1), sig_st(nx+1), Qst(nx+1);
  for(unsigned int n=nx; n>=1; --n) {
    col_st[n]=activeColSt(n);
    sig_st[n]=activeSigSt(n);
    Qst[n] = ComputeQst(n);
    colst_gsl[n-1] = col_st[n];
    sigst_gsl[n-1] = sig_st[n];
  }
  gsl_spline_init(spline_colst,mesh.x_GSL(),colst_gsl,nx);
  gsl_spline_init(spline_sigst,mesh.x_GSL(),sigst_gsl,nx);

  Qst[0]=Qst[1];

  unsigned int NN = mesh.necessaryN();
  double dn = 1.0/(((double) NN)*((double) nx));
  double yyn=0.0;
  double yynm1=0.0;
  for(unsigned int i=NN*nx; i>NN ; --i) {
    double n = ((double) i)/((double) NN);
//  for(double n=((double) nx); n-dn>=1.0; n-=dn) {  
    double xnm1 = mesh.x(((double) i-1)/((double) NN));
    double xn = mesh.x(n);
    double Qst_nm1 = sqrt(2.*(mesh.beta(xnm1)+1.0))*mesh.uu(xnm1)* gsl_spline_eval(spline_sigst,xnm1,accel_sigst)
                       / (M_PI*dim.chi()*xnm1* gsl_spline_eval(spline_colst,xnm1,accel_colst)) ;
    if(Qst_nm1 < 1.0)
      double dummy=0.0;

    if(Qst_nm1 < 0.0 || Qst_nm1!=Qst_nm1) {
      errormsg("Error computing Qst_nm1. Qst_nm1, beta, u, sigst, colst, xnm1:   "+str(Qst_nm1)+" "+str(mesh.beta(xnm1))+" "+str(mesh.uu(xnm1))+" "+str(gsl_spline_eval(spline_sigst,xnm1,accel_sigst))+" "+str(gsl_spline_eval(spline_colst,xnm1,accel_colst))+" "+str(xnm1));
    }
    if(Qst_nm1 > Qlim && !dbg.opt(10)) {
      yyn=0.0;
      yynm1=0.0;
    }
    else {
      double sigp2 = (2.0/3.0) * (mesh.psi(xn) - mesh.psi(xnm1)) 
                    +(1.0/3.0) *(pow(mesh.uu(xn),2.0) - pow(mesh.uu(xnm1),2.0)) 
                    + pow(gsl_spline_eval(spline_sigst,xn,accel_sigst),2.0);
      double yynm1 = yyn * xn* gsl_spline_eval(spline_colst,xn,accel_colst) 
                      / (xnm1*gsl_spline_eval(spline_colst,xnm1,accel_colst))
		     *(1.5 - sigp2/(2.0*pow(gsl_spline_eval(spline_sigst,xnm1,accel_sigst),2.0)))
	- max(Qlim - Qst_nm1,0.0)*mesh.uu(xnm1)*(xn-xnm1) / (2.0*M_PI*xnm1*tauHeat*Qst_nm1);
      if(yynm1>0.0) yynm1=0.0;

      if(yynm1!=yynm1 || yynm1>0.0000001 || sigp2 <0 || sigp2 > 100.0 || sigp2!=sigp2 || fabs(yynm1) > 100.0)
        errormsg("Error computing y!   n,y,sigp2,  sigp2/2sig0^2, NN, i    dPsi, sig1^2   : "+str(n)+" "+str(yynm1)+" "+str(sigp2)+"   "+str(sigp2/(2.0*pow(gsl_spline_eval(spline_sigst,xnm1,accel_sigst),2.0)))+" "+str(NN)+" "+str(i)+"   "+str(mesh.psi(xn)-mesh.psi(xnm1))+" "+str(pow(gsl_spline_eval(spline_sigst,xn,accel_sigst),2.0)));
//      if(fabs(xnm1 - mesh.x((unsigned int) (n))) < fabs(xn-xnm1)/10.0   )
      if((i-1) % NN == 0)
        yy[(i-1)/NN] = yynm1;
      yyn = yynm1;
    }
  }

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
        spsActive[k].dQdS[n] = auxpartials[k+1] * 
           (-1)*(rqp.Qsi[k])/(spsActive[k].spcol[n]);
        spsActive[k].dQds[n] = auxpartials[k+1] * 
           rqp.Qsi[k]/(spsActive[k].spsig[n]) 
          + auxpartials[k+1+spsActive.size()] * rqp.ri[k] 
               / (spsActive[k].spsig[n]);
        
        spsActive[k].dQdSerr[n] = fabs(auxpartialerrors[k+1]
                                  *(-rqp.Qsi[k])/(spsActive[k].spcol[n]));
        spsActive[k].dQdserr[n] = sqrt(pow(auxpartialerrors[k+1] 
                                  *(rqp.Qsi[k]/(spsActive[k].spsig[n])),2.) 
                                 + pow(auxpartialerrors[k+1+spsActive.size()]
                                    *rqp.ri[k]/(spsActive[k].spsig[n]),2.));
      }
    }
  }

  else {
    if(spsActive.size()>1)
      std::cerr << "WARNING: More active stellar populations than assumed: " << spsActive.size()<< std::endl;
    for(unsigned int n=1; n<=nx; ++n) {
      double col_st=activeColSt(n);
      double sig_st=activeSigSt(n);
      
      double Qst = sqrt(2.*(beta[n]+1.))*uu[n]*sig_st
                     /(M_PI*dim.chi()*x[n]*col_st);
      double Qg = sqrt(2.*(beta[n]+1.))*uu[n]*sig[n]
                     /(M_PI*dim.chi()*x[n]*col[n]);
      double rs = sig_st/sig[n];
      double W = 2./(rs + 1./rs);
      // Q_RW = 1./(W/Qst + 1./Qg) if Qst>Qg   or  1./(1./Qst + W/Qg) otherwise
      
      if(Qst>Qg) {
        dQdS[n] = -(col[n]*pow(sig[n]*sig[n]+sig_st*sig_st,2.)
                    /pow((col[n]+2*col_st)*sig[n]*sig[n] 
                         + col[n]*sig_st*sig_st,2.))*Qg;
        dQds[n] = col[n]*(2*col_st*sig[n]*sig[n]*(sig[n]*sig[n]-sig_st*sig_st) 
                  + col[n]*pow(sig[n]*sig[n]+sig_st*sig_st,2.)) / 
                  (sig[n] *pow((col[n]+2*col_st)*sig[n]*sig[n] 
                    + col[n]*sig_st*sig_st,2.)) * Qg;
        spsActive[0].dQdS[n] = -(2.*col_st*sig[n]*sig[n]*sig[n]*
                                (sig[n]*sig[n]+sig_st*sig_st) /
                                 (sig_st*pow((col[n]+2*col_st)*sig[n]*sig[n] 
                                   + col[n]*sig_st*sig_st,2.))) * Qst;
        spsActive[0].dQds[n] = 4*col_st*col_st*sig[n]*sig[n]*sig[n]
                               /pow((col[n]+2.*col_st)*sig[n]*sig[n]
                                      +col[n]*sig_st*sig_st,2.) * Qst;
      }
      else {
        dQdS[n] = -2.*col[n]*sig_st*sig_st*sig_st*
                  (sig[n]*sig[n]+sig_st*sig_st) * 
                  Qg/(sig[n]*pow(col_st*sig[n]*sig[n] 
                    + (2.*col[n]+col_st)*sig_st*sig_st,2.));
        dQds[n] = (4.*col[n]*col[n]*sig_st*sig_st*sig_st
             /pow(col_st*sig[n]*sig[n] 
                +(2.*col[n]+col_st)*sig_st*sig_st,2.))  * Qg;
        spsActive[0].dQdS[n] = -(col_st*pow(sig[n]*sig[n]+sig_st*sig_st,2.)
             /pow(col_st*sig[n]*sig[n] 
                + (2.*col[n]+col_st)*sig_st*sig_st,2.)) * Qst;
        spsActive[0].dQds[n] = (col_st*(2.*col[n]*sig_st*sig_st
             *(sig_st*sig_st-sig[n]*sig[n]) + col_st*pow(sig[n]*sig[n]
               +sig_st*sig_st,2.)) / 
           (sig_st*pow(col_st*sig[n]*sig[n] 
           + (2.*col[n]+col_st)*sig_st*sig_st,2.)))*Qst;
      }
      dQdSerr[n]=dQdserr[n]=spsActive[0].dQdSerr[n]=spsActive[0].dQdserr[n]=0.;
      
      if(dQdS[n]!=dQdS[n] || dQds[n]!=dQds[n] 
         || spsActive[0].dQdS[n]!=spsActive[0].dQdS[n] 
         || spsActive[0].dQds[n]!=spsActive[0].dQds[n]) {
        std::string spc(" ");
        errormsg(std::string("Error computing partials:  dQdS,dQds,dQdSst,dQdsst  ")
                 +std::string("Qst,Qg   W,rs  ")+str(dQdS[n])+spc+str(dQds[n])+spc
                 +str(spsActive[0].dQdS[n])+spc+str(spsActive[0].dQds[n])
                 +spc+spc+str(Qst)+spc+str(Qg)+spc+spc+str(W)+spc+str(rs));
      }

    }
  }

}

double DiskContents::activeColSt(unsigned int n)
{
  double val = 0.0;
  bool alert=false;
  for(unsigned int i=0; i!=spsActive.size(); ++i) {
	double tval = spsActive[i].spcol[n];
	val += tval;
	if(tval < 0.0) alert=true;
  }
  if(alert) {
    std::string msg=("Negative column density! n= " + str(n)+"; ");
    for(unsigned int i=0; i!=spsActive.size(); ++i) {
      msg+=str(spsActive[i].spcol[n]) + ", ";
    }
    errormsg(msg,true);
  }
  return val;
}

double DiskContents::activeSigSt(unsigned int n)
{
  double val = 0.0;
  for(unsigned int i=0; i!=spsActive.size(); ++i ) {
    val += spsActive[i].spcol[n] * spsActive[i].spsig[n]*spsActive[i].spsig[n];
  }
  val/=activeColSt(n);
  val = sqrt(val);
  return val;
}




