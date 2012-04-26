#include <vector>
#include <string>
#include <math.h>
#include "StellarPop.h"

class Cosmology;
class Dimensions;
struct RafikovQParams;
struct Initializer;


// Main container for the physical quantities which make up the disk.
// Basic structure is a set of arrays of nx elements (indexed from 1).
class DiskContents {
 public:
  DiskContents(unsigned int nx,double xm,double tH,double eta,
               double sigth,double epsff,double ql,double tol,
               bool aq, double mlf,Cosmology&,Dimensions&,
	       double thk,bool migratePassive,double Qinit,
	       double km);
  // Sum up a quantity over the entire disk, weighting 
  // by the area of each annulus
  double TotalWeightedByArea(const std::vector<double>&);

  // Fill in a RafikovQParams struct with the parameters needed to calculate
  // Q at cell n.
  void ComputeRafikovQParams(RafikovQParams*,unsigned int n);

  // Simultaneously adjust sigma and sigma_* (keeping their ratio fixed) until
  // Q is exactly its pre-adjudicated fixed value, fixedQ, at every cell 
  // in the disk.
  void EnforceFixedQ(bool fixedPhi0);

  // Fill in vectors of the partial derivatives of Q wrt 
  // each state variable
  void ComputePartials();

  // Fill in vectors of the coefficients for the torque eq, 
  // namely h2,h1,h0,H in h2 tau'' + h1 tau' + h0 tau=H
  void UpdateCoeffs(double redshift);

  // Diffuse metals in such a way that the total mass in 
  // metals is conserved. This diffusion is not meant to 
  // model a physical process and is only used to maintain 
  // numerical stability.
  void DiffuseMetals(double dt);
//  void ComputeMetalFluxes();
  void DiffuseMetallicity(double dt, double km);
  void DiffuseMetalsUnstable(double dt, double km);

  double activeColSt(unsigned int n);
  double activeSigSt(unsigned int n);
  double ComputeQst(unsigned int n);

  // At a given cell, compute the fraction of gas which is 
  // in H2, i.e. what fraction of the gas is available to 
  // form stars. This is an analytic approximation.
  double ComputeH2Fraction(unsigned int n);

  // Compute the star formation rate in the nth cell
  double dSSFdt(unsigned int n);

  // Compute the loss in column density experienced by 
  // cell n due to outflows. (Simple mass loading factor 
  // prescription)
  double dSdtOutflows(unsigned int n);

  // Compute the time rate of change of the velocity 
  // dispersion of stellar population sp.
  double dSigstdt(unsigned int n, unsigned int sp, 
		  double reshift,std::vector<StellarPop>&);

  // Append the properties of each StellarPop in the 
  // given vector to an output file.
  void WriteOutStarsFile(std::string filename, 
			 std::vector<StellarPop>&,
			 unsigned int,unsigned int step);

  // Append the radially-dependent properties of the disk to an output file, 
  // and the purely time dependent properties to a different file
  void WriteOutStepFile(std::string filename, 
                        double t, double z, double dt, 
                        unsigned int step,double **tauvec);

  // A few self-explanatory functions...
  double GetDlnx() {return dlnx;};
  std::vector<double>& GetX() {return x;};
  std::vector<double>& GetUu() { return uu;};
  std::vector<double>& GetPsi() { return psi;};
  std::vector<double>& GetBeta() {return beta;};
  std::vector<double>& GetSig() { return sig;};
  std::vector<double>& GetCol() { return col;};
  std::vector<double>& GetColSFR() { return colSFR;}
  std::vector<StellarPop>& active() { return spsActive;}
  std::vector<StellarPop>& passive() { return spsPassive;}
  std::vector<double>& GetYy() {return yy;}
  Dimensions& GetDim() { return dim;}
  Cosmology& GetCos() { return cos;}

  // Compute the dimensionless inward velocity of stars 
  // as a function of radius- computed such that 
  // dQ_*/dt = max((Qlim - Q_*),0)/(T_mig (2 pi Omega)^-1), 
  // i.e.  the stars will attempt to reach Qlim if they are 
  // less stable than Qlim.
  void ComputeY();
  void ComputeY2();
  // Compute the time derivatives of all state variables 
  // at all radii.
  void ComputeDerivs(double **tauvec);

  // Given the state variables and their derivatives, 
  // compute a time step such that no quantity is 
  // changing too quickly. The pointers record which
  // variable and which cell is limiting the time step.
  double ComputeTimeStep(const double z,int*,int*);

  // Given a time step, state variables, and their time 
  // derivatives, do a forward Euler step
  void UpdateStateVars(const double dt, 
		       const double redshift,double **);

  // Using parameters which specify the initial conditions, 
  // fill in the initial values for the state variables
  // and fixed quantities (x, beta, u,... )
  // This method assumes constant ratios sigst/sig, colst/col 
  // as functions of radius (i.e. constant fg and constant phi)
  void Initialize(double phi0,double fg0,
		  unsigned int NActive,unsigned int NPassive,
		  double velCurveTurnoverRadius);

  // Similar to the above, except put in an exponential scale 
  // length and constant velocity dispersion for the stars
  void Initialize(double Z_Init, double fcool, double fg0,
		  double sigst0, double Mh0,
		  unsigned int NActive, unsigned int NPassive,
		  double BulgeRadius, double stScaleLength);

  // Is one of the current stellar populations 'currently forming'
  //, i.e. since stars are binned by age, is the age of stars 
  // formed at the current time step such that they will belong 
  // to an existing stellar population, or does a new one need to 
  // be created? In the latter case, create a new population.  The 
  // switch 'active' lets the method know whether we're dealing 
  // with the active or the passive stellar population. In the 
  // former case, the new stellar population, if it's created, 
  // will be taken from an older population. In the latter case, 
  // the new stellar population will be just the stars formed in 
  // this time step.
  bool CheckStellarPops(const double dt, const double redshift,
			std::vector<StellarPop>&,
			unsigned int numberOfAgeBins,
			bool active);

  // Fill tauvec with the torque and its first derivative, i.e. 
  // solve the torque equation given an inner and outer boundary 
  // condition. These are such that tau(x=xmin)=IBC and tau'(x=1)=OBC
  void ComputeTorques(double **tauvec, 
		      const double IBC, 
		      const double OBC);

  // Do the same thing as ComputeTorques, except instead of 
  // solving the torque equation which enforces dQ/dt=0, this 
  // equation just computes the torques for a given alpha viscosity.
  // Also, if the torque as computed here is larger than that 
  // computed by ComputeTorques, replace tau and tau' with 
  // the values computed here. The idea is that if GI shuts down
  // and MRI still operates, let the gas be transported by MRI.
  void ComputeMRItorque(double **tauvec, const double alphaMRI);

  // Store enough information to initialize a simulation in 
  // the Initializer object in.
  void store(Initializer& in);

  // Given the information stored in the Initializer object, 
  // initialize the simulation
  void Initialize(Initializer& in, bool fixedPhi0);
 private:
  std::vector<double> col,sig; // gas column density and velocity dispersion
  std::vector<double> dQdS,dQds; // partial derivatives dQ/dS and dQ/ds
  std::vector<double> dQdSerr,dQdserr; //.. and their errors
  std::vector<double> dcoldt,dsigdt,dZDiskdt,colSFR; // time derivatives

  // store the cells where we have turned off forcing in the
  // torque equation.
  std::vector<int> keepTorqueOff;

  // store (dS/dt)*dt from artificial diffusion (not currently used)
  std::vector<double> diffused_dcoldt;

  // store the metal fluxes for metal diffusion, defined as the
  // net flux of metal mass from bin n+1 to bin n
//  std::vector<double> ZFlux;

  // a vector of stellar populations which affect the gravitational 
  // dynamics of the disk
  std::vector<StellarPop> spsActive; 
  // stellar populations which evolve passively, i.e. do not 
  // affect the gas dynamics.
  std::vector<StellarPop> spsPassive; 
  
  std::vector<double> ZDisk; // metallicity at each cell
  unsigned int nx; // number of cells
  
  //  Dimensionless values of:
  std::vector<double> 
    x,     // position of each cell
    beta,  // power law index of rotation curve
    uu,    // local circular velocity
    yy,    // inward velocity of stars
    betap, //d(beta)/dx 
    psi;

  std::vector<double>
    h2,h1,h0,H; // coefficients of the torque equation

  // inner truncation radius, logarithmic width of a cell,
  // and factor by which to reduce the timestep TOL
  const double XMIN, dlnx,TOL; 

  Dimensions& dim; // store dimensional quantities
  Cosmology& cos; // store cosmology
  
  const bool migratePassive;

  // Physical parameters of the simulation:
  const double 
    // timescale in local orbital times for Q_* to approach Q_lim
    tauHeat, 

    // dissipation rate parameter - how much of the gas's 
    // non-thermal kinetic energy Sigma sigma_nt^2 is 
    // dissipated in a scale height crossing time? 
    // eta=1.5 corresponds to all KE per crossing time
    ETA, 
    sigth, // dimensionless thermal velocity dispersion (set by T_gas)
    EPS_ff, // star formation efficiency per free fall time
    Qlim, // Q below which transient spirals heat the stellar disk
    minsigstF, // the minimum sig_st = minsigstF * sigth.
    thickness, // correction to Q owing to finite thickness

    // ratio of rate at which mass is ejected from a given cell 
    // to the star formation rate in that cell
    MassLoadingFactor,

    // the value of Q which we fix at the beginning of the simulation.
    // Reasonable values are between 1 and 2 or so, depending on the
    // thickness correction used (if any). See Elmegreen (2011) for
    // details on why one should probably choose a number >1.
    fixedQ; 

  // properties of the "bulge", i.e. the region inside the inner
  // truncation radius of the disk
  double ZBulge,MBulge;

  // parameters controlling the instantaneous recycling approx.
  double yREC, RfREC, zetaREC; 

  bool analyticQ;  // use analytic (Romeo-Wiegert 2011) or numerical (Rafikov 2001) Q
  std::vector<double> CumulativeSF; // Total number of cells ever formed in each cell
  std::vector<double> CumulativeTorqueErr2; // Cumulative error in the torque equation..
  std::vector<double> CumulativeTorqueErr;  // measured in various ways
  std::vector<double> d2taudx2;             // tau''
  std::vector<double> CuStarsOut, CuGasOut; // Cumulative stars or gas which leave a cell due to migration
  double initialStellarMass, initialGasMass; // initial total masses in stars and gas
  double cumulativeMassAccreted;             // total mass accreted at the outer edge fo the disk
  double cumulativeStarFormationMass;        // total mass of stars ever formed
  double cumulativeGasMassThroughIB;         // total mass of gas which ever flows through the inner boundary
  double cumulativeStellarMassThroughIB;     // total mass of stars which flow through the inner boundary
  double CumulativeTorque; // integral of tau(x=1) dt.

  double kappaMetals;
};
