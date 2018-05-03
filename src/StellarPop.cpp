#include "StellarPop.h"
#include "Cosmology.h"
#include "DiskContents.h"
#include "DiskUtils.h"
#include "FixedMesh.h"
#include "Debug.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <gsl/gsl_integration.h>

// typical constructor
StellarPop::StellarPop(FixedMesh & m) :
  spcol(std::vector<double>(m.nx()+1,0.)), 
  spsigR(std::vector<double>(m.nx()+1,0.)),
  spsigZ(std::vector<double>(m.nx()+1,0.)),
  dSigRdr(std::vector<double>(m.nx()+1,0.)),
  dSigZdr(std::vector<double>(m.nx()+1,0.)),
  dColdr(std::vector<double>(m.nx()+1,0.)),
  spZFe(std::vector<double>(m.nx()+1,0.)), 
  spZO(std::vector<double>(m.nx()+1,0.)), 
  //spZV(std::vector<double>(m.nx()+1,0.)),
  dQdS(std::vector<double>(m.nx()+1,0.)),
  dQdsR(std::vector<double>(m.nx()+1,0.)), 
  dQdsZ(std::vector<double>(m.nx()+1,0.)),
  dQdSerr(std::vector<double>(m.nx()+1,0.)),
  dQdserr(std::vector<double>(m.nx()+1,0.)),
  dcoldtREC(std::vector<double>(m.nx()+1,0.)),
  dcoldtIA(std::vector<double>(m.nx()+1,0.)),
  ageAtz0(-1.0),startingAge(-1.0),endingAge(-1.0),
  mesh(m),
  allocated(false),
  interpTypeK(gsl_interp2d_bilinear),
  interpTypeL(gsl_interp2d_bilinear)
{ }

StellarPop::~StellarPop()
{
  if(allocated) {
    gsl_spline2d_free( splineK );
    gsl_spline2d_free( splineL );
    gsl_interp_accel_free( accelKx );
    gsl_interp_accel_free( accelKy );
    gsl_interp_accel_free( accelLx );
    gsl_interp_accel_free( accelLy );
    free(zaK);
    free(zaL);
  }
}
void StellarPop::InitializeGSLObjs()
{
    allocated = true; // keep track that we have actually initialized the GSL objects that we're constructing here.
    const size_t nx = 300;
    const size_t ny = 300;
    zaK = (double*) malloc(nx*ny*sizeof(double));
    zaL = (double*) malloc(nx*ny*sizeof(double));
    splineK = gsl_spline2d_alloc(interpTypeK, nx, ny);
    splineL = gsl_spline2d_alloc(interpTypeL, nx, ny);
    // is there a better way? Maybe!?
    const double xa[] = { 0.001, 0.0076856187291, 0.0143712374582, 0.0210568561873, 0.0277424749164, 0.0344280936455, 0.0411137123746, 0.0477993311037, 0.0544849498328, 0.0611705685619, 0.067856187291, 0.0745418060201, 0.0812274247492, 0.0879130434783, 0.0945986622074, 0.101284280936, 0.107969899666, 0.114655518395, 0.121341137124, 0.128026755853, 0.134712374582, 0.141397993311, 0.14808361204, 0.154769230769, 0.161454849498, 0.168140468227, 0.174826086957, 0.181511705686, 0.188197324415, 0.194882943144, 0.201568561873, 0.208254180602, 0.214939799331, 0.22162541806, 0.228311036789, 0.234996655518, 0.241682274247, 0.248367892977, 0.255053511706, 0.261739130435, 0.268424749164, 0.275110367893, 0.281795986622, 0.288481605351, 0.29516722408, 0.301852842809, 0.308538461538, 0.315224080268, 0.321909698997, 0.328595317726, 0.335280936455, 0.341966555184, 0.348652173913, 0.355337792642, 0.362023411371, 0.3687090301, 0.375394648829, 0.382080267559, 0.388765886288, 0.395451505017, 0.402137123746, 0.408822742475, 0.415508361204, 0.422193979933, 0.428879598662, 0.435565217391, 0.44225083612, 0.448936454849, 0.455622073579, 0.462307692308, 0.468993311037, 0.475678929766, 0.482364548495, 0.489050167224, 0.495735785953, 0.502421404682, 0.509107023411, 0.51579264214, 0.52247826087, 0.529163879599, 0.535849498328, 0.542535117057, 0.549220735786, 0.555906354515, 0.562591973244, 0.569277591973, 0.575963210702, 0.582648829431, 0.589334448161, 0.59602006689, 0.602705685619, 0.609391304348, 0.616076923077, 0.622762541806, 0.629448160535, 0.636133779264, 0.642819397993, 0.649505016722, 0.656190635452, 0.662876254181, 0.66956187291, 0.676247491639, 0.682933110368, 0.689618729097, 0.696304347826, 0.702989966555, 0.709675585284, 0.716361204013, 0.723046822742, 0.729732441472, 0.736418060201, 0.74310367893, 0.749789297659, 0.756474916388, 0.763160535117, 0.769846153846, 0.776531772575, 0.783217391304, 0.789903010033, 0.796588628763, 0.803274247492, 0.809959866221, 0.81664548495, 0.823331103679, 0.830016722408, 0.836702341137, 0.843387959866, 0.850073578595, 0.856759197324, 0.863444816054, 0.870130434783, 0.876816053512, 0.883501672241, 0.89018729097, 0.896872909699, 0.903558528428, 0.910244147157, 0.916929765886, 0.923615384615, 0.930301003344, 0.936986622074, 0.943672240803, 0.950357859532, 0.957043478261, 0.96372909699, 0.970414715719, 0.977100334448, 0.983785953177, 0.990471571906, 0.997157190635, 1.00384280936, 1.01052842809, 1.01721404682, 1.02389966555, 1.03058528428, 1.03727090301, 1.04395652174, 1.05064214047, 1.0573277592, 1.06401337793, 1.07069899666, 1.07738461538, 1.08407023411, 1.09075585284, 1.09744147157, 1.1041270903, 1.11081270903, 1.11749832776, 1.12418394649, 1.13086956522, 1.13755518395, 1.14424080268, 1.1509264214, 1.15761204013, 1.16429765886, 1.17098327759, 1.17766889632, 1.18435451505, 1.19104013378, 1.19772575251, 1.20441137124, 1.21109698997, 1.2177826087, 1.22446822742, 1.23115384615, 1.23783946488, 1.24452508361, 1.25121070234, 1.25789632107, 1.2645819398, 1.27126755853, 1.27795317726, 1.28463879599, 1.29132441472, 1.29801003344, 1.30469565217, 1.3113812709, 1.31806688963, 1.32475250836, 1.33143812709, 1.33812374582, 1.34480936455, 1.35149498328, 1.35818060201, 1.36486622074, 1.37155183946, 1.37823745819, 1.38492307692, 1.39160869565, 1.39829431438, 1.40497993311, 1.41166555184, 1.41835117057, 1.4250367893, 1.43172240803, 1.43840802676, 1.44509364548, 1.45177926421, 1.45846488294, 1.46515050167, 1.4718361204, 1.47852173913, 1.48520735786, 1.49189297659, 1.49857859532, 1.50526421405, 1.51194983278, 1.51863545151, 1.52532107023, 1.53200668896, 1.53869230769, 1.54537792642, 1.55206354515, 1.55874916388, 1.56543478261, 1.57212040134, 1.57880602007, 1.5854916388, 1.59217725753, 1.59886287625, 1.60554849498, 1.61223411371, 1.61891973244, 1.62560535117, 1.6322909699, 1.63897658863, 1.64566220736, 1.65234782609, 1.65903344482, 1.66571906355, 1.67240468227, 1.679090301, 1.68577591973, 1.69246153846, 1.69914715719, 1.70583277592, 1.71251839465, 1.71920401338, 1.72588963211, 1.73257525084, 1.73926086957, 1.74594648829, 1.75263210702, 1.75931772575, 1.76600334448, 1.77268896321, 1.77937458194, 1.78606020067, 1.7927458194, 1.79943143813, 1.80611705686, 1.81280267559, 1.81948829431, 1.82617391304, 1.83285953177, 1.8395451505, 1.84623076923, 1.85291638796, 1.85960200669, 1.86628762542, 1.87297324415, 1.87965886288, 1.88634448161, 1.89303010033, 1.89971571906, 1.90640133779, 1.91308695652, 1.91977257525, 1.92645819398, 1.93314381271, 1.93982943144, 1.94651505017, 1.9532006689, 1.95988628763, 1.96657190635, 1.97325752508, 1.97994314381, 1.98662876254, 1.99331438127, 2.0 };     
    const double ya[] = { 0.5, 0.506688963211, 0.513377926421, 0.520066889632, 0.526755852843, 0.533444816054, 0.540133779264, 0.546822742475, 0.553511705686, 0.560200668896, 0.566889632107, 0.573578595318, 0.580267558528, 0.586956521739, 0.59364548495, 0.600334448161, 0.607023411371, 0.613712374582, 0.620401337793, 0.627090301003, 0.633779264214, 0.640468227425, 0.647157190635, 0.653846153846, 0.660535117057, 0.667224080268, 0.673913043478, 0.680602006689, 0.6872909699, 0.69397993311, 0.700668896321, 0.707357859532, 0.714046822742, 0.720735785953, 0.727424749164, 0.734113712375, 0.740802675585, 0.747491638796, 0.754180602007, 0.760869565217, 0.767558528428, 0.774247491639, 0.780936454849, 0.78762541806, 0.794314381271, 0.801003344482, 0.807692307692, 0.814381270903, 0.821070234114, 0.827759197324, 0.834448160535, 0.841137123746, 0.847826086957, 0.854515050167, 0.861204013378, 0.867892976589, 0.874581939799, 0.88127090301, 0.887959866221, 0.894648829431, 0.901337792642, 0.908026755853, 0.914715719064, 0.921404682274, 0.928093645485, 0.934782608696, 0.941471571906, 0.948160535117, 0.954849498328, 0.961538461538, 0.968227424749, 0.97491638796, 0.981605351171, 0.988294314381, 0.994983277592, 1.0016722408, 1.00836120401, 1.01505016722, 1.02173913043, 1.02842809365, 1.03511705686, 1.04180602007, 1.04849498328, 1.05518394649, 1.0618729097, 1.06856187291, 1.07525083612, 1.08193979933, 1.08862876254, 1.09531772575, 1.10200668896, 1.10869565217, 1.11538461538, 1.1220735786, 1.12876254181, 1.13545150502, 1.14214046823, 1.14882943144, 1.15551839465, 1.16220735786, 1.16889632107, 1.17558528428, 1.18227424749, 1.1889632107, 1.19565217391, 1.20234113712, 1.20903010033, 1.21571906355, 1.22240802676, 1.22909698997, 1.23578595318, 1.24247491639, 1.2491638796, 1.25585284281, 1.26254180602, 1.26923076923, 1.27591973244, 1.28260869565, 1.28929765886, 1.29598662207, 1.30267558528, 1.30936454849, 1.31605351171, 1.32274247492, 1.32943143813, 1.33612040134, 1.34280936455, 1.34949832776, 1.35618729097, 1.36287625418, 1.36956521739, 1.3762541806, 1.38294314381, 1.38963210702, 1.39632107023, 1.40301003344, 1.40969899666, 1.41638795987, 1.42307692308, 1.42976588629, 1.4364548495, 1.44314381271, 1.44983277592, 1.45652173913, 1.46321070234, 1.46989966555, 1.47658862876, 1.48327759197, 1.48996655518, 1.49665551839, 1.50334448161, 1.51003344482, 1.51672240803, 1.52341137124, 1.53010033445, 1.53678929766, 1.54347826087, 1.55016722408, 1.55685618729, 1.5635451505, 1.57023411371, 1.57692307692, 1.58361204013, 1.59030100334, 1.59698996656, 1.60367892977, 1.61036789298, 1.61705685619, 1.6237458194, 1.63043478261, 1.63712374582, 1.64381270903, 1.65050167224, 1.65719063545, 1.66387959866, 1.67056856187, 1.67725752508, 1.68394648829, 1.69063545151, 1.69732441472, 1.70401337793, 1.71070234114, 1.71739130435, 1.72408026756, 1.73076923077, 1.73745819398, 1.74414715719, 1.7508361204, 1.75752508361, 1.76421404682, 1.77090301003, 1.77759197324, 1.78428093645, 1.79096989967, 1.79765886288, 1.80434782609, 1.8110367893, 1.81772575251, 1.82441471572, 1.83110367893, 1.83779264214, 1.84448160535, 1.85117056856, 1.85785953177, 1.86454849498, 1.87123745819, 1.8779264214, 1.88461538462, 1.89130434783, 1.89799331104, 1.90468227425, 1.91137123746, 1.91806020067, 1.92474916388, 1.93143812709, 1.9381270903, 1.94481605351, 1.95150501672, 1.95819397993, 1.96488294314, 1.97157190635, 1.97826086957, 1.98494983278, 1.99163879599, 1.9983277592, 2.00501672241, 2.01170568562, 2.01839464883, 2.02508361204, 2.03177257525, 2.03846153846, 2.04515050167, 2.05183946488, 2.05852842809, 2.0652173913, 2.07190635452, 2.07859531773, 2.08528428094, 2.09197324415, 2.09866220736, 2.10535117057, 2.11204013378, 2.11872909699, 2.1254180602, 2.13210702341, 2.13879598662, 2.14548494983, 2.15217391304, 2.15886287625, 2.16555183946, 2.17224080268, 2.17892976589, 2.1856187291, 2.19230769231, 2.19899665552, 2.20568561873, 2.21237458194, 2.21906354515, 2.22575250836, 2.23244147157, 2.23913043478, 2.24581939799, 2.2525083612, 2.25919732441, 2.26588628763, 2.27257525084, 2.27926421405, 2.28595317726, 2.29264214047, 2.29933110368, 2.30602006689, 2.3127090301, 2.31939799331, 2.32608695652, 2.33277591973, 2.33946488294, 2.34615384615, 2.35284280936, 2.35953177258, 2.36622073579, 2.372909699, 2.37959866221, 2.38628762542, 2.39297658863, 2.39966555184, 2.40635451505, 2.41304347826, 2.41973244147, 2.42642140468, 2.43311036789, 2.4397993311, 2.44648829431, 2.45317725753, 2.45986622074, 2.46655518395, 2.47324414716, 2.47993311037, 2.48662207358, 2.49331103679, 2.5};

    accelKx = gsl_interp_accel_alloc();
    accelKy = gsl_interp_accel_alloc();
    accelLx = gsl_interp_accel_alloc();
    accelLy = gsl_interp_accel_alloc();

    std::string kfn = "Lacey84_table_K.txt";
    std::string lfn = "Lacey84_table_L.txt";
    std::ifstream kf;
    kf.open("Lacey84_table_K.txt");
    if(!kf) {
        errormsg("Couldn't find Lacey84_table_K.txt");
    }
    int i,j;
    double val;
    while( !kf.eof() ) {
        kf >> i >> j >> val;
	gsl_spline2d_set( splineK, zaK, i, j, val );

    }
    gsl_spline2d_init( splineK, xa, ya, zaK, nx, ny );
    kf.close();

    std::ifstream lf;
    lf.open("Lacey84_table_L.txt");
    if(!lf) {
        errormsg("Couldn't find Lacey84_table_L.txt");
    }
    while( !lf.eof() ) {
        lf >> i >> j >> val;
	gsl_spline2d_set( splineL, zaL, i, j, val );

    }
    gsl_spline2d_init( splineL, xa, ya, zaL, nx, ny );
    lf.close();
}

void StellarPop::ComputeSpatialDerivs()
{
  std::vector<double> & x = mesh.x();
  unsigned int nx = mesh.nx(); 
  for(unsigned int n=1; n<=nx; ++n) {
    dSigRdr[n] = ddx(spsigR,n,x,false,true);
    dSigZdr[n] = ddx(spsigZ,n,x,false,true);
    dColdr[n] = ddx(spcol,n,x,false,true);
  }
}


// Merge sp2 into sp1. sp2 should be unaffected by the procedure.
void StellarPop::MergeStellarPops(const StellarPop& sp2,DiskContents& disk)
{
  for(unsigned int i=1; i<=spcol.size()-1; ++i) {
    if(sp2.spcol[i]>0.0) { // if there are actually stars to add.
      // otherwise do nothing.
      
      // sig3 ^2 = (col1*sig1^2 + col2*sig2^2)/(col1+col2):
      spsigR[i]= sqrt(((spcol[i])*((*this).spsigR[i])*((*this).spsigR[i]) 
                              + (sp2.spcol[i])*(sp2.spsigR[i])*(sp2.spsigR[i]))
                             /((*this).spcol[i]+sp2.spcol[i]));
      spsigZ[i]= sqrt((spcol[i]*spsigZ[i]*spsigZ[i] 
		       + sp2.spcol[i]*sp2.spsigZ[i]*sp2.spsigZ[i])
		      /(spcol[i]+sp2.spcol[i]));
      //    (*this).spZ[i] = ((*this).spZ[i] * (*this).spcol[i] + sp2.spZ[i] * sp2.spcol[i]) / ((*this).spcol[i] + sp2.spcol[i]);
      // explicitly state the moments of each metallicity distribution:
      double wtdAvgFe = (spZFe[i]*spcol[i] + sp2.spZFe[i]*sp2.spcol[i])/(spcol[i] + sp2.spcol[i]);
      double wtdAvgO = (spZO[i]*spcol[i] + sp2.spZO[i]*sp2.spcol[i])/(spcol[i] + sp2.spcol[i]);
      //double avg1 = spZ[i];
      //double avg2 = sp2.spZ[i];
      //double var1 = spZV[i];
      //double var2 = sp2.spZV[i];
      double wt1 = spcol[i];
      double wt2 = sp2.spcol[i];
      
      // Merge the two distributions:
      (*this).spZFe[i] = wtdAvgFe;
      (*this).spZO[i] = wtdAvgO;
 //     (*this).spZV[i] = wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
 //       + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);
 //     spZV[i] = ComputeVariance(spcol[i],0.0, sp2.spcol[i],
//				spZ[i],sp2.spZ[i],spZV[i],sp2.spZV[i]);
      
      if((*this).spsigR[i]!=(*this).spsigR[i] || spsigZ[i]!=spsigZ[i])
        errormsg("Error merging populations:  this spcol,spsig  sp2 spcol,spsig  "+str((*this).spcol[i])+" "+str((*this).spsigR[i])+"  "+str(sp2.spcol[i])+" "+str(sp2.spsigR[i]));
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
    spsigR[n] = sp2.spsigR[n];
    spsigZ[n] = sp2.spsigZ[n];
    spcol[n] = frac*sp2.spcol[n];
    spZFe[n] = sp2.spZFe[n];
    spZO[n] = sp2.spZO[n];
    // spZV[n] = sp2.spZV[n];
    sp2.spcol[n] -= spcol[n];;
  }
  ageAtz0 = sp2.ageAtz0;

}

void StellarPop::ComputeSNIArate(DiskContents& disk, double z)
{
    unsigned int nx=spcol.size()-1;
    double C0=0.046;
    double lambda = 2.76e5 * speryear;
    double RfREC = 0.77;
    double current_age  = (ageAtz0 - disk.GetCos().lbt(z));// in seconds
    double fml = C0 * log(1.0 + current_age/lambda) ;
    if(fml<0.23) {
       fml=0.23;
    }
    double rateAdjustIA = 3.0e-3/0.0013;
    for( unsigned int n=1; n<=nx; ++n) {
        double col_orig_est = spcol[n]/(1-fml); // estimate of the column density of stars formed originally.
	if (current_age>0.1*speryear*1.0e9 && current_age<10.0*speryear*1.0e9)  {
	    dcoldtIA[n] = col_orig_est * rateAdjustIA * 0.0013* 0.14476/(current_age * disk.GetDim().vphiR/(2.0*M_PI*disk.GetDim().Radius));  //// this is the surface density of SNIA explosions per time, in code units.

	    
	    //std::cout << "IA rate: "<<dcoldtIA[n]<<" "<<col_orig_est<<" "<<disk.GetColSFR()[n] << std::endl;
	}

        else
	    dcoldtIA[n] = 0.0;
    }

}
void StellarPop::ComputeRecycling(DiskContents& disk, double z)
{
    unsigned int nx=spcol.size()-1;
    double C0 = 0.046;
    double lambda = 2.76e5 * speryear;
    //double RfREC = disk.GetRfRECinst();
    double RfREC = 0.77; // 1- fractional mass lost at 40 Myr
    // double RfREC = 0.9; //// override the user.
    // At this time step, frac will be ~constant over the whole disk.
    double fracDot = C0/((ageAtz0 - disk.GetCos().lbt(z) + lambda)*disk.GetDim().vphiR/(2.0*M_PI*disk.GetDim().Radius));
    // RfREC is an input parameter, namely the remnant fraction.
    // For an instantaneous recycling approximation, this should equal 1 - return fraction,
    // where return fraction is the fraction of mass returned after some specified time.
    // If we set that time to be very large, e.g. a Hubble time, return fraction = .497
    // using the formula in Leitner11, for a Chabrier IMF.
    // What we want to do here is require that the user-specified return fraction = 1-RfREC
    // is returned immediately, and anything after that is returned slowly.
    // If this stellar pop has not yet reached an age where the returned fraction is greater than this 
    // instantaneous amount, then return nothing, i.e. set frac=0
    // As an example, if the user sets RfREC=0.9, only 10% of star-forming gas is immediately returned
    // and the remaining .397 is returned according to the rate specified by frac above.
    // Whereas if the user specifies a low remnant fraction, an artificially large fraction of the gas
    // is immediately returned, which is similar to just artificially slowing the star formation rate 
    // at fixed gas column density.
    double frac = C0 * log(1.0 + (ageAtz0 - disk.GetCos().lbt(z))/lambda);
    if( frac < 1.0 - RfREC)  {
        fracDot=0.0;
	frac = 1.0 - RfREC;
    }
    for(unsigned int n=1; n<=nx; ++n) {
        dcoldtREC[n] = spcol[n]/(1.0-frac) * fracDot;
    }
}

double tanfac(double s)
{
    if (s>=0) {
        return atan(sqrt(s))/sqrt(s);
    }
    else {
        return atanh(sqrt(-s))/sqrt(-s);
    }
    return 0;
};


double Lintegrand(double chi, void * params)
{
    std::vector<double> * par =  (std::vector<double> *) params; // yikes
    double alpha = par[0][0];
    double beta = par[0][1];
    double sinchi = sin(chi);
    double coschi = cos(chi);
    double betasq = beta*beta;
    double b = 2.0 - (betasq*sinchi*sinchi + (1.0/betasq)*coschi*coschi); 
    double a = sinchi*sinchi + (1.0/betasq)*coschi*coschi;
    double s = a/(alpha*alpha) - 1.0;
    double ret = (-3.0 + (s+3.0)*tanfac(s))/(alpha*s);
    return ret;
};

//double L(double alpha, double beta)
//{
//    int nintervals = 1000; // nothing varying rapidly with chi I think...
//    gsl_integration_workspace *w = gsl_integration_workspace_alloc( nintervals ); 
//    gsl_function F;
//    double result, error;
//
//    std::vector<double> par(2);
//    par[0] = alpha;
//    par[1] = beta;
//    F.function = &Lintegrand;
//    F.params = & par;
//
//    gsl_integration_qags( &F, 0, M_PI/2.0, 1.0e-2, 1.0e-2, nintervals, w, &result, &error );
//
//    gsl_integration_workspace_free(w);
//    return result*2.0/M_PI;
//};


double Kintegrand(double chi, void * params)
{
    std::vector<double> * par =  (std::vector<double> *) params; // yikes
    double alpha = par[0][0];
    double beta = par[0][1];
    //printf("chi, alpha, beta: % .18f % .18f % .18f \n", chi,alpha,beta);
    double sinchi = sin(chi);
    double coschi = cos(chi);
    double betasq = beta*beta;
    double b = 2.0 - (betasq*sinchi*sinchi + (1.0/betasq)*coschi*coschi); 
    double a = sinchi*sinchi + (1.0/betasq)*coschi*coschi;
    double s = a/(alpha*alpha) - 1.0;
    double ret = (3.0 - (b*s+3.0)*tanfac(s))/(alpha*alpha*alpha*s*(s+1.0));
    return ret;
};

//double K(double alpha, double beta)
//{
//    int nintervals = 1000; // nothing varying rapidly with chi I think...
//    gsl_integration_workspace *w = gsl_integration_workspace_alloc( nintervals ); 
//    gsl_function F;
//    double result, error;
//
//    std::vector<double> par(2);
//    par[0] = alpha;
//    par[1] = beta;
//    F.function = &Kintegrand;
//    F.params = & par;
//
//    gsl_integration_qags( &F, 0, M_PI/2.0, 1.0e-2, 1.0e-2, nintervals, w, &result, &error );
//
//    gsl_integration_workspace_free(w);
//    return result*2.0/M_PI;
//};

double StellarPop::L(const double alpha,const double beta)
{
    if(!allocated) {
        (*this).InitializeGSLObjs();
    }
    // ugh hack for now:
    double alphaEval, betaEval;
    if (alpha>1.9999) {
	alphaEval = 1.9999;
    }
    else if( alpha<0.0 ) {
	alphaEval = 0.0;
    }
    else {
        alphaEval = alpha;
    }
    if (beta>2.4999) {
	betaEval = 2.4999;
    }
    else if (beta<0.5) {
        betaEval = 0.5;
    }
    else {
        betaEval = beta;
    }
    const double alphaEvalReal = alphaEval;
    const double betaEvalReal = betaEval;
    //std::cout << "Evaluating L spline for alpha, beta "<<alpha<<" "<<beta<<std::endl;
    return gsl_spline2d_eval( splineL, alphaEvalReal, betaEvalReal, accelLx, accelLy );

}

double StellarPop::K(const double alpha,const double beta)
{
    if(!allocated) {
        (*this).InitializeGSLObjs();
    }
    //std::cout << "Evaluating K spline for alpha, beta "<<alpha<<" "<<beta<<std::endl;
    // ugh hack for now:
    double alphaEval, betaEval;
    if (alpha>1.9999) {
	alphaEval = 1.9999;
    }
    else if( alpha<0.01 ) {
	alphaEval = 0.01;
    }
    else {
        alphaEval = alpha;
    }
    if (beta>2.4999) {
	betaEval = 2.4999;
    }
    else if (beta<0.501) {
        betaEval = 0.501;
    }
    else {
        betaEval = beta;
    }
    const double alphaEvalReal = alphaEval;
    const double betaEvalReal = betaEval;

    return gsl_spline2d_eval( splineK, alphaEvalReal, betaEvalReal, accelKx, accelKy );

}


void StellarPop::CloudHeatStellarPop(double dt, DiskContents& disk, double heatingRate)
{
  std::vector<double> dsigRdt(spcol.size());
  std::vector<double> dsigZdt(spcol.size());
  unsigned int nx=spcol.size()-1;
  std::vector<double>& uu = disk.GetUu();
  std::vector<double>& beta = disk.GetBeta();
  std::vector<double>& x = disk.GetX();
  for(unsigned int n=1; n<=nx; ++n) {
    double hc = disk.hGas(n); // in cm!
    double alpha = spsigZ[n]/spsigR[n];
    double q = disk.hStars(n)/hc ;  
    double gamma = sqrt(2.0/(1.0+beta[n])); // 2\Omega/\kappa. This is \beta in Lacey 1984's notation.
    double Mcloud = 1.0e7 * 2.0e33 /( disk.GetDim().MdotExt0 * (2.0*M_PI*disk.GetDim().Radius / disk.GetDim().vphiR));
    double ac = 10.0 * 3.086e18; // 10 pc in cm
    // this formula is probably `correct' but there's no guarantee Lambda>1 For now neglect the case that causes that.
    //double lnLambda = log(sqrt(2.0)* spsigR[n]*x[n]/(uu[n]*sqrt(2.0*(beta[n]+1.0))) * min( disk.GetDim().Radius/ac, spsigR[n]*spsigR[n] / (2.0*M_PI*Mcloud * disk.GetDim().chi())));
    double lnLambda = log(sqrt(2.0)* spsigR[n]*x[n]/(uu[n]*sqrt(2.0*(beta[n]+1.0))) *disk.GetDim().Radius/ac);
    if (lnLambda<0.3) {
	//errormsg("Unphysical coulomb logarithm in stellar heating");
	lnLambda=0.3;
    }
    //if (lnLambda<1)
	//std::cerr << "WARNING: the coulomb logarithm is kind of small "<<lnLambda<<" "<<spsigR[n]<<std::endl;
    // the factor of 0.1 comes from assuming the surface mass density of "clouds" is 0.1 of the total gas disk column density.
    double C = heatingRate* 2.0 * (4.0*M_PI*M_PI) * disk.GetDim().chi() * disk.GetDim().chi() * 0.1 * disk.GetCol()[n] * Mcloud * lnLambda * disk.GetDim().Radius/hc;
    dsigRdt[n] = 0.5/spsigR[n] * C*K(alpha,gamma)/(spsigR[n]*sqrt(1.0+q*q));
    dsigZdt[n] = 0.5/spsigZ[n] * C*L(alpha,gamma)/(spsigR[n]*sqrt(1.0+q*q));
    //std::cout << "dbgCH q, alpha, Mcloud, lnLambda, C, disgRdtt, dsigZdt, sigR, sigZ, n, dt: " <<q<<" "<<alpha<<" "<<Mcloud<<" "<<lnLambda<<" "<<C<<" "<<dsigRdt[n]<<" "<<dsigZdt[n]<<" "<<spsigR[n]<<" "<<spsigZ[n]<<" "<<n<<" "<<dt<<std::endl;
  }
  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    spsigR[n] += dsigRdt[n]*dt;
    spsigZ[n] += dsigZdt[n]*dt;
  }


}

void StellarPop::MigrateStellarPop(double dt, double ** tauvecStar, DiskContents& disk, std::vector<double>& MdotiPlusHalf)
{
  // A few convenience vectors to store data before updating the state variables.
  // These hold derivatives:
  std::vector<double> dcoldt(spcol.size());
  std::vector<double> dsigRdt(spcol.size());
  std::vector<double> dZOdt(spcol.size());
  std::vector<double> dZFedt(spcol.size());
  //std::vector<double> MdotiPlusHalf(spcol.size());
  // These refer back to mesh variables referenced by the disk object:
  std::vector<double>& uu = disk.GetUu();
  std::vector<double>& beta = disk.GetBeta();
  std::vector<double>& x = disk.GetX();
  // FixedMesh & mesh = disk.GetMesh();
  // These store information about the metal fluxes to calculate the new variance of Z, spZV.
  std::vector<double> incomingMass(spcol.size(),0.0);
  std::vector<double> outgoingMass(spcol.size(),0.0);
  std::vector<double> incomingZFe(spcol.size());
  std::vector<double> incomingZO(spcol.size());
  std::vector<double> cellMass(spcol.size());
  Debug& dbg = disk.GetDbg();
  double spsigZp1,spsigRp1,spcolp1;
  unsigned int nx=spcol.size()-1;

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
  //  tauvecStar[2][n] = (tauvecStar[1][n+1]-tauvecStar[1][n-1])/(mesh.x(n+1)-mesh.x(n-1));
  //  MdotiPlusHalf[n] = -1.0/mesh.u1pbPlusHalf(n) * (tauvecStar[1][n+1]-tauvecStar[1][n])/(mesh.x(n+1)-x[n]);
  }
//  MdotiPlusHalf[0]= -1.0/mesh.u1pbPlusHalf(0) * (tauvecStar[1][1]-tauvecStar[1][0])/(x[1]-mesh.x(0.0));
  
  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    double f = (spcol[n]/disk.activeColSt(n));
    dcoldt[n] = f*(MdotiPlusHalf[n]-MdotiPlusHalf[n-1])*2.0*M_PI/mesh.area(n);
    double MdotCentered = (-tauvecStar[2][n]*f
		          /(uu[n]*(1+beta[n]))); 
    if(n<nx) {
      spsigZp1=spsigZ[n+1];
      spsigRp1=spsigR[n+1];
      spcolp1 =spcol[n+1];
    }
    else {
      spsigZp1=spsigZ[nx];
      spsigRp1=spsigR[nx];
      spcolp1 =spcol[nx];
    }

//    dsigRdt[n] = MdotCentered* 
//      (1.0/(x[n]*spcol[n]*(spsigR[n] + spsigZ[n]))) *
//      (2.0*spsigZ[n]* dSigZdr[n] //ddx(spsigZ,n,x,false)
//       + 3.0* spsigR[n]* dSigRdr[n] //ddx(spsigR,n,x,false) 
//       + spsigR[n]*spsigR[n]/spcol[n]* dColdr[n] //ddx(spcol,n,x,false)
//       + (spsigR[n]*spsigR[n] - spsigZ[n]*spsigZ[n])/x[n]);
    //dsigRdt[n] = 1.0/(x[n]*spcol[n]*(spsigR[n]+spsigZ[n])) * ((beta[n]-1.)*uu[n]*tauvecStar[1][n]/(x[n]*x[n]) + (2.0*spsigR[n]*dSigRdr[n] - uu[n]*uu[n]*(1.+beta[n]))/x[n] * MdotCentered + spsigR[n]*spsigR[n]*(MdotiPlusHalf[n]-MdotiPlusHalf[n])/mesh.dx(n));
    //
    //Compute some derivatives. Check that it makes sense to do so, i.e. there's material there.
    if(spcol[n] > 0.0) {
        dsigRdt[n] =  1.0/(x[n]*spcol[n]*(spsigR[n]+spsigZ[n])) * ((beta[n]-1.)*uu[n]*f*tauvecStar[1][n]/(x[n]*x[n]) + (3.0*spsigR[n]*dSigRdr[n] + 2.0*spsigZ[n]*dSigZdr[n]) *(-tauvecStar[2][n]*f/(uu[n]*(1.+beta[n]))) + spsigR[n]*spsigR[n]*f*(MdotiPlusHalf[n]-MdotiPlusHalf[n])*2.0*M_PI/mesh.area(n));

        dZOdt[n] =  MdotCentered*ddx(spZO,n,x,true,true)/(x[n]*spcol[n]);
        dZFedt[n] =  MdotCentered*ddx(spZFe,n,x,true,true)/(x[n]*spcol[n]);
    }
    else { 
        dsigRdt[n] = 0.0;
        dZOdt[n] = 0.0;
        dZFedt[n] = 0.0;
    }


    if(dsigRdt[n]!=dsigRdt[n] || dZOdt[n]!=dZOdt[n] ||  dZFedt[n]!=dZFedt[n] ||  dcoldt[n]!=dcoldt[n]) {
        errormsg("Something has gone wrong in calculating time derivs (MigrateStellarPop, StellarPop.cpp)");
    }

    // Now we proceed to do what looks like a ridiculous amount of work to compute spZV.
    cellMass[n] = spcol[n]*mesh.area(n)/(2.0*M_PI);
    bool fromRight = MdotiPlusHalf[n]>0.0;
    bool fromLeft = MdotiPlusHalf[n-1]<0.0;
    unsigned nx = spcol.size()-1;
    double spZFep1, spZOp1, spZFem1, spZOm1;
    if(n<nx) { 
      spZFep1=spZFe[n+1];
      spZOp1=spZO[n+1];
    }
    else { // these values shouldn't matter in theory.
      spZFep1 = 0.0;
      spZOp1 = 0.0;
    }
    if(n>1) {
      spZOm1=spZO[n-1];
      spZFem1=spZFe[n-1];
      //spZVm1=spZV[n-1];
    }
    else { // again, these values should not affect the calculation.
      spZOm1= 0.0;
      spZFem1= 0.0;
      //spZVm1= 0.0;
    }
    if(fromRight) {
      incomingMass[n] += MdotiPlusHalf[n]*dt;
      if( !fromLeft) {
        incomingZO[n] = spZO[n+1];
        incomingZFe[n] = spZFe[n+1];
        //incomingZV[n] = spZV[n+1];
      }
    }
    if(!fromRight) {
      outgoingMass[n] -= MdotiPlusHalf[n]*dt; 
    }
    if(fromLeft) {
      incomingMass[n] -= MdotiPlusHalf[n-1]*dt;
      if(n>1 && !fromRight) {
        incomingZO[n] = spZO[n-1];
        incomingZFe[n] = spZFe[n-1];
        //incomingZV[n] = spZV[n-1];
      }
    }
    if(!fromLeft) {
      outgoingMass[n] += MdotiPlusHalf[n-1]*dt;
    }
    if(fromLeft && fromRight) {
      incomingZFe[n] = (MdotiPlusHalf[n]*dt*spZFep1 - MdotiPlusHalf[n-1]*dt*spZFem1)/incomingMass[n];
      incomingZO[n] = (MdotiPlusHalf[n]*dt*spZOp1 - MdotiPlusHalf[n-1]*dt*spZOm1)/incomingMass[n];
      //incomingZV[n] = ComputeVariance(MdotiPlusHalf[n]*dt,0.0,-MdotiPlusHalf[n-1]*dt,
      //                                spZp1, spZm1, spZVp1, spZVm1);
    }

  }

  for(unsigned int n=1; n<=spcol.size()-1; ++n) {
    spcol[n] += dcoldt[n]*dt;
    spsigR[n] += dsigRdt[n]*dt;
    spsigZ[n] += .5*dsigRdt[n]*dt;
    spZO[n] += dZOdt[n]*dt;
    spZFe[n] += dZFedt[n]*dt;
    //spZV[n] = ComputeVariance(cellMass[n],outgoingMass[n],incomingMass[n],
    //                          spZ[n],incomingZ[n],spZV[n],incomingZV[n]);
    if(spcol[n]<0 || spsigR[n]<0 || spsigZ[n]<0 || spZO[n]<0  || spZFe[n]<0 ||
       spcol[n]!=spcol[n] || spsigR[n]!=spsigR[n] || spsigZ[n]!=spsigZ[n] ||
       spZFe[n]!=spZFe[n] || spZO[n]!=spZO[n] )
       errormsg("Migrating the populations has produced nonsensical results! "+str(n)+" "+str(spcol[n]) +" "+str(spsigR[n])+" "+str(spsigZ[n])+" "+str(spZO[n])+" "+str(spZFe[n])+" "+str(dZOdt[n])+" "+str(dZFedt[n]));
  }
}




double ComputeVariance(double cellMass, double outgoingMassINPUT, double incomingMass, 
                       double Z, double incomingZ, double ZV, double incomingZV)
{
    double outgoingMass = outgoingMassINPUT;
    double wt1 = cellMass - outgoingMass;
    if(wt1 < 0.0) wt1=0.0;
    double wt2 = incomingMass;
    double avg1 = Z;
    double avg2 = incomingZ;
    double var1 = ZV;
    double var2 = incomingZV;
    double wtdAvg = (wt1*avg1 + wt2*avg2)/(wt1+wt2);
    double val;
    val =  wt1/(wt1+wt2) * (avg1*avg1 + var1 - 2*wtdAvg*avg1 + wtdAvg*wtdAvg)
         + wt2/(wt1+wt2) * (avg2*avg2 + var2 - 2*wtdAvg*avg2 + wtdAvg*wtdAvg);

    if(val >= 0.0) return val;
    else if(val <-1.0e-10)
      errormsg("Something has gone wrong in computing the variance of metallicity. We were given the following: cellMass, outgoingMass, incomingMass, Z, incomingZ, ZV, incoming ZV:  "+str(cellMass)+" "+str(outgoingMass)+" "+str(incomingMass)+" "+str(Z)+" "+str(incomingZ)+" "+str(ZV)+" "+str(incomingZV));
    else return 0.0;


    errormsg("Something strange has happened in ComputeVariance in StellarPop.cpp");
    return 0; // never get here.

}


