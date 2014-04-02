#include <string>
#include <vector>
#include <math.h>

#include "FixedMesh.h"
#include "Debug.h"
#include "AccretionProfile.h"
#include "Max.h"
#include "Str.h"
#include "Errors.h"

AccretionProfile::AccretionProfile(FixedMesh & theMesh,
        int wp,
        double a,
        Debug& db,
        double rs,
        double theWidth):
    mesh(theMesh),
    whichProfile(wp),
    alpha(a),
    dbg(db),
    nx(theMesh.nx()),
    profile(std::vector<double>(theMesh.nx()+1,-1.0)),
    sigma(rs*theWidth),
    radialScale(rs),
    currentRadialScale(rs),
    normalization(1.0),
    width(theWidth)
{
    return;
}

void AccretionProfile::UpdateProfile(double r200)
{
    currentRadialScale = radialScale*r200;
    if(alpha != 0.0 || profile[1]<0.0) {
        // The following few lines are only relevant for Gaussian profiles,
        // although since they only need to be evaluated once (not for every cell),
        // they're up here.
        sigma = width*currentRadialScale;
        // normalization = Integral (2 pi r * NormalPDF dr) from 0 to inf.
        // Here NormalPDF = 1/sqrt(2pi sigma^2) exp(-(r-rs)^2/(2 sigma^2))
//        normalization = currentRadialScale*M_PI + exp(-currentRadialScale*currentRadialScale / (2.0*sigma*sigma)) * sqrt(2.0*M_PI)*sigma + currentRadialScale*M_PI*erf(currentRadialScale/(sqrt(2.0)*sigma));
        normalization = exp(-currentRadialScale*currentRadialScale/(2.0*sigma*sigma))*sqrt(2.0*M_PI)*sigma + currentRadialScale*M_PI*(1.0 + erf(currentRadialScale / (sqrt(2.0)*sigma)));

        double currentFracOuter = fOuter();

	if(dbg.opt(1))
            currentFracOuter = 0.0;

        for(unsigned int n=1; n<=nx; ++n) {
            double nD = ((double) n);
            double xlo = mesh.x(nD-0.5);
            double xhi = mesh.x(nD+0.5);

            if(whichProfile == 0) { // exponential
                profile[n] = (( (currentRadialScale +xlo)*exp(-xlo/currentRadialScale)
                            - (currentRadialScale+xhi)*exp(-xhi/currentRadialScale))
                        * 2.0/(currentRadialScale*(xhi-xlo)*(xhi+xlo)))*
                     1.0 / (1.0 - currentFracOuter);
            }

            if(whichProfile == 1) { // flat
                if(xhi < currentRadialScale) { // this cell is contained within the accretion radius.
                    profile[n] = 2.0 /(currentRadialScale*currentRadialScale) * 1.0/(1.0-currentFracOuter);
                }
                else if( xlo > currentRadialScale) { // this cell is completely oustide the accretion radius.
                    profile[n] = 0.0;
                }
                else { // this cell contains the accretion radius.
                    profile[n] = (2.0 / (currentRadialScale*currentRadialScale)) *
                        (currentRadialScale*currentRadialScale - xlo*xlo)/(xhi*xhi-xlo*xlo)*
                        1.0/(1.0-currentFracOuter);
                }

            }

	    double fraction;
            if(whichProfile == 2) { // gaussian
                // this is the fraction of the accreted mass which will be going into this cell,
                // namely Integral (2 pi r * NormalPDF dr) from xlo to xhi, all over normalization.
                fraction = max(0.0,  ((exp(-(currentRadialScale-xlo)*(currentRadialScale-xlo)/(2.0*sigma*sigma)) - exp(-(currentRadialScale-xhi)*(currentRadialScale-xhi)/(2.0*sigma*sigma)))*sqrt(2.0*M_PI)*sigma + currentRadialScale*M_PI*(erf((currentRadialScale-xlo)/(sigma*sqrt(2.0))) - erf((currentRadialScale-xhi)/(sigma*sqrt(2.0)))))) / normalization;
                // the fraction of accreted mass which nominally falls outside the disk boundary.
                // To properly conserve mass, it is added back into the computational domain
                // in proportion to the already-computed accretion there.

                profile[n] = 2.0/((xlo+xhi)*(xhi-xlo)) * // ~ 1/area of this cell
                    fraction * 1.0/(1.0-currentFracOuter);
            }
            if(profile[n]<0.0) errormsg("Negative accretion profile! This could mean you've selected an invalid/unsupported accretion profile "+str(n)+" "+str(profile[n])+" "+str(fraction)+" "+str(currentFracOuter) +" "+ str(currentRadialScale)+" "+str(normalization)+"   "+str(radialScale)+" "+str(r200)+" "+str(alpha));
        }
    }

}

// second two arguments only relevant for Gaussian column density profile.
double AccretionProfile::fOuter()
{
    double xb = mesh.x(0.5 + ((double) nx));
    double fracOuter=-1.0;
    if(whichProfile == 0) {
        fracOuter = (1.0 + xb/currentRadialScale) * exp(-xb/currentRadialScale);
    }
    if(whichProfile == 1) {
        fracOuter = max((currentRadialScale*currentRadialScale-xb*xb)/(currentRadialScale*currentRadialScale), 0.0);
    }
    if(whichProfile == 2) {
        //fracOuter = (currentRadialScale*M_PI + exp(-(currentRadialScale-xb)*(currentRadialScale-xb)/(2.0*sigma*sigma)) * sqrt(2.0*M_PI)*sigma + currentRadialScale*M_PI*erf((currentRadialScale-xb)/(sqrt(2.0)*sigma)))/ normalization;
        fracOuter = (exp(-(currentRadialScale-xb)*(currentRadialScale-xb)/(2.0*sigma*sigma))*sqrt(2.*M_PI)*sigma + currentRadialScale*M_PI*(1. + erf((currentRadialScale - xb)/(sqrt(2.0)*sigma)))) / normalization;
    }
    if(fracOuter > 1.0 || fracOuter <0.0) errormsg("A nonphysical fraction of the accretion has occurred outside the computational domain");
    return fracOuter;
}
double AccretionProfile::fInner()
{
    double xinner = mesh.x(0.5);
    double fracInner = -1.0;
    if(whichProfile == 0)
        fracInner = 1.0 - exp(-xinner/currentRadialScale) * (1.0 + xinner/currentRadialScale);
    if(whichProfile == 1)
        fracInner = xinner*xinner / (currentRadialScale*currentRadialScale);
    if(whichProfile == 2) {
        fracInner = ((exp(-currentRadialScale*currentRadialScale / (2.0*sigma*sigma)) - exp(-(currentRadialScale - xinner)*(currentRadialScale - xinner)/(2.0*sigma*sigma))) * sqrt(2.0*M_PI)*sigma + currentRadialScale*M_PI*(erf(currentRadialScale/(sqrt(2.0)*sigma)) - erf((currentRadialScale - xinner)/(sqrt(2.0)*sigma))))/normalization;
        // The above alone yields negative fracInner's. If this is the result of catastrophic cancellation, 
        // the absolute value of fracInner should be less than about 1.0e-15 * (either of the terms).
        // Let's check whether this is the case:
        if(fabs(fracInner)<1.0e-11 && fracInner<0.0)
            fracInner = 0.0;
    }

   if(fracInner > 1.0 || fracInner < 0.0) errormsg("A nonphysical fraction of the accretion has occurred interior to the innermost radius.   RS cRS sigma "+str(radialScale)+" "+str(currentRadialScale)+" "+str(sigma)+" "+str(normalization)+" "+str(xinner)+" "+str(fracInner)+" "+str(currentRadialScale*M_PI*(erf(currentRadialScale/(sqrt(2.0)*sigma)) - erf((currentRadialScale - xinner)/(sqrt(2.0)*sigma)))/normalization)+" "+str(exp(-currentRadialScale*currentRadialScale/(2.0)))); 
   return fracInner;
}




