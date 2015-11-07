import numpy as np
import pdb
import struct
import glob
import os
import random
import copy
import verticalProfile
import halo
from behroozi import *
from scipy.interpolate import interp1d
from cosmolopy import *
from math import pi,log,sinh,sin,cos,sqrt,log10
#from bitstring import Bits
import pyfits
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
import matplotlib.cm as cmx
from makeGidgetMovies import makeMovies
import subprocess

speryear = 31557600.0
cmperkpc = 3.08567758e21
cmperpc = 3.08567758e18
pcperkpc = 1.0e3
kmperkpc = 1.0e-5*cmperkpc
gpermsun = 1.9891e33
Gcgs = 6.67384e-8
kB = 1.3806488e-16
gperH = 1.008*1.66053892e-24

#RedBlueCM = cm = plt.get_cmap('RdBu')

#cm = plt.get_cmap('gist_rainbow')
cm = plt.get_cmap('gnuplot_r')
cmPer = plt.get_cmap('spring')

outward = [0,1,-1,2,-2,3,-3,4,-4]

gidgetdir = os.environ['GIDGETDIR']+'/'



def filledPlots(qvecs, perc, thisCM, ax, xx):

    if(perc is not None):
        if(len(perc) % 2 == 0):
            # a unique fill between each quantile
            for k in range(len(perc)-1):
                ax.fill_between(xx, qvecs[k,:], qvecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(perc))),alpha=0.3)
        else:
            # symmetric fills around the central quantile.
            v2 = (len(perc)-1)/2
            for k in range(v2):
                col = thisCM(float(k+.5)/float(v2+.5))
                #pdb.set_trace()
                ax.fill_between(xx, qvecs[k,:], qvecs[(k+1),:], facecolor=col, alpha=0.3)
                ax.fill_between(xx, qvecs[(-k-2),:], qvecs[(-k-1),:], facecolor=col, alpha=0.3)
                ax.plot( xx, qvecs[k,:], lw=k+1, color='k' )
                ax.plot( xx, qvecs[(-k-1),:], lw=k+1, color='k' )





class RadialFunction:
    def __init__(self, arr, name, cgsConv=1.0, sensibleConv=1.0,
            texString='', inner=None, outer=None, log=True, theRange=None):
        self.arr=arr 
        self.name=name
        self.cgsConv=cgsConv
        self.sensibleConv = sensibleConv
        self.texString = texString
        if(inner is not None):
            self.innerVal = inner
        else:
            self.innerVal = arr[:,0]
        if(outer is not None):
            self.outerVal = outer
        else:
            self.outerVal = arr[:,-1]
        self.log = log
        self.theRange = theRange
    def inner(self,timeIndex=None,cgs=False):
        if(cgs):
            if(timeIndex is not None and hasattr(self.innerVal,"__len__")):
                return self.innerVal[timeIndex]*self.cgsConv
            return self.innerVal*self.cgsConv
        else:
            if(timeIndex is not None and hasattr(self.innerVal,"__len__")):
                return self.innerVal[timeIndex]*self.sensibleConv
            return self.innerVal*self.sensibleConv
    def outer(self,timeIndex=None,cgs=False):
        if(cgs):
            if(timeIndex is not None and hasattr(self.outerVal,"__len__")):
                return self.outerVal[timeIndex]*self.cgsConv
            return self.outerVal*self.cgsConv
        else:
            if(timeIndex is not None and hasattr(self.outerVal,"__len__")):
                return self.outerVal[timeIndex]*self.sensibleConv
            return self.outerVal*self.sensibleConv
    def cgs(self,timeIndex=None,locIndex=None):
        if(timeIndex is None and locIndex is None):
            return self.arr[:,:]*self.cgsConv
        elif(timeIndex is not None and locIndex is not None):
            return self.arr[timeIndex,locIndex]*self.cgsConv
        elif(timeIndex is not None and locIndex is None):
            return self.arr[timeIndex,:]*self.cgsConv
        elif(timeIndex is None and locIndex is not None):
            return self.arr[:,locIndex]*self.cgsConv
        else:
            print "Something has gone wrong in cgs"
    def sensible(self,timeIndex=None,locIndex=None):
        if(timeIndex is None and locIndex is None):
            return self.arr[:,:]*self.sensibleConv
        elif(timeIndex is not None and locIndex is not None):
            return self.arr[timeIndex,locIndex]*self.sensibleConv
        elif(timeIndex is not None and locIndex is None):
            return self.arr[timeIndex,:]*self.sensibleConv
        elif(timeIndex is None and locIndex is not None):
            return self.arr[:,locIndex]*self.sensibleConv
        else:
            print "Something has gone wrong in sensible"
    def atR(self, rNew, rVec, timeIndex, sensible=True):
        #if(len(rVec) != len(self.sensible(timeIndex=timeIndex))):
        #    pdb.set_trace()
        assert len(rVec) == len(self.sensible(timeIndex=timeIndex))
        if(sensible):
            f = interp1d(rVec,self.sensible(timeIndex=timeIndex),kind='linear',bounds_error=False)
        else:
            f = interp1d(rVec,self.cgs(timeIndex=timeIndex),kind='linear',bounds_error=False)
        rMin = np.min(rVec)
        rMax = np.max(rVec)
        onesVec = np.ones(np.shape(rNew))
        ret = np.where(rNew < rMin*onesVec, self.inner(timeIndex)*onesVec, f(rNew))
        ret = np.where(rNew > rMax*onesVec, self.outer(timeIndex)*onesVec, ret)
        return ret

    def range(self):
        if(self.theRange is None):
            return self.theRange
        else:
            return np.percentile(self.arr,(.2,99.8))*sensibleConv

def constructQuantiles(rVecs, varVecs, quantiles, ti):
    ''' Given a vector of radius vectors, a vector of variable vectors, and a list of quantiles:
        Set up a radial grid, interpolate the data from each model to these radii, and find the requested quantiles.
        Return a vector of vectors to be plotted.'''
    rmin = np.min(rVecs)
    rmax = np.max(rVecs)
    r = np.linspace(rmin,rmax,num=500)
    newVars = np.zeros((len(varVecs),500))
    for k in range(len(varVecs)):
        newVars[k,:] = varVecs[k].atR(r, rVecs[k,ti,:], ti)
    perc = np.percentile(newVars,quantiles,axis=0,overwrite_input=False)
    return r,perc 



        


class TimeFunction:
    def __init__(self,arr,name,cgsConv=1.0,sensibleConv=1.0,texString='No label provided',log=True,theRange=None):
        self.arr = np.array(arr)
        self.name = name
        self.cgsConv = cgsConv
        self.sensibleConv = sensibleConv
        self.texString = texString
        self.log = log
        self.theRange=theRange
    def cgs(self,timeIndex=None,locIndex=None):
        if(timeIndex is not None):
            return self.arr[timeIndex]*self.cgsConv
        return self.arr * self.cgsConv
    def sensible(self,timeIndex=None,locIndex=None):
        if(timeIndex is not None):
            return self.arr[timeIndex]*self.sensibleConv
        return self.arr * self.sensibleConv
    def range(self):
        if(self.theRange is None):
            return self.theRange
        else:
            return np.percentile(self.arr,(.2,99.8))*sensibleConv


class SingleModel:
    def __init__(self,path):
        self.path=path
        rfs = path.rfind('/') # location of last '/' in the path.
        self.name = path[rfs+1:] # the run name follows the '/'
        self.dirname = path[:rfs] # the directory name preceds the '/'
        self.p={}
    def plotSequenceOfTimes(self,rfx, rfy, times=None):
        figure,ax = plt.subplots(1,1)
        if times is None:
            times = range(1,self.nTimeSteps())
        colorby = 't'
        for ti in times:
            colors,fail,log,_ = self.constructQuantity(colorby,timeIndex=ti)
            cNorm = mpcolors.Normalize(vmin=np.min(colors),vmax=np.max(colors))
            scalarMap = cmx.ScalarMappable(cmap=cm,norm=cNorm)
            scalarMap.set_array(colors)
            scalarMap.autoscale()
            cbar = plt.colorbar(scalarMap,ax=ax)
            cbar.set_label(colorby)
            plt.plot(self.var[rfx].sensible(ti), self.var[rfy].sensible(ti),c=rgb[ti])
        plt.savefig(self.name+'_'+rfy+'_vs_'+rfx+'.png')
        plt.close(figure)
    def TwoDimSFR(self,filename,inclination,orientation,timeIndex):
        cosmo = {'omega_M_0' : .258, 'omega_lambda_0':1-.258 , 'omega_b_0':.17*.258 , 'omega_k_0':0, 'h':.7 }
        d_a = cd.angular_diameter_distance(self.var['z'].sensible(timeIndex), **cosmo) # ang diam dist in Mpc
        dx = (d_a *1.0e3) *.005 /(60.0 * 60.0 * 180/pi) # .005 arsec
        Nx = 2000 # 10 arsec / .005 arcsec = 10^3 'pixels'
        x0 = -dx*Nx/2.0*(cos(orientation) - sin(orientation))
        y0 = -dx*Nx/2.0*(sin(orientation) + cos(orientation))
        pixelVals = np.zeros((Nx,Nx))
        pixelCentersX = np.zeros((Nx,Nx))
        pixelCentersY = np.zeros((Nx,Nx))
        annuliA = self.var['rb'].sensible(timeIndex) 
        annuliB = annuliA*cos(inclination)
        # set up our grid.
        for i in range(Nx):
            for j in range(Nx):
                pixelCentersX[i,j] = x0 + i*dx*cos(orientation) - j*dx*sin(orientation)
                pixelCentersY[i,j] = y0 + i*dx*sin(orientation) + j*dx*cos(orientation)
        for i in range(Nx):
            for j in range(Nx):
                # Are we within any annulus?
                if ((pixelCentersX[i,j]/annuliA[-1])**2.0 + (pixelCentersY[i,j]/annuliB[-1])**2.0 <= 1.0):
                    failsAt = -1
                    for k in range(len(annuliA)-1):
                        if((pixelCentersX[i,j]/annuliA[-2-k])**2.0 + (pixelCentersY[i,j]/annuliB[-2-k])**2.0 > 1.0):
                            failsAt = k
                            pixelVals[i,j] = self.var['colsfr'].sensible(timeIndex,-1-failsAt)
                            break
                    # if now that we've gone through all the annuli and we didn't find any matches:
                    if(failsAt == -1):
                        pixelVals[i,j] = (self.var['mdotBulgeG'].sensible(timeIndex)*self.p['RfREC']/(self.p['RfREC']+self.var['MassLoadingFactor'].inner()))/ (pi * annuliA[0]**2.0 )#annuliB[0]) # SFR (Msun/yr) per area (kpc^2)
                
        pixelVals = pixelVals / cos(inclination) # projected surface density increases with inclination.
        hdu = pyfits.PrimaryHDU(pixelVals)
        hdulist = pyfits.HDUList([hdu])
        if(os.path.exists(filename+'.fits')):
            os.remove(filename+'.fits')
        hdulist.writeto(filename+'.fits')
    def nTimeSteps(self):
        return self.nt
    def TotalWithinR(self,radialQuantity,centralValue,cutoffRadius):
        nr= len(self.var['r'].sensible(timeIndex=0))
        sums = np.zeros(self.nt)
        for ti in range(self.nt):
            theSum = centralValue[ti]
            for ri in range(nr):
                if(self.var['r'].sensible(timeIndex=ti,locIndex=ri) > cutoffRadius):
                    break
                theSum += radialQuantity[ti,ri]
            sums[ti] = theSum
        return sums

    def read(self, keepOnly=[], keepStars=False, paramsOnly=False):
        with open(self.path+'_comment.txt','r') as comment:
            lines = comment.readlines()
            # paramnames - copied from exper.py's experiment class.
            paramnames = ['nx','eta','epsff','tauHeat','analyticQ', \
                    'cosmologyOn','xmin','NActive','NPassive','vphiR', \
                    'R','gasTemp','Qlim','fg0','phi0', \
                    'zstart','tmax','stepmax','TOL','muNorm', \
                    'muColScaling','b','innerPowerLaw','softening', \
                    'diskScaleLength','whichAccretionHistory','alphaMRI' \
                    ,'thickness','migratePassive','fixedQ','kappaMetals', \
                    'Mh0','minSigSt','NChanges','dbg','accScaleLength', \
                    'zquench','zrelax','xiREC','RfREC','deltaOmega', \
                    'Noutputs','accNorm','accAlphaZ','accAlphaMh', \
                    'accCeiling','fscatter','invMassRatio','fcool', \
                    'whichAccretionProfile','alphaAccretionProfile', \
                    'widthAccretionProfile','fH2Min','tDepH2SC','ZIGM','yREC', \
                    'concentrationRandomFactor','muFgScaling']
            params=[]
            line = lines[-1] # get the last line
            tmp = line.split() # split the line into a list of strings
            for k,pname in enumerate(paramnames):
                self.p[pname] = float(tmp[k+2]) # k=0 is the run name which we already have
            #for k,line in enumerate(lines):
            #    cloc = line.find(':')
            #    if(cloc != -1):
            #        params.append(float(line[cloc+1:-1]))
            for k,p in enumerate(params):
                self.p[paramnames[k]] = p
        auxparams=[]
        with open(self.path+'_aux.txt','r') as aux:
            lines = aux.readlines()
            for k,line in enumerate(lines):
                cloc = line.find(':')
                if(cloc!=-1):
                    auxparams.append(float(line[cloc+1:-1]))
        self.p['md0'] = auxparams[0]
        self.p['ND08attempts'] = auxparams[1]
        irfs=[]
        try:
            with open(self.path+'_inputRandomFactors.txt','r') as irf:
                # Record input random factors, if they exist, as x0, x1, ...
                lines = irf.readlines()
                for k,line in enumerate(lines):
                    irfs.append('x'+str(k))
                    self.p[irfs[-1]] = float(line)
        except:
            print "Failed to read in <path>_inputRandomFactors.txt"
        irfYs=[]
        try:
            with open(self.path+'_inputRandomFactorsY.txt','r') as irfY:
                # Record input random factors, if they exist, as y0, y1, ...
                lines = irfY.readlines()
                for k,line in enumerate(lines):
                    irfYs.append('y'+str(k))
                    self.p[irfYs[-1]] = float(line)
        except:
            print "Failed to read in <path>_inputRandomFactorsY.txt"
        self.pLog={}
        for par in self.p.keys():
            self.pLog[par] = True # set all parameters to be logarithmic by default
        # except for the following
        for par in ['xiREC','b','innerPowerLaw','zrelax','zstart','zquench','concentrationRandomFactor','accAlphZ','accAlphaMh','alphaAccretionProfile','whichAccretionHistory','whichAccretionProfile','concentrationRandomFactor','muFgScaling','ND08attempts','muColScaling']+irfs+irfYs:
            self.pLog[par] = False
        if paramsOnly:
            return
        with open(self.path+'_evolution.dat','r') as evolution:
            ncolev, =struct.unpack('i',evolution.read(4))
            self.ncolev=ncolev
            arr = evolution.read()
            evarray = np.fromstring(arr)
            self.nsteps = len(evarray)/ncolev
            try:
                self.evarray = np.reshape(evarray,(self.nsteps,self.ncolev))
            except:
                print "Failed to reshape array. Off by one error?", np.shape(evarray),self.nsteps,self.ncolev,self.nsteps*self.ncolev,self.path
                
        with open(self.path+'_radial.dat','r') as radial:
            dataCube=[]
            for i in range(self.nsteps):
                ncolstep, = struct.unpack('i',radial.read(4))
                nrowstep, = struct.unpack('i',radial.read(4))
                strSqArr = radial.read(8*ncolstep*nrowstep)
                sqArr = np.fromstring(strSqArr)
                dataCube.append(np.reshape(sqArr,(nrowstep,ncolstep)))
            self.dataCube = np.array(dataCube)
        if keepStars:
            with open(self.path+'_stars.dat','r') as stars:
                for i in range(self.nsteps):
                    NABp1, = struct.unpack('i', stars.read(4)) # NAgeBins+1
                    sz, = struct.unpack('i', stars.read(4)) # sps.size()
                    nnx, = struct.unpack('i', stars.read(4)) # nx
                    x = np.fromstring(stars.read(8*nnx))
                    if i==0:
                        age = np.zeros((NABp1, self.nsteps))
                        startingAge = np.zeros((NABp1, self.nsteps))
                        endingAge = np.zeros((NABp1, self.nsteps))
                        col = np.zeros((NABp1, self.nsteps, self.p['nx']))
                        sigR = np.zeros((NABp1, self.nsteps, self.p['nx']))
                        sigZ = np.zeros((NABp1, self.nsteps, self.p['nx']))
                        Zst = np.zeros((NABp1, self.nsteps, self.p['nx']))
                        VarZst = np.zeros((NABp1, self.nsteps, self.p['nx']))
                    for j in range(sz):
                        temp = np.fromstring( stars.read(8*1))
                        age[j,i] = temp[0]
                        startingAge[j,i] = np.fromstring(stars.read(8*1))[0]
                        endingAge[j,i] =np.fromstring( stars.read(8*1))[0]
                        col[j,i,:] = np.fromstring(stars.read(8*nnx))
                        sigR[j,i,:] = np.fromstring(stars.read(8*nnx))
                        sigZ[j,i,:] = np.fromstring(stars.read(8*nnx))
                        Zst[j,i,:] = np.fromstring(stars.read(8*nnx))
                        VarZst[j,i,:] = np.fromstring(stars.read(8*nnx))

        # Alright, at this point we've read in the critical data.
        # Let's try to store it in a more comprehensible manner.
        # Keep in mind that self.dataCube ~ (timestep, nx, var)
        # And self.evarray ~ (timestep, var)
        self.var={}
        starList=[]
        if keepStars:
            for j in range(NABp1):
                stj=str(j).zfill(2)
                self.var['colst'+stj] = RadialFunction( copy.deepcopy(col[j]), 'colst'+stj, \
                    self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                    self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                    r'$\Sigma_{*,'+stj+'} (M_\odot\ pc^{-2})$',theRange=[0.5,3000])
                self.var['sigstR'+stj] = RadialFunction( \
                        np.copy(sigR[j]*self.p['vphiR']), 'sigstR'+stj, \
                        1.0e5,1.0, r'$\sigma_{r,*,'+stj+'}$ (km s$^{-1}$)',log=True)
                self.var['sigstZ'+stj] = RadialFunction( \
                        np.copy(sigZ[j]*self.p['vphiR']), 'sigstZ'+stj, \
                        1.0e5,1.0, r'$\sigma_{z,*,'+stj+'}$ (km s$^{-1}$)',log=True)
                self.var['Zst'+stj] = RadialFunction( \
                        np.copy(Zst[j]),'Zst'+stj,cgsConv=1.0,sensibleConv=1.0/.02, \
                        texString=r'$Z_{*,'+stj+'} (Z_\odot)$')
                self.var['VarZst'+stj] = RadialFunction( \
                        np.copy(np.sqrt(VarZst[j])),'VarZst'+stj, \
                        cgsConv=1.0,sensibleConv=1.0/.02,texString=r'$Z_{*,'+stj+'} (Z_\odot)$')
                self.var['ageSt'+stj] = TimeFunction( \
                        np.copy(age[j]), 'ageSt'+stj, \
                        cgsConv = speryear, texString='Age (Gyr)',log=False)
                self.var['startingAgeSt'+stj] = TimeFunction( \
                        np.copy(startingAge[j]), 'startingAgeSt'+stj, \
                        cgsConv = speryear, texString='Age (Gyr)',log=False)
                self.var['endingAgeSt'+stj] = TimeFunction( \
                        np.copy(endingAge[j]), 'endingAgeSt'+stj, \
                        cgsConv = speryear, texString='Age (Gyr)',log=False)
                starList = starList + [ 'colst'+stj, 'sigstR'+stj, 'sigstZ'+stj, 'Zst'+stj, 'VarZst'+stj, 'ageSt'+stj, 'startingAgeSt'+stj, 'endingAgeSt'+stj ]


        self.var['step'] = TimeFunction( \
                np.copy(self.evarray[:,0]), \
                'step',1,1,'Number of Steps',log=False)
        self.var['t'] = TimeFunction(\
                np.copy(self.evarray[:,1]),'t', \
                2.0*pi*cmperkpc*self.p['R']/(1.0e5*self.p['vphiR']) ,\
                2.0*pi*cmperkpc*self.p['R']/(1.0e5*self.p['vphiR']*speryear*1.0e9), \
                'Time since zstart (Gyr)',log=False)
        self.nt = len(self.evarray[:,1])
        self.var['z'] = TimeFunction(np.copy(self.evarray[:,9]),'z',1,1,'z',log=False)
        self.var['r'] = RadialFunction( \
                np.copy(self.dataCube[:,:,0]),'r', \
                self.p['R']*cmperkpc,self.p['R'],'r (kpc)',log=False)
        dlnx=-log(self.p['xmin'])/(self.p['nx']-1.0)
        nxI = int(self.p['nx'])
        # This assumes the grid is fixed and logarithmic. If the code is modified so that
        # this changes, dx could easily be printed by the code.
        self.var['dr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,0]) * dlnx, 'dr', \
                self.p['R']*cmperkpc,self.p['R'],'$\Delta$r (kpc)')
        # r on the boundaries between cells (code units, i.e. x=r/R)
        internalR = np.sqrt(self.dataCube[:,0:-1,0]*self.dataCube[:,1:,0])
        self.var['rb'] = RadialFunction( \
            np.column_stack((internalR[:,0]-self.dataCube[:,0,0]*dlnx, \
                             internalR, \
                             internalR[:,-1]+self.dataCube[:,-1,0]*dlnx)),
            'rb',self.p['R']*cmperkpc,self.p['R'],'r border (kpc)',log=False)
#        self.var['dA'] = RadialFunction( \
#                2.0*pi* self.dataCube[:,:,0] *self.dataCube[:,:,0] * sinh(dlnx), 'dA', \
#                (self.p['R']*cmperkpc)**2.0,self.p['R']**2.0,'$\Delta$A (kpc$^2$)')
        self.var['MassLoadingFactor'] = RadialFunction( \
                np.copy(self.dataCube[:,:,54]),'MassLoadingFactor', \
                1.0,1.0, r'$\dot{\Sigma}_{out}/\dot{\Sigma}_*^{SF}$', theRange=[0.1,100.0])
        self.var['dA'] = RadialFunction( \
                pi*(np.power(self.getData('rb',locIndex=range(1,nxI+1),cgs=True),2.0) \
                - np.power(self.getData('rb',locIndex=range(nxI),cgs=True),2.0)), \
                'dA',1.0,1.0/cmperkpc**2.0, r'dA (kpc$^{2}$)',inner=pi*np.power(self.getData('rb',locIndex=0,cgs=True),2.0))
        self.var['col'] =RadialFunction( \
                np.copy(self.dataCube[:,:,3]),'col', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma (M_\odot\ pc^{-2})$',theRange=[0.5,3000])
        # Conversion from code units to Msun/pc^2 for column density-like units
        colSensibleConv = self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc)
        self.var['colst'] =RadialFunction( \
                np.copy(self.dataCube[:,:,5]),'colst', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma_* (M_\odot\ pc^{-2})$', \
                inner=self.evarray[:,3]*2.0*self.p['md0']*self.p['R']*kmperkpc/(speryear*self.p['vphiR']*self.var['rb'].sensible(None,0)**2.0 * pcperkpc**2.0 *colSensibleConv), theRange=[0.1,3.0e4])
        self.var['sigstR'] = RadialFunction( \
                np.copy(self.dataCube[:,:,6]*self.p['vphiR']), 'sigstR', \
                1.0e5,1.0, r'$\sigma_{r,*}$ (km s$^{-1}$)',log=True)
        self.var['sigstZ'] = RadialFunction( \
                np.copy(self.dataCube[:,:,27]*self.p['vphiR']), 'sigstZ', \
                1.0e5,1.0, r'$\sigma_{z,*}$ (km s$^{-1}$)',log=True)
        self.var['hStars'] = RadialFunction( \
                self.var['sigstZ'].cgs()*self.var['sigstZ'].cgs() / (np.pi*Gcgs*(self.var['col'].cgs()+self.var['colst'].cgs())), \
                'hStars', cgsConv=1.0, sensibleConv=1.0/cmperpc, texString=r'$h_* (pc)$',log=True, theRange=[1.0,3000.0])
        self.var['sig'] =RadialFunction( \
                np.copy(self.dataCube[:,:,4]*self.p['vphiR']),'sig', \
                1.0e5,1.0, r'$\sigma$ (km s$^{-1}$)',log=True, theRange=[5.0,110.0])
        self.var['hGas'] = RadialFunction( \
                self.var['sig'].cgs()*self.var['sig'].cgs() / (np.pi*Gcgs*(self.var['col'].cgs()+ self.var['sig'].cgs()/self.var['sigstZ'].cgs()*self.var['colst'].cgs())), \
                'hGas', cgsConv=1.0, sensibleConv=1.0/cmperpc, texString=r'$h_g (pc)$',log=True, theRange=[1.0,3000.0])
        self.var['maxsig'] = TimeFunction(np.amax(self.dataCube[:,:,4]*self.p['vphiR'],axis=1), \
                'maxsig',1.0e5,1.0,r'$\max(\sigma)$',log=True, theRange=[5,110])
        self.var['mdotBulgeG'] = TimeFunction(np.copy(self.evarray[:,8]),'mdotBulgeG', \
                self.p['md0']*gpermsun/speryear,self.p['md0'],r'$\dot{M}_{\mathrm{Bulge}\ \mathrm{(gas)}} (M_\odot/yr)$', theRange=[1.0e-7,100])
        self.var['mdotAccr'] = TimeFunction( \
                self.evarray[:,19],'mdotAccr',gpermsun/speryear,1.0,r'$\dot{M}_{ext}$')
        self.var['feedingEfficiency'] = TimeFunction( \
                self.var['mdotBulgeG'].sensible()/self.evarray[:,19],'feedingEfficiency', \
                1.0,1.0,r'$\dot{M}_\mathrm{bulge}/\dot{M}_{ext}$', theRange=[1.0e-7,10.0])
        self.var['dcoldt'] = RadialFunction( \
                np.copy(self.dataCube[:,:,7]),'dcoldt',self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']*cmperkpc)**2.0), \
                self.p['md0']/(2.0*pi*(self.p['R'])**2.0),r'$\partial \Sigma/\partial t$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)')
        self.var['dcolstdt'] = RadialFunction( \
                np.copy(self.dataCube[:,:,9]),'dcoldt',self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']*cmperkpc)**2.0), \
                self.p['md0']/(2.0*pi*(self.p['R'])**2.0),r'$\partial \Sigma_*/\partial t$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)')
        self.var['colAccr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,29]),'colAccr', \
                self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']**2.0)*cmperkpc*cmperkpc),\
                self.p['md0']/(2.0*pi*self.p['R']**2.0), \
                r'$\dot{\Sigma}_{cos} (M_\odot\ yr^{-1}\ kpc^{-2})$', \
                inner=2.0*(self.evarray[:,20])*self.p['RfREC']/((self.p['RfREC']+self.var['MassLoadingFactor'].inner())*np.power(internalR[:,0]-self.dataCube[:,0,0]*dlnx,2.0)), \
                theRange = [1.0e-5,1.0])
        mdotCodeUnits = np.column_stack((self.evarray[:,8],np.copy(self.dataCube[:,:,38])/self.p['md0']))
        self.var['Mdot'] = RadialFunction( \
                mdotCodeUnits, 'Mdot', \
                self.p['md0']*gpermsun/speryear, \
                self.p['md0'],r'$\dot{M}$ (M$_\odot$ yr$^{-1}$)',log=False, theRange=[-5,5])
        self.var['colTr'] = RadialFunction( \
                (self.var['Mdot'].cgs(locIndex=range(1,nxI+1))-self.var['Mdot'].cgs(locIndex=range(nxI)))/self.var['dA'].cgs(), \
                'colTr',1.0, speryear*cmperkpc**2.0/gpermsun, r'$\dot{\Sigma}_{tr}$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)',log=False, theRange=[-10,10])
        self.var['colsfr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,28]),'colsfr', \
                self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']**2.0)*cmperkpc*cmperkpc),\
                self.p['md0']/(2.0*pi*self.p['R']**2.0), \
                r'$\dot{\Sigma}_*^{SF} (M_\odot\ yr^{-1}\ kpc^{-2})$', \
                inner=2.0*(self.evarray[:,8]+self.evarray[:,20])*self.p['RfREC']/((self.p['RfREC']+self.var['MassLoadingFactor'].inner())*np.power(internalR[:,0]-self.dataCube[:,0,0]*dlnx,2.0)) , theRange = [1.0e-5,10.0])
        self.var['colOut'] = RadialFunction( \
                self.var['colsfr'].sensible()*self.var['MassLoadingFactor'].sensible(), 'colOut', \
                cgsConv = gpermsun/speryear/cmperkpc**2, sensibleConv=1.0,
                inner = self.var['colsfr'].inner()*self.var['MassLoadingFactor'].inner(), texString=r'$\dot{\Sigma}_\mathrm{out}$', theRange=[1.0e-5, 10.0])
        self.var['colTrPerAccr'] = RadialFunction( np.abs(self.var['colTr'].cgs())/self.var['colAccr'].cgs(), 'colTrPerAccr', \
                texString=r'$|\dot{\Sigma}_{tr}|/\dot{\Sigma}_{accr}$', theRange=[1.0e-5,10.0])
        self.var['colOutPerAccr'] = RadialFunction( self.var['colOut'].cgs()/self.var['colAccr'].cgs(), 'colOutPerAccr', \
                texString=r'$\dot{\Sigma}_{out}/\dot{\Sigma}_{accr}$', theRange=[1.0e-5, 10.0])
        self.var['colSfrPerAccr'] = RadialFunction( self.var['colsfr'].cgs()/self.var['colAccr'].cgs(), 'colSfrPerAccr', \
                texString=r'$\dot{\Sigma}_{out}/\dot{\Sigma}_{accr}$', theRange=[1.0e-5, 10.0])

        self.var['Q'] = RadialFunction( \
                np.copy(self.dataCube[:,:,11]),'Q',1.0,1.0,r'Q',theRange=[1.0,10.0])
        self.var['Qg'] = RadialFunction( \
                np.copy(self.dataCube[:,:,23]),'Qg',1.0,1.0,r'$Q_g$',theRange=[.5,10.0])
        self.var['Qst'] = RadialFunction( \
                np.copy(self.dataCube[:,:,22]),'Qst',1.0,1.0,r'$Q_*$',theRange=[.5,10.0])
        self.var['tDepRadial'] = RadialFunction( \
                self.var['col'].cgs()/self.var['colsfr'].cgs(), 'tDepRadial',\
                1.0, 1.0/speryear, r'$t_\mathrm{dep} = \Sigma/\dot{\Sigma}_*^{SF} (yr)$',theRange=[3.0e8,1.0e12])
        self.var['fH2']= RadialFunction(np.copy(self.dataCube[:,:,47]),'fH2',1.0,1.0,r'$f_{\mathrm{H}_2}$',log=False,theRange=[0.0,1.0])
        self.var['colH2']= RadialFunction(np.copy(self.dataCube[:,:,47]) * self.var['col'].sensible(),'colH2', \
                cgsConv = gpermsun/cmperpc**2, texString=r'$\Sigma_{\mathrm{H}_2}$',log=True,theRange=[0.1,3000.0])
        self.var['colHI']= RadialFunction((1.0-np.copy(self.dataCube[:,:,47])) * self.var['col'].sensible(),'colHI', \
                cgsConv = gpermsun/cmperpc**2, sensibleConv=1, texString=r'$\Sigma_{\mathrm{HI}}$',log=True,theRange=[0.1,3000.0])
        self.var['Z'] = RadialFunction(np.copy(self.dataCube[:,:,21]),'Z',cgsConv=1.0,sensibleConv=1.0/.02,texString=r'$Z_g (Z_\odot)$')
        self.var['vPhi'] = RadialFunction(np.copy(self.dataCube[:,:,15]),'vPhi',self.p['vphiR']*1.0e5,self.p['vphiR'], \
                 r'$v_\phi$ (km/s)',log=False)
        self.var['Mh'] = TimeFunction(self.evarray[:,18],'Mh',gpermsun,1.0,r'$M_h (M_\odot)$')
        self.var['mCentral'] = TimeFunction( \
                self.evarray[:,3], 'mCentral',
                self.p['md0']*2.0*pi*self.p['R']/self.p['vphiR'] *gpermsun* kmperkpc/speryear, \
                self.p['md0']*2.0*pi*self.p['R']/self.p['vphiR'] * kmperkpc/speryear, \
                r'$M_\mathrm{center}\ (M_\odot)$', theRange=[3.0e7, 5.0e10])
        self.var['mstar'] = TimeFunction( \
                self.var['mCentral'].sensible() + 
                np.sum( self.var['dA'].sensible()*self.var['colst'].sensible()*1.0e6, 1 ), \
                'mstar',gpermsun,1.0,r'$M_*$ (M$_\odot$)')


        rAcc = -(self.var['r'].sensible(locIndex=1) - self.var['r'].sensible(locIndex=2)) \
                /np.log(self.var['colAccr'].sensible(locIndex=1)/self.var['colAccr'].sensible(locIndex=2))
        if(rAcc[0] < 0): # we're not using an exponential profile apparently!
            rAccInd = np.argmax(self.var['colAccr'].sensible()*self.var['r'].sensible(), axis=1)
            rAcc = self.var['r'].sensible(timeIndex=range(self.nt),locIndex=rAccInd)
#        self.var['scaleRadius'] = TimeFunction( \
#                self.p['accScaleLength']*np.power(self.var['Mh'].sensible()/self.p['Mh0'],self.p['alphaAccretionProfile']), \
#                'accScaleLength',cmperkpc,1.0,r'$r_\mathrm{acc}$ (kpc)',log=False)
        self.var['accretionRadius'] = TimeFunction( \
                rAcc, 'accretionRadius', cgsConv=cmperkpc, texString=r'$r_\mathrm{acc}$',log=False)



        mbulge,_,fcentral,scaleLengths = self.computeMBulge()
        self.var['mBulge'] = TimeFunction(mbulge,'mBulge',gpermsun,1.0,r'$M_B (M_\odot)$',theRange=[3.0e7,5.0e10])
        self.var['scaleLength'] = TimeFunction(scaleLengths, 'scaleLength', cmperkpc, 1.0, r'$r_\mathrm{scale}$')
        self.var['BT'] = TimeFunction(self.var['mBulge'].cgs()/(self.var['mBulge'].cgs()+self.var['mstar'].cgs()), \
                'BT',1,1,r'Bulge to Total Ratio',log=False)
        self.var['fCentral'] = TimeFunction(fcentral,'fCentral',1,1,r'Central Excess/$M_B$',log=False)


        cos = halo.Cosmology()
        vPhiDM = np.zeros((len(self.var['Mh'].sensible()), len(self.var['r'].sensible(timeIndex=0))))
        vPhiBulge = np.zeros((len(self.var['Mh'].sensible()), len(self.var['r'].sensible(timeIndex=0))))
        for i in range(len(self.var['Mh'].sensible())):
            thisHalo = halo.halo(self.var['Mh'].sensible(timeIndex=i),self.var['z'].sensible(timeIndex=i),cos,0) 
            for j in range(len(self.var['r'].sensible(timeIndex=0))):
                vPhiDM[i,j] = np.sqrt( Gcgs * thisHalo.mInterior(self.var['r'].sensible(timeIndex=0,locIndex=j)) / self.var['r'].cgs(timeIndex=0,locIndex=j) )
                vPhiBulge[i,j] = np.sqrt( Gcgs * self.var['mBulge'].cgs(timeIndex=i)*self.var['fCentral'].sensible(timeIndex=i) / self.var['r'].cgs(timeIndex=0,locIndex=j) )

        self.var['vPhiDM'] = RadialFunction(np.copy(vPhiDM),'vPhiDM', cgsConv=1.0, sensibleConv=1.0e-5, \
                 texString=r'$v_{\phi,\mathrm{DM}}$ (km/s)',log=False)
        self.var['vPhiBulge'] = RadialFunction(np.copy(vPhiBulge), 'vPhiBulge', cgsConv=1.0, sensibleConv=1.0e-5, \
                 texString=r'$v_{\phi, \mathrm{bulge}}$ (km/s)',log=False)
        self.var['vPhiDisk'] = RadialFunction( np.sqrt( np.power(self.var['vPhi'].cgs(),2.0) - np.power(self.var['vPhiDM'].cgs(),2.0) - np.power(self.var['vPhiBulge'].cgs(),2.0)), 'vPhiDisk', cgsConv=1.0, sensibleConv=1.0e-5, texString=r'$v_{\phi, \mathrm{disk}} (km/s)$', log=False)


        self.var['vOverSigGas'] = RadialFunction(self.var['vPhi'].sensible()/self.var['sig'].sensible(),'vOverSigGas',texString=r'$v_\phi/\sigma$',log=True)
        self.var['vOverSigStR'] = RadialFunction(self.var['vPhi'].sensible()/self.var['sigstR'].sensible(),'vOverSigStR',texString=r'$v_\phi/\sigma_{*,r}$',log=True)
        self.var['vOverSigStZ'] = RadialFunction(self.var['vPhi'].sensible()/self.var['sigstZ'].sensible(),'vOverSigStZ',texString=r'$v_\phi/\sigma_{*,z}$',log=True)
        self.var['spEnergy'] = RadialFunction(np.sqrt(np.power(self.var['vPhi'].sensible(),2.0) + 1.5*np.power(self.var['sig'].sensible(),2.0)), 'spEnergy', sensibleConv=1.0, cgsConv=1.0e5, texString=r'$\sqrt{v_\phi^2 + 1.5 \sigma^2}$')

        self.var['vrst'] = RadialFunction(np.copy(self.dataCube[:,:,34]),'vrst',self.p['vphiR']*1.0e5,self.p['vphiR'], \
                 r'$v_{r,*}$',log=False)
        self.var['vrg'] = RadialFunction(np.copy(self.dataCube[:,:,36]),'vrg',self.p['vphiR']*1.0e5,self.p['vphiR'], \
                 r'$v_{r,g}$',log=False)
        self.var['NHI'] = RadialFunction(self.getData('col',cgs=True)*(1.0-self.getData('fH2'))*(1.0-self.getData('Z',cgs=True))/gperH,\
                'NHI', 1.0,1.0,r'$N_{\mathrm{HI}}$ (cm$^{-2}$)',theRange=[1.0e19,3.0e21])
        self.var['rx'] = RadialFunction( \
                np.array([self.var['r'].sensible(ti)/self.var['accretionRadius'].sensible(timeIndex=ti) for ti in range(self.nt)]), \
                'rx',1.0,1.0,r'r/r$_{acc}$',log=False)
        self.var['sfr'] = TimeFunction( \
                self.var['mdotBulgeG'].sensible()*self.p['RfREC']/(self.p['RfREC']+self.var['MassLoadingFactor'].inner()) \
                +np.sum( self.var['dA'].sensible()*self.var['colsfr'].sensible(), 1 ), \
                'sfr',gpermsun/speryear, 1.0, r'SFR (M$_\odot$ yr$^{-1}$)')
        self.var['sfrPerAccr'] = TimeFunction( \
                self.var['sfr'].cgs()/self.var['mdotAccr'].cgs(),'sfrPerAccr',1.0,1.0,r'SFR / $\dot{M}_{ext}$', theRange=[1.0e-3,10.0])
        def areaWeightedWithin( varName, maxRadius):
            radiusInd = np.searchsorted( self.var['r'].sensible(timeIndex=-1), maxRadius)
            mInt = np.sum(self.var['dA'].cgs(locIndex=range(radiusInd))*self.var[varName].cgs(locIndex=range(radiusInd)), 1)
            return mInt
        eightkpc = np.searchsorted( self.var['r'].sensible(timeIndex=-1), 8.0 )
        m1kpc = self.var['mCentral'].cgs() + areaWeightedWithin('colst', 1.0) # mass within 1 kpc 
        self.var['Sigma1'] = TimeFunction( \
                m1kpc/(1.0*cmperkpc)**2.0, 'Sigma1', cgsConv=1.0, sensibleConv = cmperkpc**2/gpermsun, \
                texString=r'$\Sigma_1$')
        self.var['mgas']=TimeFunction( \
                np.sum(self.var['dA'].cgs()*self.var['col'].cgs(), 1),'mgas', \
                1.0,1.0/gpermsun,r'$M_g$ (M$_\odot$)', theRange=[1.0e8,1.0e11])
        self.var['tdep']=TimeFunction( \
                self.var['mgas'].cgs()/self.var['sfr'].cgs(),'tdep',1.0,1.0/speryear,r't$_{dep}$ (yr)',theRange=[1.0e8,3.0e10])
        self.var['efficiency'] = TimeFunction( \
                self.var['mstar'].cgs()/self.var['Mh'].cgs(), 'efficiency',1.0,1.0,r'$M_*/M_h$',log=True)
        self.var['fg'] = TimeFunction( \
                self.var['mgas'].cgs() / (self.var['mgas'].cgs()+self.var['mstar'].cgs()), \
                'fg',1.0,1.0,r'$f_g$', log=False)
        self.var['fgRadial'] = RadialFunction( \
                self.var['col'].cgs()/(self.var['col'].cgs()+self.var['colst'].cgs()), \
                'fgRadial',1.0,1.0,r'$f_g = \Sigma/(\Sigma_*+\Sigma)$', log=False)

        self.var['tDepH2Radial'] = RadialFunction( \
                self.var['fH2'].cgs()*self.var['col'].cgs()/self.var['colsfr'].cgs(), 'tDepH2Radial',\
                1.0, 1.0/speryear, r'$t_{\mathrm{dep},H_2}\ \mathrm{(yr)}$', theRange=[3.0e8,1.0e11])
        self.var['sSFR'] = TimeFunction(self.getData('sfr',cgs=True)/self.getData('mstar',cgs=True) , \
                'sSFR',1.0,1.0e9*speryear,r'sSFR (Gyr$^{-1}$)', theRange=[1.0e-4,1.0e2])
        mg = self.getData('col',cgs=True)*self.getData('dA',cgs=True)
        self.var['integratedZ'] = TimeFunction(np.sum(mg*self.getData('Z',cgs=True),axis=1)/np.sum(mg,axis=1), \
                'integratedZ', cgsConv=1.0, sensibleConv=1.0/.02, texString=r'$Z_g (Z_\odot)$')
        sfRad = self.getData('colsfr',cgs=True)*self.getData('dA',cgs=True)
        self.var['sfZ'] = TimeFunction(np.sum(sfRad*self.getData('Z',cgs=True),axis=1)/np.sum(sfRad,axis=1), \
                'sfZ',cgsConv=1.0,sensibleConv=1.0/.02,texString=r'$Z_g (Z_\odot)$ weighted by SFR')
        self.var['fH2Integrated'] = TimeFunction( \
                np.sum(mg*self.var['fH2'].cgs(),axis=1)/np.sum(mg,axis=1), \
                'fH2Integrated',1.0,1.0,r'$f_{\mathrm{H}_2}$',theRange=[0,1], log=False)
        self.var['vPhiGas'] = TimeFunction( \
                np.sum(self.var['vPhi'].cgs() * self.var['dA'].cgs() * self.var['col'].cgs(),axis=1)/self.var['mgas'].cgs(),
                'vPhiGas',cgsConv=1.0,sensibleConv=1.0e-5,texString=r'Gas mass weighted $v_\phi$', log=True)
        self.var['vPhiStars'] = TimeFunction( \
                np.sum(self.var['vPhi'].cgs() * self.var['dA'].cgs() * self.var['colst'].cgs(),axis=1)/self.var['mstar'].cgs(),
                'vPhiStars',cgsConv=1.0,sensibleConv=1.0e-5,texString=r'Stellar mass weighted $v_\phi$', log=True)
        self.var['vPhiOuter'] = TimeFunction( \
                self.var['vPhi'].cgs(locIndex=-1), 'vPhiOuter',cgsConv=1.0,sensibleConv=1.0e-5,texString=r'$v_\phi$ at outer radius', log=True)
        self.var['integratedMLF'] = TimeFunction( \
                np.sum(self.var['MassLoadingFactor'].cgs() * self.var['dA'].cgs() * self.var['colsfr'].cgs(),axis=1)/self.var['sfr'].cgs(),
                'integratedMLF',cgsConv=1.0,sensibleConv=1.0,texString=r'Mass loading factor', theRange=[0.03, 50.0])
        mJeansMask = np.zeros(np.shape(self.var['sig'].sensible()),dtype=float)
        mJeansMask[self.var['fH2'].sensible()>0.1] = 1.0
        self.var['MJeans'] = RadialFunction( \
                np.power(self.var['sig'].cgs(),4.0)/(Gcgs**2.0 *self.var['col'].cgs())*mJeansMask, \
                'MJeans', 1.0, 1.0/gpermsun, r'2D Jeans Mass where $f_{\mathrm{H}_2}>0.1$ ($M_\odot$)', \
                theRange=[1.0e5,1.0e9])

        self.var['MJeansUnmasked'] = RadialFunction( \
                np.power(self.var['sig'].cgs(),4.0)/(Gcgs**2.0 *self.var['col'].cgs()), \
                'MJeansUnmasked', 1.0, 1.0/gpermsun, r'2D Jeans Mass ($M_\odot$)', \
                theRange=[1.0e5,1.0e9])
        mja = np.sum(self.var['MJeansUnmasked'].sensible()*self.var['colsfr'].cgs()*self.var['dA'].cgs(),axis=1)/np.sum(self.var['colsfr'].cgs()*self.var['dA'].cgs(), axis=1)
        self.var['MJeansAvg'] = TimeFunction( mja, \
                'MJeansAvg', cgsConv = gpermsun, texString=r'$M_{Jeans}\ (M_\odot)$', theRange=[1.0e5, 1.0e9])
                
        self.var['centralDensity'] = TimeFunction(  \
            self.TotalWithinR(self.var['dA'].cgs()*self.var['colst'].cgs(),self.var['mCentral'].cgs(),1.0)/(gpermsun*np.pi),
            'centralDensity',gpermsun/cmperkpc**2.0,1.0,r'Average $\Sigma_*$ with $r<1$ kpc')
        mst = self.var['mstar'].cgs()
        mj = self.var['MJeans'].cgs()
        for ti in range(len(mst)):
            mj[ti,:]/=mst[ti]
        self.var['ClumpMassPerDisk'] = RadialFunction( \
                mj,'ClumpMassPerDisk',1.0,1.0,r'$M_J/M_*$')
        colAccr = self.getData('colAccr',cgs=True)
        colTr = self.getData('colTr',cgs=True)
        colSFR = self.getData('colsfr',cgs=True)*(self.var['MassLoadingFactor'].sensible()+self.p['RfREC'])
        dcoldt = self.getData('dcoldt',cgs=True)
        shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)
        self.var['equilibrium'] = RadialFunction( \
                dcoldt/shareNorm, 'equilibrium', 1.0, 1.0, r'Equilibrium',log=False)

        eqnorm = shareNorm*self.var['dA'].cgs()
        self.var['integratedEquilibrium'] = TimeFunction( \
                np.sum(self.var['equilibrium'].cgs()*eqnorm,axis=1)/np.sum(eqnorm,axis=1), \
                'integratedEquilibrium',1.0,1.0,r'Equilibrium',log=False)
        self.var['Sigma1p5'] = TimeFunction( \
                self.var['mstar'].sensible()/np.power(self.var['scaleLength'].sensible(),1.5), 'Sigma1p5', sensibleConv=1.0, cgsConv=gpermsun/cmperkpc**1.5,texString=r'$\Sigma_{1.5}\ (M_\odot/\mathrm{kpc}^{1.5})$', theRange=[1.0e8,1.0e12])
        self.var['specificJRadial'] = RadialFunction( self.var['r'].cgs() * self.var['vPhi'].cgs(), 'specificJRadial', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JColGas'] = RadialFunction( self.var['r'].cgs() *self.var['r'].cgs() * self.var['vPhi'].cgs() * self.var['col'].cgs(), 'JColGas', sensibleConv=1.0e-5/cmperkpc**2/gpermsun*cmperpc**2, cgsConv=1.0, texString=r'$r^2 v_\phi\Sigma\ (\mathrm{kpc}^2\ \mathrm{km}/\mathrm{s}\ M_\odot/\mathrm{pc}^2$')
        self.var['JColStars'] = RadialFunction( self.var['r'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs() * self.var['colst'].cgs(), 'JColSt', sensibleConv=1.0e-5/cmperkpc**2/gpermsun*cmperpc**2, cgsConv=1.0, texString=r'$r^2v_\phi\Sigma_*\ (\mathrm{kpc}^2\ \mathrm{km}/\mathrm{s}\ M_\odot/\mathrm{pc}^2$')

#        try:
#            maxGIRadius = np.zeros(len(self.var['z'].cgs()))
#            for i in range(len(self.var['z'].cgs())):
#                unstable = self.var['Q'].sensible(timeIndex=i) < self.p['fixedQ']  # locations in space/time where Q<fixedQ
#                maxGIRadius[i] = np.max( self.var['r'].sensible(timeIndex=i)[unstable] )
#            self.var['maxGIRadius'] = TimeFunction( maxGIRadius, 'maxGIRadius', cgsConv=cmperkpc, sensibleConv=1.0, texString=r'Largest GI radius (kpc)')
#        except:
#            print "WARNING: failed to calculate max GI radius"
        indRange = range(eightkpc)

        self.var['specificJStars'] = TimeFunction( np.sum(self.var['colst'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colst'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJStars', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ M_*\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJGas'] = TimeFunction( np.sum(self.var['col'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['col'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJGas', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \Sigma_\mathrm{gas}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJH2'] = TimeFunction( np.sum(self.var['colH2'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colH2'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJH2', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \Sigma_{\mathrm{H}_2}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJHI'] = TimeFunction( np.sum(self.var['colHI'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colHI'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJHI', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \Sigma_\mathrm{HI}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJOut'] = TimeFunction( np.sum(self.var['colOut'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colOut'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJOut', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \Sigma_\mathrm{out}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJSFR'] = TimeFunction( np.sum(self.var['colsfr'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colsfr'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJSFR', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \mathrm{SFR}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])
        self.var['specificJAccr'] = TimeFunction( np.sum(self.var['colAccr'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1)/np.sum(self.var['colAccr'].cgs()*self.var['dA'].cgs(),axis=1), 'specificJAccr', sensibleConv=1.0e-5/cmperkpc, cgsConv=1.0, texString=r'$j\ \mathrm{weighted}\ \mathrm{by}\ \mathrm{Accr}\ (\mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])

        self.var['JStars'] = TimeFunction( np.sum(self.var['colst'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JStars', sensibleConv=1.0e-5/cmperkpc/gpermsun, cgsConv=1.0, texString=r'$J_* (M_\odot \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JGas'] = TimeFunction( np.sum(self.var['col'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JGas', sensibleConv=1.0e-5/cmperkpc/gpermsun, cgsConv=1.0, texString=r'$J_g\ (M_\odot \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JH2'] = TimeFunction( np.sum(self.var['colH2'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JH2', sensibleConv=1.0e-5/cmperkpc/gpermsun, cgsConv=1.0, texString=r'$J_{H_2}\ (M_\odot \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JHI'] = TimeFunction( np.sum(self.var['colHI'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JHI', sensibleConv=1.0e-5/cmperkpc/gpermsun, cgsConv=1.0, texString=r'$J_{HI}\ (M_\odot\ \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JOut'] = TimeFunction( np.sum(self.var['colOut'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JOut', sensibleConv=1.0e-5/cmperkpc/gpermsun*speryear, cgsConv=1.0, texString=r'$\dot{J}_\mathrm{out}\ (M_\odot/\mathrm{yr}\ \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JSFR'] = TimeFunction( np.sum(self.var['colsfr'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JSFR', sensibleConv=1.0e-5/cmperkpc/gpermsun*speryear, cgsConv=1.0, texString=r'$\dot{J}_\mathrm{SFR}\ (M_\odot/\mathrm{yr}\ \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$')
        self.var['JAccr'] = TimeFunction( np.sum(self.var['colAccr'].cgs()*self.var['dA'].cgs()*self.var['r'].cgs() * self.var['vPhi'].cgs(), axis=1), 'JAccr', sensibleConv=1.0e-5/cmperkpc/gpermsun * speryear, cgsConv=1.0, texString=r'$\dot{J}_\mathrm{accr}\  (M_\odot/\mathrm{yr}\ \mathrm{kpc}\ \mathrm{km}/\mathrm{s})$', theRange=[10,5000])



        rhoCrit = cos.rhocrit( self.var['z'].sensible() )
        self.var['Rvir'] = TimeFunction(np.power(self.var['Mh'].cgs()*3.0/(4.0*np.pi*200.0*rhoCrit),1./3.), 'Rvir', sensibleConv=1/cmperkpc, cgsConv=1,texString=r'$R_\mathrm{vir}$') # r200
        self.var['Vvir'] = TimeFunction(np.power(Gcgs*self.var['Mh'].cgs()/self.var['Rvir'].cgs(),1./2.), 'Vvir', sensibleConv=1.0e-5, cgsConv=1,texString=r'$V_\mathrm{vir}$') # V
        lambdaRange = [8.0e-3, 1.8e-1]
        self.var['dimensionlessSpinStars'] = TimeFunction(self.var['specificJStars'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinStellarMass', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_*/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinSFR'] = TimeFunction(self.var['specificJSFR'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinSFR', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{SFR}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinGas'] = TimeFunction(self.var['specificJGas'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinGas', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{gas}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinH2'] = TimeFunction(self.var['specificJH2'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinH2', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{H2}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinHI'] = TimeFunction(self.var['specificJHI'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinHI', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{HI}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinOut'] = TimeFunction(self.var['specificJOut'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinOut', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{out}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)
        self.var['dimensionlessSpinAccr'] = TimeFunction(self.var['specificJAccr'].cgs()/(np.sqrt(2.0)*self.var['Rvir'].cgs()*self.var['Vvir'].cgs()), 'dimensionlessSpinAccr', sensibleConv=1.0, cgsConv=1.0, texString=r'$j_\mathrm{accr}/\sqrt{2} R_\mathrm{vir} V_\mathrm{vir}$', theRange=lambdaRange)



        theKeys = self.var.keys()
        whitelist = ['rb','r','dA','rx','dr'] + keepOnly + starList
        if 'vPhi' in whitelist:
            if 'vPhiDM' not in whitelist:
                whitelist += ['vPhiDM']
            if 'vPhiBulge' not in whitelist:
                whitelist += ['vPhiBulge']
        for key in theKeys:
            if(isinstance(self.var[key],RadialFunction) and not (key in whitelist)):
                del self.var[key]

        #npd = (np.diff(np.sign(self.getData('fH2')-0.5),axis=1) != 0)*1
        #LI = []
        #for i in range(len(np.shape(npd)[0])):
        #    ind = np.argmax(npd[i,:]==1) 
        #    if(ind==0):
        #        ind=-1

        #self.var['rHI'] = TimeFunction(self.getData('r',locIndex=npd,cgs=True),'rHI', \
        #        1.0,1.0/cmperkpc,r'$r_{H\mathrm{I}}$ (kpc)')
        del self.dataCube
        del self.evarray
    def getData(self,name,timeIndex=None,locIndex=None,cgs=False):
        '''Get the data associated with RadialFunction or TimeFunction named name. '''
        if name in self.var.keys():
            if(not cgs):
                return self.var[name].sensible(timeIndex=timeIndex,locIndex=locIndex)
            else:
                return self.var[name].cgs(timeIndex=timeIndex,locIndex=locIndex)
        raise ValueError('key not found: '+name)
    def get(self,name):
        ''' Get the RadialFunction or TimeFunction named name. '''
        if name in self.p.keys():
            return self.p[name]
        if name in self.var.keys():
            return self.var[name]
        raise ValueError('key not found: '+name)
    def getRadialFunctions(self):
        rf=[]
        blacklist=['rb','r','dA','rx','dr']
        for key in self.var.keys():
            if(isinstance(self.var[key],RadialFunction) and not key in blacklist):
                rf.append(key)
        return rf
    def getTimeFunctions(self):
        tf=[]
        blacklist=['t','z','step']
        for key in self.var.keys():
            if(isinstance(self.var[key],TimeFunction) and not key in blacklist and not 'ageSt' in key and not 'startingAgeSt' in key and not 'endingAgeSt' in key):
                tf.append(key)
        return tf
    def computeMBulge(self):
        ''' Compute the mass of stars an observer might interpret as being part of a bulge. 
        To do this, we project an exponential profile from 1.5x the scale radius inwards
        We make this profile shallower until the projection lies below the true stellar column
        density at all radii.'''
        maxTries=200
        r = self.var['r'].sensible(0) # Only valid if the grid is fixed!
        mb=[]
        fc=[]
        rs=[]
        failures=[]
        for timeIndex in range(len(self.var['z'].sensible())):
            doneFlag = 0
            scaleRadius =  self.var['accretionRadius'].sensible(timeIndex)
            for i in range(maxTries):
                f=2.5
                ind,theMin = Nearest(r,scaleRadius*f)
                if ind>(len(r)-1):
                    ind = len(r)-2
                elif ind<3:
                    ind = len(r)/4
                logcol = np.log10(self.var['colst'].cgs(timeIndex))
                slope = float(maxTries/2.0-i)/float(maxTries/2.0) * (logcol[ind]-logcol[ind-1])/(r[ind]-r[ind-1])
                yguess = logcol[ind] + (r[:]-r[ind])*slope
                fail = (yguess[0:ind] - logcol[0:ind]).clip(min=0)
                if len(fail)==0:
                    print "Something odd has happened in computeMBulge."
                    print "ind: ",ind
                    print "r,scaleRadius*f: ",r,scaleRadius*f
                    failflag = True
                else:
                    failflag = np.max(fail) > 0
                if(not failflag):
                    theFit = np.power(10.0,logcol[ind] + (r[:]-r[ind])*slope)
                    excess = self.var['dA'].cgs(timeIndex)*(np.power(10.0,logcol) - theFit).clip(min=0)
                    excessL = np.sum(excess[0:ind])
                    excessR = np.sum(excess[ind:])
                    centralExcess =self.var['mCentral'].cgs(timeIndex) - np.power(10.0,logcol[ind] - r[ind]*slope) *pi*self.var['rb'].cgs(timeIndex,0)**2.0
                    if(centralExcess < 0):
                        pass
                    else:
                        doneFlag =1
                        whichTry = i
                        break
            # At this point we've either succeeded(doneFlag=1) or failed (doneFlag=0)
            if(doneFlag==1):
                mb.append(excessL+centralExcess)
                rs.append(-1.0/slope/np.log(10.0))
                fc.append(centralExcess/(excessL+centralExcess))
                failures.append(0)
            else:
                mb.append(0)
                rs.append(0)
                fc.append(0)
                failures.append(1)
        # End loop over time.
        return np.array(mb)/gpermsun,failures,fc,rs
    def quantityAt(self,radius,key,timeIndex):
        '''Given a radius (in kpc), find which annulus we're in. We need to handle a couple edge cases,
        namely the radius is within the inner edge, or the radius is beyond the outer edge of the 
        computational domain.'''
        # Beyond the outer edge
        rb = self.var['rb'].sensible(timeIndex)
        if(radius > rb[-1]):
            return self.var[key].outer(timeIndex)
        elif(radius < rb[0]):
            return self.var[key].inner(timeIndex)
        else:
            guessIndex = int(self.p['nx'] * log(radius / rb[0]) / log(rb[-1]/rb[0]))
            for k in outward:
                if(guessIndex+k>0 and guessIndex+k<len(rb)):
                    if(rb[guessIndex+k]>radius and radius>rb[guessIndex+k-1]):
                        # found it!
                        return self.var[key].sensible(timeIndex,guessIndex+k-1)
                else:
                    pass # went outside the boundary - hopefully not a big deal.
            # Hopefully we don't get here.
            print "Something has gone wrong in quantityAt!"
            pdb.set_trace()
            return 1.0e20

def Nearest(arr,val):
    index = (np.abs(arr-val)).argmin()
    return index,arr[index]

def reachedRedshiftZero(fn):
    ''' Read in the last few lines of the standard output and check whether the code successfully reached redshift zero.
        This can be used to check whether a model has finshed successfully, or whether it's either still running or exited early.'''
    with open(fn+'_stdo.txt','r') as f:
        try:
            f.seek(-200, 2)
            lastFewLines = f.read(200)
        except:
            print "WARNING: failed to read last 200 elements of ",fn+'_stdo.txt'
            lastFewLines=''
        return 'Reached redshift zero' in lastFewLines


class Experiment:
    def __init__(self,name):
        '''Collect every model in every experiment matching name. '''
        self.name=name
        # Search for comment files
        fnameKey = gidgetdir+'analysis/*'+name+'*/*_comment.txt'
        fnames = sorted(glob.glob(fnameKey))
        self.models=[]
        for i,fn in enumerate(fnames):
            fnames[i] = fn[:-12] #
            # Only include a model if it has not left an error message AND if it has reached redshift zero.
            # The former should be a subset of the latter.
            # I recently added the latter to allow the analysis of large experiments which are still running,
            #    excluding the models that are running as all models are read in.
            if(os.stat(fnames[i]+'_stde_aux.txt')[6]==0) and reachedRedshiftZero(fnames[i]):
                self.models.append(SingleModel(fnames[i]))
            else:
                pass
        self.fn = fnames
        if(len(fnames)==0):
            raise ValueError("No models found in experiment "+name+" using key "+fnameKey)
    def merge(self,expt):
        '''Merge two experiments - essentially just append two lists of models. '''
        self.fn = self.fn + expt.fn # concatenate the lists of models
    def assignIndices(self):
        ''' For each model, assign it a parameter equal to its index in the Experiment objects self.models list '''
        for i,model in enumerate(self.models):
            model.p['experIndex'] = float(i)
            model.pLog['experIndex'] = False
    def read(self, keepOnly=[],paramsOnly=False,keepStars=False):
        ''' Read in every model in the experiment. '''
        n=0
        for model in self.models:
            try:
                model.read(keepOnly=keepOnly,paramsOnly=paramsOnly,keepStars=keepStars)
            except:
                print "Failed to read in model n=",n
            n+=1
            if(n % 50 == 0):
                print "Reading in model ",n," of ",len(self.models)
    def rankBy(self,var=None,timeIndex=None,locIndex=None,keepMagnitude=False):
        ''' For each model in the experiment, add a TimeSeries object called rankBy<var> to the model's var structure.
        This variable will store the model's rank as a function of time.'''
        if(var is None):
            var = 'Mh0'
        qu,_,_,_ = self.constructQuantity(var,timeIndex,locIndex) 
        # possible results: 
        #  qu ~ nModels x times x nx   --- Can't deal with this one
        #   --- var is a radial variable with neither timeIndex nor locIndex specified
        #  qu ~ nModels x times --- This is ok
        #   --- one of: var is a timeVar w/ no timeIndex, or var is a radial w/ no ti but a loc specified
        #  qu ~ nModels x nx  --- Can't deal with this one
        #   --- one of: var is a radial w/ only ti specified
        #  qu ~ nModels   --- typical case
        #   --- one of: var is a parameter, var is a timeVariable but we've specified a timeIndex, var is radial w/ ti, loc
        # Which can we deal with? we're making plots or movies vs time. How about this? rank has to be a timeVariable.
        # So we can't deal with nModels x nx or nModels x times x nx:
        nmodels = len(self.models)
        nt = self.models[0].nTimeSteps()
        ranks = np.ones((len(self.models),nt))
        if(qu.ndim == 3 or (qu.ndim==2 and isinstance(self.models[0].var[var],RadialFunction) and locIndex is None)):
            print "WARNING: you asked me to rank models by ",var," with timeIndex and locIndex",timeIndex,locIndex," but I can't do that."
        else:
            if(qu.ndim==2):
                for i in range(nt):
                    indices = np.argsort(qu[:,i])
                    for j,ind in enumerate(indices):
                        ranks[j,i] = ind
            else:
                indices = np.argsort(qu)
                for i,ind in enumerate(indices):
                    ranks[i,:] = ind
        for i,model in enumerate(self.models):
            model.var['rankBy'+var]=TimeFunction(np.copy(ranks[i,:]), 'rankBy'+var, 1.0, 1.0, 'Ranked by '+var, log=False, theRange=[0,nmodels])
            


    def radialPlot(self,timeIndex=None,variables=None,colorby=None,percentiles=None,logR=False,scaleR=False,movie=True):
        ''' Plot quantities vs radius. '''
        if(variables is None):
            variables = self.models[0].getRadialFunctions()
        if(timeIndex is None):
            timeIndex=range(1,self.models[0].nTimeSteps(),3)
        if(colorby is None):
            colorby='Mh0'

        for i,v in enumerate(variables):

            overallVar,overallFail,overallLog,overallRange = self.constructQuantity(v)
            indVar = 'r'
            if(scaleR):
                indVar='rx'
            allR,_,_,_ = self.constructQuantity(indVar)
            rRange = [np.min(allR),np.max(allR)]
            allRb,_,_,_ = self.constructQuantity('rb')
            if(scaleR):
                rRange=[0,4]

            dirname = 'movie_'+self.name+'_'+v+'_cb'+colorby+'_vs'+indVar
            #dirname = self.name+'_'+v+'_cb'+colorby+'_vs'+indVar
            if(percentiles is not None and len(self.models) > 10):
                dirname +='_perc'
            if(not os.path.exists(dirname)):
                os.makedirs(dirname)
            print "Making movie: ",v
# from SO
#            rate=29.97
#            cmdstring = ('ffmpeg','-r','%d' % rate,
#                    '-y','-an','-f','image2pipe',
#                    #'-loglvel','quiet'
#                    '-vcodec','png',
#                    '-i','pipe:',
#                    dirname+'.avi')
#            pipe = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)


            _,_,logColor,overallColorRange = self.constructQuantity(colorby)
            if(logColor):
                overallColorRange = np.log10(overallColorRange)

            counter = 0 
            plotModels = range(len(self.models))
            if(percentiles is not None and len(self.models) > 10):
                #plotModels = np.random.choice(plotModels,size=10,replace=False)
                #plotModels = [random.choice(plotModels) for qq in range(3)]
                plotModels = []

            if v=='vPhi':
                rotCurve = np.loadtxt('../Bhattacharjee2014.txt',skiprows=15)
                rc = rotCurve[51:101, 2:5]

            if not movie:
                fig,ax = plt.subplots(1,1)
            for ti in timeIndex:
                if movie:
                    fig,ax = plt.subplots(1,1)
                theVar,varFail,varLog,varRange = self.constructQuantity(v,ti)
                r,rFail,rLog,_ = self.constructQuantity(indVar,ti)
                #rx,rxFail,rxLog,_ = self.constructQuantity('rx',ti)

                colors,fail,log,_ = self.constructQuantity(colorby,timeIndex=ti)
                #rgb, scaledColors,origColors = getRGB(colors,log)
                if(log):
                    colors=np.log10(colors)
                
                model = self.models[0]
                ax.set_xlabel(model.get(indVar).texString)
                ax.set_ylabel(model.get(v).texString)
                sc = ax.scatter(r[:,0],theVar[:,0],c=colors,cmap=cm,vmin=overallColorRange[0],vmax=overallColorRange[1],lw=0,s=4)
                if v=='vPhi':
                    ax.errorbar(rc[:,0], rc[:,1], yerr=rc[:,2],color='k',lw=3)
                if v=='colst':
                    ax.errorbar([8.3], [38.0], yerr=[4.0],color='k')
                ## The following is from Nakanishi & Sofue (2006), errorbars included.
                if v=='colH2':
                    vx = np.arange(11)+0.5
                    vxe = np.zeros(11)+0.5
                    vy = [26.7, 6.0, 3.5, 3.5, 4.6, 3.9, 2.8, 1.9, 0.9, 0.5, 0.3]
                    vye = [20.3, 5.4, 1.8, 1.6, 2.0, 2.0, 1.4, 1.0, 0.7, 0.4, 0.2]
                    ax.errorbar(vx,vy,xerr=vxe,yerr=vye,color='k',lw=3,label='Nakanishi06')

                    ## From the Herschel paper - Pineda et al (2013) A&A 554 103
                    vx = np.array([1.2264620876811554, 1.7322710031600081, 2.2387544688954497, 2.7385962090275573, \
                            3.2361029675022435, 3.7372419196662534, 4.222606773522346, 4.721514520991485, \
                            5.227271547989061, 5.7285661655969, 6.228978679023045, 6.736707468309111, \
                            7.245733469627077, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, \
                            12.25, 12.75, 13.25, 13.75, 14.25, 14.75] )
                    vy =np.array( [6.548974943052395, 4.498861047835993, 2.1526195899772205, 2.7220956719817764, \
                            4.316628701594535, 4.316628701594535, 11.241457858769932, 12.220956719817767, \
                            10.19362186788155, 10.125284738041003, 10.444191343963555, 7.551252847380409, \
                            4.088838268792713, 4.7038724373576315, 3.6332574031890665, 3.5193621867881575, \
                            2.5626423690205016, 3.041002277904326, 2.927107061503417, 1.9476082004555835, \
                            2.175398633257405, 1.6059225512528492, 1.287015945330296, 2.1070615034168583, \
                            1.150341685649206, 1.7425968109339394, 1.6514806378132114, 2.835990888382689] )
                    vyp = np.array( [  9.09185803757829, 6.628392484342381, 4.1022964509394555, 4.77035490605428, \
                            6.816283924843425, 7.0459290187891455, 14.749478079331944, 16.29436325678497, \
                            14.373695198329855, 14.436325678496871, 14.937369519832988, 11.784968684759917, \
                            7.839248434237998, 8.298538622129438, 7.3382045929018815, 8.653444676409187, \
                            6.9624217118997915, 7.171189979123174, 7.171189979123174, 6.064718162839251, \
                            6.127348643006265, 4.937369519832988, 4.54070981210856, 5.041753653444678, \
                            3.810020876826723, 4.039665970772443, 3.6638830897703585, 4.979123173277662 ] )
                    ax.errorbar(vx,vy,xerr=vx*0+0.25,yerr=vyp-vy,color='r',lw=3, label=r'Pineda13')
                ## The following HI data are from Nakanishi & Sofue (2003) -- errorbars are not specified
                if v=='colHI':
                    vy = [1.85, .54, 1.84, 1.81, 1.86, 1.5, 2.6, 4.5, 4.25, 4.0, 4.25, 4.31, 3.6, 2.8, 2.0, 1.4, 1.05, .8, .7, .6, .5, .38, .25, .16, .1]
                    vx = np.arange(25)+0.5
                    vxe = vx*0+0.5
                    vy = np.array([1.79665, 0.591, 1.775, 1.746, 1.818, 1.488, 2.608, 4.510, 4.301, 4.114, 4.309, 4.388, 3.691, 2.780, 2.055, 1.431, 1.072, 0.821, 0.641, 0.498, 0.390, 0.311, 0.246, 0.189, 0.160])
                    vye = vy*0.1
                    ax.errorbar(vx,vy,xerr=vxe,yerr=vye,color='k',lw=3, label=r'Nakanishi03')

                    # Following from Pineda+ 2013 A&A (as above). Figure 21
                    vx = np.array([.75+i*.5 for i in range(39)])
                    vxe = vx*0+0.25
                    vy = np.array([ 2.175398633257405, 1.6514806378132114, 1.5375854214123024, 1.5831435079726646, \
                            1.7198177676537618, 2.1070615034168583, 2.357630979498861, 2.9498861047836016, \
                            3.3826879271070602, 3.610478359908882, 3.7243735763097945, 3.7927107061503413, \
                            3.5193621867881575, 3.2915717539863323, 3.0865603644646917, 3.1093394077448764, \
                            4.56719817767654, 3.906605922551254, 3.5876993166287043, 3.5193621867881575, \
                            3.6332574031890665, 3.4738041002277917, 2.995444191343964, 2.927107061503417, \
                            2.6082004555808673, 2.3120728929384953, 1.9476082004555835, 1.6287015945330303, \
                            1.5148063781321177, 1.3097949886104807, 1.3781321184510276, 1.4009111617312051, \
                            1.1958997722095681, 1.3097949886104807, 1.1731207289293835, 1.1958997722095681, \
                            1.3781321184510276, 1.264236902050115, 1.287015945330296 ])
                    vye = vy*0+0.1
                    ax.errorbar(vx,vy,xerr=vxe,yerr=vye,color='r',lw=3, label=r'Pineda13')
                if v=='hGas':
                    vx = np.array([1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0])
                    vy = np.array([103.448, 142.241, 161.638, 183.190, 316.810, 349.138, 441.810, 497.844, 685.345, 862.069, 928.879, 935.345, 939.655])
                    vxe = vx*0+1.0
                    vye = vy*0.1
                    ax.errorbar(vx,vy,xerr=vxe,yerr=vye,color='k',lw=3)
                #if counter==0:
                if movie:
                    cbar = plt.colorbar(sc,ax=ax)
                    if(log):
                        cbar.set_label(r'$\log_{10}$'+colorby)
                    else:
                        cbar.set_label(colorby)
                normColors = (colors - float(overallColorRange[0]))/ float(overallColorRange[1]-overallColorRange[0])
                theRGB = cm(normColors)
                lwnorm = 1.0 + log10(float(len(self.models)))

                if movie:
                    thisCM = cmPer
                else:
                    thisCM = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]
                    thisCM = thisCM[counter]
                # fill in some percentiles
                if(percentiles is not None and len(percentiles) < len(self.models)):
                    try:
                        rr, qVecs = constructQuantiles(allR, [model.var[v] for model in self.models], percentiles, ti)
                    except:
                        rr, qVecs = constructQuantiles(allRb, [model.var[v] for model in self.models], percentiles, ti)
                    qVecs=np.array(qVecs)
                    np.clip(qVecs, overallRange[0], overallRange[1], out=qVecs)
                    if(len(percentiles) % 2 == 0):
                        # a unique fill between each quantile
                        for k in range(len(percentiles)-1):
                            #ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=cmPer(percentiles[k]/100.0))
                            ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(percentiles))),alpha=0.3)
                    else:
                        # symmetric fills around the central quantile.
                        v2 = (len(percentiles)-1)/2
                        for k in range(v2):
                            col = thisCM(float(k+.5)/float(v2+.5))
                            ax.fill_between(rr, qVecs[k,:], qVecs[(k+1),:], facecolor=col, alpha=0.2)
                            ax.fill_between(rr, qVecs[(-k-2),:], qVecs[(-k-1),:], facecolor=col, alpha=0.2)
                            ax.plot( rr, qVecs[k,:], lw=2*k+1, color=thisCM(1.0) )
                            ax.plot( rr, qVecs[(-k-1),:], lw=2*k+1, color=thisCM(1.0) )
                            #print "dbg perc: ",v,k,ti
                            #print "dbg perc: ",rr[0],rr[-1],qVecs[k,0],qVecs[k,-1],qVecs[(k+1),0],qVecs[(k+1),-1]
                            #print "dbg perc: ",rr[0],rr[-1],qVecs[(-k-2),0],qVecs[(-k-2),-1],qVecs[(-k-1),0],qVecs[(-k-1),-1]

                for k in plotModels:
                    model = self.models[k]
                    try:
                        ax.plot(r[k],theVar[k],c=theRGB[k],lw=2.0/lwnorm)
                        if v=='colst':
                            pass
                            #ax.plot(r[k], theVar[k][len(theVar[k])/2] *  np.exp( -r[k] / model.getData('scaleLength',timeIndex=ti)),ls='--')
                        if v=='vPhi':
                            ax.plot(r[k], model.var['vPhiDM'].sensible( timeIndex=ti), ls='--')
                            ax.plot(r[k], model.var['vPhiBulge'].sensible( timeIndex=ti), ls=':')
                    except ValueError:
                        rr = model.getData('rb',ti)
                        if(scaleR):
                            rr/=model.getData('accretionRadius',ti)
                        ax.plot(rr,theVar[k],c=theRGB[k],lw=2.0/lwnorm)


                ax.set_xlim(rRange[0],rRange[1])
                ax.set_ylim(overallRange[0],overallRange[1])
                if(overallLog):
                    ax.set_yscale('log')
                if(logR):
                    ax.set_xscale('log')
                if movie:
                    z=model.getData('z')
                    dispz = "%.3f" % z[ti]
                    plt.text(0.7,0.9,'z='+dispz,transform=ax.transAxes)
                    plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
                    #plt.savefig(pipe.stdin, format='png')
                    plt.close(fig)
                counter = counter+1
            if not movie:
                plt.savefig(dirname[6:]+'.png')
            plt.close(fig)
            ##makeMovies(self.name+'_'+v+'_cb'+colorby+'_vs'+indVar)
            if movie:
                makeMovies(dirname[6:])
            #pipe.terminate()

    def plotAgeFuncs(self):
        ''' Plot properties of the stellar population'''
        dirname = self.name+'_age-velocity-dispersion'

        # First we need to identify the fields we're interested in.
        vars = sorted(self.models[0].var.keys()) # Get everything
        cols = []
        sigRs = []
        sigZs = []
        Zs = []
        VZs = []
        ages =[]
        startingAges=[]
        endingAges=[]
        N=100 # meh.
        for i in range(N):
            sti = str(i).zfill(2)
            if 'colst'+sti in vars:
                cols.append('colst'+sti)
            if 'sigstR'+sti in vars:
                sigRs.append('sigstR'+sti)
            if 'sigstZ'+sti in vars:
                sigZs.append('sigstZ'+sti)
            if 'Zst'+sti in vars:
                Zs.append('Zst'+sti)
            if 'VarZst'+sti in vars:
                VZs.append('VarZst'+sti)
            if 'ageSt'+sti in vars:
                ages.append('ageSt'+sti)
            if 'stargingAgeSt'+sti in vars:
                startingAges.append('startingAgeSt'+sti)
            if 'endingAgeSt'+sti in vars:
                endingAges.append('endingAgeSt'+sti)
        # Now we'll make some plots.

        # First up, the age-velocity dispersion correlation at 8.3 kpc.
        fig,ax = plt.subplots()
        for i in range(len(sigRs)):
            theSigRs = []
            theSigZs = []
            theAges = []
            for j,model in enumerate(self.models):
                r= model.var['r'].sensible(timeIndex=-1)
                theSigRs.append( model.var[ sigRs[i] ].atR([8.3], r, -1)  )
                theSigZs.append( model.var[ sigZs[i] ].atR([8.3], r, -1)  )
                theAges.append( model.var[ ages[i] ].sensible(timeIndex=-1) )
            ax.errorbar( [np.mean(theAges)], [np.mean(theSigRs)], xerr=[np.std(theAges)], yerr=[np.std(theSigRs)], ecolor='blue')
            ax.errorbar( [np.mean(theAges)], [np.mean(theSigZs)], xerr=[np.std(theAges)], yerr=[np.std(theSigZs)], ecolor='red')



        tW = np.array([ 1.1545454471140513, 1.4562262504912922, 1.6174725466389603, 1.7283910560022944, 1.7965734639435877, 1.8777958523969331, 1.9626902734826572, 2.040115430494593, 2.108906312357657, 2.180016761710267, 2.253524992309725, 2.316671738324658, 2.394787853464255, 2.4892586052022443, 2.5590109183314973, 2.6452984955626317, 2.7648914130211915, 2.8898910798240056, 3.003892903558393, 3.2098854713099376, 3.5066815356965195, 3.830920293793729, 4.374347699826242, 5.337385567099643, 6.128288264102687, 7.075387313922233, 8.034520783900259, 9.225089476610762, 10.592078624278916, 12.571709892264574 ])
        sW = np.array([ 10.15947094711839, 10.119366364989414, 11.807187869366663, 9.80417679462689, 11.394270225407617, 12.678456671422186, 11.039371014586864, 12.043018110237995, 10.909152427547689, 12.043018110237995, 13.560276180165056, 15.884716272874831, 13.400321414471529, 15.884716272874831, 13.614017609237196, 14.793111980906472, 12.779149095990569, 13.189979561162446, 15.029019138748827, 15.512179549310623, 14.618614888860135, 16.72285843631423, 17.12447052897269, 20.379765023842037, 20.86920049387745, 24.83635342658302, 24.44650273214034, 27.526432797204478, 30.26749576398466, 30.99439256572452])
        sWp = np.array([10.952387051506143, 10.823194512848927,12.678456671422186,10.48608269127704,12.235069260658351,13.667972023735317,11.900960668746963,12.931689153515185,11.714153946296314,12.931689153515185,14.676550674587839,17.12447052897269    ,14.389149754895975,17.056871562617136,14.676550674587839,15.947669810219324,13.776523188728843,14.219417734176508,16.201988811116035,16.72285843631423,15.759553745770539,18.028027163560893,18.460983870288672,21.970344294667754,22.49797872882931,    26.774755997810324,26.458925422884747,29.674787960072617,32.75909940704369,33.545833954643506])
        
        tU = np.array([ 1.1545454471140513, 1.4481996267037256, 1.6174725466389603, 1.7283910560022944, 1.8065309442770299, 1.8777958523969331, 1.9626902734826572, 2.040115430494593, 2.108906312357657, 2.1920994705262253, 2.253524992309725, 2.316671738324658, 2.394787853464255, 2.5030552821005823, 2.5449058176283663, 2.6452984955626317, 2.7496515125055185, 2.8898910798240056, 3.003892903558393, 3.192192792668354, 3.4873529365723437, 3.830920293793729, 4.350236581591234, 5.337385567099643, 6.128288264102687, 7.036388242108341, 8.079051989720647, 9.17424138716694, 10.650785054440707, 12.571709892264574 ])
        sU = np.array([ 23.03828470297126,26.25044414355292,24.83635342658302,21.118308554196396,26.04360557809312,26.35447863100062,24.640657091469063,26.774755997810324,27.745047990670184,21.455084220864002,30.38745043564325,31.738746343902054,28.523964850804024,23.12958891586919,31.488663079506047,32.24488713788331,33.41341154210631,33.94625838061393,25.94079826418102 ,32.24488713788331,30.38745043564325,36.740695743681265,34.62428217248994,37.326603098484604,39.45184274341137,41.533482242701545,45.13064926230107,48.46090259220439,53.497788341709885,48.26960284262397] )
        sUp = np.array([ 24.83635342658302, 28.299212235621507, 26.774755997810324, 22.676657891822398, 27.965399428100582, 28.411366301434906, 26.458925422884747, 28.864445050736656, 29.910465410861065, 22.947340914087935, 32.75909940704369, 34.080792545905844, 30.50788050584042, 24.83635342658302, 33.94625838061393, 34.761503448799075, 35.87903361647077, 36.4511999637219, 27.855005819617247, 34.761503448799075, 32.75909940704369, 39.45184274341137, 37.326603098484604, 40.239832032642134, 42.53094023006707, 44.775045421823734, 48.65296049167699, 52.243130062368444, 57.673129572180216, 51.831483919616105 ])
        
        tV = np.array([ 1.154781984689458, 1.4489643232788298, 1.6096474622839103, 1.7201919979837987, 1.8080545883896824, 1.879483258284322, 1.9645777937496076, 2.042190040422458, 2.0994976101348772, 2.1824399740999354, 2.2561365673928373, 2.3194478660758224, 2.4110795441388793, 2.5063312063142997, 2.548296747979347, 2.648969287610528, 2.753618977957703, 2.8942661247167507, 3.0253054571024536, 3.1974789370674914, 3.4935682400446213, 3.8382619912089933, 4.3593634327528035, 5.350142932297119,  6.144146161448629, 7.095149298539764, 8.013941241504929, 9.254354972049763, 10.627773930104956, 12.687188271180933 ])
        sV = np.array([ 16.855671538516322,14.618614888860135,14.332348486294974,14.618614888860135,15.268688330380503,15.759553745770546,15.884716272874831,16.330665161755658,13.996218708819354,17.192337400274095,16.656844836483952,14.851739325048687,15.697342797897065,16.460363461207066,19.667049629525607,20.299315759857507,17.60522438531374,20.139368700347838,19.823245526441635,19.901808051959897,18.099474966706556,20.70476306688229,21.118308554196396,24.158135846053955,23.96778347962331, 25.736399532087084,27.309540167088848,28.864445050736656,28.97883921983696,38.223029701586235] )
        sVp = np.array([ 18.099474966706556,15.822011245056265,15.389952399644592,15.697342797897065,16.395386062114245,16.92247311869177,17.056871562617147,17.53572763977108,15.08858141870079,18.460983870288672,17.956861401114423,15.947669810219333, 16.92247311869177,17.60522438531374,21.20200361071367,21.797230480984602,18.979259116958353,21.71118576226263,21.455084220864002,21.455084220864002,19.512084467414365,22.320707455986533,22.676657891822398,26.04360557809312,25.838396782854662,27.745047990670184,29.557646528486323,31.11722803954335,31.364361363287067,41.20622203183994 ])

        ax.errorbar( tU*1.0e9, sU, yerr=sUp-sU, label=r'U from Holmberg09', c='b', ls='--' )
        ax.errorbar( tW*1.0e9, sW, yerr=sWp-sW, label=r'W from Holmberg09', c='r', ls='--' )
        ax.errorbar( tV*1.0e9, sV, yerr=sVp-sV, label=r'V from Holmberg09', c='g', ls='--' )
        ax.set_xlabel( 'Age (yr)' )
        ax.set_ylabel( r'$\sigma_*$ (km/s)' )
        plt.savefig( self.name+'_age-vel-disp.png' )
        plt.close(fig)


        # Next up, the age-metallicity correlation at 8.3 kpc.
        fig,ax = plt.subplots()
        for i in range(len(sigRs)):
            theZs = []
            theAges = []
            for j,model in enumerate(self.models):
                r= model.var['r'].sensible(timeIndex=-1)
                theZs.append( model.var[ Zs[i] ].atR([8.3], r, -1)  )
                theAges.append( model.var[ ages[i] ].sensible(timeIndex=-1) )
            ax.errorbar( [np.mean(theAges)], [np.mean(theZs)], xerr=[np.std(theAges)], yerr=[np.std(theZs)], ecolor='red')
        ax.set_yscale('log')
        ax.set_xlabel( 'Age (yr)' )
        ax.set_ylabel( r'$Z_*$' )
        plt.savefig( self.name+'_age-Z.png' )
        plt.close(fig)

        # simple stuff: radial profiles of each quantity.
        fig,ax = plt.subplots()
        allR,_,_,_ = self.constructQuantity('r')
        ti = -1
        for i in range(len(sigRs)):
            sti = str(i).zfill(2)
            v = 'colst'+sti
            thisCM = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]*20
            thisCM = thisCM[i]
            percentiles = [2.5, 16, 50, 84, 97.5]

            # fill in some percentiles
            rr, qVecs = constructQuantiles(allR, [model.var[v] for model in self.models], percentiles, ti)
            qVecs=np.array(qVecs)
            if(len(percentiles) % 2 == 0):
                # a unique fill between each quantile
                for k in range(len(percentiles)-1):
                    ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(percentiles))),alpha=0.3)
            else:
                # symmetric fills around the central quantile.
                v2 = (len(percentiles)-1)/2
                for k in range(v2):
                    col = thisCM(float(k+.5)/float(v2+.5))
                    ax.fill_between(rr, qVecs[k,:], qVecs[(k+1),:], facecolor=col,alpha=0.3)
                    ax.fill_between(rr, qVecs[(-k-2),:], qVecs[(-k-1),:], facecolor=col,alpha=0.3)
                    ax.plot( rr, qVecs[k+1,:], lw=k, color=thisCM(1.0) )
                    ax.plot( rr, qVecs[(-k-1),:], lw=k, color=thisCM(1.0) )
        ax.set_yscale('log')
        ax.set_ylim(0.5, 3.0e3)
        ax.set_xlabel(r'r (kpc)')
        ax.set_ylabel(r'$\Sigma_{*,i}\ (M_\odot/\mathrm{pc}^2)$')
        plt.savefig( self.name+'_stcols_vsr.png' )



        fig,ax = plt.subplots()
        allR,_,_,_ = self.constructQuantity('r')
        ti = -1
        for i in range(len(sigRs)):
            sti = str(i).zfill(2)
            v = 'sigstR'+sti
            thisCM = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]*20
            thisCM = thisCM[i]
            percentiles = [2.5, 16, 50, 84, 97.5]

            # fill in some percentiles
            rr, qVecs = constructQuantiles(allR, [model.var[v] for model in self.models], percentiles, ti)
            qVecs=np.array(qVecs)
            if(len(percentiles) % 2 == 0):
                # a unique fill between each quantile
                for k in range(len(percentiles)-1):
                    ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(percentiles))),alpha=0.3)
            else:
                # symmetric fills around the central quantile.
                v2 = (len(percentiles)-1)/2
                for k in range(v2):
                    col = thisCM(float(k+.5)/float(v2+.5))
                    ax.fill_between(rr, qVecs[k,:], qVecs[(k+1),:], facecolor=col, alpha=0.3)
                    ax.fill_between(rr, qVecs[(-k-2),:], qVecs[(-k-1),:], facecolor=col, alpha=0.3)
                    ax.plot( rr, qVecs[k+1,:], lw=k, color=thisCM(1.0) )
                    ax.plot( rr, qVecs[(-k-1),:], lw=k, color=thisCM(1.0) )
        ax.set_xlabel(r'r (kpc)')
        ax.set_ylabel(r'$\sigma_{*,R,i}\ (\mathrm{km}/\mathrm{s})$')
        plt.savefig( self.name+'_sigstRs_vsr.png' )



        fig,ax = plt.subplots()
        allR,_,_,_ = self.constructQuantity('r')
        ti = -1
        for i in range(len(sigRs)):
            sti = str(i).zfill(2)
            v = 'sigstZ'+sti
            thisCM = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]*20
            thisCM = thisCM[i]
            percentiles = [2.5, 16, 50, 84, 97.5]

            # fill in some percentiles
            rr, qVecs = constructQuantiles(allR, [model.var[v] for model in self.models], percentiles, ti)
            qVecs=np.array(qVecs)
            if(len(percentiles) % 2 == 0):
                # a unique fill between each quantile
                for k in range(len(percentiles)-1):
                    ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(percentiles))),alpha=0.3)
            else:
                # symmetric fills around the central quantile.
                v2 = (len(percentiles)-1)/2
                for k in range(v2):
                    col = thisCM(float(k+.5)/float(v2+.5))
                    ax.fill_between(rr, qVecs[k,:], qVecs[(k+1),:], facecolor=col, alpha=0.3)
                    ax.fill_between(rr, qVecs[(-k-2),:], qVecs[(-k-1),:], facecolor=col, alpha=0.3)
                    ax.plot( rr, qVecs[k+1,:], lw=k, color=thisCM(1.0) )
                    ax.plot( rr, qVecs[(-k-1),:], lw=k, color=thisCM(1.0) )
        ax.set_xlabel(r'r (kpc)')
        ax.set_ylabel(r'$\sigma_{*,Z,i}\ (\mathrm{km}/\mathrm{s})$')
        plt.savefig( self.name+'_sigstZs_vsr.png' )


        fig,ax = plt.subplots()
        allR,_,_,_ = self.constructQuantity('r')
        ti = -1
        for i in range(len(sigRs)):
            sti = str(i).zfill(2)
            v = 'Zst'+sti
            thisCM = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]*20
            thisCM = thisCM[i]
            percentiles = [2.5, 16, 50, 84, 97.5]

            # fill in some percentiles
            rr, qVecs = constructQuantiles(allR, [model.var[v] for model in self.models], percentiles, ti)
            qVecs=np.array(qVecs)
            if(len(percentiles) % 2 == 0):
                # a unique fill between each quantile
                for k in range(len(percentiles)-1):
                    ax.fill_between(rr, qVecs[k,:], qVecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(percentiles))),alpha=0.3)
            else:
                # symmetric fills around the central quantile.
                v2 = (len(percentiles)-1)/2
                for k in range(v2):
                    col = thisCM(float(k+.5)/float(v2+.5))
                    ax.fill_between(rr, qVecs[k,:], qVecs[(k+1),:], facecolor=col, alpha=0.3)
                    ax.fill_between(rr, qVecs[(-k-2),:], qVecs[(-k-1),:], facecolor=col, alpha=0.3)
                    ax.plot( rr, qVecs[k+1,:], lw=k, color=thisCM(1.0) )
                    ax.plot( rr, qVecs[(-k-1),:], lw=k, color=thisCM(1.0) )
        ax.set_xlabel(r'r (kpc)')
        ax.set_ylabel(r'$Z_{*,i}$')
        ax.set_yscale('log')
        plt.savefig( self.name+'_Zst_vsr.png' )




#def plotPercentiles(

    def hist1d(self,vars=None,timeIndex=None,movie=True):
        model = self.models[0]
        if timeIndex is None:
            timeIndex=[1]
            #timeIndex = range(1,model.nTimeSteps(),3)
        if vars is None:
            vars =  model.p.keys()

        for i,var in enumerate(vars):

            print "Making histogram for ",var
            overallX,_,overallXLog,overallXRange = self.constructQuantity(var)
            if(len(np.shape(overallX))==1):
                overallX = np.tile(overallX, (self.models[0].nTimeSteps(),1)).T
            dirname = 'movie_'+self.name+'_hist_'+var
            if not movie:
                fig,ax = plt.subplots()
            counter=0
            colors=['b','g','r','orange','purple','k']*10
            for ti in timeIndex:
                if movie:
                    fig,ax = plt.subplots()
                try:
                    if overallXLog:
                        h,e = np.histogram(np.log10(overallX[:,ti]), bins=80, density=True, range=overallXRange)
                    else:
                        h,e = np.histogram(overallX[:,ti], bins=80, density=True, range=overallXRange)
                except:
                    pdb.set_trace()

                width = 0.3 * (e[1] - e[0])
                center = (e[:-1] + e[1:]) / 2
                ax.bar(center, h, align='center', width=width, color=colors[counter])

                counter+=1
                #ax.set_xlim(e[0],e[-1])
                ax.set_ylim(0,np.max(h))
                try:
                    ax.set_xlabel(model.get(var).texString)
                except:
                    ax.set_xlabel(var)
                if(overallXLog):
                    try:
                        ax.set_xlabel(r'$\log_{10}$' + model.get(var).texString)
                    except:
                        ax.set_xlabel(r'$\log_{10}$' +var)
                if movie:
                    dispz = "%.3f" % z[ti]
                    plt.text(0.7,0.9,'z='+dispz,transform=ax.transAxes)
                    plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
                    plt.close(fig)
            if movie: 
                makeMovies(dirname[7:])
            else:
                plt.savefig(dirname[6:]+'.png')
                plt.close(fig)


    def angularMomentumAnalysis(self):
        # Barro movie
        nts = int(self.models[0].p['Noutputs']+1)
        stepsize=1
        cb = 'z'
        self.ptMovie(timeIndex=range(1,nts+1,stepsize)+[nts],xvar='Sigma1',yvar=['sSFR'],prev=0,colorby=cb,movie=False)
        self.ptMovie(timeIndex=range(1,nts+1,stepsize)+[nts],xvar='specificJAccr',yvar=['specificJSFR','specificJStars','specificJOut','specificJGas'],prev=0,colorby=cb,movie=False)
        self.ptMovie(timeIndex=range(1,nts+1,stepsize)+[nts],xvar='dimensionlessSpinAccr',yvar=['dimensionlessSpinSFR','dimensionlessSpinStars','dimensionlessSpinOut','dimensionlessSpinGas'],prev=0,colorby=cb,movie=False)
        cb = 'dimensionlessSpinAccr'
        self.ptMovie(timeIndex=range(1,nts/2+1,stepsize),xvar='Sigma1',yvar=['sSFR'],prev=0,colorby=cb,movie=True)
        self.ptMovie(timeIndex=range(1,nts/2+1,stepsize),xvar='specificJAccr',yvar=['specificJSFR','specificJStars','specificJOut','specificJGas'],prev=0,colorby=cb,movie=True)
        self.ptMovie(timeIndex=range(1,nts/2+1,stepsize),xvar='dimensionlessSpinAccr',yvar=['dimensionlessSpinSFR','dimensionlessSpinStars','dimensionlessSpinOut','dimensionlessSpinGas'],prev=0,colorby=cb,movie=True)
        
        # Just plot all components of the angular momentum normalized to the instantaneous accretion ang. mom. on the same plt vs. z.
        fig,ax = plt.subplots()
        ls=['-','--','-.']*30
        for i in range(len(self.models)):
            if i==0:
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinAccr'].sensible(), c='k', lw=3, ls=ls[i], label=r'Accretion')
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinGas'].sensible(), c='grey', lw=3, ls=ls[i], label=r'Gas')
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinStars'].sensible(), c='r', lw=3, ls=ls[i], label=r'Stars')
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinSFR'].sensible(), c='blue', lw=3, ls=ls[i], label=r'SFR')
            else:
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinAccr'].sensible(), c='k', lw=3, ls=ls[i])
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinGas'].sensible(), c='grey', lw=3, ls=ls[i])
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinStars'].sensible(), c='r', lw=3, ls=ls[i])
                ax.plot( self.models[i].var['z'].sensible(),  self.models[i].var['dimensionlessSpinSFR'].sensible(), c='blue', lw=3, ls=ls[i])
        ax.set_xlabel(r'z') 
        ax.set_ylabel(r'$\lambda$')
        ax.set_yscale('log')
        ax.legend()
        ax.set_ylim(self.models[0].var['dimensionlessSpinAccr'].theRange[0], self.models[0].var['dimensionlessSpinAccr'].theRange[1])
        plt.savefig( self.name+'_angularMomentumEvolution.png' )

        import balanceplot
        balanceplot.budgetAM(self.models, name=self.name, sortby='Mh0', ncols=5, nrows=3)
        balanceplot.balanceAM(self.models, name=self.name, sortby='Mh0', ncols=5, nrows=3)


    def ptMovie(self,xvar='mstar',yvar=None,colorby=None,timeIndex=None,prev=0,movie=True):
        if(yvar is None):
            variables = self.models[0].getTimeFunctions()
        else:
            variables = yvar
        if(timeIndex is None):
            timeIndex=range(1,self.models[0].nTimeSteps(),3)
        if(colorby is None):
            colorby='Mh0'


        overallX,_,overallXLog,overallXRange = self.constructQuantity(xvar)
        colors,_,log,overallColorRange = self.constructQuantity(colorby)
        # Make sure that the X variable is a 2-d array (model x time). If it's e.g. a 
        # parameter, it will only have one dimension indexing the number of models.
        if(len(np.shape(overallX))==1):
            overallX = np.tile(overallX, (self.models[0].nTimeSteps(),1)).T
        if(log):
            colors=np.log10(colors)
            overallColorRange = np.log10(overallColorRange)
        for i,v in enumerate(variables):
            if((xvar=='Mh' and v=='efficiency') or (xvar=='Mh' and v=='mstar')):
                b = behroozi()
                b.readSmmr()
            dirname = 'movie_'+self.name+'_'+v+'_cb'+colorby+'_vs'+xvar
            if(not os.path.exists(dirname) and movie):
                os.makedirs(dirname)
            print "Making movie : ",v," vs ",xvar

            overallVar,_,overallLog,overallRange = self.constructQuantity(v)


            cmaps = [plt.cm.Blues, plt.cm.Greens, plt.cm.Reds, plt.cm.Purples, plt.cm.Greys]*10
            counter=0
            if not movie:
                fig,ax = plt.subplots(1,1)

            # As we did above for the X variable, make sure that the y-variable conforms to later
            # expectations, i.e. a 2-d array (model x time)
            if(len(np.shape(overallVar))==1):
                overallVar = np.tile(overallVar, (self.models[0].nTimeSteps(),1)).T
                
            for ti in timeIndex:
                if movie:
                    fig,ax = plt.subplots(1,1)
                model = self.models[0]
                try:
                    ax.set_xlabel(model.get(xvar).texString)
                except:
                    ax.set_xlabel(xvar)
                try:
                    ax.set_ylabel(model.get(v).texString)
                except:
                    ax.set_ylabel(v)
                
                colorLoc,_,_,_ = self.constructQuantity(colorby,timeIndex=ti)
                if(log):
                    colorLoc = np.log10(colorLoc)

                if movie:
                    lw=1
                    s=30
                else:
                    lw=0
                    s=1
                nx = int(model.p['nx'])
                if v == 'Mdot':
                    nx+=1
                if(len(np.shape(overallX[:,ti]))==2 and len(np.shape(overallVar[:,ti]))==2):
                    # X and Y are 2D, so we have to flatten them, and adjust colorLoc.
                    colorLoc = np.tile(colorLoc, (nx, 1)).T.flatten()
                    xx = overallX[:,ti].flatten()
                    yy = overallVar[:,ti].flatten()
                elif(len(np.shape(overallX[:,ti]))==1 and len(np.shape(overallVar[:,ti]))==2):
                    # Y is 2D but X is not, so we have to flatten them, and adjust colorLoc.
                    colorLoc = np.tile(colorLoc, (nx, 1)).T.flatten()
                    xx = np.tile(overallX[:,ti], (nx,1)).T.flatten()
                    yy = overallVar[:,ti].flatten()
                elif(len(np.shape(overallX[:,ti]))==2 and len(np.shape(overallVar[:,ti]))==1):
                    # X is 2D but Y is not, so we have to flatten them, and adjust colorLoc.
                    colorLoc = np.tile(colorLoc, (nx, 1)).T.flatten()
                    xx = overallX[:,ti].flatten()
                    yy = np.tile(overallVar[:,ti], (nx,1)).T.flatten()
                else:
                    xx=overallX[:,ti]
                    yy=overallVar[:,ti]

                if np.max(colorLoc)-np.min(colorLoc)<(overallColorRange[1]-overallColorRange[0])*0.1 and len(xx)>200:
                    if(overallLog):
                        yscale = 'log'
                    else:
                        yscale = 'linear'
                    if(overallXLog):
                        xscale = 'log'
                    else:
                        xscale = 'linear'

                    histFlag = True
                    if overallLog:
                        if counter==0:
                            try:
                                assert overallRange[0]>0 and overallRange[1]>0
                            except:
                                pdb.set_trace()
                            overallRange[0] = np.log10(overallRange[0])
                            overallRange[1] = np.log10(overallRange[1])
                        yy= np.log10(yy)
                    if overallXLog:
                        if counter==0 and i==0:
                            try:
                                assert overallXRange[0]>0 and overallXRange[1]>0
                            except:
                                pdb.set_trace()
                            overallXRange[0] = np.log10(overallXRange[0])
                            overallXRange[1] = np.log10(overallXRange[1])
                        xx = np.log10(xx)
                    weights=None
                    if xvar=='r' and not overallXLog:
                        weights=copy.deepcopy(xx)
                    try:
                        h, xe, ye = np.histogram2d(xx,yy, bins=30, range=[[overallXRange[0],overallXRange[1]],[overallRange[0],overallRange[1]]] , weights=weights, normed=True) 
                    except:
                        pdb.set_trace()
                    print "hist min, max: ",np.min(h),np.max(h)
                    #ax.imshow(np.log10(h.T+1.0), cmap=cmaps[counter], origin='lower', extent=( overallXRange[0],overallXRange[1],overallRange[0],overallRange[1] ), interpolation='gaussian', alpha=.8, aspect='auto')
                    ax.contour(np.log10(1.0+h.T), extent=( overallXRange[0],overallXRange[1],overallRange[0],overallRange[1] ), cmap=cmaps[counter]  )
                    #plt.hexbin( xx, yy, extent=, alpha = 1.0/float(len(timeIndex)), xscale='linear', yscale='linear', cmap=cmaps[counter] )
                else:
                    histFlag=False
                    sc = ax.scatter(xx,yy,c=colorLoc,vmin=overallColorRange[0],vmax=overallColorRange[1],cmap=cm,lw=lw,s=s)
                    if movie or (not movie and counter==0):
                        cbar = plt.colorbar(sc,ax=ax)
                        if(log):
                            cbar.set_label(r'$\log_{10}$'+colorby)
                        else:
                            cbar.set_label(colorby)
                    if(prev>0 and movie):
                        for k in range(max(ti-prev,0),ti):
                            #ax.scatter(overallX[:,k],overallVar[:,k],c=colorLoc,cmap=cm,s=6,lw=0,vmin=overallColorRange[0],vmax=overallColorRange[1])
                            ax.scatter(overallX[:,k],overallVar[:,k],c=colorLoc,s=6,lw=0,vmin=overallColorRange[0],vmax=overallColorRange[1],cmap=cm)
                z=model.getData('z')
                if(xvar=='Mh' and v=='efficiency'):
                    b.plotSmmr(z[ti], ax, minMh=overallXRange[0], maxMh=overallXRange[1])
                elif (xvar=='Mh' and v=='mstar'):
                    b.plotSmmr(z[ti], ax, minMh=overallXRange[0], maxMh=overallXRange[1],eff=False)
                if('specificJ' in xvar and 'specficJ' in v or 'dimensionlessSpin' in xvar and 'dimensionlessSpin' in v):
                    ax.plot([overallXRange[0],overallXRange[0]], [overallRange[0],overallRange[1]], c='k',ls='--')

                # Moster
                if not histFlag:
                    def Moster(Mh, mparams):
                        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
                        logM1z = M10 + M11*z[ti]/(z[ti]+1.0)
                        Nz = N10 + N11*z[ti]/(z[ti]+1.0)
                        betaz = beta10 + beta11*z[ti]/(z[ti]+1.0)
                        gammaz = gamma10 + gamma11*z[ti]/(z[ti]+1.0)
                        M1 = np.power(10.0, logM1z)
                        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
                        return eff
                    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329])
                    unc =     np.array([0.236, 0.353, 0.0058, 0.0069, 0.153, 0.225, 0.059, 0.173])
                    logMhs = np.linspace(np.log10(overallXRange[0]), np.log10(overallXRange[1]), num=100)
                    Mhs = np.power(10.0, logMhs)
                    eff = Moster(Mhs,central)
                    for i in range(len(unc)):
                        theseParams = copy.copy(central)
                        theseParams[i] = theseParams[i]+unc[i]
                        eff = np.vstack([eff, Moster(Mhs, theseParams)])
                        theseParams = copy.copy(central)
                        theseParams[i] = theseParams[i]-unc[i]
                        eff = np.vstack([eff, Moster(Mhs, theseParams)])
                    effM = np.min(eff, axis=0)
                    effP = np.max(eff, axis=0)
                    if(xvar=='Mh' and (v=='efficiency' or v=='mstar')):
                        if v=='efficiency':
                            ax.fill_between(Mhs, effM, effP, alpha=.1)
                        if v=='mstar':
                            ax.fill_between(Mhs, effM*Mhs, effP*Mhs,alpha=.1)


                ax.set_xlim(overallXRange[0],overallXRange[1])
                ax.set_ylim(overallRange[0],overallRange[1])
                if(overallLog and not histFlag):
                    ax.set_yscale('log')
                if(overallXLog and not histFlag):
                    ax.set_xscale('log')
                if movie:
                    dispz = "%.3f" % z[ti]
                    plt.text(0.7,0.9,'z='+dispz,transform=ax.transAxes)
                    plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
                    plt.close(fig)
                counter+=1
            if movie: 
                makeMovies(dirname[7:])
            else:
                plt.savefig(dirname[6:]+'.png')

    def timePlot(self,variables=None,colorby=None,perc=None,vsz=False):
        if(variables is None):
            variables = self.models[0].getTimeFunctions()
        time,fail,log,_ = self.constructQuantity('t')
        z,fail,log,_ = self.constructQuantity('z')
        if(colorby is None):
            colorby = 'Mh0'

        for i,v in enumerate(variables):
            fig,ax=plt.subplots(1,1)
            model = self.models[0]

            theVar,varFail,varLog,varRange = self.constructQuantity(v)
            #sc=plt.scatter(model.getData('t'),model.getData(v),c=colors,cmap=cm)
            plt.cla()

            colors,fail,log,overallColorRange = self.constructQuantity(colorby)
            if(colors.ndim==1):
                colorsFlat = np.copy(colors[:])
                for j in range(self.models[0].nTimeSteps()-1):
                    colors = np.vstack((colors,colorsFlat))
                colors = colors.T
            if(log):
                colors=np.log10(colors)
                overallColorRange=np.log10(overallColorRange)

            if vsz:
                ax.set_xlabel(model.get('z').texString)
            else:
                ax.set_xlabel(model.get('t').texString)
            ax.set_ylabel(model.get(v).texString)
            z = model.getData('z')
            t = model.getData('t')
            zz=0.0
            zs = [zz]
            while zz<np.max(z):
                if zz>1.5:
                    zz=zz+1.0
                else:
                    zz=zz+0.5
                zs.append(zz)
            #zs = [0.0,0.5,1.0,1.5,2.0]
            zind = [Nearest(z,zl)[0] for zl in zs]
            correspondingTs = [t[zi] for zi in zind]
            lwnorm = 1.0 + log10(float(len(self.models)))

            if vsz:
                xc = [.02]
            else:
                xc = [t[-1]*.99]
            if v=='BT':
                ax.errorbar(xc, [0.150 + (0.028-0.019)/2.0], yerr=[.028],color='k',lw=3)
            if v=='sfr':
                ax.errorbar(xc, [1.65], yerr=[0.19],color='k',lw=3)
            if v=='mstar':
                ax.errorbar(xc, [6.08e10], yerr=[1.14e10],color='k',lw=3)
            if v=='scaleLength':
                ax.errorbar(xc, [2.15], yerr=[0.14],color='k',lw=3)

            allData = np.zeros( (len(self.models), self.models[0].nt) )
            if perc is not None:
                qvecs = np.percentile(np.array(theVar), perc, axis=0)
                qvecs = np.array(qvecs)

            if perc is None:
                if len(theVar) != len(self.models):
                    pdb.set_trace() # Something has gone wrong and we don't have all the data we expected to have!
                for j,model in enumerate(self.models):
                    #ax.plot(model.getData('t'),theVar[j],c=scalarMap.to_rgba(colors[j]),lw=2.0/lwnorm)
                    if(len(model.getData('t')) != len(theVar[j])):
                        # pdb.set_trace() -- looks like this situation occurs when bulge fraction calculation fails.
                        pass
                    else:
                        if vsz:
                            xx = model.getData('z')
                        else:
                            xx = model.getData('t')
                        sc = ax.scatter(xx,
                                    theVar[j],
                                    c=(colors[j]),
                                    lw=0.1/lwnorm,
                                    s=10.0/lwnorm,
                                    vmin=overallColorRange[0],
                                    vmax=overallColorRange[1],
                                    cmap=cm)
                cbar = plt.colorbar(sc,ax=ax)
                if(log):
                    cbar.set_label(r'$\log_{10}$'+colorby)
                else:
                    cbar.set_label(colorby)

            if vsz:
                xx = model.getData('z')
            else:
                xx = model.getData('t')

            # fill in some percentiles
            if(perc is not None):
                np.clip(qvecs,  varRange[0], varRange[1], out=qvecs )
                if(len(perc) < len(self.models)):
                    if(len(perc) % 2 == 0):
                        # a unique fill between each quantile
                        for k in range(len(perc)-1):
                            ax.fill_between(xx, qvecs[k,:], qvecs[k+1,:], facecolor=cmPer(float(k+.5)/float(len(perc))),alpha=0.3)
                    else:
                        # symmetric fills around the central quantile.
                        v2 = (len(perc)-1)/2
                        for k in range(v2):
                            col = cmPer(float(k+.5)/float(v2+.5))
                            #pdb.set_trace()
                            ax.fill_between(xx, qvecs[k,:], qvecs[(k+1),:], facecolor=col, alpha=0.3)
                            ax.fill_between(xx, qvecs[(-k-2),:], qvecs[(-k-1),:], facecolor=col, alpha=0.3)
                            ax.plot( xx, qvecs[k,:], lw=k+1, color='k' )
                            ax.plot( xx, qvecs[(-k-1),:], lw=k+1, color='k' )

                    

            if vsz:
                ax.set_xlim(zs[0],zs[-1])
                ax.invert_xaxis()
            else:
                ax2 = ax.twiny()
                ax2.set_xticks(correspondingTs)
                ax2.set_xticklabels([str(zl) for zl in zs] , rotation=50, ha='center')
                ax2.set_xlabel(model.get('z').texString)

                ax.set_xlim(correspondingTs[-1],correspondingTs[0])

            ax.set_ylim(varRange[0],varRange[1])
            #### Take this out for now. Hopefully it's fine, and will let us see the errorbars when it's available.
            #try:
            #    ax.set_ylim(varRange[0],varRange[1])
            #except:
            #    pdb.set_trace()

            #if model.get(v).log:
            #    ax.set_yscale('log')
            if varLog:
                ax.set_yscale('log')
            ## ax2.set_xlim(zs[-1],zs[0])

            #scalarMap.set_array(colors)
            #scalarMap.autoscale()
            print "Making plot v: ",v

            #plt.colorbar(sc,ax=ax)
            try:
                if vsz:
                    plt.savefig(self.name+'_vsz_'+v+'_cb'+colorby+'.png')
                else:
                    plt.savefig(self.name+'_vst_'+v+'_cb'+colorby+'.png')

            except ValueError:
                print "Caught ValueError while trying to plot variable ",v," colored by ",colorby
            plt.close(fig)
        #### End of timePlot(self,variables=None,colorby=None,perc=None):
    def customPlotPPD(self, expensive=False):
        ''' Make a custom plot where we assume the experiment is composed of a sample from the posterior
            predictive distribution.'''
        # hg vs MLF
        # r vs mInterior
        # r vs rhoDM
        # vcirc decomp??
        # column density decomp
        # Mh vs other measurements. Ex. Guglielmo+ MNRASS 444, 1759-1774 (2014)
        fig,ax = plt.subplots()
        muNorms,_,_,_ = self.constructQuantity('muNorm')
        muFgScaling,_,_,_ = self.constructQuantity('muFgScaling')
        muColScaling,_,_,_ = self.constructQuantity('muColScaling')
        fgs = np.power(10.0, np.linspace(-2,0,50))
        perc=[2.5, 16.0, 50.0, 86.0, 97.5]
        cms = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')] * 5
        for counter, col in enumerate([1.0, 10.0, 100.0, 1000.0]):
            mus = []
            for i in range(len(muNorms)):
                mus.append( muNorms[i] * np.power(fgs, muFgScaling[i]) * np.power(col, muColScaling[i]) )
            # mus ~ models x fgs
            # qvecs ~ percentiles x fgs
            qvecs = np.percentile(np.array(mus), perc, axis=0)
            qvecs = np.array(qvecs)
            xx = fgs
            # fill in some percentiles
            #np.clip(qvecs,  varRange[0], varRange[1], out=qvecs )
            thisCM = cms[counter]
            if(len(perc) < len(self.models)):
                filledPlots(qvecs, perc, thisCM, ax, xx)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'Gas Fraction')
        ax.set_ylabel(r'Local mass loading factor ($\dot{\Sigma}_\mathrm{out}/\dot{\Sigma}_*$)')
        plt.savefig(self.name+'_fg_vs_mlf.png')
        plt.close(fig)

        # rho(r) and mInterior(r) for the DM halo
        concentrationRandomFactor,_,_,_ = self.constructQuantity('concentrationRandomFactor')
        Mh,_,_,_ = self.constructQuantity('Mh',timeIndex=-1)
        halos = []
        cos = halo.Cosmology()
        rs = np.power(10.0, np.linspace(-1,2.5,50))
        rhos = np.zeros((len(Mh), len(rs)))
        mInteriors = np.zeros((len(Mh), len(rs)))
        for i in range(len(Mh)):
            thisHalo = halo.halo(Mh[i],0,cos,concentrationRandomFactor[i] ) 
            for j in range(len(rs)):
                rhos[i,j] = thisHalo.rho(rs[j])
                mInteriors[i,j] = thisHalo.mInterior(rs[j])
        fig,ax = plt.subplots()
        qvecs = np.percentile( rhos, perc, axis=0 )
        qvecs = np.array(qvecs)
        thisCM = cms[0]
        xx = rs
        if(len(perc) < len(self.models)):
            filledPlots(qvecs, perc, thisCM, ax, xx)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel(r'$\rho_{\rm DM}\ ({\rm g}\ {\rm cm}^{-3})$')
        ax.set_xlabel(r'$r\ ({\rm kpc})$')
        plt.savefig(self.name+'_r_vs_rhodm.png')
        plt.close(fig)

        fig,ax = plt.subplots()
        qvecs = np.percentile( mInteriors, perc, axis=0 )
        qvecs = np.array(qvecs)
        if(len(perc) < len(self.models)):
            filledPlots(qvecs, perc, thisCM, ax, xx)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel(r'$M_{\rm DM}(r)\ (M_\odot)$')
        ax.set_xlabel(r'$r\ ({\rm kpc})$')
        plt.savefig(self.name+'_r_vs_rhodm.png')
        plt.close(fig)




        fig,ax = plt.subplots()
        cols = np.power(10.0, np.linspace(-1,3,50))
        perc=[2.5, 16.0, 50.0, 86.0, 97.5]
        cms = [plt.get_cmap('Reds'),plt.get_cmap('Greens'),plt.get_cmap('Blues')]*5
        for counter, fg in enumerate([.05,.1,.2]):
            mus = []
            for i in range(len(muNorms)):
                mus.append( muNorms[i] * np.power(fg, muFgScaling[i]) * np.power(cols, muColScaling[i]) )
            qvecs = np.percentile(np.array(mus), perc, axis=0)
            qvecs = np.array(qvecs)
            xx = cols
            # fill in some percentiles
            #np.clip(qvecs,  varRange[0], varRange[1], out=qvecs )
            thisCM = cms[counter]
            if(perc is not None):
                if(len(perc) < len(self.models)):
                    if(len(perc) % 2 == 0):
                        # a unique fill between each quantile
                        for k in range(len(perc)-1):
                            ax.fill_between(xx, qvecs[k,:], qvecs[k+1,:], facecolor=thisCM(float(k+.5)/float(len(perc))),alpha=0.3)
                    else:
                        # symmetric fills around the central quantile.
                        v2 = (len(perc)-1)/2
                        for k in range(v2):
                            col = thisCM(float(k+.5)/float(v2+.5))
                            #pdb.set_trace()
                            ax.fill_between(xx, qvecs[k,:], qvecs[(k+1),:], facecolor=col, alpha=0.3)
                            ax.fill_between(xx, qvecs[(-k-2),:], qvecs[(-k-1),:], facecolor=col, alpha=0.3)
                            ax.plot( xx, qvecs[k,:], lw=k+1, color='k' )
                            ax.plot( xx, qvecs[(-k-1),:], lw=k+1, color='k' )

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'Gas Column Density ($M_\odot/\mathrm{pc}^2$)')
        ax.set_ylabel(r'Local mass loading factor ($\dot{\Sigma}_\mathrm{out}/\dot{\Sigma}_*$)')
        plt.savefig(self.name+'_col_vs_mlf.png')
        plt.close(fig)


        if expensive:
            # Alright, here we're going to do something pretty expensive: reproduce the Z_* vs. r @ different heights from JYC's paper.
            fig,ax = plt.subplots(4)
            fig.subplots_adjust(hspace=0.04)
            rs = [6,9,12,15,18]
            whichModels = []
            # This is such an expensive operation that we only want to includ a random subset of the full PPD sample.
            for i in range(len(self.models)):
                # randomly select ~ 1% of the models
                if np.random.uniform(0,1) < 0.02: 
                    whichModels.append(i)
            theseModels = [self.models[wi] for wi in whichModels]
            allZ = np.zeros( (len(rs), 4 ,len(theseModels))  )  ## r x z x models
            slopes = np.zeros( (4, len(theseModels)) ) # z x models

            for i,model in enumerate(theseModels):
                rr = model.var['r'].sensible(timeIndex=-1)
                # For each of these we need to collect the Sigst's, sigst's, and Zst's.
                # With the first two we compute the weights for the different Z bins.
                colData = []
                sigData = []
                zData = []
                # search for passive stellar populations (up to 100)
                for k in range(100):
                    stk=str(k).zfill(2)
                    if 'colst'+stk in model.var: 
                        # for each stellar population found in this run of the model, 
                        # collect surface density, velocity dispersion, and metallicity at the specified radii (rs)
                        colData.append( model.var['colst'+stk].atR(rs, rr, -1) )
                        sigData.append( model.var['sigstZ'+stk].atR(rs, rr, -1) )
                        zData.append( model.var['Zst'+stk].atR(rs,rr,-1) )
                # Convert these to numpy arrays
                colData = np.array(colData) #  nStellarPops  x len(rs) 
                sigData = np.array(sigData) 
                zData = np.array(zData)

                # For each of our pre-specified radii...
                for j,r in enumerate(rs):
                    try:
                        col = colData[:,j]
                    except:
                        pdb.set_trace()
                    assert not np.any(col<0)
                    sigsq = np.power(sigData[:,j], 2.0)
                    # Use the Narayan and Jog (2002) formalism to compute vertical profiles assuming equilibrium.
                    y,zs = verticalProfile.run( col, sigsq, r, model.getData('Mh',timeIndex=-1), 0.0, model.get('concentrationRandomFactor') )
                    # From these vertical profiles extract a weight matrix. This matrix has a shape
                    #   components  x z
                    #  For each z window (second index), the weight of each component (first index) must sum to 1.
                    weightMatrix = verticalProfile.weights(y,zs)

                    # allZ  ~ r x z x models
                    for m in range(4):
                        # allZ is the weighted average of the metallicity at this radius and height.
                        allZ[j,m,i] =  np.sum(weightMatrix[:,m] * zData[:,j])   # one per height bin.
            ## compute the slope of each ppd model.
            for i in range(len(theseModels)):
                for j in range(4):
                    polyfitResult = np.polyfit(rs, np.log10(allZ[:,j,i]), 1)
                    # slopes ~ z x model
                    slopes[j,i] = polyfitResult[0]
            # Now, at each r and z, identify quantiles of the [average] metallicity distribution across our ensemble of models
            # qvec ~ r x z x quantiles
            qvec = np.percentile( allZ, perc, axis=2 )
            shp2 = np.shape(qvec)
            assert shp2[0] == len(perc)
            assert shp2[1] == len(rs)
            assert shp2[2] == 4
            # Produce the plots
            for m in range(4):
                # one for each window in z

                if(len(perc) < len(self.models)):
                    if(len(perc) % 2 == 0):
                        # a unique fill between each quantile
                        for kk in range(len(perc)-1):
                            ax[3-m].fill_between(rs, qvec[kk,:,m], qvec[kk+1,:,m], facecolor=thisCM(float(kk+.5)/float(len(perc))),alpha=0.3)
                    else:
                        # symmetric fills around the central quantile.
                        v2 = (len(perc)-1)/2
                        for kk in range(v2):
                            col = thisCM(float(kk+.5)/float(v2+.5))
                            #pdb.set_trace()
                            ax[3-m].fill_between(rs, qvec[kk,:,m], qvec[(kk+1),:,m], facecolor=col, alpha=0.3)
                            ax[3-m].fill_between(rs, qvec[(-kk-2),:,m], qvec[(-kk-1),:,m,], facecolor=col, alpha=0.3)
                            ax[3-m].plot( rs, qvec[kk,:,m], lw=kk+1, color='k' )
                            ax[3-m].plot( rs, qvec[(-kk-1),:,m], lw=kk+1, color='k' )
                ax[3-m].set_yscale('log')
                ax[3-m].set_ylabel(r'$\log_{10}(Z_*/Z_\odot)$')
                #ax[3-m].set_ylim(-1.8,0.5)

            # Now overplot data extracted from Judy's paper
            ax[3].set_xlabel(r'r (kpc)')
            jycR = [8.073, 9.103, 9.485]
            jycZ = np.power(10.0, [-0.158385, -0.2329193, -0.24534])
            ax[3].plot(jycR, jycZ, lw=2, c='r')
            jycR = [ 7.27346, 8.302845, 8.967, 9.3, 9.98, 10.4446, 11.125]
            jycZ = np.power(10.0, [-0.175, -0.2375, -0.2875, -0.3, -0.35, -0.375, -0.425])
            ax[2].plot(jycR, jycZ, lw=2, c='r')
            jycR = [6.495, 7.542, 8.7375, 9.0864, 9.6179, 10.199, 10.98, 11.5615, 12.19269, 13.35548]
            jycZ = np.power(10.0, [-0.35124, -0.36505, -0.37905, -0.379485, -0.3801495, -0.3933762, -0.40685, -0.42, -0.4333679, -0.43482])
            ax[1].plot(jycR, jycZ, lw=2, c='r')
            jycR = [6.196, 7.973, 9.40, 10.465, 11.7275, 12.392, 12.9734, 14.0033, 16.1462]
            jycZ = np.power(10.0, [-0.5375, -0.5375, -0.5375, -0.5375, -0.525, -0.525, -0.525, -0.525, -0.5125])
            ax[0].plot(jycR, jycZ, lw=2, c='r')

            ax[0].text(14.0,1.3, r'$1.0\ {\rm kpc} < |z| < 1.5\ {\rm kpc}$')
            ax[1].text(14.0,1.3, r'$0.50\ {\rm kpc} < |z| < 1.0\ {\rm kpc}$')
            ax[2].text(14.0,1.3, r'$0.25\ {\rm kpc} < |z| < 0.50\ {\rm kpc}$')
            ax[3].text(14.0,1.3, r'$0.15\ {\rm kpc} < |z| < 0.25\ {\rm kpc}$')
            for i in range(4):
                ax[i].set_xlim(6,18)
                ax[i].set_ylim(10.0**-1.2, 10.0**0.4)
                if i!=3:
                    ax[i].get_xaxis().set_ticks([])

            plt.savefig(self.name+'_Zsti_r_z.png')
            plt.close(fig)



            fig,ax = plt.subplots()
            ## slopes ~ z x models
            qvec = np.percentile( slopes, perc, axis=1 )
            # Produce the plots
            zs = [0.2, .75/2, .75, 1.25]
            shp1 = np.shape(qvec)
            assert shp1[0] == len(perc)
            assert shp1[1] == len(zs)
            if(len(perc) < len(self.models)):
                if(len(perc) % 2 == 0):
                    # a unique fill between each quantile
                    for kk in range(len(perc)-1):
                        ax.fill_between(zs, qvec[:,kk], qvec[:,kk+1], facecolor=thisCM(float(kk+.5)/float(len(perc))),alpha=0.3)
                else:
                    # symmetric fills around the central quantile.
                    v2 = (len(perc)-1)/2
                    for kk in range(v2):
                        col = thisCM(float(kk+.5)/float(v2+.5))
                        ax.fill_between(zs, qvec[kk,:], qvec[(kk+1),:], facecolor=col, alpha=0.3)
                        ax.fill_between(zs, qvec[(-kk-2),:], qvec[(-kk-1),:], facecolor=col, alpha=0.3)
                        ax.plot( zs, qvec[kk,:], lw=kk+1, color='k' )
                        ax.plot( zs, qvec[(-kk-1),:], lw=kk+1, color='k' )
            ax.set_ylabel(r'$\partial \log_{10}Z/\partial r\ ({\rm dex}/{\rm kpc})$')
            ax.errorbar(zs, [-0.066, -0.064, -0.013, 0.0028], yerr=[[0.044,0.0041,0.0,0.0032],[0.030,0.015,0.0093,0.0071]], color='r', lw=2)
            ax.set_xlabel(r'$z\ ({\rm kpc})$')
            #ax.set_ylim(-1.8,0.5)
            ax.set_xlim(0.19, 1.26)
            
            plt.savefig(self.name+'_Zsti_slope_z.png')
            plt.close(fig)


                
    def constructQuantity(self,name,timeIndex=None,locIndex=None,flatten=False):
        '''Loop over our list of models and construct the corresponding list of the quantity name.
        name can be a key in model[i].p (parameters) or model[i].var (variables), including either
        time or radial functions.'''
        construction = []
        #construction = np.zeros(())
        failures=[]
        log = True
        nameIsParam=False
        for i,model in enumerate(self.models):
            params = model.p.keys()
            varNames = model.var.keys()
            failure=False
            if(name in params):
                construction.append(model.p[name])
                nameIsParam=True
                if(not model.pLog[name]):
                    log=False
            elif(name in varNames):
                if(isinstance(model.var[name],RadialFunction)):
                    construction.append(model.var[name].sensible(timeIndex,locIndex))
                    log = model.var[name].log
                if(isinstance(model.var[name],TimeFunction)):
                    construction.append(model.var[name].sensible(timeIndex))
                    log = model.var[name].log
            else:
                failure=True
                print "Couldn't find the requested variable! ",name
            failures.append(failure)
        construction = np.array(construction)
        if flatten:
            construction = construction.flatten()
        #theRange = np.percentile(construction,(.2,99.8))
        if log:
            try:
                sgn = np.sign(construction)
                npmax = np.max(construction)
                if(npmax>0):
                    minAllowed = np.min(construction[sgn==1])
                else:
                    minAllowed = npmax*1.0e-7
                construction[sgn!=1] = minAllowed
            except:
                pdb.set_trace() # Something is going wrong in the above and I don't know what!
        if nameIsParam:
            theRange = [np.min(construction),np.max(construction)]
        elif self.models[0].get(name).theRange is None:
            try:
                theRange = [np.min(construction),np.max(construction)]
            except:
                pdb.set_trace()
        else:
            theRange= copy.deepcopy(self.models[0].get(name).theRange)

        # If the arrays weren't all the same size from model to model, we have a problem.
        if(hasattr(theRange[0],"__len__")):
            print "WARNING: arrays of different sizes for variable ",name
            theRange[0] = np.min(theRange[0])
            theRange[1] = np.max(theRange[1])

        if(log):
            try:
                assert theRange[0]>0 and theRange[1]>0
            except:
                pdb.set_trace()
            if(theRange[1]/theRange[0] < 10.0):
                log=False

        if theRange[0]==theRange[1]:
            theRange[1]=theRange[0]+1.0
            if theRange[0]>0:
                theRange[0]=theRange[0]*0.9
            else:
                theRange[0]=theRange[0]-1.0

        return construction,failures,log,theRange
        ##### End of constructQuantity(self,name,timeIndex=None,locIndex=None):

    def modelClosestTo(self,mstar, sfr, z):
        '''Locate the model closest to these mstar, sfr, z coordinates as follows:
            - Find the nearest redshift to z in the grid of outputs - call it zCoarse
            - Compare this mstar and sfr to every mstar and sfr by computing a euclidean
            distance in log space'''
        minDist = 1000.0
        minmodel=-1
        for i,model in enumerate(self.models):
            timeIndex,zCoarse = Nearest(model.var['z'].sensible(),z)
            logInputMstar = log(mstar)
            logModelMstar = log(model.var['mstar'].sensible(timeIndex))
            logInputSFR = log(sfr)
            logModelSFR = log(model.var['sfr'].sensible(timeIndex))
            dist = sqrt((logModelMstar- logInputMstar)**2.0 + (logModelSFR - logInputSFR)**2.0)
            if(dist<minDist):
                minDist=dist
                minmodel = i
        return self.models[minmodel],timeIndex
    def msCoords(self,z):
        coords = []
        for model in self.models:
            timeIndex,zCoarse = Nearest(model.var['z'].sensible(),z)
            coords.append((model.var['mstar'].sensible(timeIndex), \
                            model.var['sfr'].sensible(timeIndex), \
                            zCoarse))
        return coords

    def fitMS(self,mscoords):
        coords = np.array(mscoords)
        mst = np.log10(coords[:,0])
        sfr = np.log10(coords[:,1])
        params = np.polyfit(mst,sfr,1)
        residuals = sfr - (mst*params[0] + params[1])
        scatter = np.std(residuals)
        return [params[0],params[1],scatter],residuals
    def storeMS(self):
        allResiduals = []
        msParams = []
        for z in self.models[0].var['z'].sensible():
            params,residuals = self.fitMS(self.msCoords(z))
            allResiduals.append(residuals)
            msParams.append(params)
        npres = np.array(allResiduals)
        self.msParams = np.array(msParams)
        rng = [-0.7*np.max(self.msParams[:,2]),0.7*np.max(self.msParams[:,2])]
        for i, model in enumerate(self.models):
            model.var['deltaMS'] = TimeFunction(npres[:,i], 'deltaMS', 1.0, 1.0, r'$\Delta$ MS (dex)',log=False,theRange=rng)
        # Save data as time(Gyr), slope, zp, scatter
        np.savetxt(self.name+'_msParams.dat',np.vstack((self.models[0].var['t'].sensible(),self.msParams.T)).T)
    def setParam(self,keyname,value,models=None):
        if(models is None):
            models=range(len(self.models))
        for i in models:
            self.models[i].p[keyname] = value




    
if __name__=='__main__':
    pass

