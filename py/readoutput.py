import numpy as np
import pdb
import struct
import glob
import os
import random
from cosmolopy import *
from math import pi,log,sinh,sin,cos,sqrt,log10
#from bitstring import Bits
import pyfits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mpcolors
import matplotlib.cm as cmx

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

cm = plt.get_cmap('jet')

outward = [0,1,-1,2,-2,3,-3,4,-4]

gidgetdir = '../'

class RadialFunction:
    def __init__(self,arr,name,cgsConv,sensibleConv,texString,inner=0,outer=0,log=True,theRange=None):
        self.arr=arr
        self.name=name
        self.cgsConv=cgsConv
        self.sensibleConv = sensibleConv
        self.texString = texString
        self.innerVal = inner
        self.outerVal = outer
        self.log = log
        self.theRange = theRange
    def inner(self,timeIndex=None,cgs=False):
        if(cgs):
            if(timeIndex is not None):
                return self.innerVal[timeIndex]*self.cgsConv
            return self.innerVal*self.cgsConv
        else:
            if(timeIndex is not None):
                return self.innerVal[timeIndex]*self.sensibleConv
            return self.innerVal*self.sensibleConv
    def outer(self,timeIndex=None,cgs=False):
        if(cgs):
            if(timeIndex is not None):
                return self.outerVal[timeIndex]*self.cgsConv
            return self.outerVal*self.cgsConv
        else:
            if(timeIndex is not None):
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
    def range(self):
        if(self.theRange is None):
            return self.theRange
        else:
            return np.percentile(self.arr,(.2,99.8))*sensibleConv

class TimeFunction:
    def __init__(self,arr,name,cgsConv,sensibleConv,texString,log=True,theRange=None):
        self.arr = arr
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
            plt.plot(self.vars[rfx].sensible(ti), self.vars[rfy].sensible(ti),c=rgb[ti])
        plt.savefig(self.name+'_'+rfy+'_vs_'+rfx+'.png')
        plt.clf()
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
                        pixelVals[i,j] = (self.var['mdotBulgeG'].sensible(timeIndex)*self.p['RfREC']/(self.p['RfREC']+self.p['mu']))/ (pi * annuliA[0]**2.0 )#annuliB[0]) # SFR (Msun/yr) per area (kpc^2)
                
        pixelVals = pixelVals / cos(inclination) # projected surface density increases with inclination.
        hdu = pyfits.PrimaryHDU(pixelVals)
        hdulist = pyfits.HDUList([hdu])
        if(os.path.exists(filename+'.fits')):
            os.remove(filename+'.fits')
        hdulist.writeto(filename+'.fits')
    def nTimeSteps(self):
        return self.nt
    def read(self):
        with open(self.path+'_comment.txt','r') as comment:
            lines = comment.readlines()
            # paramnames - copied from exper.py's experiment class.
            paramnames = ['nx','eta','epsff','tauHeat','analyticQ', \
                    'cosmologyOn','xmin','NActive','NPassive','vphiR', \
                    'R','gasTemp','Qlim','fg0','phi0', \
                    'zstart','tmax','stepmax','TOL','mu', \
                    'b','innerPowerLaw','softening','diskScaleLength','whichAccretionHistory', \
                    'alphaMRI','thickness','migratePassive','fixedQ','kappaMetals', \
                    'Mh0','minSigSt','NChanges','dbg','accScaleLength', \
                    'zquench','zrelax','xiREC','RfREC','deltaOmega', \
                    'Noutputs','accNorm','accAlphaZ','accAlphaMh','accCeiling', \
                    'fscatter','invMassRatio','fcool','whichAccretionProfile','alphaAccretionProfile', \
                    'widthAccretionProfile','fH2Min','tDepH2SC','ZIGM','yREC']
            params=[]
            for k,line in enumerate(lines):
                cloc = line.find(':')
                if(cloc != -1):
                    params.append(float(line[cloc+1:-1]))
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
        with open(self.path+'_evolution.dat','r') as evolution:
            ncolev, =struct.unpack('i',evolution.read(4))
            self.ncolev=ncolev
            arr = evolution.read()
            evarray = np.fromstring(arr)
            self.nsteps = len(evarray)/ncolev
            self.evarray = np.reshape(evarray,(self.nsteps,self.ncolev))
        with open(self.path+'_radial.dat','r') as radial:
            dataCube=[]
            for i in range(self.nsteps):
                ncolstep, = struct.unpack('i',radial.read(4))
                nrowstep, = struct.unpack('i',radial.read(4))
                strSqArr = radial.read(8*ncolstep*nrowstep)
                sqArr = np.fromstring(strSqArr)
                dataCube.append(np.reshape(sqArr,(nrowstep,ncolstep)))
            self.dataCube = np.array(dataCube)
        # Alright, at this point we've read in the critical data.
        # Let's try to store it in a more comprehensible manner.
        # Keep in mind that self.dataCube ~ (timestep, nx, var)
        # And self.evarray ~ (timestep, var)
        self.var={}
        self.var['step'] = TimeFunction( \
                np.copy(self.evarray[:,0]), \
                'step',1,1,'Number of Steps')
        self.var['t'] = TimeFunction(\
                np.copy(self.evarray[:,1]),'t', \
                2.0*pi*cmperkpc*self.p['R']/(1.0e5*self.p['vphiR']) ,\
                2.0*pi*cmperkpc*self.p['R']/(1.0e5*self.p['vphiR']*speryear*1.0e9), \
                'Time since z=2 (Gyr)')
        self.nt = len(self.evarray[:,1])
        self.var['z'] = TimeFunction(np.copy(self.evarray[:,9]),'z',1,1,'z')
        self.var['r'] = RadialFunction( \
                np.copy(self.dataCube[:,:,0]),'r', \
                self.p['R']*cmperkpc,self.p['R'],'r (kpc)')
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
        self.var['dA'] = RadialFunction( \
                pi*(np.power(self.getData('rb',locIndex=range(1,nxI+1),cgs=True),2.0) \
                - np.power(self.getData('rb',locIndex=range(nxI),cgs=True),2.0)), \
                'dA',1.0,1.0/cmperkpc**2.0, r'dA (kpc$^{2}$)')
        self.var['col'] =RadialFunction( \
                np.copy(self.dataCube[:,:,3]),'col', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma (M_\odot\ pc^{-2})$',theRange=[0.1,1000])
        # Conversion from code units to Msun/pc^2 for column density-like units
        colSensibleConv = self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc)
        self.var['colst'] =RadialFunction( \
                np.copy(self.dataCube[:,:,5]),'colst', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma_* (M_\odot\ pc^{-2})$', \
                inner=self.evarray[:,3]*2.0*self.p['md0']*self.p['R']*kmperkpc/(speryear*self.p['vphiR']*self.var['rb'].sensible(None,0)**2.0 * pcperkpc**2.0 *colSensibleConv), theRange=[0.1,3.0e4])
        self.var['sig'] =RadialFunction( \
                np.copy(self.dataCube[:,:,4]*self.p['vphiR']),'sig', \
                1.0e5,1.0, r'$\sigma$ (km s$^{-1}$',log=True)
        self.var['maxsig'] = TimeFunction(np.amax(self.dataCube[:,:,4]*self.p['vphiR'],axis=1), \
                'maxsig',1.0e5,1.0,r'$\max(\sigma)$',log=True)
        self.var['mdotBulgeG'] = TimeFunction(np.copy(self.evarray[:,8]),'mdotBulgeG', \
                self.p['md0']*gpermsun/speryear,self.p['md0'],r'$\dot{M}_\mathrm{Gas to Bulge}$')
        self.var['dcoldt'] = RadialFunction( \
                np.copy(self.dataCube[:,:,7]),'dcoldt',self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']*cmperkpc)**2.0), \
                self.p['md0']/(2.0*pi*(self.p['R'])**2.0),r'$\partial \Sigma/\partial t$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)')
        self.var['colAccr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,29]),'colAccr', \
                self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']**2.0)*cmperkpc*cmperkpc),\
                self.p['md0']/(2.0*pi*self.p['R']**2.0), \
                r'$\dot{\Sigma}_{cos} (M_\odot\ yr^{-1}\ kpc^{-2})$', \
                inner=2.0*(self.evarray[:,20])*self.p['RfREC']/((self.p['RfREC']+self.p['mu'])*np.power(internalR[:,0]-self.dataCube[:,0,0]*dlnx,2.0)) )
        mdotCodeUnits = np.column_stack((self.evarray[:,8],np.copy(self.dataCube[:,:,38])/self.p['md0']))
        self.var['Mdot'] = RadialFunction( \
                mdotCodeUnits, 'Mdot', \
                self.p['md0']*gpermsun/speryear, \
                self.p['md0'],r'$\dot{M}$ (M$_\odot$ yr$^{-1}$)',log=False)
        self.var['colTr'] = RadialFunction( \
                (self.var['Mdot'].cgs(locIndex=range(1,nxI+1))-self.var['Mdot'].cgs(locIndex=range(nxI)))/self.var['dA'].cgs(), \
                'colTr',1.0, speryear*cmperkpc**2.0/gpermsun, r'$\dot{\Sigma}_{tr}$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)',log=False)
        self.var['colsfr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,28]),'colsfr', \
                self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']**2.0)*cmperkpc*cmperkpc),\
                self.p['md0']/(2.0*pi*self.p['R']**2.0), \
                r'$\dot{\Sigma}_*^{SF} (M_\odot\ yr^{-1}\ kpc^{-2})$', \
                inner=2.0*(self.evarray[:,8]+self.evarray[:,20])*self.p['RfREC']/((self.p['RfREC']+self.p['mu'])*np.power(internalR[:,0]-self.dataCube[:,0,0]*dlnx,2.0)) )
        self.var['fH2']= RadialFunction(np.copy(self.dataCube[:,:,47]),'fH2',1.0,1.0,r'$f_{\mathrm{H_2}}$',log=False)
        self.var['Z'] = RadialFunction(np.copy(self.dataCube[:,:,21]),'Z',1.0,1.0,'Z')
        self.var['NHI'] = RadialFunction(self.getData('col',cgs=True)*(1.0-self.getData('fH2'))*(1.0-self.getData('Z'))/gperH,\
                'NHI', 1.0,1.0,r'$N_{H\mathrm{I}}$ cm$^{-2}$')
        self.var['Mh'] = TimeFunction(self.evarray[:,18],'Mh',gpermsun,1.0,r'$M_h (M_\odot)$')
        self.var['scaleRadius'] = TimeFunction( \
                self.p['accScaleLength']*np.power(self.var['Mh'].sensible()/self.p['Mh0'],self.p['alphaAccretionProfile']), \
                'accScaleLength',cmperkpc,1.0,r'$r_\mathrm{acc}$ (kpc)',log=False)
        self.var['lambda'] = TimeFunction( self.getData('scaleRadius')/np.power(self.getData('Mh'),1.0/3.0), \
                'lambda', 1.0,1.0,r'$\lambda$')
        self.var['rx'] = RadialFunction( \
                np.array([self.var['r'].sensible(ti)/self.var['scaleRadius'].sensible(ti) for ti in range(self.nt)]), \
                'rx',1.0,1.0,r'r/r$_{acc}$',log=False)
        self.var['sfr'] = TimeFunction( \
                self.var['mdotBulgeG'].sensible()*self.p['RfREC']/(self.p['RfREC']+self.p['mu']) \
                +np.sum( self.var['dA'].sensible()*self.var['colsfr'].sensible(), 1 ), \
                'sfr',gpermsun/speryear, 1.0, r'SFR (M$_\odot$ yr$^{-1}$)')
        self.var['mCentral'] = TimeFunction( \
                self.evarray[:,3], 'mCentral',
                self.p['md0']*2.0*pi*self.p['R']/self.p['vphiR'] *gpermsun* kmperkpc/speryear, \
                self.p['md0']*2.0*pi*self.p['R']/self.p['vphiR'] * kmperkpc/speryear, \
                r'$M_\mathrm{center}\ (M_\odot)$')
        self.var['mstar'] = TimeFunction( \
                self.var['mCentral'].sensible() + 
                np.sum( self.var['dA'].sensible()*self.var['colst'].sensible()*1.0e6, 1 ), \
                'mstar',gpermsun,1.0,r'$M_*$ (M$_\odot$)',theRange=[1.0e9,1.0e11])
        self.var['mBulge'] = TimeFunction(self.computeMBulge()[0],'mBulge',gpermsun,1.0,r'$M_B (M_\odot)$')
        self.var['BT'] = TimeFunction(self.var['mBulge'].cgs()/(self.var['mBulge'].cgs()+self.var['mstar'].cgs()), \
                'BT',1,1,r'Bulge to Total Ratio',log=False)
        self.var['sSFR'] = TimeFunction(self.getData('sfr',cgs=True)/self.getData('mstar',cgs=True) , \
                'sSFR',1.0,1.0e9*speryear,r'sSFR (Gyr$^{-1}$)')
        mg = self.getData('col')
        self.var['integratedZ'] = TimeFunction(np.sum(mg*self.getData('Z'),axis=1)/np.sum(mg,axis=1), \
                'integratedZ', 1.0, 1.0, r'Z_g')
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
            if(isinstance(self.var[key],TimeFunction) and not key in blacklist):
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
        failures=[]
        for timeIndex in range(len(self.var['z'].sensible())):
            doneFlag = 0
            scaleRadius = self.var['scaleRadius'].sensible(timeIndex)
            for i in range(maxTries):
                f=1.5
                ind,theMin = Nearest(r,scaleRadius*f)
                logcol = np.log10(self.var['colst'].cgs(timeIndex))
                slope = float(maxTries/2.0-i)/float(maxTries/2.0) * (logcol[ind]-logcol[ind-1])/(r[ind]-r[ind-1])
                yguess = logcol[ind] + (r[:]-r[ind])*slope
                fail = (yguess[0:ind] - logcol[0:ind]).clip(min=0)
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
                fc.append(centralExcess/(excessL+centralExcess))
                failures.append(0)
            else:
                mb.append(0)
                failures.append(1)
        # End loop over time.
        return np.array(mb)/gpermsun,failures
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
            if(os.stat(fnames[i]+'_stde_aux.txt')[6]==0):
                self.models.append(SingleModel(fnames[i]))
            else:
                pass
        self.fn = fnames
    def merge(self,expt):
        '''Merge two experiments - essentially just append two lists of models. '''
        self.fn = self.fn + expt.fn # concatenate the lists of models
    def read(self):
        ''' Read in every model in the experiment. '''
        n=0
        for model in self.models:
            model.read()
            n+=1
            if(n % 50 == 0):
                print "Reading in model ",n," of ",len(self.models)
    def radialPlot(self,timeIndex=None,variables=None,colorby=None,percentiles=None,logR=False,scaleR=False):
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
            #logR=False
            #if(rRange[0]/rRange[1] < .001 ):
            #    logR=True


            dirname = 'movie_'+self.name+'_'+v+'_cb'+colorby+'_vs'+indVar
            if(not os.path.exists(dirname)):
                os.makedirs(dirname)
            print "Making movie: ",v

            counter = 0 
            for ti in timeIndex:
                theVar,varFail,varLog,varRange = self.constructQuantity(v,ti)
                r,rFail,rLog,_ = self.constructQuantity(indVar,ti)
                #rx,rxFail,rxLog,_ = self.constructQuantity('rx',ti)

                colors,fail,log,_ = self.constructQuantity(colorby,timeIndex=ti)
                #rgb, scaledColors,origColors = getRGB(colors,log)
                if(log):
                    colors=np.log10(colors)
                
                fig,ax = plt.subplots(1,1)
                model = self.models[0]
                ax.set_xlabel(model.get(indVar).texString)
                ax.set_ylabel(model.get(v).texString)
                cNorm = mpcolors.Normalize(vmin=np.min(colors),vmax=np.max(colors))
                scalarMap = cmx.ScalarMappable(cmap=cm,norm=cNorm)
                scalarMap.set_array(colors)
                scalarMap.autoscale()
                cbar = plt.colorbar(scalarMap,ax=ax)
                cbar.set_label(colorby)
                #ax2 = ax.twin()
                #s = ax.scatter(r,model.getData(v,ti),c=colors,cmap=cm,alpha=0)
                lwnorm = 1.0 + log10(float(len(self.models)))
                for k,model in enumerate(self.models):
                    try:
                        ax.plot(r[k],theVar[k],c=scalarMap.to_rgba(colors[k]),lw=2.0/lwnorm)
                    except ValueError:
                        ax.plot(model.getData('rb',ti),theVar[k],c=scalarMap.to_rgba(colors[k]),lw=2.0/lwnorm)
                ax.set_xlim(rRange[0],rRange[1])
                ax.set_ylim(overallRange[0],overallRange[1])
                if(overallLog):
                    ax.set_yscale('log')
                if(logR):
                    ax.set_xscale('log')
                z=model.getData('z')
                dispz = "%.3f" % z[ti]
                plt.text(0.7,0.9,'z='+dispz,transform=ax.transAxes)
                plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
                counter = counter+1
                plt.clf()
                del fig
                del ax
    def ptMovie(self,xvar='mstar',yvar=None,colorby=None,timeIndex=None,prev=0):
        if(yvar is None):
            variables = self.models[0].getTimeFunctions()
        if(timeIndex is None):
            timeIndex=range(1,self.models[0].nTimeSteps(),3)
        if(colorby is None):
            colorby='Mh0'

        overallX,_,overallXLog,overallXRange = self.constructQuantity(xvar)
        for i,v in enumerate(variables):
            dirname = 'movie_'+self.name+'_'+v+'_cb'+colorby+'_vs'+xvar
            if(not os.path.exists(dirname)):
                os.makedirs(dirname)
            print "Making movie : ",v," vs ",xvar

            overallVar,_,overallLog,overallRange = self.constructQuantity(v)
            counter=0
            for ti in timeIndex:
                fig,ax = plt.subplots(1,1)
                model = self.models[0]
                ax.set_xlabel(model.get(xvar).texString)
                ax.set_ylabel(model.get(v).texString)

                colors,fail,log,_ = self.constructQuantity(colorby,timeIndex=ti)
                if(log):
                    colors=np.log10(colors)
                cNorm = mpcolors.Normalize(vmin=np.min(colors),vmax=np.max(colors))
                scalarMap = cmx.ScalarMappable(cmap=cm,norm=cNorm)
                scalarMap.set_array(colors)
                scalarMap.autoscale()
                cbar=plt.colorbar(scalarMap,ax=ax)
                cbar.set_label(colorby)
                #pdb.set_trace()
                ax.scatter(overallX[:,ti],overallVar[:,ti],c=colors,cmap=cm,norm=cNorm)
                if(prev>0):
                    for k in range(max(ti-prev,0),ti):
                        ax.scatter(overallX[:,k],overallVar[:,k],c=colors,cmap=cm,s=6,lw=0)
                ax.set_xlim(overallXRange[0],overallXRange[1])
                ax.set_ylim(overallRange[0],overallRange[1])
                if(overallLog):
                    ax.set_yscale('log')
                if(overallXLog):
                    ax.set_xscale('log')
                z=model.getData('z')
                dispz = "%.3f" % z[ti]
                plt.text(0.7,0.9,'z='+dispz,transform=ax.transAxes)
                plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
                counter+=1
                plt.clf()
                del fig
                del ax

    def timePlot(self,variables=None,colorby=None):
        if(variables is None):
            variables = self.models[0].getTimeFunctions()
        time,fail,log,_ = self.constructQuantity('t')
        z,fail,log,_ = self.constructQuantity('z')
        if(colorby is None):
            colorby = 'Mh0'
        colors,fail,log,_ = self.constructQuantity(colorby)
        #scalarMap = cmx.ScalarMappable(cmap=cm,norm=plt.normalize(vmin=0,vmax=1))
        if(log):
            colors=np.log10(colors)
        #scalarMap.set_array(colors)

        #scalarMap.autoscale()
        #pdb.set_trace()
        #rgb, colors = getRGB(colors,log)
        for i,v in enumerate(variables):
            fig,ax=plt.subplots(1,1)
            model = self.models[0]

            theVar,varFail,varLog,varRange = self.constructQuantity(v)
            #sc=plt.scatter(model.getData('t'),model.getData(v),c=colors,cmap=cm)
            plt.cla()

            cNorm=mpcolors.Normalize(vmin=np.min(colors),vmax=np.max(colors))
            scalarMap = cmx.ScalarMappable(cmap=cm,norm=cNorm)
            scalarMap.set_array(colors)
            scalarMap.autoscale()
            
            ax.set_xlabel(model.get('t').texString)
            ax.set_ylabel(model.get(v).texString)
            cbar = plt.colorbar(scalarMap,ax=ax)
            cbar.set_label(colorby)
            ax2 = ax.twiny()
            z = model.getData('z')
            t = model.getData('t')
            zs = [0.0,0.5,1.0,1.5,2.0]
            zind = [Nearest(z,zl)[0] for zl in zs]
            correspondingTs = [t[zi] for zi in zind]
            #print "z,ind,ts: ", zs, zind, correspondingTs
            ax2.set_xticks(correspondingTs)
            ax2.set_xticklabels([str(zl) for zl in zs] )
            ax2.set_xlabel(model.get('z').texString)
            lwnorm = 1.0 + log10(float(len(self.models)))
            for j,model in enumerate(self.models):
                ax.plot(model.getData('t'),theVar[j],c=scalarMap.to_rgba(colors[j]),lw=2.0/lwnorm)
                    
            ax.set_xlim(correspondingTs[-1],correspondingTs[0])
            ax.set_ylim(varRange[0],varRange[1])
            if model.get(v).log:
                ax.set_yscale('log')
            ## ax2.set_xlim(zs[-1],zs[0])

            #scalarMap.set_array(colors)
            #scalarMap.autoscale()
            print "Making plot v: ",v

            #plt.colorbar(sc,ax=ax)
            try:
                plt.savefig(self.name+'_'+v+'_cb'+colorby+'.png')
            except ValueError:
                print "Caught ValueError! "
                pdb.set_trace()
                abc=114
            plt.clf()
            plt.close(fig)
                
    def constructQuantity(self,name,timeIndex=None,locIndex=None):
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
        #theRange = np.percentile(construction,(.2,99.8))
        if log:
            sgn = np.sign(construction)
            npmax = np.max(construction)
            if(npmax>0):
                minAllowed = np.min(construction[sgn==1])
            else:
                minAllowed = npmax*1.0e-7
            construction[sgn!=1] = minAllowed
        if nameIsParam:
            theRange = [np.min(construction),np.max(construction)]
        elif self.models[0].get(name).theRange is None:
            theRange = [np.min(construction),np.max(construction)]
        else:
            theRange= self.models[0].get(name).theRange

        return construction,failures,log,theRange
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
        self.msParams = msParams
        for i, model in enumerate(self.models):
            model.var['deltaMS'] = TimeFunction(npres[:,i], 'deltaMS', 1.0, 1.0, r'$\Delta$ MS (dex)',log=False)
    def setParam(self,keyname,value,models=None):
        if(models is None):
            models=range(len(self.models))
        for i in models:
            self.models[i].p[keyname] = value




    
if __name__=='__main__':
    pass

