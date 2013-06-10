import numpy as np
import pdb
import struct
import glob
import os
import random
from cosmolopy import *
from math import pi,log,sinh,sin,cos,sqrt
#from bitstring import Bits
import pyfits
import matplotlib.pyplot as plt 

speryear = 31557600.0
cmperkpc = 3.08567758e21
cmperpc = 3.08567758e18
gpermsun = 1.9891e33
Gcgs = 6.67384e-8
kB = 1.3806488e-16

outward = [0,1,-1,2,-2,3,-3,4,-4]

gidgetdir = '../'

class RadialFunction:
    def __init__(self,arr,name,cgsConv,sensibleConv,texString,inner=0,outer=0):
        self.arr=arr
        self.name=name
        self.cgsConv=cgsConv
        self.sensibleConv = sensibleConv
        self.texString = texString
        self.innerVal = inner
        self.outerVal = outer
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

class TimeFunction:
    def __init__(self,arr,name,cgsConv,sensibleConv,texString):
        self.arr = arr
        self.name = name
        self.cgsConv = cgsConv
        self.sensibleConv = sensibleConv
        self.texString = texString
    def cgs(self,timeIndex=False):
        if(timeIndex):
            return self.arr[timeIndex]*self.cgsConv
        return self.arr * self.cgsConv
    def sensible(self,timeIndex=False):
        if(timeIndex):
            return self.arr[timeIndex]*self.sensibleConv
        return self.arr * self.sensibleConv

class SingleModel:
    def __init__(self,path):
        self.path=path
        rfs = path.rfind('/') # location of last '/' in the path.
        self.name = path[rfs+1:] # the run name follows the '/'
        self.dirname = path[:rfs] # the directory name preceds the '/'
        self.p={}
    def plotSequenceOfTimes(self,rfx, rfy, N):
        figure = plt.figure()
        plt.plot(self.vars[rfx].sensible(ti), self.vars[rfy].sensible(ti))
    def plotSequenceOfRedshifts(self,rf1, rf2, N):
        pass
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
        self.var['z'] = TimeFunction(np.copy(self.evarray[:,9]),'z',1,1,'z')
        self.var['r'] = RadialFunction( \
                np.copy(self.dataCube[:,:,0]),'r', \
                self.p['R']*cmperkpc,self.p['R'],'r (kpc)')
        dlnx=-log(self.p['xmin'])/(self.p['nx']-1.0)
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
            'rb',self.p['R']*cmperkpc,self.p['R'],'r border (kpc)')
        self.var['dA'] = RadialFunction( \
                2.0*pi* self.dataCube[:,:,0] *self.dataCube[:,:,0] * sinh(dlnx), 'dA', \
                (self.p['R']*cmperkpc)**2.0,self.p['R']**2.0,'$\Delta$A (kpc$^2$)')
        self.var['col'] =RadialFunction( \
                np.copy(self.dataCube[:,:,3]),'col', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma (M_\odot\ yr^{-1}$')
        self.var['colst'] =RadialFunction( \
                np.copy(self.dataCube[:,:,3]),'col', \
                self.p['md0']*gpermsun/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                self.p['md0']*cmperpc*cmperpc/(self.p['vphiR']*self.p['R']*speryear*1.0e5*cmperkpc), \
                r'$\Sigma (M_\odot\ pc^{-1}$')
        self.var['sig'] =RadialFunction( \
                np.copy(self.dataCube[:,:,4]),'sig', \
                1.0e5*self.p['vphiR'],self.p['vphiR'],\
                r'$\sigma$ (km s$^{-1}$')
        self.var['mdotBulgeG'] = TimeFunction(np.copy(self.evarray[:,8]),'mdotBulgeG', \
                self.p['md0']*gpermsun/speryear,self.p['md0'],r'$\dot{M}_\mathrm{Gas to Bulge}$')
        self.var['colsfr'] = RadialFunction( \
                np.copy(self.dataCube[:,:,28]),'colsfr', \
                self.p['md0']*gpermsun/(speryear*2.0*pi*(self.p['R']**2.0)*cmperkpc*cmperkpc),\
                self.p['md0']/(2.0*pi*self.p['R']**2.0), \
                r'$\dot{\Sigma}_*^{SF} (M_\odot\ yr^{-1}\ kpc^{-2})$', \
                inner=2.0*(self.evarray[:,8]+self.evarray[:,20])*self.p['RfREC']/((self.p['RfREC']+self.p['mu'])*np.power(internalR[:,0]-self.dataCube[:,0,0]*dlnx,2.0)) )
        self.var['fH2']= RadialFunction(np.copy(self.dataCube[:,:,47]),'fH2',1.0,1.0,r'$f_{\mathrm{H_2}}$')
        self.var['Z'] = RadialFunction(np.copy(self.dataCube[:,:,21]),'Z',1.0,1.0,'Z')

        self.var['sfr'] = TimeFunction( \
                self.var['mdotBulgeG'].sensible()*self.p['RfREC']/(self.p['RfREC']+self.p['mu']) \
                +np.sum( self.var['dA'].sensible()*self.var['colsfr'].sensible(), 1 ), \
                'sfr',gpermsun/speryear, 1.0, r'SFR (M$_\odot$ yr$^{-1}$)')
        self.var['mBulge'] = TimeFunction( \
                self.evarray[:,3], 'mBulge',
                self.p['md0']*gpermsun/(speryear*2.*pi*self.p['R']*cmperkpc/(self.p['vphiR']*1.0e5)),
                self.p['md0']/(speryear*2.*pi*self.p['R']*cmperkpc/(self.p['vphiR']*1.0e5)),
                r'$M_\mathrm{center}\ (M_\odot)$')
        self.var['mstar'] = TimeFunction( \
                self.var['mBulge'].sensible() + 
                np.sum( self.var['dA'].sensible()*self.var['colst'].sensible()*1.0e6, 1 ), \
                'mstar',gpermsun,1.0,r'$M_*$ (M_\odot)')
        del self.dataCube
        del self.evarray
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
        self.name=name
        fnameKey = gidgetdir+'analysis/*'+name+'*/*_comment.txt'
        fnames = glob.glob(fnameKey)
        self.models=[]
        for i,fn in enumerate(fnames):
            fnames[i] = fn[:-12]
            self.models.append(SingleModel(fnames[i]))
        self.fn = fnames
    def merge(self,expt):
        self.fn = self.fn + expt.fn
    def read(self):
        n=0
        for model in self.models:
            model.read()
            n+=1
            if(n % 50 == 0):
                print "Reading in model ",n," of ",len(self.models)
    def modelClosestTo(self,mstar, sfr, z):
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

def fitMS(mscoords):
    coords = np.array(mscoords)
    mst = np.log10(coords[:,0])
    sfr = np.log10(coords[:,1])
    params = np.polyfit(mst,sfr,1)
    residuals = sfr - (mst*params[0] + params[1])
    scatter = np.std(residuals)
    return [params[0],params[1],scatter]

def generateMockSampleCoords(N,MSparams,z):
    masses = drawFromSchechter(N,-1.3,10.0**10.9,10.0**9.0)
    coords = []
    locs=[]
    for mass in masses:
        rnd = random.gauss(0,1)*MSparams[2]
        coord = (mass, 10.0**((MSparams[0]*np.log10(mass)+MSparams[1])+rnd), z)
        locs.append(rnd/MSparams[2])
        coords.append(coord)
    return coords,locs

def imposeCuts(coords,sigmaMS,zrange,massrange,deltamsRange):
    newCoords = []
    newSigmaMS = []
    for i,coord in enumerate(coords):
        if(massrange[0]<coord[0] and coord[0] < massrange[1] and \
           zrange[0]<coord[2] and coord[2] < zrange[1] and \
           deltamsRange[0] < sigmaMS[i] and sigmaMS[i] < deltamsRange[1]):
            newCoords.append(coord)
            newSigmaMS.append(sigmaMS[i])
    return newCoords,newSigmaMS

def drawFromSchechter(N,alpha,Mcutoff, Mmin):
    vals = []
    attempts=0
    while len(vals)<N:
        draw = np.random.gamma(alpha+2,Mcutoff)
        if draw < Mmin:
            pass
        else:
            if(random.random() > Mmin/draw):
                pass
            else:
                vals.append(draw)
        attempts+=1
    print "To generate N= ",N," draws from a Schechter fn required ",attempts," attempts."
    return vals

def makeStackFrames(expt,coords,stackname):
#    avgProfR = 10.0 * np.power((200.0 - np.arange(200)[::-1])/200.0) # 10 kpc, logarithmic (untested!)
    # N is the number of cells in the radial direction.
    N=500
    avgProfR = 10.0*np.arange(N)/(N-1.0) # 10 kpc linear
    avgProfColSFR = np.zeros(N)
    # For every given coordinate in M* - SFR - z space, find the closest model galaxy
    for i in range(len(coords)):
        if(i % 50 == 0):
            print "Adding stack frame "+str(i)+" of "+str(len(coords))
        mstar,sfr,z = coords[i]
        model,timeIndex = expt.modelClosestTo(mstar,sfr,z)
        model.TwoDimSFR(expt.name+'_'+stackname+'_frame_'+str(i).zfill(4), random.random()*pi/2.0, random.random()*2.0*pi,timeIndex)
        for k in range(N):
            qa = model.quantityAt(avgProfR[k],'colsfr',timeIndex)
            if(type(qa) != type(np.float64(43))):
                pdb.set_trace()
            avgProfColSFR[k] += model.quantityAt(avgProfR[k],'colsfr',timeIndex) / float(len(coords))
    with open(expt.name+'_'+stackname+'_prof.dat','w') as plotfile:
        for k in range(N):
            plotfile.write(str(avgProfR[k])+' '+str(avgProfColSFR[k])+' \n')

def testStackAnalysis():
    test2 = Experiment('rg78')
    test2.read()
    z=1.0
    msCoords = test2.msCoords(z)
    msparams = fitMS(msCoords)
    sampleCoords, locs = generateMockSampleCoords(1000,msparams,z)
    massbounds = [1.0e9,3.14e9,1.0e10,3.14e10,1.0e11]
    masscuts = [(massbounds[i],massbounds[i+1]) for i in range(len(massbounds)-1)]
    for i,cut in enumerate(masscuts):
        sampleCoordsStack, locsStack = imposeCuts(sampleCoords,locs, (-.1,3.0), cut, (-5,5))
        makeStackFrames(test2, sampleCoordsStack, "massStack_"+str(i))
    mscuts = [(-50,-1),(-1,0),(0,1),(1,50)]
    for i,cut in enumerate(mscuts):
        sampleCoordsStack, locsStack = imposeCuts(sampleCoords,locs, (-.1,3.0), (1.0e10,5.0e10), mscuts)
        makeStackFrames(test2, sampleCoordsStack, "msStack_"+str(i)+"_")
    
if __name__=='__main__':
    #test = Experiment('rg78')
    #test.read()


