import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import glob
import os
import pdb


gidgetdir = os.environ['GIDGETDIR']+'/'

class behroozi:
    def __init__(self):
        #self.dir = '/Users/jforbes/behroozi/release-sfh_z0_z8_052913/'
        self.dir = gidgetdir + '../behroozi/release-sfh_z0_z8_052913/'
    def readSmmr(self):
        zStr = ['0.10','1.00','2.00','3.00','4.00','5.00','6.00','7.00','8.00']
        self.smmr_z = np.array([0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        self.smmr_Mh = []
        self.smmr_eff_hi=[]
        self.smmr_eff_lo=[]
        self.smmr_eff_me=[]
        self.smmr_mst_hi=[]
        self.smmr_mst_lo=[]
        self.smmr_mst_me=[]
        self.smmr_fn_eff_hi=[]
        self.smmr_fn_eff_lo=[]
        self.smmr_fn_eff_me=[]
        self.smmr_fn_mst_hi=[]
        self.smmr_fn_mst_lo=[]
        self.smmr_fn_mst_me=[]
        for k, z in enumerate(zStr):
            arr = np.loadtxt(glob.glob(self.dir+'smmr/c_smmr_z'+z+'*.dat')[0])
            self.smmr_Mh.append(arr[:,0])
            self.smmr_eff_me.append(arr[:,1])
            self.smmr_eff_hi.append(arr[:,1]+arr[:,2])
            self.smmr_eff_lo.append(arr[:,1]-arr[:,3])
            self.smmr_mst_me.append(arr[:,0]+arr[:,1])
            self.smmr_mst_hi.append(arr[:,0]+arr[:,1]+arr[:,2])
            self.smmr_mst_lo.append(arr[:,0]+arr[:,1]-arr[:,3])

            self.smmr_fn_eff_me.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_eff_me[k]))
            self.smmr_fn_eff_hi.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_eff_hi[k]))
            self.smmr_fn_eff_lo.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_eff_lo[k]))
            self.smmr_fn_mst_me.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_mst_me[k]))
            self.smmr_fn_mst_hi.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_mst_hi[k]))
            self.smmr_fn_mst_lo.append(interpolate.interp1d(self.smmr_Mh[k],self.smmr_mst_lo[k]))
    def readSfh(self):
        pass
    def getSmmr(self, Mh, z, eff=True):
        # First find the nearest z, and get the indices of the surround outputs: zi and zi2.
        zi = np.argmin(np.abs(self.smmr_z-z))
        if self.smmr_z[zi]<z:
            zi2 = zi+1
        else:
            zi = zi-1
            zi2 = zi+1
        if zi<0:
            zi=0
        if zi==len(self.smmr_z)-1:
            zi2=zi
        Mh1 = self.smmr_Mh[zi]
        Mh2 = self.smmr_Mh[zi2]

        Flag=False

        if(Mh < min(Mh1) and Mh < min(Mh2)):
            raise ValueError
        elif Mh>max(Mh1) and Mh>max(Mh2):
            raise ValueError
        elif Mh < min(Mh1) or Mh>max(Mh1):
            zi=zi2
            Mh1=Mh2[:]
            Flag = True
        elif Mh < min(Mh2) or Mh>max(Mh2):
            zi2=zi
            Mh2=Mh1[:]
            Flag=True
        if eff:
            left = np.array([self.smmr_fn_eff_lo[zi](Mh), self.smmr_fn_eff_me[zi](Mh), self.smmr_fn_eff_hi[zi](Mh)])
            right= np.array([self.smmr_fn_eff_lo[zi2](Mh), self.smmr_fn_eff_me[zi2](Mh), self.smmr_fn_eff_hi[zi2](Mh)])
        else:
            left = np.array([self.smmr_fn_mst_lo[zi](Mh), self.smmr_fn_mst_me[zi](Mh), self.smmr_fn_mst_hi[zi](Mh)])
            right= np.array([self.smmr_fn_mst_lo[zi2](Mh), self.smmr_fn_mst_me[zi2](Mh), self.smmr_fn_mst_hi[zi2](Mh)])
        
        if zi==zi2:
            ret = left
        else:
            ret = (left + (z-self.smmr_z[zi])*(right-left)/(self.smmr_z[zi2]-self.smmr_z[zi]))
        return ret,Flag

    def plotSmmr(self, z, ax, minMh=None, maxMh=None, eff=True):
        if minMh is None:
            minMh = 1.0e10
        if maxMh is None:
            maxMh = 1.0e13


        mh=[]
        me=[]
        hi=[]
        lo=[]
        mhFlagged=[]
        meFlagged=[]
        hiFlagged=[]
        loFlagged=[]
        for i in range(100):
            mmh = minMh * (maxMh/minMh)**(float(i)/(99.0)) 
            try:
                ret,flag = self.getSmmr(np.log10(mmh),z,eff)
                if not flag:
                    me.append(10.0**ret[1])
                    hi.append(10.0**ret[2])
                    lo.append(10.0**ret[0])
                    mh.append(mmh)
                if flag:
                    meFlagged.append(10.0**ret[1])
                    hiFlagged.append(10.0**ret[2])
                    loFlagged.append(10.0**ret[0])
                    mhFlagged.append(mmh)
            except ValueError:
                pass
        ax.plot(mh,me,color='r',lw=3,ls='--')
        ax.fill_between(mh,lo,hi,color='r',alpha=0.2)
        ax.plot(mhFlagged,meFlagged,color='red',lw=2,ls='--')
        ax.fill_between(mhFlagged,loFlagged,hiFlagged,color='red',alpha=0.1)
        #print "Plotting SMMR!"
        #pdb.set_trace()


if __name__=='__main__':
    b = behroozi()
    b.readSmmr()
    fig,ax = plt.subplots()
    b.plotSmmr(.3,ax)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.savefig('testsmmr.png')
    plt.close(fig)


