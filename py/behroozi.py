import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import glob
import os
import pdb


class behroozi:
    def __init__(self):
        self.dir = '/Users/jforbes/behroozi/release-sfh_z0_z8_052913/'
    def readSmmr(self):
        zStr = ['0.10','0.10','0.10','0.10','1.00','1.00','2.00','2.00','3.00','3.00','4.00','4.00','5.00','5.00','6.00','6.00','7.00','7.00','8.00','8.00']
        zs = [0.0,0.01,0.1,0.11, 1.0,1.01, 2.0,2.01, 3.0,3.01, 4.0,4.01, 5.0,5.01, 6.0,6.01, 7.0,7.01, 8.0,8.01]
        self.smmr_z = []
        self.smmr_Mh = []
        self.smmr_eff_hi=[]
        self.smmr_eff_lo=[]
        self.smmr_eff_me=[]
        self.smmr_mst_hi=[]
        self.smmr_mst_lo=[]
        self.smmr_mst_me=[]
        self.smmr_all_z=[]
        for k, z in enumerate(zStr):
            arr = np.loadtxt(glob.glob(self.dir+'smmr/c_smmr_z'+z+'*.dat')[0])
            for i in range(len(arr[:,0])):
                self.smmr_Mh.append(arr[i,0])
                self.smmr_eff_me.append(arr[i,1])
                self.smmr_eff_hi.append(arr[i,1]+arr[i,2])
                self.smmr_eff_lo.append(arr[i,1]-arr[i,3])
                self.smmr_mst_me.append(arr[i,0]+arr[i,1])
                self.smmr_mst_hi.append(arr[i,0]+arr[i,1]+arr[i,2])
                self.smmr_mst_lo.append(arr[i,0]+arr[i,1]-arr[i,3])
                self.smmr_z.append(zs[k])
        self.smmr_fn_eff_me = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_eff_me, kind='linear',bounds_error=True) 
        self.smmr_fn_eff_hi = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_eff_hi, kind='linear',bounds_error=True) 
        self.smmr_fn_eff_lo = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_eff_lo, kind='linear',bounds_error=True) 
        self.smmr_fn_mst_me = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_mst_me, kind='linear',bounds_error=True) 
        self.smmr_fn_mst_hi = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_mst_hi, kind='linear',bounds_error=True) 
        self.smmr_fn_mst_lo = interpolate.interp2d(self.smmr_Mh,self.smmr_z,self.smmr_mst_lo, kind='linear',bounds_error=True) 
    def readSfh(self):
        pass
    def getSmmr(self, which, Mh, z): 
        pass

    def plotSmmr(self, z, ax, minMh=None, maxMh=None, eff=True):
        if minMh is None:
            minMh = 1.0e10
        if maxMh is None:
            maxMh = 1.0e13
        mh=[]
        me=[]
        hi=[]
        lo=[]
        for i in range(100):
            mmh = minMh * (maxMh/minMh)**(float(i)/(99.0)) 
            try:
                if eff:
                    me.append(10.0**self.smmr_fn_eff_me(np.log10(mmh),z)[0])
                    hi.append(10.0**self.smmr_fn_eff_hi(np.log10(mmh),z)[0])
                    lo.append(10.0**self.smmr_fn_eff_lo(np.log10(mmh),z)[0])
                    mh.append(mmh)
                else:
                    me.append(10.0**self.smmr_fn_mst_me(np.log10(mmh),z)[0])
                    hi.append(10.0**self.smmr_fn_mst_hi(np.log10(mmh),z)[0])
                    lo.append(10.0**self.smmr_fn_mst_lo(np.log10(mmh),z)[0])
                    mh.append(mmh)
            except:
                pass
        ax.plot(mh,me,color='r',lw=3,ls='--')
        ax.fill_between(mh,lo,hi,color='r',alpha=0.2)
        

if __name__=='__main__':
    b = behroozi()
    b.readSmmr()
    fig,ax = plt.subplots()
    b.plotSmmr(.5,ax)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.savefig('testsmmr.png')
    plt.close(fig)


