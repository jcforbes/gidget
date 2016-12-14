import numpy as np
import scipy.interpolate.interp1d as interp1d

# This is the simplest thing that will currently do well in my setup. 
# 
class DataSet:
    def __init__(self, xvar, yvar, xval, yval, yLower=None, yUpper=None, xerr=None, yerr=None, zmin=-1, zmax=0.5, label=None, logx=True, logy=True):
        self.xvar = xvar
        self.yvar = yvar
        self.xval = xval
        self.yval = yval
        self.yLower=yLower
        self.yUpper=yUpper
        self.xerr = xerr
        if yerr is not None:
            if yUpper is None:
                yUpper = yval+yerr
            if yLower is None:
                yLower = yval-yerr
        self.zmin = zmin
        self.zmax = zmax
        self.label = label
        self.logx = logx
        self.logy = logy

    def plot(self, z, axIn=None, color='k', lw=2):
        ''' Plot the requested variables at redshift z. If z is outside the z
            range (specified by zmin, zmax in __init__), plot a dotted line instead'''
        pass
    def interior(self, x,y, sigma=1.0):
        ''' Return true if the given point is within this dataset's uncertainties'''
        xv = self.xval
        yv = self.yval
        yu = self.yUpper
        yl = self.yLower
        xEval = x
        yEval = y
        if self.logx:
            xv = np.log10(xv)
            xEval = np.log10(xEval)
        if self.logy:
            yv = np.log10(yv)
            yu = np.log10(yu)
            yl = np.log10(yl)
            yEval = np.log10(yEval)
        f = interp1d(xv,yv,kind='linear',bounds_error=False)
        fu = interp1d(xv,yv+sigma*(yu-yv),kind='linear',bounds_error=False)
        fl = interp1d(xv,yv-sigma*(yv-yl),kind='linear',bounds_error=False)
        #yf = f(xEval) # don't need this
        yfu = fu(xEval)
        yfl = fl(xEval)
        inx = np.logical_and( np.min(xv) < xEval, xEval<np.max(xv) )
        ret = np.logical_and( np.logical_and( yfl<yEval, yEval < yfu ), inx)
        return ret




def defineMoster(z):
    def Moster(Mh,z, mparams):
        M10, M11, N10, N11, beta10, beta11, gamma10, gamma11 = mparams
        logM1z = M10 + M11*z/(z+1.0)
        Nz = N10 + N11*z/(z+1.0)
        betaz = beta10 + beta11*z/(z+1.0)
        gammaz = gamma10 + gamma11*z/(z+1.0)
        M1 = np.power(10.0, logM1z)
        eff = 2.0*Nz / (np.power(Mh/M1,-betaz) + np.power(Mh/M1,gammaz))
        return eff
    central = np.array([11.590, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329])
    unc =     np.array([0.236, 0.353, 0.0058, 0.0069, 0.153, 0.225, 0.059, 0.173])
    logMhs = np.linspace(10, 14, num=100)
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

    mosterDS = DataSet('Mh', 'efficiency', Mhs, eff, yLower=effM, yUpper=effP, zmin=z-0.5, zmax=z+0.5, label='Moster10')
    mosterDS2 = DataSet('Mh', 'mstar', Mhs, eff*Mhs, yLower=effM*Mhs, yUpper=effP*Mhs, zmin=z-0.5, zmax=z+0.5, label='Moster10')
    datasets['Moster10eff'+str(z)] = mosterDS
    datasets['Moster10'+str(z)] = mosterDS2


def defineBroeils():
    mhis = np.power(10.0, np.linspace(8.0, 10.5, 100))
    DHI = np.power(10.0, (np.log10(mhis) - 6.52)/1.96)/2.0
    datasets.append( DataSet('MHI', 'broeilsHI', mhis, DHI, yLower=DHI/10**0.13, yUpper=DHI*10**0.13, label='Broeils97') )

def defineMetalRelations(z):
    thisMst = np.power(10.0, np.linspace(10.0, 11.5))
    ZHayward = -8.69 + 9.09*np.power(1.0+z[ti],-0.017) - 0.0864*np.power(np.log10(thisMst) - 11.07*np.power(1.0+z[ti],0.094),2.0)
    ZHayward = np.power(10.0, ZHayward) * 0.02
            
    b = 10.4 + 4.46*np.log10(1.0+z)-1.78*np.log10(1.0+z)**2.0
    ZGenzel15 = np.power(10.0,  8.74 - 0.087 * np.power(np.log10(thisMst) -b,2.0) -8.69) # Equation 12a normalized to solar
    #bb = np.log10(mst) - 0.32*np.log10(sfr) - 10
    #ZMannucci10 = np.power(10.0, 0.21 + 0.39*bb - 0.2*bb*bb - 0.077*bb*bb*bb + 0.064*bb*bb*bb*bb
    datasets['Hayward16Z'+str(z)] = DataSet( 'mstar', 'sfZ', thisMst, ZHayward/0.02, yLower=ZHayward/0.02/10.0**0.15, ZHayward/0.02*10.0**0.15, label='Hayward16', zmin=z-0.5, zmax=z+0.5 )
    datasets['Genzel15Z'+str(z)] = DataSet('mstar', 'sfZ', thisMst, ZGenzel15, yLower=ZGenzel15/10.0**0.15, yUpper=ZGenzel15*10.0**0.15, label='Genzel15', zmin=z-0.5, zmax=z+0.5 )
    mLee = np.power(10.0, np.linspace(5.89,9.29, 100))
    ZLee = np.power(10.0, 5.65 + 0.298*np.log10(mLee) - 8.7) # Z in Zsun
    datasets['Lee06'] = DataSet( 'mstar', 'sfZ', mLee, ZLee, yLower=ZLee/10.0**0.117, yUpper=ZLee*10.0**0.117, label='Lee06', zmax=0.5 )

    mTremonti = np.power(10.0, np.array([8.57, 8.67, 8.76, 8.86, 8.96, 9.06, 9.16, 9.26, 9.36, 9.46, 9.57, 9.66, 9.76, 9.86, 9.96, 10.06, 10.16, 10.26, 10.36, 10.46, 10.56, 10.66, 10.76, 10.86, 10.95, 11.05, 11.15, 11.25]))
    Z16Tremonti = np.power(10.0, np.array([8.25, 8.28, 8.32, 8.37, 8.46, 8.56, 8.59, 8.60,8.63, 8.66, 8.69, 8.72, 8.76, 8.80, 8.83, 8.85, 8.88, 8.92, 8.94, 8.96, 8.98, 9.00, 9.01, 9.02, 9.03, 9.03, 9.04, 9.03 ]) )
    Z50Tremonti = np.power(10.0, np.array([8.44, 8.48, 8.57, 8.61, 8.63, 8.66, 8.68, 8.71, 8.74, 8.78, 8.82, 8.84, 8.87, 8.90, 8.94, 8.97, 8.99, 9.01, 9.03, 9.05, 9.07, 9.08, 9.09, 9.10, 9.11, 9.11, 9.12, 9.12 ]))
    Z84Tremonti = np.power(10.0, np.array([8.64, 8.65, 8.70, 8.73, 8.75, 8.82, 8.82, 8.86, 8.88, 8.92, 8.94, 8.96, 8.99, 9.01, 9.05, 9.06, 9.09, 9.10, 9.11, 9.12, 9.14, 9.15, 9.15, 9.16, 9.17, 9.17, 9.18, 9.18 ]) )
    datasets['Tremonti04'] = DataSet( 'mstar', 'sfZ', mTremonti, Z50Tremonti, yLower=Z16Tremonti, yUpper=Z84Tremonti, label='Tremonti04', zmin=-0.5, zmax=0.5 )


def defineStructureRelations()
    mstThis = np.power(10.0, np.linspace(9.75, 11.25, 100))
    fang = np.power(10.0, 9.29 + 0.64*(np.log10(mstThis)-10.25))  # Fang et al 2013
    datasets['Fang13'] = DataSet( 'mstar', 'Sigma1', mstThis, fang, yLower=fang/10.0**0.16, yUpper=fang*10.0**0.16, label='Fang13 Red and Green', zmin=-0.5, zmax=3.0)


    from readoutput import nToC82
    mstThis = np.power(10.0, np.linspace(8.75,11.5, 30))
    gamma, M0, logn1, logn2 = (1.95, 10.0**10.16, 0.12, 0.62) # dutton09 all galaxies
    lognAll = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    gamma, M0, logn1, logn2 = (1.70, 10.0**9.90, 0.21, 0.62) # dutton09 red galaxies
    lognRed = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    gamma, M0, logn1, logn2 = (1.70, 10.0**10.41, 0.12, 0.61) # dutton09 blue galaxies
    lognBlue = logn2 + (logn1-logn2)/(1.0+np.power(10.0, gamma*np.log10(mstThis/M0)))
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognAll]
    datasets['Dutton09All'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 All')
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognRed]
    datasets['Dutton09Red'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 Red')
    c82s = [ nToC82( np.power(10.0,lognThis) ) for lognThis in lognBlue]
    datasets['Dutton09Blue'] = DataSet( 'mstar', 'c82', mstThis, c82s, label='Dutton09 Blue')


    mstBarro = np.power(10.0, np.linspace(10.0, 11.0, 10))
    er = 10.0**0.14
    barrohq = np.power(10.0, 0.65*(np.log10(mstBarro)-10.5) + 9.53)
    datasets['Barro15HQ'] = DataSet( 'mstar', 'Sigma1', mstBarro, barrohq, yLower=barrohq/er, yUpper=barrohq*er, label='Barro15 Red', zmin=0.5, zmax=1.0)
    barro1q = np.power(10.0, 0.65*(np.log10(mstBarro)-10.5) + 9.64)
    datasets['Barro151Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro1q, yLower=barro1q/er, yUpper=barro1q*er, label='Barro15 Red', zmin=1.0, zmax=1.5)
    barro2q = np.power(10.0, 0.64*(np.log10(mstBarro)-10.5) + 9.76)
    datasets['Barro152Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro2q, yLower=barro2q/er, yUpper=barro2q*er, label='Barro15 Red', zmin=1.5, zmax=2.2)
    barro3q = np.power(10.0, 0.67*(np.log10(mstBarro)-10.5) + 9.80)
    datasets['Barro153Q'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro3q, yLower=barro3q/er, yUpper=barro3q*er, label='Barro15 Red', zmin=2.2, zmax=3.05)
    mstBarro = np.power(10.0, np.linspace(9.0, 10.5, 10))
    er = 10.0**0.25
    barrohs = np.power(10.0, 0.89*(np.log10(mstBarro)-10.5) + 9.12)
    datasets['Barro15HS'] = DataSet( 'mstar', 'Sigma1', mstBarro, barrohs, yLower=barrohs/er, yUpper=barrohs*er, label='Barro15 Blue', zmin=0.5, zmax=1.0)
    barro1s = np.power(10.0, 0.88*(np.log10(mstBarro)-10.5) + 9.16)
    datasets['Barro151S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro1s, yLower=barro1s/er, yUpper=barro1s*er, label='Barro15 Blue', zmin=1.0, zmax=1.5)
    barro2s = np.power(10.0, 0.86*(np.log10(mstBarro)-10.5) + 9.25)
    datasets['Barro152S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro2s, yLower=barro2s/er, yUpper=barro2s*er, label='Barro15 Blue', zmin=1.5, zmax=2.2)
    barro3s = np.power(10.0, 0.89*(np.log10(mstBarro)-10.5) + 9.33)
    datasets['Barro153S'] = DataSet( 'mstar', 'Sigma1', mstBarro, barro3s, yLower=barro3s/er, yUpper=barro3s*er, label='Barro15 Blue', zmin=2.2, zmax=3.05)


    obsmst = np.power(10.0, np.linspace(9.0, 11.5, 100))
    obsmstETG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([0.03, 0.04, 0.13, 0.42, 0.65]))
    reffETGvdW50 = np.power(10.0, np.array([0.27, 0.27, 0.38, 0.67, 0.76]))
    reffETGvdW84 = np.power(10.0, np.array([0.46, 0.46, 0.58, 0.92, 1.08]))
    datasets['vdW14ETG0']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmax=0.5)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75 ]))
    reffLTGvdW16 = np.power(10.0, np.array([0.24, 0.36, 0.42, 0.61]))
    reffLTGvdW50 = np.power(10.0, np.array([0.49, 0.61, 0.66, 0.83]))
    reffLTGvdW84 = np.power(10.0, np.array([0.70, 0.80, 0.85, 1.01]))
    datasets['vdW14LTG0']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmax=0.5)
                
    obsmstETG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.02, -0.14, 0.02, 0.26, 0.62]))
    reffETGvdW50 = np.power(10.0, np.array([0.23, 0.21, 0.23, 0.45, 0.81]))
    reffETGvdW84 = np.power(10.0, np.array([0.43, 0.44, 0.42, 0.64, 0.97]))
    datasets['vdW14ETG1']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=0.5, zmax=1.0)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.18, 0.32, 0.39, 0.51, 0.77]))
    reffLTGvdW50 = np.power(10.0, np.array([0.43, 0.56, 0.64, 0.75, 0.90]))
    reffLTGvdW84 = np.power(10.0, np.array([0.65, 0.76, 0.83, 0.90, 1.12]))
    datasets['vdW14LTG1']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=0.5, zmax=1.0)

    obsmstETG = np.power(10.0, np.array([9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.15, -0.15, 0.07, 0.41]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.18, 0.09, 0.30, 0.58]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.42, 0.36, 0.54, 0.81]))
    datasets['vdW14ETG1h']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=1.0, zmax=1.5)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.11, 0.23, 0.33, 0.47, 0.62]))
    reffLTGvdW50 = np.power(10.0, np.array([0.37, 0.48, 0.57, 0.67, 0.82]))
    reffLTGvdW84 = np.power(10.0, np.array([0.60, 0.69, 0.77, 0.83, 0.96]))
    datasets['vdW14LTG1h']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=1.0, zmax=1.5)

    obsmstETG = np.power(10.0, np.array([9.75, 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.02, -0.27, -0.04, 0.28]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.22, 0.02, 0.19, 0.45]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.48, 0.35, 0.50, 0.74]))
    datasets['vdW14ETG2']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=1.5, zmax=2.0)

    obsmstLTG = np.power(10.0, np.array([9.25, 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([0.07, 0.16, 0.28, 0.35, 0.53]))
    reffLTGvdW50 = np.power(10.0, np.array([0.33, 0.42, 0.52, 0.61, 0.70]))
    reffLTGvdW84 = np.power(10.0, np.array([0.57, 0.65, 0.72, 0.80, 0.87]))
    datasets['vdW14LTG2']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=1.5, zmax=2.0)


    obsmstETG = np.power(10.0, np.array([ 10.25, 10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.37, -0.20, 0.16]))
    reffETGvdW50 = np.power(10.0, np.array([-0.04, 0.08, 0.36]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.36, 0.54, 0.55]))
    datasets['vdW14ETG2h']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=2.0, zmax=2.5)

    obsmstLTG = np.power(10.0, np.array([ 9.75, 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([ 0.10, 0.17, 0.26, 0.40]))
    reffLTGvdW50 = np.power(10.0, np.array([ 0.35, 0.44, 0.53, 0.64]))
    reffLTGvdW84 = np.power(10.0, np.array([ 0.57, 0.64, 0.70, 0.84]))
    datasets['vdW14LTG2h']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=2.0, zmax=2.5)

    obsmstETG = np.power(10.0, np.array([  10.75, 11.25]))
    reffETGvdW16 = np.power(10.0, np.array([-0.22, 0.07]))
    reffETGvdW50 = np.power(10.0, np.array([ 0.10, 0.39]))
    reffETGvdW84 = np.power(10.0, np.array([ 0.50, 0.68]))
    datasets['vdW14ETG3']  = DataSet( 'mstar', 'halfMassStars', obsmstETG, reffETGvdW50, yLower=reffETGvdW16, yUpper=reffETGvdW84, label='van der Wel 14 ETG', zmin=2.5, zmax=3.0)

    obsmstLTG = np.power(10.0, np.array([ 10.25, 10.75, 11.25]))
    reffLTGvdW16 = np.power(10.0, np.array([  0.16, 0.19, 0.33]))
    reffLTGvdW50 = np.power(10.0, np.array([  0.43, 0.47, 0.55]))
    reffLTGvdW84 = np.power(10.0, np.array([  0.65, 0.71, 0.76]))
    datasets['vdW14LTG3']  = DataSet( 'mstar', 'halfMassStars', obsmstLTG, reffLTGvdW50, yLower=reffLTGvdW16, yUpper=reffLTGvdW84, label='van der Wel 14 LTG', zmin=2.5, zmax=3.0)




def defineAngularMomenta()
    datasets['Fall13Disks'] = DataSet('mstar', 'specificJStars', [1.0e9,1.0e11], [10**2.3, 10**3.45], label='Fall13 Disks')
    datasets['Fall13Ellipticals'] = DataSet('mstar', 'specificJStars', [1.0e10, 5.0e11], [10**2.1, 10**2.1*50**0.6], label='Fall13 Ellipticals')
            if z[ti]>0.8 and z[ti]<2.6 and v=='specificJStars' or movie:
                thisLabel='Burkert16'
                thisLS='-'
    datasets['Burkert16'] = DataSet('mstar', 'specificJStars', [10**9.8, 10**11.4],[10.0**(3.33+2.0/3.0*(9.8-11.0)), 10.0**(3.33+2.0/3.0*(11.4-11.0))],  yLower=[10.0**(3.33+2.0/3.0*(9.8-11.0)-0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)-0.17)], yUpper=[10.0**(3.33+2.0/3.0*(9.8-11.0)+0.17), 10.0**(3.33+2.0/3.0*(11.4-11.0)+0.17)], zmin=0.8, zmax=2.6, label='Burkert16' )


#### Work in progress!

#
#        if v=='sfr' or v=='sSFR':
#            cos = halo.Cosmology()
#            mstSpeagle = np.power(10.0, np.linspace(9.7,11.1, 10)) # just linear - no need for a large number of samples
#            psiSpeagle = (0.84 - 0.026*cos.age(z[ti]))*np.log10(mstSpeagle) - (6.51 - 0.11*cos.age(z[ti]))
#            fac = 1.0
#            if v=='sSFR':
#                fac = mstSpeagle*1.0e-9
#            ax.fill_between(mstSpeagle, np.power(10.0,psiSpeagle-0.28)/fac, np.power(10.0, psiSpeagle+0.28)/fac, facecolor='gray', alpha=0.3)
#            theLabel=None
#            if z[ti]<0.5 or movie:
#                theLabel='Speagle14'
#            ax.plot(mstSpeagle, np.power(10.0, psiSpeagle)/fac, c='gray', label=theLabel)
#            
#            def plotBrinchmann(theAx, theLS, theColor, theLabel=None, fill=False, specific=True):
#                brinchmannMode = np.loadtxt('brinchmann04_msmode.csv', delimiter=',')
#                brinchmann975 = np.loadtxt('brinchmann04_ms975.csv', delimiter=',')
#                fbmode = interp1d( brinchmannMode[:,0], brinchmannMode[:,1], kind='linear' )
#                fb975 = interp1d( brinchmann975[:,0], brinchmann975[:,1], kind='linear' )
#                brinchmannm = np.linspace( 6.51, 11.88, 1000 ) 
#                fac = 1.0
#                if specific:
#                    fac = np.power(10.0, brinchmannm)*1.0e-9
#                theAx.plot( np.power(10.0, brinchmannm), np.power(10.0, fbmode(brinchmannm))/fac, c=theColor, ls=theLS, label=theLabel)
#                delt = fb975(brinchmannm) - fbmode(brinchmannm)
#                if fill:
#                    theAx.fill_between( np.power(10.0, brinchmannm), np.power(10.0, fbmode(brinchmannm)-delt/2.0)/fac, np.power(10.0, fbmode(brinchmannm)+delt/2.0)/fac , facecolor=theColor, alpha=0.3)
#
#
#            theLabel=None
#            theLS='--'
#            if z[ti]<0.5 or movie:
#                theLabel='Brinchmann04'
#                theLS = '-'
#            plotBrinchmann(ax, theLS, 'orange', theLabel=theLabel, fill=(not theLabel is None), specific=(v=='sSFR'))
#
#
#            thisMst = np.power(10.0, np.linspace(10,12,100))
#            whitaker12 = np.power(10.0, -1.12 + 1.14*z[ti] - 0.19*z[ti]*z[ti] - (0.3+0.13*z[ti])*(np.log10(thisMst)-10.5)) # Gyr^-1 -- Genzel+15 eq 1
#            mmin=None
#            if 2.0<z[ti] and z[ti]<2.5:
#                a = -19.99
#                sa = 1.87
#                b = 3.44
#                sb = 0.36
#                c = -0.13
#                sc = 0.02
#                mmin= 9.3
#                mmax = 11.4
#            if 1.5<z[ti] and z[ti]<2.0:
#                a = -24.04
#                sa = 2.08
#                b = 4.17
#                sb = 0.40
#                c = -0.16
#                sc = 0.02
#                mmin= 9.2
#                mmax = 11.5
#            if 1.0<z[ti] and z[ti]<1.5:
#                a = -26.03
#                sa = 1.69 
#                b = 4.62
#                sb = 0.34
#                c = -0.19
#                sc = 0.02
#                mmin= 8.8
#                mmax = 11.3
#            # Technically the following should only apply for z[ti]>0.5, but I want something to be plotted at z<0.5
#            if z[ti]<1.0:
#                a = -27.40
#                sa = 1.91
#                b = 5.02
#                sb = 0.39
#                c = -0.22
#                sc = 0.02
#                mmin= 8.4
#                mmax = 11.2
#
#            if not mmin is None and (movie or iti==0 ):
#                mwhitaker = np.power(10.0, np.linspace(mmin,mmax,100))
#                thisLabel=None
#                thisLS='-'
#                if movie or iti==0 and z[ti]<0.5:
#                    thisLabel='Whitaker14'
#                if z[ti]<0.5:
#                    thisLS='--'
#                if v=='sfr':
#                    ax.plot( mwhitaker, np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0), c='b', label=thisLabel, lw=1, ls=thisLS )
#                    if thisLS!='--':
#                        ax.fill_between( mwhitaker, np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 - 0.34), np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 + 0.34), facecolor='b', label=thisLabel )
#                if v=='sSFR':
#                    ax.plot( mwhitaker, (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0))/mwhitaker * 1.0e9, c='b', label=thisLabel, lw=1, ls=thisLS, alpha=0.3 )
#                    if thisLS!='--':
#                        ax.fill_between( mwhitaker, (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 - 0.34))/mwhitaker * 1.0e9, (np.power(10.0, a + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0 + 0.34))/mwhitaker * 1.0e9, facecolor='b', alpha=0.3 )
#                #ax.plot( mwhitaker, np.power(10.0, a+sa + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#                #ax.plot( mwhitaker, np.power(10.0, a-sa + b*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#                #ax.plot( mwhitaker, np.power(10.0, a + (b+sb)*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#                #ax.plot( mwhitaker, np.power(10.0, a + (b-sb)*np.log10(mwhitaker) + c*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#                #ax.plot( mwhitaker, np.power(10.0, a + b*np.log10(mwhitaker) + (c+sb)*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#                #ax.plot( mwhitaker, np.power(10.0, a + b*np.log10(mwhitaker) + (c-sb)*(np.log10(mwhitaker))**2.0), c='b', lw=1, ls=':' )
#            if z[ti]<2:
#                lilly13 = 0.117 * np.power(thisMst/3.16e10, -0.1) * np.power(1.0+z[ti],3.0)
#            else:
#                lilly13 = 0.5 * np.power(thisMst/3.16e10, -0.1) * np.power(1.0+z[ti],1.667)
#
#            thisLabel=None
#            if movie or iti==0 and z[ti]<0.5:
#                thisLabel='Whitaker12'
#            if v=='sfr':
#                ax.plot(thisMst, whitaker12*thisMst*1.0e-9, c='k',label=thisLabel)
#            if v=='sSFR':
#                ax.plot(thisMst, whitaker12, c='k',label=thisLabel)
#            if movie or iti==0 and z[ti]<0.5:
#                thisLabel='Lilly13'
#            if v=='sfr':
#                ax.plot(thisMst, lilly13*thisMst*1.0e-9, c='r', label=thisLabel)
#            if v=='sSFR':
#                ax.plot(thisMst, lilly13, c='r', label=thisLabel)
#            labelled=True
#        if v=='stZ':
#            MstKirby = np.power(10.0, np.linspace(3,9,20))
#            ZKirby = np.power(10.0, -1.69 + 0.30* np.log10(MstKirby/1.0e6))
#
#            # These are actually from Gallazzi05. (Tremonti is on the paper, so I don't feel that bad)
#            MstTremonti = np.power(10.0, np.array([8.91, 9.11, 9.31, 9.51, 9.72, 9.91, 10.11, 10.31, 10.51, 10.72, 10.91, 11.11, 11.31, 11.51, 11.72, 11.91]))
#            ZTremontiMed = np.power(10.0, np.array([-0.60, -0.61, -0.65, -0.61, -0.52, -0.41, -0.23, -0.11, -0.01, 0.04, 0.07, 0.10, 0.12, 0.13, 0.14, 0.15]))
#            ZTremonti16 = np.power(10.0, np.array([-1.11, -1.07, -1.10, -1.03, -0.97, -0.90, -0.80, -0.65, -0.41, -0.24, -0.14, -0.09, -0.06, -0.04, -0.03, -0.03]))
#            ZTremonti84 = np.power(10.0, np.array([-0.00, -0.00, -0.05, -0.01, 0.05, 0.09, 0.14, 0.17, 0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.30]))
#            thisLS = '--'
#            thisLabelKirby = None
#            thisLabelTremonti = None
#            if z[ti]<0.3:
#                thisLS = '-'
#                thisLabelKirby = 'Kirby13'
#                thisLabelTremonti = 'Gallazzi05'
#                ax.fill_between( MstKirby, ZKirby*np.power(10.0, -0.17), ZKirby*np.power(10.0, 0.17), facecolor='b', alpha=0.3)
#                ax.fill_between( MstTremonti, ZTremonti16, ZTremonti84, facecolor='r', alpha=0.3 )
#            ax.plot( MstKirby, ZKirby, c='b', ls=thisLS, label=thisLabelKirby)
#            ax.plot( MstTremonti, ZTremontiMed, c='r', ls=thisLS, label=thisLabelTremonti)
#            labelled=True
#
#        if v=='sSFR' and False:
#            whitaker12 = np.power(10.0, -1.12 + 1.14*z[ti] - 0.19*z[ti]*z[ti] - (0.3+0.13*z[ti])*(np.log10(mst)-10.5)) # Gyr^-1 -- Genzel+15 eq 1
#            if z[ti]<2:
#                lilly13 = 0.117 * np.power(mst/3.16e10, -0.1) * np.power(1.0+z[ti],3.0)
#            else:
#                lilly13 = 0.5 * np.power(mst/3.16e10, -0.1) * np.power(1.0+z[ti],1.667)
#            thisLabel=None
#            if movie or iti==0 and z[ti]<0.5:
#                thisLabel='Whitaker12'
#            ax.plot(mst, whitaker12, c='k', label=thisLabel)
#            if movie or iti==0 and z[ti]<0.5:
#                thisLabel='Lilly13'
#            ax.plot(mst, lilly13, c='r', label=thisLabel)
#            labelled=True
#        f0 = 1.0/(1.0 + np.power(mst/10.0**9.15,0.4)) # from Hayward & Hopkins (2015) eq. B2
#        tau4 = (12.27-t[ti])/(12.27+1.60) # fractional lookback time at z=4
#        fgz4 = f0*np.power(1.0 - tau4*(1.0-np.power(f0,1.5)), -2.0/3.0)
#        if v=='fg' or v=='fghm' or v=='fgsf' or v=='fgh2' or v=='fg2hm':
#            #thisLabel = 'Hayward16'
#            #if not movie and iti>0:
#                thisLabel=None
#            #ax.plot(mst,fgz4,c='k', label=thisLabel)
#            #labelled=True
#        if v=='fgh2':
#            thisMst = np.power(10.0, np.linspace(10.0,11.5, 100))
#            def molToStarGenzel( af2, xif2, xig2, xih2):
#                return np.power(10.0, af2 + xif2*np.log10(1.0+z[ti]) + xig2*0.0 + xih2*np.log10(thisMst/10.0**10.77))
#            def ratioToFraction(rat):
#                ## f_g = Mg/(Mg+M*) = 1/(1+M*/Mg) = 1/(1+1/rat)
#                return 1.0/(1.0 + 1.0/rat)
#                
#            thisLabel1 = 'Genzel15 (Lilly)'
#            thisLabel2 = 'Genzel15 (FMR)'
#            if not movie and z[ti]>0.5:
#                thisLabel1=None
#                thisLabel2=None
#            ax.plot(thisMst, ratioToFraction(molToStarGenzel(-0.98,2.65,0.5,-0.25)), c='green', label=thisLabel1)
#            ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
#            #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
#            #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
#            labelled=True
#        if v=='fghi':
#            thisMst = np.power(10.0, np.linspace(7.0,11.5))
#            thisMHI = np.power(10.0, -0.43*np.log10(thisMst) + 3.75) * thisMst
#            fgThis = thisMHI/(thisMHI+thisMst)
#            thisLS = '-'
#            thisLabel = 'Papastergis12'
#            if not movie and z[ti]>0.5:
#                thisLS = '--'
#                thisLabel = None
#            ax.plot(thisMst, fgThis, c='red', label=thisLabel, ls=thisLS)
#            labelled=True
#        if v=='gasToStellarRatio':
#            thisMst = np.power(10.0, np.linspace(7.0,11.5))
#            thisMHI = np.power(10.0, -0.43*np.log10(thisMst) + 3.75) * thisMst
#            thisLS = '-'
#            thisLabel = 'Papastergis12'
#            if not movie and z[ti]>0.5:
#                thisLS = '--'
#                thisLabel = None
#            ax.plot(thisMst, thisMHI/thisMst, c='red', label=thisLabel, ls=thisLS)
#            labelled=True
#        if v=='gasToStellarRatioHI':
#            def plotPeeples(theAx, theLS, theColor, theLabel=None, fill=False):
#                peeplesMode = np.loadtxt('peeples11_fgmode.csv', delimiter=',')
#                peeples84 = np.loadtxt('peeples11_fg68.csv', delimiter=',')
#                peeplesMass = peeplesMode[:,0]
#                theAx.plot( np.power(10.0, peeplesMass), np.power(10.0, peeplesMode[:,1]), c=theColor, ls=theLS, label=theLabel)
#                delt = peeples84[:,1] - peeplesMode[:,1]
#                if fill:
#                    theAx.fill_between( np.power(10.0, peeplesMass), np.power(10.0, peeplesMode[:,1]-delt), np.power(10.0, peeplesMode[:,1]+delt), facecolor=theColor, alpha=0.3 )
#
#            theLabel=None
#            theLS='--'
#            if z[ti]<0.5 or movie:
#                theLabel='Peeples11'
#                theLS = '-'
#            plotPeeples(ax, theLS, 'blue', theLabel=theLabel, fill=(not theLabel is None))
#
#            thisMst = np.power(10.0, np.linspace(7.0,11.5))
#            thisMHI = np.power(10.0, -0.43*np.log10(thisMst) + 3.75) * thisMst
#            fgThis = thisMHI/thisMst
#            thisLS = '-'
#            thisLabel = 'Papastergis12'
#            if not movie and z[ti]>0.5:
#                thisLS = '--'
#                thisLabel = None
#            ax.plot(thisMst, fgThis, c='red', label=thisLabel, ls=thisLS)
#            labelled=True
#        if v=='gasToStellarRatioH2':
#            thisMst = np.power(10.0, np.linspace(10.0,11.5, 100))
#            def molToStarGenzel( af2, xif2, xig2, xih2):
#                return np.power(10.0, af2 + xif2*np.log10(1.0+z[ti]) + xig2*0.0 + xih2*np.log10(thisMst/10.0**10.77))
#            def ratioToFraction(rat):
#                ## f_g = Mg/(Mg+M*) = 1/(1+M*/Mg) = 1/(1+1/rat)
#                return 1.0/(1.0 + 1.0/rat)
#                
#            thisLabel1 = 'Genzel15 (Lilly)'
#            thisLabel2 = 'Genzel15 (FMR)'
#            if not movie and z[ti]>0.5:
#                thisLabel1=None
#                thisLabel2=None
#            ax.plot(thisMst, molToStarGenzel(-0.98,2.65,0.5,-0.25), c='green', label=thisLabel1)
#            ax.plot(thisMst, molToStarGenzel(-1.05,2.6,0.54,-0.41), c='orange', label=thisLabel2)
#            #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
#            #ax.plot(thisMst, ratioToFraction(molToStarGenzel(-1.05,2.6,0.54,-0.41)), c='orange', label=thisLabel2)
#
#            thisLabel3=None
#            thisLS = '--'
#            if z[ti]<0.5:
#                thisLS = '-'
#                thisLabel3 = 'Saintonge11'
#            saintongeM = np.power(10.0, np.linspace(10.3, 11.4, 20))
#            saintongeFr = np.power(10.0, -1.607 - 0.455*(np.log10(saintongeM) - 10.70))
#            uncertainty = 10.0**.331 # uncertainty on b in Table 3 (undetected@upper limits) of Saintonge11.
#            if z[ti]<0.5:
#                ax.fill_between(saintongeM, saintongeFr/uncertainty, saintongeFr*uncertainty, facecolor='r', alpha=0.3)
#            ax.plot(saintongeM, saintongeFr, ls=thisLS, label=thisLabel3, c='r')
#            labelled=True
#        #if v=='gasToStellarRatio' or v=='gasToStellarRatioH2':
#        #    thisLabel = 'Hayward16'
#        #    if not movie and iti>0:
#        #        thisLabel=None
#        #    ax.plot(mst, fgz4/(1-fgz4), c='k', label=thisLabel)
#        #    labelled=True
#        if v=='vPhiOuter':
#            thisLabel = 'Hayward16'
#            thisLS = '-'
#            if not movie and z[ti]>0:
#                thisLabel=None
#                thisLS = '--'
#            ax.plot(mst,  147*np.power(mst/1.0e10, 0.23),c='k', label=thisLabel, ls=thisLS) # Hayward & Hopkins eq. B1
#            labelled=True
#        if v=='vPhi22':
#            thisLabel = 'Miller11'
#            thisLS = '-'
#            if not movie and z[ti]>1.25:
#                thisLabel=None
#                thisLS = '--'
#            mstMiller = np.power(10.0, np.linspace(9.0, 11.5, 20))
#            vMiller = np.power(10.0, (np.log10(mstMiller) - 1.718)/3.869)
#            ax.plot(mstMiller, vMiller,c='k', label=thisLabel, ls=thisLS) 
#            if thisLS=='-':
#                ax.fill_between(mstMiller, vMiller/10**0.058, vMiller*10**0.058,facecolor='k', alpha=0.3) 
#            labelled=True
#
#
#if xvar=='gbar' and v=='gtot':
#    gbars = np.power(10.0, np.linspace(-12,-8, 100))
#    gdagger = 1.20e-10
#    gobss = gbars/(1.0 - np.exp(-np.sqrt(gbars/gdagger)))
#    thisLabel=None
#    thisLS = '--'
#    if z[ti]<0.5:
#        thisLabel='McGaugh16'
#        thisLS = '-'
#    if histFlag:
#        ax.plot( np.log10(gbars), np.log10(gobss), label=thisLabel, ls=thisLS, c='b', lw=3 )
#        ax.plot( np.log10(gbars), np.log10(gbars), ls='--', c='gray', lw=3)
#    else:
#        ax.plot( gbars, gobss, label=thisLabel, ls=thisLS, c='b', lw=3 )
#        ax.plot( gbars, gbars, ls='--', c='gray', lw=3)
#    labelled=True


# Now let's put in some data.
datasets = []
defineMoster(0)
