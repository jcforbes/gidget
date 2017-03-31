from readoutput import *
from balanceplot import balance
import argparse
import cleanup
import pdb



def makeThePlots(args):
    if args.time is False and args.radial is False and args.scaled is False and args.mass==0 and args.balance is False:
        print "Warning: you did not ask me to produce any plots! Use python plots.py -h to look at the options available."
    print "plots.py in makeThePlots: args.models: ",args.models
    balanceArgs=[]
    AMargs=[]
    MONDargs=[]
    if(args.balance):
        balanceArgs=['colTr','colAccr','colsfr','dcoldt','MassLoadingFactor','Mdot','colREC']
    if(args.angularMomentum):
        AMargs = ['colTr','r','vPhi','dA']
    if(args.quick):
        MONDargs = ['gbar', 'gtot', 'hGas', 'sSFRRadial', 'rxl', 'colstNormalizedKravtsov', 'colNormalizedKravtsov', 'colHI', 'colH2', 'colst', 'fH2', 'vPhi', 'sigstR', 'sigstZ', 'ageRadial', 'colsfr', 'Z']
    for modelName in args.models:
        print "Beginning to analyze experiment ",modelName
        theExp = Experiment(modelName)
        print "Reading in the experiment keeping: ", args.vsr + balanceArgs + AMargs +MONDargs
        theExp.read(args.vsr+balanceArgs+AMargs+MONDargs, keepStars=(args.stellarPops or args.quick), computeFit=args.fit)
        nts = int(theExp.models[0].p['Noutputs']+1)
        tis = [nts/5,nts/2,nts]
        if args.scalings or args.genzel:
            theExp.storeScalingRelation('MS', 'mstar','sfr')
            theExp.storeScalingRelation('MZR', 'mstar','integratedZ')
            theExp.storeScalingRelation('MFG', 'mstar','gasToStellarRatio')
            theExp.storeScalingRelation('TF', 'mstar','vPhiOuter')
            theExp.storeScalingRelation('MsTd', 'mstar','tdep')
            theExp.storeScalingRelation('MsTdH2', 'mstar','tDepH2')
            theExp.storeScalingRelation('Rg', 'mstar','halfMassGas')
            theExp.storeScalingRelation('mRho', 'mstar','rho1')
        theExp.assignIndices()
        for i,rankby in enumerate(args.rankby):
            tti=None
            if(args.rbz[i]>=0):
                tti,_ = Nearest(theExp.models[0].var['z'].sensible(),args.rbz[i])
            theExp.rankBy(var=rankby,timeIndex=tti)
        stepsize = args.step
        for cb in args.colorby:
            if(args.time):
                if(args.percentiles):
                    per = [2.5, 16, 50, 84, 97.5]
                    theExp.timePlot(colorby=cb, perc=per,vsz=False)
                    theExp.timePlot(colorby=cb, perc=per,vsz=True)
                else:
                    theExp.timePlot(colorby=cb,vsz=False)
                    theExp.timePlot(colorby=cb,vsz=True)
            if(args.radial):
                theExp.radialPlot(timeIndex=range(1,nts+1,stepsize)+[nts],variables=args.vsr,colorby=cb,logR=args.logR)
            if(args.scaled):
                theExp.radialPlot(timeIndex=range(1,nts+1,stepsize)+[nts],variables=args.vsr,scaleR=True,colorby=cb,logR=args.logR)
            #if(args.mass):
            #    theExp.ptMovie(timeIndex=range(1,202,stepsize)+[201],prev=args.prev,colorby=cb)
            #    theExp.ptMovie(timeIndex=range(1,202,stepsize)+[201],xvar='Mh',prev=args.prev,colorby=cb)
            if len(args.mass)!=0:
                for xv in args.mass:
                    theExp.ptMovie(timeIndex=range(1,nts+1,stepsize)+[nts],xvar=xv,prev=args.prev,colorby=cb,movie=True)
                    theExp.ptMovie(timeIndex=range(1,nts+1,stepsize)+[nts],xvar=xv,prev=0,colorby=cb,movie=False)
        if len(args.mass)!=0 and args.snapshot:
            for xv in args.mass:
                theExp.ptMovie(timeIndex=tis,xvar=xv,prev=0,colorby='t',movie=False)
                theExp.ptMovie(timeIndex=tis,xvar=xv,yvar=args.vsr,prev=0,colorby='t',movie=False)
        if(args.radial and args.snapshot):
            theExp.radialPlot(timeIndex=tis,variables=args.vsr,colorby='t',logR=args.logR,movie=False)
        if(args.scaled and args.snapshot):
            theExp.radialPlot(timeIndex=tis,variables=args.vsr,scaleR=True,colorby='t',logR=args.logR,movie=False)
        if args.snapshot:
            pass
            #theExp.hist1d(timeIndex=tis, vars=None, movie=False)

        # END loop over for cb in colorby
        # OK, this is just a test..
        #theExp.ptMovie(timeIndex=tis,xvar='colsfr',yvar=['MassLoadingFactor'],prev=0,colorby='t',movie=False)
        #theExp.ptMovie(timeIndex=tis,xvar='hGas',yvar=['MassLoadingFactor'],prev=0,colorby='t',movie=False)
        #theExp.ptMovie(timeIndex=tis,xvar='Mh',yvar=['MassLoadingFactor','mstar','fg','integratedMLF','integratedZ'],prev=0,colorby='t',movie=False)
        #theExp.ptMovie(timeIndex=tis,xvar='x3',yvar=['MassLoadingFactor'],prev=0,colorby='t',movie=False)
        #theExp.ptMovie(timeIndex=tis,xvar='colsfr',yvar=['colTr','colAccr','colOut'],prev=0,colorby='t',movie=False)
        #theExp.ptMovie(timeIndex=tis,xvar='r',yvar=['col','colst','Z','vPhi','Q','MassLoadingFactor','hGas','colsfr','fH2','fgRadial','colTr','colAccr','colOut','equilibrium'],prev=0,colorby='t',movie=False)
        theExp.krumholzAnalysis()
        if args.genzel:
            theExp.globalGenzelAnalysis()
        if args.quick:
            theExp.quickCheck()
        if args.fraction:
            theExp.globalFractionAnalysis()
            theExp.globalFractionAnalysis(funcs=[np.std])
        if args.stellarPops:
            theExp.plotAgeFuncs()
        if(args.percentiles):
            per = [2.5, 16, 50, 84, 97.5]
            if(args.radial):
                theExp.radialPlot(timeIndex=range(1,nts+1,stepsize)+[nts],variables=args.vsr,colorby=args.colorby[0],logR=args.logR,percentiles=per)
                if args.snapshot:
                    theExp.radialPlot(timeIndex=tis,variables=args.vsr,colorby=args.colorby[0],logR=args.logR,percentiles=per,movie=False)
            if(args.scaled):
                theExp.radialPlot(timeIndex=range(1,nts+1,stepsize)+[nts],variables=args.vsr,colorby=args.colorby[0],logR=args.logR,percentiles=per,scaleR=True)
        if(args.balance):
            #balance(theExp.models,timeIndex=range(1,nts+1,stepsize)+[nts],name=modelName,sortby=args.colorby[0],logR=args.logR, ncols=5, nrows=3)
            balance(theExp.models,timeIndex=range(1,nts+1,stepsize)+[nts],name=modelName,sortby=args.colorby[0],logR=args.logR, ncols=1, nrows=1)
        #theExp.customPlotPPD( expensive=True)
       # theExp.customPlotPPD( expensive=False)

        if args.angularMomentum:
            theExp.angularMomentumAnalysis()



if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Analyze data from GIDGET experiments.')
    parser.add_argument('models', metavar='experiment',type=str,nargs='+',
            help='Each experiment specified will be constructed and analyzed. For instance, rh01 rh02 rh03 will analyze as a single Experiment any experiment matching *rh01*, then any experiment matching *rh02*, etc. All of them could be analyzed together (e.g. all the models would appear simultaneously on the same plots) by giving rh0 as an argument (this would also include rh04, rh05,...).')
    parser.add_argument('--step',type=int,default=3,help='Stepsize among snapshots to make movies (default: every 3rd)')
    parser.add_argument('--time',dest='time',action='store_true',help="Make plots vs t.")
    parser.set_defaults(time=False)
    parser.add_argument('--radial',dest='radial',action='store_true',help="Make plots vs r.")
    parser.set_defaults(radial=False)
    parser.add_argument('--scaled',dest='scaled',action='store_true',help="Make plots vs r/racc.")
    parser.set_defaults(scaled=False)
    parser.add_argument('--balance',dest='balance',action='store_true',help="Make balance plots.")
    parser.set_defaults(balance=False)
    parser.add_argument('--angularMomentum',dest='angularMomentum',action='store_true',help="Make some simple plots to try to understand angular momentum.")
    parser.set_defaults(angularMomentum=False)
    #parser.add_argument('--mass',dest='mass',action='store_true',help="Make plots vs mstar")
    #parser.set_defaults(mass=False)
    parser.add_argument('--mass',type=str,nargs='+',default=[],help='List of x- variables to use in point movies. Eg. mstar, Mh, sSFR,...')
    parser.add_argument('--colorby',type=str,nargs='+',default=['Mh0'],help='List of variables to color points by in vs mstar plots. Default is halo mass at z=0')
    oldVsrDefaults = ['colsfr','colst','NHI','sig','col','Z','fH2','Mdot','MJeans','tDepRadial','tDepH2Radial','Q','Qg','Qst','fgRadial','equilibrium','colHI','colH2','vPhi','colTrPerAccr','hGas','hStars','sigstR','sigstZ', 'JColGas', 'JColStars','MassLoadingFactor']
    parser.add_argument('--vsr',type=str,nargs='+',default=oldVsrDefaults,help='List of variables to plot vs mstar. Default is '+repr(oldVsrDefaults))
    parser.add_argument('--prev',type=int,default=5,help='Number of previous points to plot in vsmstar movies.')
    parser.add_argument('--rankby',type=str,nargs='+',default=[],help='Sort the models according to these arguments.')
    parser.add_argument('--rbz',type=float,nargs='+',default=[],help='Sort at a particular redshift (use -1 to keep the full time information)')
    parser.add_argument('--logR',dest='logR',action='store_true',help="Use logarithmic radial coordinate in plots vs. r")
    parser.set_defaults(logR=False)
    parser.add_argument('--percentiles',dest='percentiles',action='store_true',help="In radial plots overplot percentiles.")
    parser.add_argument('--snapshot',dest='snapshot',action='store_true',help="Produce 1d histograms of all parameters and time variables.")
    parser.add_argument('--stellarPops',dest='stellarPops',action='store_true',help="Produce plots relating to the passive stellar pops")
    parser.add_argument('--genzel', dest='genzel', action='store_true',help="Produce plots relating to Genzel et al 2015")
    parser.add_argument('--fraction', dest='fraction', action='store_true',help="Produce heatmaps of quantites as fn of Mh and z")
    parser.add_argument('--fit', dest='fit', action='store_true',help="Fit the stellar column density profiles. May be time-consuming")
    parser.add_argument('--scalings', dest='scalings', action='store_true',help="Fit galaxy scaling relations - only really makes sense in runs where you have a decent range of masses and some source of variability between galaxies at a given mass.")
    parser.add_argument('--quick', dest='quick', action='store_true',help="Quickly check fit to scaling relations.")
    args = parser.parse_args()

    weNeed = len(args.rankby) - len(args.rbz) 
    if(weNeed>0):
        for i in weNeed:
            args.rbz.append(-1)

    makeThePlots(args)
    cleanup.moveFiles()



