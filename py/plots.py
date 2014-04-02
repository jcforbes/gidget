from readoutput import *
from balanceplot import balance
import argparse
import cleanup



def makeThePlots(args):
    if args.time is False and args.radial is False and args.scaled is False and args.mass==0 and args.balance is False:
        print "Warning: you did not ask me to produce any plots! Use python plots.py -h to look at the options available."
    print "plots.py in makeThePlots: args.models: ",args.models
    for modelName in args.models:
        print "Beginning to analyze experiment ",modelName
        theExp = Experiment(modelName)
        theExp.read()
        theExp.storeMS()
        theExp.assignIndices()
        for i,rankby in enumerate(args.rankby):
            tti=None
            if(args.rbz[i]>=0):
                tti,_ = Nearest(theExp.models[0].var['z'].sensible(),args.rbz[i])
            theExp.rankBy(var=rankby,timeIndex=tti)
        stepsize = args.step
        for cb in args.colorby:
            if(args.time):
                theExp.timePlot(colorby=cb)
            if(args.radial):
                theExp.radialPlot(timeIndex=range(1,202,stepsize)+[201],variables=args.vsr,colorby=cb,logR=args.logR)
            if(args.scaled):
                theExp.radialPlot(timeIndex=range(1,202,stepsize)+[201],variables=args.vsr,scaleR=True,colorby=cb,logR=args.logR)
            #if(args.mass):
            #    theExp.ptMovie(timeIndex=range(1,202,stepsize)+[201],prev=args.prev,colorby=cb)
            #    theExp.ptMovie(timeIndex=range(1,202,stepsize)+[201],xvar='Mh',prev=args.prev,colorby=cb)
            if len(args.mass)!=0:
                for xv in args.mass:
                    theExp.ptMovie(timeIndex=range(1,202,stepsize)+[201],xvar=xv,prev=args.prev,colorby=cb)
        if(args.balance):
            balance(theExp.models,timeIndex=range(1,202,stepsize)+[201],name=modelName,sortby=args.colorby[0])



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
    #parser.add_argument('--mass',dest='mass',action='store_true',help="Make plots vs mstar")
    #parser.set_defaults(mass=False)
    parser.add_argument('--mass',type=str,nargs='+',default=[],help='List of x- variables to use in point movies. Eg. mstar, Mh, sSFR,...')
    parser.add_argument('--colorby',type=str,nargs='+',default=['deltaMS','lambda','integratedZ'],help='List of variables to color points by in vs mstar plots. Default is deltaMS,lambda,integratedZ')
    parser.add_argument('--vsr',type=str,nargs='+',default=['colsfr','colst','NHI','sig','col','Z','fH2','Mdot','MJeans','ClumpMassPerDisk','tDepRadial','tDepH2Radial','Q','Qg','Qst','fgRadial','equilibrium'],help='List of variables to plot vs mstar. Default is colsfr,colst,NHI,sig,col,Z,fH2,Mdot,MJeans,ClumpMassPerDisk,tDepRadial,tDepH2Radial,Q,Qg,Qst,fgRadial,equilibrium')
    parser.add_argument('--prev',type=int,default=5,help='Number of previous points to plot in vsmstar movies.')
    parser.add_argument('--rankby',type=str,nargs='+',default=[],help='Sort the models according to these arguments.')
    parser.add_argument('--rbz',type=float,nargs='+',default=[],help='Sort at a particular redshift (use -1 to keep the full time information)')
    parser.add_argument('--logR',dest='logR',action='store_true',help="Use logarithmic radial coordinate in plots vs. r")
    parser.set_defaults(logR=False)
    args = parser.parse_args()

    weNeed = len(args.rankby) - len(args.rbz) 
    if(weNeed>0):
        for i in weNeed:
            args.rbz.append(-1)

    makeThePlots(args)
    cleanup.moveFiles()



