from readoutput import *
from balanceplot import balance
from guppy import *

hp = hpy()

def fiducialAnalysis():
    test = Experiment('rg78')
    test.read()
    test.storeMS()
    #test.timePlot()
    #test.radialPlot(variables=['colsfr','colst','NHI','sig','col','Z','fH2','Mdot'])
    #test.ptMovie(timeIndex=range(1,202),prev=5)
    #test.ptMovie(timeIndex=range(1,202),prev=5,colorby='integratedZ')
    balance(test.models[0:16])
    #test.ptMovie(timeIndex=range(1,202),prev=5,colorby='deltaMS')
    test.ptMovie(timeIndex=range(1,202),prev=5,colorby='lambda')

def massStudy():
    massStudy = Experiment('rh45')
    massStudy.read()
    balance(massStudy.models,name='rh45')
    massStudy.timePlot()
    massStudy.radialPlot(timeIndex=range(1,202,2)+[201],variables=['colsfr','colst','NHI','sig','col','Z','fH2','Mdot'],scaleR=True)
    massStudy.radialPlot(timeIndex=range(1,202,2)+[201],variables=['colsfr','colst','NHI','sig','col','Z','fH2','Mdot'])

    print hp.heap()


if __name__=='__main__':
    massStudy()
    #fiducialAnalysis()


