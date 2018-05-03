## This script aims to plot the results of a "suite" of GIDGET experiments, where each experiment contains N different ~evenly spaced masses, and each has one parameter that has been varied.
## The plot of any given quantity should find either the M most important physical parameters affecting it, or maybe parameters above some threshold distance in the space. 
## The plot will show the fiducial result (i.e. the one best matching observations), and arrows indicating the directions the points may move under the different values of the parameters.
## Expected pitfalls: Each run will have a somewhat different accretion history. Some models may have failed.

import readoutput as ro
import matplotlib.pyplot as plt
import numpy as np
import glob
import pdb

def runner(basename, fidname, fh=0.0):
    fid = ro.Experiment(fidname)
    MONDargs = ['gbar', 'gtot', 'hGas', 'sSFRRadial', 'rxl', 'colstNormalizedKravtsov', 'colNormalizedKravtsov', 'colHI', 'colH2', 'colst', 'fH2', 'vPhi', 'sigstR', 'sigstZ', 'ageRadial', 'colsfr', 'Z', 'sig']
    fid.read(keepOnly= MONDargs, fh=fh, keepStars=True)
    experimentNameList = sorted(glob.glob( '../analysis/'+basename+'*'))
    experimentList = []
    for exp in experimentNameList:
        experimentList.append( ro.Experiment(exp[12:]) ) # 12 gets rid of the "../analysis/" part of the string used to find all the relevant models.
        experimentList[-1].read(keepOnly= MONDargs, keepStars=True)
    if len(experimentList)<1:
        pdb.set_trace() # we didn't find any experiments.. what's going on?


    zinds = [ ro.Nearest( fid.models[0].var['z'].sensible(), z)[0] for z in [0,1,2,3,4]]

    ###  raccRvir, rstarRed, rgasRed, fg0mult, muColScaling,
    ###   muFgScaling, muNorm, muMhScaling, ZIGMfac, zmix,
    ###   eta, Qf, alphaMRI, epsquench, accCeiling, 
    ###   conRF, kZ, xiREC, epsff, scaleAdjust, 
    ###   mquench, enInjFac, fh = emceeparams
    featureNames = [r'$\alpha_r$', r'$\alpha_{r,*,0}$',  r'$\alpha_{r,g,0}$',  r'$\chi_{f_{g,0}}$', r'$\alpha_\Sigma$', 
            r'$\alpha_{f_g}$', r'$\mu_0$', r'$\alpha_{M_h}$', r'$\chi_{Z_\mathrm{IGM}}$', r'$\xi_\mathrm{acc}$', 
            r'$\eta$', r'$Q_f$', r'$\alpha_\mathrm{MRI}$', r'$\epsilon_\mathrm{quench}$', r'$\epsilon_\mathrm{ceil}$', 
            r'$\alpha_\mathrm{con}$', r'$k_Z$', r'$\xi$', r'$\epsilon_\mathrm{ff}$', r'$\Delta\beta$', r'$M_Q$', r'$E_\mathrm{inj}$']
    colorNames = ['k', 'olive', 'b', 'pink', 'lightblue', 'darkgreen', 'lightgreen', 'gray', 'orange', 'yellow', 'magenta', 'cyan', 'maroon', 'purple', 'lightslategray', 'moccasin', 'fuchsia', 'mediumorchid', 'cadetblue', 'darkkhaki', 'chocolate', 'red', 'darkorange']*2
    

    fig,ax = plt.subplots(5,4, figsize=(12,12))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    colorby = 'Mh0'
    for j in range(4):
        fid.ptMovie(xvar='mstar', yvar=['sSFR'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6)
        expand(ax[0,j])
        singlePanel('mstar', 'sSFR', zinds[j], ax[0,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['sfZ'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6)
        expand(ax[1,j])
        singlePanel('mstar', 'sfZ', zinds[j], ax[1,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['stZ'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6)
        expand(ax[2,j])
        singlePanel('mstar', 'stZ', zinds[j], ax[2,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['gasToStellarRatioH2'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6)
        expand(ax[3,j])
        singlePanel('mstar', 'gasToStellarRatioH2', zinds[j], ax[3,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['gasToStellarRatioHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[4,j], textsize=6)
        expand(ax[4,j])
        singlePanel('mstar', 'gasToStellarRatioHI', zinds[j], ax[4,j], featureNames, colorNames, fid, experimentList, 1.0, None)

        ax[2,j].text(1.0e10, 2.0e-2, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[4,j].get_xaxis().set_ticks([1.0e5, 1.0e7, 1.0e9, 1.0e11])
    plt.savefig(fidname+'_directionalCalibrationLo1.pdf')
    plt.close(fig)



    fig,ax = plt.subplots(4,4, figsize=(12,11))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    for j in range(4):
        fid.ptMovie(xvar='mstar', yvar=['halfMassStars'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6)
        expand(ax[0,j], bottom=1.0)
        singlePanel('mstar', 'halfMassStars', zinds[j], ax[0,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['vPhi22'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6)
        expand(ax[1,j], bottom=1.0)
        singlePanel('mstar', 'vPhi22', zinds[j], ax[1,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['c82'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6)
        expand(ax[2,j], bottom=1.0)
        singlePanel('mstar', 'c82', zinds[j], ax[2,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        fid.ptMovie(xvar='mstar', yvar=['Sigma1'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6)
        expand(ax[3,j], bottom=1.0)
        singlePanel('mstar', 'Sigma1', zinds[j], ax[3,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        ax[1,j].text(1.0e10, 80.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(4):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<3:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[3,j].get_xaxis().set_ticks([1.0e5, 1.0e7, 1.0e9, 1.0e11])
    plt.savefig(fidname+'_directionalCalibrationLo2.pdf')
    plt.close(fig)




    #percentiles=[16,50,84]
    percentiles=None
    fig,ax = plt.subplots(5,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        fid.radialPlot(timeIndex=[zinds[j]],variables=['colst'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[0,j])
        singlePanelR( 'colst', zinds[j], ax[0,j], featureNames, colorNames, fid, experimentList, 1.0, None)

        fid.radialPlot(timeIndex=[zinds[j]],variables=['colsfr'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[1,j])
        singlePanelR( 'colsfr', zinds[j], ax[1,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        
        fid.radialPlot(timeIndex=[zinds[j]],variables=['sSFRRadial'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[2,j])
        singlePanelR( 'sSFRRadial', zinds[j], ax[2,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        
        fid.radialPlot(timeIndex=[zinds[j]],variables=['colH2'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[3,j])
        singlePanelR( 'colH2', zinds[j], ax[3,j], featureNames, colorNames, fid, experimentList, 1.0, None)

        fid.radialPlot(timeIndex=[zinds[j]],variables=['colHI'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[4,j])
        singlePanelR( 'colHI', zinds[j], ax[4,j], featureNames, colorNames, fid, experimentList, 1.0, None)
        #ax[0,j].text(1.0e12, 1.0e7, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
    plt.savefig(fidname+'_directionalCalibration3.pdf')
    plt.close(fig)


def expand(axIn, bottom=10.0, top=3.2, right=10.0, left=1.0):
    ylim = axIn.get_ylim()
    xlim = axIn.get_xlim()
    axIn.set_xlim(xlim[0]/left, xlim[1]*right)
    axIn.set_ylim(ylim[0]/bottom, ylim[1]*top)

def getJindices(experiment, referenceMh0):
    indicesJ=[]
    for j, Mh0 in enumerate(referenceMh0):
        for model in experiment.models:
            #print np.log10(model.p['Mh0']), np.log10(Mh0), j, 
            if np.abs( np.log10( model.p['Mh0'] / Mh0)) < 0.0001  and j not in indicesJ:
                indicesJ.append(j)
    return indicesJ

def singlePanelR( yvar, zind, axIn, featureNames, colorNames, fidExper, experimentList, thresholdDiff, labelledAlready ):
    #### overallVar,_,overallLog,overallRange = self.constructQuantity(v)
    xvar = 'rx'
    maxnm = 0 # maximum number of models in given experiment. In principle they should all be the same, but some number of models may have failed in any given experiment.
    referenceMh0s = []
    for experiment in experimentList:
        nm = len(experiment.models)
        if nm>maxnm:
            maxnm = nm
            referenceMh0s = [model.p['Mh0'] for model in experiment.models]

    fidPts = np.zeros((2, maxnm)) # x&y for each model
    newPts = np.zeros((2, maxnm, len(experimentList))) # x&y for each model for each experiment
    diffs = np.zeros((maxnm, len(experimentList))) - 1.0 # one for each experiment and model; -1 if there was no corresponding model between experiments)
    jFid = getJindices(fidExper, referenceMh0s)
    rxrefs = [1.0, 2.0, 3.0]
    for k, rxref in enumerate(rxrefs):
        for i, experiment in enumerate(experimentList):
            jIndices = getJindices(experiment, referenceMh0s)
            counterFid = 0
            counterExper = 0
            for j in range(maxnm):
                if j in jIndices and j in jFid:
                    ### if the experiment in question actually has a model at the appropriate Mh0, grab it.

                    # Find our reference pt:
                    fidLoc =  ro.Nearest( fidExper.models[counterFid].var['rx'].sensible(timeIndex=zind), rxref)[0] # the index closest to rx=1 in the fiducial model
                    experLoc =  ro.Nearest( experiment.models[counterExper].var['rx'].sensible(timeIndex=zind), rxref)[0] # the index closest to rx=1 in the experimental model

                    try:
                        fidPt = ( rxref,  fidExper.models[counterFid].var[yvar].sensible(timeIndex=zind, locIndex=fidLoc) )
                        fidPts[:,j] = np.array([fidPt[0], fidPt[1]])

                        newPt = ( rxref,  experiment.models[counterExper].var[yvar].sensible(timeIndex=zind, locIndex=experLoc) )
                        newPts[:,j,i] = np.array([newPt[0],newPt[1]])
                    except:
                        pdb.set_trace()

                    if experiment.models[0].var[xvar].log:
                        dx = np.log10(newPt[0]/fidPt[0])
                    else:
                        dx = newPt[0] - fidPt[0]
                    if experiment.models[0].var[yvar].log:
                        dy = np.log10(newPt[1]/fidPt[1])
                    else:
                        dy = newPt[1] - fidPt[1]

                    thisDiff = np.sqrt( dx*dx + dy*dy )
                    diffs[j, i] = thisDiff

                if j in jIndices:
                    counterExper+=1
                if j in jFid:
                    counterFid+=1

        meanDiffs = np.zeros(len(experimentList)) ## the avg diff for each experiment
        for i in range(len(experimentList)):
            counter=0
            for j in range(maxnm):
                if diffs[j,i]>-0.5: 
                    # only include diffs where there's a model for both the fiducial and experimental case
                    meanDiffs[i] += diffs[j,i]
                    counter = counter+1
            if counter>0:
                meanDiffs[i] = meanDiffs[i]/float(counter)

        zipped = zip(meanDiffs, range(len(meanDiffs)))
        zsorted = sorted(zipped, key=lambda x: -x[0])
        highDiffInds = [zsorted[i][1] for i in range(3)]
        for i in range(len(highDiffInds)):
            jIndices = getJindices( experimentList[highDiffInds[i]], referenceMh0s )
            offset = np.random.normal()*0.2
            for j in range(maxnm):
                if j in jIndices and j in jFid and j%5==0:
                    axIn.plot( [fidPts[0,j]*1.003**j+offset, newPts[0,j, highDiffInds[i]]*1.003**j+offset], [fidPts[1,j], newPts[1,j,highDiffInds[i]]], lw=1, alpha=0.6, c=colorNames[highDiffInds[i]] )
                if j in jIndices and j in jFid and j==highDiffInds[i] or j==highDiffInds[i]+1:
                    #axIn.text(newPts[0,j,highDiffInds[i]], newPts[1,j,highDiffInds[i]], featureNames[highDiffInds[i]], color=colorNames[highDiffInds[i]])
                    bottomLeft = [axIn.get_xlim()[0], axIn.get_ylim()[0]]
                    bottomRight= [axIn.get_xlim()[1], axIn.get_ylim()[0]]
                    upperRight = [axIn.get_xlim()[1], axIn.get_ylim()[1]]
                    dy = axIn.get_ylim()[1] / axIn.get_ylim()[0]
                    dx = axIn.get_xlim()[1] / axIn.get_xlim()[0]
                    if k==0:
                        axIn.text( bottomRight[0]/dx**0.35, bottomRight[1]*dy**(0.04+0.075*i), featureNames[highDiffInds[i]] + ' '+str( np.round(meanDiffs[highDiffInds[i]],2) ), color=colorNames[highDiffInds[i]])


#                if thisDiff>thresholdDiff:
#                    axIn.annotate('', xy=(newPt[0],newPt[1]), xytext=(fidPt[0],fidPt[1]), arrowprops=dict(color=colorNames[i], width=1, headwidth=4) )
#                    if counterFid==0:


    return None 

def singlePanel(xvar, yvar, zind, axIn, featureNames, colorNames, fidExper, experimentList, thresholdDiff, labelledAlready ):

    maxnm = 0 # maximum number of models in given experiment. In principle they should all be the same, but some number of models may have failed in any given experiment.
    referenceMh0s = []
    for experiment in experimentList:
        nm = len(experiment.models)
        if nm>maxnm:
            maxnm = nm
            referenceMh0s = [model.p['Mh0'] for model in experiment.models]

    fidPts = np.zeros((2, maxnm)) # x&y for each model
    newPts = np.zeros((2, maxnm, len(experimentList))) # x&y for each model for each experiment
    diffs = np.zeros((maxnm, len(experimentList))) - 1.0 # one for each experiment and model; -1 if there was no corresponding model between experiments)
    jFid = getJindices(fidExper, referenceMh0s)
    for i, experiment in enumerate(experimentList):
        jIndices = getJindices(experiment, referenceMh0s)
        counterFid = 0
        counterExper = 0
        for j in range(maxnm):
            if j in jIndices and j in jFid:
                # if the experiment in question actually has a model at the appropriate Mh0, grab it.

                try:
                    fidPt = ( fidExper.models[counterFid].var[xvar].sensible(timeIndex=zind),  fidExper.models[counterFid].var[yvar].sensible(timeIndex=zind) )
                    fidPts[:,j] = np.array([fidPt[0], fidPt[1]])


                    newPt = ( experiment.models[counterExper].var[xvar].sensible(timeIndex=zind),  experiment.models[counterExper].var[yvar].sensible(timeIndex=zind) )
                    newPts[:,j,i] = np.array([newPt[0],newPt[1]])
                except:
                    pdb.set_trace()

                if experiment.models[0].var[xvar].log:
                    dx = np.log10(newPt[0]/fidPt[0])
                else:
                    dx = newPt[0] - fidPt[0]
                if experiment.models[0].var[yvar].log:
                    dy = np.log10(newPt[1]/fidPt[1])
                else:
                    dy = newPt[1] - fidPt[1]

                thisDiff = np.sqrt( dx*dx + dy*dy )
                diffs[j, i] = thisDiff

            if j in jIndices:
                counterExper+=1
            if j in jFid:
                counterFid+=1

    meanDiffs = np.zeros(len(experimentList)) ## the avg diff for each experiment
    for i in range(len(experimentList)):
        #if yvar=='sfZ':
        #    pdb.set_trace()
        counter=0
        for j in range(4): #### only include the first 6 galaxies in the average!
            if diffs[j,i]>-0.5: 
                # only include diffs where there's a model for both the fiducial and experimental case
                meanDiffs[i] += diffs[j,i]
                counter = counter+1
            else:
                print "Ignoring i,j = ",i,j
        if counter>0:
            meanDiffs[i] = meanDiffs[i]/float(counter)

    zipped = zip(meanDiffs, range(len(meanDiffs)))
    zsorted = sorted(zipped, key=lambda x: -x[0])
    highDiffInds = [zsorted[i][1] for i in range(3)]
    for i in range(len(highDiffInds)):
        jIndices = getJindices( experimentList[highDiffInds[i]], referenceMh0s )
        for j in range(maxnm):
            if j in jIndices and j in jFid and j%2==0:
                #axIn.annotate('', xy=(newPts[0,j, highDiffInds[i]],newPts[1,j,highDiffInds[i]]), xytext=(fidPts[0,j],fidPts[1,j]), arrowprops=dict(color=colorNames[highDiffInds[i]], width=1, frac=0.1, headwidth=3, alpha=0.6, shrink=0.3) )
                axIn.plot( [fidPts[0,j], newPts[0,j, highDiffInds[i]]], [fidPts[1,j], newPts[1,j,highDiffInds[i]]], lw=2, alpha=1.0, c=colorNames[highDiffInds[i]] )
            if j in jIndices and j in jFid and j==highDiffInds[i] or j==highDiffInds[i]+1:
                #axIn.text(newPts[0,j,highDiffInds[i]], newPts[1,j,highDiffInds[i]], featureNames[highDiffInds[i]], color=colorNames[highDiffInds[i]])
                bottomLeft = [axIn.get_xlim()[0], axIn.get_ylim()[0]]
                bottomRight= [axIn.get_xlim()[1], axIn.get_ylim()[0]]
                upperRight = [axIn.get_xlim()[1], axIn.get_ylim()[1]]
                dy = axIn.get_ylim()[1] / axIn.get_ylim()[0]
                dx = axIn.get_xlim()[1] / axIn.get_xlim()[0]
                axIn.text( bottomRight[0]/dx**0.35, bottomRight[1]*dy**(0.04+0.075*i), featureNames[highDiffInds[i]] + ' '+str( np.round(meanDiffs[highDiffInds[i]],2) ), color=colorNames[highDiffInds[i]])


#                if thisDiff>thresholdDiff:
#                    axIn.annotate('', xy=(newPt[0],newPt[1]), xytext=(fidPt[0],fidPt[1]), arrowprops=dict(color=colorNames[i], width=1, headwidth=4) )
#                    if counterFid==0:


    return None 




if __name__ =='__main__':
    #runner('rf100', 'rf100fid', fh=0.4)
    #runner('rf109', 'rf107', fh=0.533)
    runner('rf220_', 'rf220fid')

