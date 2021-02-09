from readoutput import *
import matplotlib.colors as pltcolors

def heatingExperiment( hopkins_experiments_no_heating, hopkins_experiments_with_heating, map_no_heating, map_with_heating):
    # the goal here is to look at the differences between these two sets of sims, with paritcular emphasis on 
    pass


def quickCheck(experiment_list):
    ''' Take a look at Mh vs M_*, and a variety of other easy diagnostics:
            M_* vs. ReSF, ReMst, Z_g, Z_SF, fg, fgh2, vPhi, sfr, Sigma1 '''
    def markAs(ax, fit, adj=0):
        ''' For the given axis, typically a panel in a many-panel plot, mark whether data plotted here was used in a fit (fit=1) or not (fit=0). In the latter case, mark this as a prediction '''
        if fit==1:
            [i.set_linewidth(3) for i in ax.spines.values()]
            ax.text(.45+adj,.89, 'FIT', transform=ax.transAxes, fontsize=9)
        else:
            ax.text(.45+adj,.9, 'PRED', transform=ax.transAxes, fontsize=9)
    #def markAsPred(ax, adj=0):
    #    ax.text(.45+adj,.9, 'PRED', transform=ax.transAxes, fontsize=9)
    def Theta(x):
        if x>0:
            return 1
        else:
            return 0

    bigcolors = ['k','r','b','orange', 'green', 'pink', 'purple', 'tan','lightblue']*20
    #colors = [None]*10
    ncolors = len(experiment_list)
    colors = [ (1.0-float(i)/float(ncolors),0.2,float(i+1)/float(ncolors)) for i in range(ncolors) ]
    colorby = 'Mh'

    nz = len(experiment_list[0].models[0].var['z'].sensible())
    skip6 = range(1,nz,6)
    #self.ptMovie(xvar='Mh', yvar=['mstar'],colorby='Mh0',prev=1,timeIndex=skip6, movie=True)
    #self.ptMovie(xvar='mstar', yvar=['halfMassStars','halfMassSFR','integratedZ','sfZ','Z1','fg','fgh2','fghm','vPhiOuter','sfr','Sigma1','gasToStellarRatio','gasToStellarRatioH2'], colorby='Mh0',prev=1,timeIndex=skip6, movie=True)
    #self.ptMovie(xvar='MHI', yvar=['broeilsHI'], colorby='Mh0',prev=1,timeIndex=skip6, movie=True)
    #self.ptMovie(xvar='Sigma1p5', yvar=['sSFR'],colorby='Mh0',prev=1,timeIndex=skip6, movie=True)

    # Same as the GIFs above, but try to make it more suitable for a paper
    zinds = [ Nearest( experiment_list[0].models[0].var['z'].sensible(), z)[0] for z in [0,1,2,3,4]]
    #self.ptMovie(xvar='Mh', yvar=['mstar'],colorby='Mh0',prev=0,timeIndex=zinds, movie=False)
    #self.ptMovie(xvar='mstar', yvar=['halfMassStars','halfMassSFR','integratedZ','sfZ','Z1','fg','fgh2','fghm','vPhiOuter','sfr','Sigma1','gasToStellarRatio','gasToStellarRatioH2','c82'], colorby='z',prev=0,timeIndex=zinds, movie=False)
    #self.ptMovie(xvar='MHI', yvar=['broeilsHI'], colorby='z',prev=0,timeIndex=zinds, movie=False)
    #self.ptMovie(xvar='Sigma1p5', yvar=['sSFR'],colorby='z',prev=0,timeIndex=zinds, movie=False)

    #colorby = 'Mh0'
    #colorby = 'experIndex'
   
    basename = ''
    for ex in experiment_list:
        basename+=ex.name+'_'
    print ("*******************************************************")
    print ("ANALYZING GALAXIES WITH THE FOLLOWING BREAKDOWN:")
    for ex in experiment_list:
        print (ex.name, len(ex.models), [np.log10(model.p['Mh0']) for model in ex.models])
    print ("*******************************************************")

    # Go through and compare some runs to the Hopkins 
    if len(experiment_list)==3:
        fig,ax = plt.subplots(3,3, figsize=(7,7))
        if len(experiment_list)==3:
            zHopkins = []
            for i,z in enumerate([0,0,2]):
                zHopkins.append( Nearest( experiment_list[i].models[0].var['z'].sensible(), z)[0] )
        else:
            print("Don't expect reasonable results in the comparison to Hopkins because something is wrong in the assumptions")
            zHopkins = [Nearest( experiment_list[0].models[0].var['z'].sensible(), z)[0] for z in [0,0,2] ]
        for i, ex in enumerate(experiment_list):
            experiment_list[i].ptMovie( xvar='mstar', yvar=['Mh'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[0,0], textsize=6, plotObs=True)
            ax[0,0].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [1.4e12, 1.6e12, 2.0e10], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['mbar'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[0,1], textsize=6, plotObs=True)
            ax[0,1].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [1.07e11, 7.13e10, 8.9e8], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['mCentral'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[0,2], textsize=6, plotObs=True)
            ax[0,2].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [7e9, 1.5e10, 1e7], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['mdisk'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[1,0], textsize=6, plotObs=True)
            ax[1,0].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [3e10,4.73e10, 1.3e8], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['naiveScaleLength'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[1,1], textsize=6, plotObs=True)
            ax[1,1].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [1.6,3.0,0.7], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['mgas'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[1,2], textsize=6, plotObs=True)
            ax[1,2].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [7e10, 0.9e10, 7.5e8], c='k', marker='*', s=200)
            
            experiment_list[i].ptMovie( xvar='mstar', yvar=['naiveGasScaleLength'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[2,0], textsize=6, plotObs=True)
            ax[2,0].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [3.2, 6.0, 2.1], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['Qavg'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[2,1], textsize=6, plotObs=True)
            ax[2,1].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [1,1,1], c='k', marker='*', s=200)

            experiment_list[i].ptMovie( xvar='mstar', yvar=['fghm'], colorby='Mh0', prev=0, timeIndex=[zHopkins[i]], movie=False, axIn=ax[2,2], textsize=6, plotObs=True)
            ax[2,2].scatter( [7e9+3e10, 1.5e10+4.73e10, 1e7+1.3e8], [.49, .09, .56], c='k', marker='*', s=200)
        plt.tight_layout()
        plt.savefig(basename+'callibrationHopkins.pdf')
        plt.close(fig)
        



    fig,ax = plt.subplots(1,1, figsize=(7,7))
    for i, ex in enumerate(experiment_list):
        for j in range(1,nz,2):
            experiment_list[i].ptMovie(xvar='mstar', yvar=['halfMassStars'], colorby='deltaBeta', prev=0, timeIndex=[j], movie=False, axIn=ax, textsize=6, plotObs=(i==0 and (j in zinds) or (j+1 in zinds)) )
    plt.savefig(basename+'calibration80.pdf')
    plt.close(fig)

    fig,ax = plt.subplots(1,1, figsize=(7,7))
    for i, ex in enumerate(experiment_list):
        for j in range(1,nz,2):
            experiment_list[i].ptMovie(xvar='mstar', yvar=['Sigma1'], colorby='deltaBeta', prev=0, timeIndex=[j], movie=False, axIn=ax, textsize=6, plotObs=(i==0 and (j in zinds) or (j+1 in zinds)) )
    plt.savefig(basename+'calibration81.pdf')
    plt.close(fig)

    fig,ax = plt.subplots(1,1, figsize=(7,7))
    for i, ex in enumerate(experiment_list):
        for j in range(1,nz,2):
            experiment_list[i].ptMovie(xvar='mstar', yvar=['sfZ'], colorby='deltaBeta', prev=0, timeIndex=[j], movie=False, axIn=ax, textsize=6, plotObs=(i==0 and (j in zinds) or (j+1 in zinds)) )
    plt.savefig(basename+'calibration82.pdf')
    plt.close(fig)

    fig,ax = plt.subplots(1,1, figsize=(7,7))
    for i, ex in enumerate(experiment_list):
        for j in range(1,nz,2):
            experiment_list[i].ptMovie(xvar='mstar', yvar=['BTExtrap'], colorby='deltaBeta', prev=0, timeIndex=[j], movie=False, axIn=ax, textsize=6, plotObs=(i==0 and (j in zinds) or (j+1 in zinds)) )
    plt.savefig(basename+'calibration83.pdf')
    plt.close(fig)

    # just try making the plots you want directly!
    fig,ax = plt.subplots(1,4, figsize=(8,3))
    fig.subplots_adjust(wspace=0.06, hspace=0.3, bottom=0.18)
    LOL = []
    for j in range(4):
        ax[j].set_title( r'$z=$'+str(j))
        for i,ex in enumerate(experiment_list):
            if j==0:
                LOL=[1,2]
            elif j==1:
                LOL=[0,3]
            else:
                LOL=[]
            experiment_list[i].ptMovie(xvar='Mh', yvar=['mstar'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[j], textsize=6, plotObs=(i==0), labelObsList=LOL)
            #self.ptMovie(xvar='Rvir', yvar=['halfMassStars'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6)
            #self.ptMovie(xvar='MHI', yvar=['broeilsHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6)
    for j in range(4):
        for i in range(1):
            if j>0:
                ax[j].set_ylabel('')
                ax[j].get_yaxis().set_ticks([])
        ax[j].get_xaxis().set_ticks([1.0e9, 1.0e11,1.0e13])
        ax[j].set_xlim(1.0e8, 2.0e13)
        ax[j].set_ylim(1.0e6, 2.0e11)
        markAs(ax[j], 1)
    plt.savefig(basename+'calibration0.pdf')
    plt.close(fig)




    fig,ax = plt.subplots(5,4, figsize=(8,8))
    fig.subplots_adjust(wspace=0.042, hspace=0.07, bottom=0.1)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))

        for i,ex in enumerate(experiment_list):
            # customize which labels to show on which panel. Gross.
            if j==0:
                LOL = [0,1]
            elif j==1:
                LOL = [2]
            elif j==2:
                LOL = [3,4]
            else:
                LOL = []
            experiment_list[i].ptMovie(xvar='mstar', yvar=['sSFR'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0), labelObsList=LOL )
            if j==0:
                LOL = [2,3]
            elif j==1:
                LOL = [0,1]
            else:
                LOL = []
            experiment_list[i].ptMovie(xvar='mstar', yvar=['sfZ'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0), labelObsList=LOL)
            if j==0:
                LOL=[0,1]
            else:
                LOL=[]
            experiment_list[i].ptMovie(xvar='mstar', yvar=['stZ'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, plotObs=(i==0), labelObsList=LOL)
            if j==0:
                LOL = [2]
            elif j==3:
                LOL = [0,1]
            else:
                LOL = []
            experiment_list[i].ptMovie(xvar='mstar', yvar=['gasToStellarRatioH2'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, plotObs=(i==0), labelObsList=LOL)
            experiment_list[i].ptMovie(xvar='mstar', yvar=['gasToStellarRatioHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[4,j], textsize=6, plotObs=(i==0), labelObsList=[0,1])

        
        markAs(ax[0,j], 1)
        markAs(ax[1,j], 1)
        markAs(ax[3,j], 1, (1-Theta(j))*-0.2)

        markAs(ax[2,j], 1-Theta(j), (1-Theta(j))*0.1)
        markAs(ax[4,j], 1-Theta(j), (1-Theta(j))*-0.25)
        #ax[2,j].text(1.0e10, 1.0e-1, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[4,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    plt.savefig(basename+'calibration1.pdf')
    plt.close(fig)
        







    fig,ax = plt.subplots(4,4, figsize=(8,7))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    LOL=[]
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i,ex in enumerate(experiment_list):
            if j==0:
                LOL=[0,1]
            elif j==3:
                LOL=[2,3]
            else:
                LOL=[]
            ex.ptMovie(xvar='mstar', yvar=['halfMassStars'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6,  plotObs=(i==0), labelObsList=LOL)
            if j==1:
                LOL=[1]
            else:
                LOL=[]
            ex.ptMovie(xvar='mstar', yvar=['vPhi22'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0), labelObsList=LOL)
            ex.ptMovie(xvar='mstar', yvar=['c82'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, plotObs=(i==0))
            if j==0:
                LOL=[2]
            elif j==3:
                LOL=[1,2]
            else:
                LOL=[]
            ex.ptMovie(xvar='mstar', yvar=['Sigma1'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, plotObs=(i==0), labelObsList=LOL)
        markAs(ax[0,j], 1)
        markAs(ax[1,j], 1-Theta(j-1), (1-Theta(j))*0.1)
        #markAs(ax[2,j], 1-Theta(j), (1-Theta(j))*-0.3)
        markAs(ax[2,j], 0)
        markAs(ax[3,j], 1)

    for j in range(4):
        for i in range(4):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<3:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[3,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    ax[1,1].get_yaxis().set_ticks([])
    ax[1,2].get_yaxis().set_ticks([])
    ax[1,3].get_yaxis().set_ticks([])
    ax[1,1].get_yaxis().set_visible(False)
    ax[1,2].get_yaxis().set_visible(False)
    ax[1,3].get_yaxis().set_visible(False)
    plt.savefig(basename+'calibration2.pdf')
    plt.close(fig)


    ### Not really callibration, but the same sort of plots
    fig,ax = plt.subplots(4,4, figsize=(8,7))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i,ex in enumerate(experiment_list):
            ex.ptMovie(xvar='mstar', yvar=['specificJStars'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, color=colors[i], plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['specificJH2'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, color=colors[i], plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['specificJHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, color=colors[i], plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['specificJAccr'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, color=colors[i], plotObs=(i==0))
        #ax[1,j].text(1.0e10, 80.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(4):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<3:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[3,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    plt.savefig(basename+'calibration3.pdf')
    plt.close(fig)


    ### Not really callibration, but the same sort of plots
    fig,ax = plt.subplots(4,4, figsize=(8,7))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.ptMovie(xvar='mstar', yvar=['metallicityGradient2kpc'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['metallicityGradientR90'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['metallicityGradient'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6,  plotObs=(i==0))
            ex.ptMovie(xvar='mstar', yvar=['massWeightedMetallicityGradient'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6,  plotObs=(i==0))
        #ax[1,j].text(1.0e10, 80.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(4):
            ax[i,j].axhline(0.0, c='gray', ls='--') # Mark zero.
            kmosMasses = np.power(10.0, np.linspace(9.88,11.13, 20))
            kmosMedian = kmosMasses*0 - 0.008
            if j==1 or j==2:
                ax[i,j].fill_between( kmosMasses, kmosMedian-0.04, kmosMedian+0.04, facecolor='r', alpha=0.3) 
                ax[i,j].plot( kmosMasses, kmosMedian, c='r', label='Wuyts16')
            if j==2:
                jonesMstars = np.power(10.0, np.array([9.4, 9.1, 9.9, 10.1]))
                J0744means = np.array([-0.06, -0.13, 0.02, 0.02, 0.10])
                J0744stds = np.array([0.04, 0.05, 0.04, 0.04, 0.07])
                J1038means = np.array([0.08, 0.15, 0.25, 0.31, 0.37])
                J1038stds = np.array([0.03, 0.07, 0.07, 0.08, 0.09])
                J1148means = np.array([-0.28, -0.51, -0.33, -0.37, -0.37])
                J1148stds = np.array([ 0.05, 0.11, 0.12, 0.11, 0.19])
                J1206means = np.array([-0.25, -0.45, -0.40, -0.46, -0.40])
                J1206stds = np.array([0.06, 0.08, 0.08, 0.08, 0.11])
                ax[i,j].errorbar([jonesMstars[0]], np.mean(J0744means), yerr=np.sqrt(np.var(J0744means)+np.mean(J0744stds*J0744stds)) , c='b', label='Jones13', fmt='o')
                ax[i,j].errorbar([jonesMstars[1]], np.mean(J1038means), yerr=np.sqrt(np.var(J1038means)+np.mean(J1038stds*J1038stds)) , c='b', fmt='o')
                ax[i,j].errorbar([jonesMstars[2]], np.mean(J1148means), yerr=np.sqrt(np.var(J1148means)+np.mean(J1148stds*J1148stds)) , c='b', fmt='o')
                ax[i,j].errorbar([jonesMstars[3]], np.mean(J1206means), yerr=np.sqrt(np.var(J1206means)+np.mean(J1206stds*J1206stds)) , c='b', fmt='o')
            if j==0:
                ax[i,j].fill_between( [10**8.4, 10**11.37], [-.041-0.03]*2,[-0.041+0.03]*2, facecolor='k', alpha=0.3)
                ax[i,j].plot( [10**8.4, 10**11.37], [-.041]*2, c='k', label='Rupke10')
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<3:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[3,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    plt.savefig(basename+'calibration4.pdf')
    plt.close(fig)

    ### Not really callibration, but the same sort of plots
    fig,ax = plt.subplots(5,4, figsize=(8,8))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.ptMovie(xvar='mstar', yvar=['MJeansAvg'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['BTcen'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['maxsig'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['mdotBulgeG'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['stellarHaloFraction'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[4,j], textsize=6, plotObs=(i==0))#, color=colors[i])
        #ax[1,j].text(1.0e10, 80.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[4,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    plt.savefig(basename+'calibration5.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(4,4, figsize=(8,6))
    fig.subplots_adjust(wspace=0.01, hspace=0.04, bottom=0.1)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.ptMovie(xvar='mstar', yvar=['fractionGI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['tdep'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['tDepH2'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='mstar', yvar=['maxsig'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, plotObs=(i==0))#, color=colors[i])
        #self.ptMovie(xvar='mstar', yvar=['v1kpc'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6)
        #self.ptMovie(xvar='mstar', yvar=['v1PerSFR'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[4,j], textsize=6)
        #self.ptMovie(xvar='mstar', yvar=['vOverSigGlobal'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[5,j], textsize=6)
        #ax[1,j].text(1.0e10, 80.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(4):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<3:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
        ax[2,j].get_xaxis().set_ticks([1.0e7, 1.0e9, 1.0e11])
    plt.savefig(basename+'calibration6.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(6,4, figsize=(8,10))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    #cax = fig.add_axes([.9, .1, .05, .12])
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i,ex in enumerate(experiment_list):
        #self.ptMovie(xvar='sfr', yvar=['LXProxy'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6)
            ex.ptMovie(xvar='MHI', yvar=['broeilsHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='sfr', yvar=['sfsig'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='sfr', yvar=['maxsig'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            #self.ptMovie(xvar='mbar', yvar=['vPhiOuter'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[2,j], textsize=6)
            ex.ptMovie(xvar='gbar', yvar=['gtot'], colorby='z', prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[3,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='stellarToGasMass', yvar=['sfZ'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[4,j], textsize=6)#, color=colors[i])
            ex.ptMovie(xvar='Mh', yvar=['integratedMLF'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[5,j], textsize=6, plotObs=(i==0))#, color=colors[i])
        #ax[1,j].text(1.0e7, 2.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
    plt.savefig(basename+'calibration7.pdf')
    plt.close(fig)





    percentiles=None
    fig,ax = plt.subplots(5,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colst'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[0,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colsfr'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[1,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sSFRRadial'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[2,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colH2'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[3,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colHI'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[4,j], color=colors[i])
        #ax[0,j].text(1.0e12, 1.0e7, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
    plt.savefig(basename+'calibration38.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(6,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['Z'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[0,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['fH2'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[1,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sig'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[2,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['hGas'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[3,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['ageRadial'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[4,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhi'],colorby='Mh0',percentiles=percentiles,logR=True,scaleR=False,movie=False, axIn=ax[5,j], color=colors[i])
        #ax[2,j].text(2, 70.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(6):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].set_yticks([])
            if i<5:
                ax[i,j].set_xlabel('')
                ax[i,j].set_xticks([])
    plt.savefig(basename+'calibration39.pdf')
    plt.close(fig)


    percentiles=None
    fig,ax = plt.subplots(3,4, figsize=(8,7))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['kappaZ'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[0,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['kappaZlimit'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[1,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['kappaZconservativeLimit'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[2,j]) #, color=colors[i])
        #ax[0,j].text(1.0e12, 1.0e7, r'$z=$'+str(j))
    for j in range(4):
        for i in range(3):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<2:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
    plt.savefig(basename+'calibration27kz.pdf')
    plt.close(fig)


    percentiles=None
    fig,ax = plt.subplots(5,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colst'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[0,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colsfr'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[1,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sSFRRadial'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[2,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colH2'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[3,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colHI'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[4,j]) #, color=colors[i])
        #ax[0,j].text(1.0e12, 1.0e7, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
    plt.savefig(basename+'calibration28.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(6,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['Z'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[0,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['fH2'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[1,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sig'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[2,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['hGas'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[3,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['ageRadial'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[4,j]) #, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhi'],colorby='Mh0',percentiles=percentiles,logR=False,scaleR=False,movie=False, axIn=ax[5,j]) #, color=colors[i])
        #ax[2,j].text(2, 70.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(6):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].set_yticks([])
            if i<5:
                ax[i,j].set_xlabel('')
                ax[i,j].set_xticks([])
    plt.savefig(basename+'calibration29.pdf')
    plt.close(fig)








    #percentiles=[16,50,84]
    percentiles=None
    fig,ax = plt.subplots(5,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colst'],colorby='Mh',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[0,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colsfr'],colorby='Mh',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[1,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sSFRRadial'],colorby='Mh',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[2,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colH2'],colorby='Mh',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[3,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['colHI'],colorby='Mh',percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[4,j])#, color=colors[i])
        #ax[0,j].text(1.0e12, 1.0e7, r'$z=$'+str(j))
    for j in range(4):
        for i in range(5):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
            if i<4:
                ax[i,j].set_xlabel('')
                ax[i,j].get_xaxis().set_ticks([])
    plt.savefig(basename+'calibration8.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(6,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    cbthis = 'Mh'
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['Z'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[0,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['fH2'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[1,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sig'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[2,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['hGas'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[3,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['ageRadial'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[4,j])#, color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhi'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[5,j])#, color=colors[i])
        #ax[2,j].text(2, 70.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(6):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].set_yticks([])
            if i<5:
                ax[i,j].set_xlabel('')
                ax[i,j].set_xticks([])
    plt.savefig(basename+'calibration9.pdf')
    plt.close(fig)



    fig,ax = plt.subplots(6,4, figsize=(8,9))
    fig.subplots_adjust(wspace=0.01, hspace=0.03)
    cbthis = 'Mh'
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i, ex in enumerate(experiment_list):
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhi'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[0,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhiGasRadial'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[1,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['vPhiStarsRadial'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[2,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sig'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[3,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sigstR'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[4,j], color=colors[i])
            ex.radialPlot(timeIndex=[zinds[j]],variables=['sigstZ'],colorby=cbthis,percentiles=percentiles,logR=False,scaleR=True,movie=False, axIn=ax[5,j], color=colors[i])
        #ax[2,j].text(2, 70.0, r'$z=$'+str(j))
    for j in range(4):
        for i in range(6):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].set_yticks([])
            if i<5:
                ax[i,j].set_xlabel('')
                ax[i,j].set_xticks([])
    plt.savefig(basename+'calibration9q.pdf')
    plt.close(fig)


    colormap_list = ['Reds', 'Blues', 'Oranges']*10
    #colormap_list = ['viridis']
    fig,ax = plt.subplots(1,1, figsize=(6,6))
    #self.ptMovie(xvar='gbar', yvar=['gtot'], colorby='rx', prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6)
    for i,ex in enumerate(experiment_list):
        xs = []
        ys = []
        weights = []
        for model in ex.models:
            xs = xs + list( np.log10(model.var['gbar'].sensible(timeIndex=zinds[0])) )
            ys = ys + list( np.log10(model.var['gtot'].sensible(timeIndex=zinds[0])) )
            weights = weights + list( model.var['dr'].sensible(timeIndex=zinds[0]) )
            #inner = model.var['rx'].sensible(timeIndex=zinds[0])<2.0 # within 1 half mass radius
        xs = np.array(xs)
        ys = np.array(ys)
        valid = np.isfinite(xs)
        xs[np.logical_not(valid)] = -13
        valid = np.isfinite(ys)
        ys[np.logical_not(valid)] = -13
        ax.hist2d( xs , ys , weights=weights, cmap=colormap_list[i], alpha=0.8, bins=60, norm=pltcolors.LogNorm())
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlim(-12, -8+.5)
    ax.set_ylim(-12, -8+.5)
    ax.set_xlabel(r'$g_\mathrm{bar} (\mathrm{m}/\mathrm{s}^2)$')
    ax.set_ylabel(r'$g_\mathrm{tot} (\mathrm{m}/\mathrm{s}^2)$')

    gbars = np.power(10.0, np.linspace(-12,-8, 100))
    gdagger = 1.20e-10
    gobss = gbars/(1.0 - np.exp(-np.sqrt(gbars/gdagger)))
    ax.plot( np.log10(gbars), np.log10(gbars), ls='--', c='gray', lw=3)
    ax.plot( np.log10(gbars), np.log10(gobss), label='McGaugh16', ls='-', c='b', lw=3 )
    ax.legend(loc=0)
    plt.savefig(basename+'calibration10.pdf')
    plt.close(fig)

    fig,ax = plt.subplots(1,1, figsize=(6,6))
    for i,ex in enumerate(experiment_list):
        ex.ptMovie(xvar='ageRadial', yvar=['colst'], colorby='mstar', prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6, color=colors[i])
    plt.savefig(basename+'calibration11.pdf')
    plt.close(fig)

    fig,ax = plt.subplots(1,1, figsize=(6,6))
    for i,ex in enumerate(experiment_list):
        ex.ptMovie(xvar='ageRadial', yvar=['sigstR'], colorby='mstar', prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6, color=colors[i])
    plt.savefig(basename+'calibration12.pdf')
    plt.close(fig)

    import observationalData as obsDat

    # Same as next, but now filter a bit.
    fig,ax = plt.subplots(1,1, figsize=(6,6))
    for i,ex in enumerate(experiment_list):
        xs = []
        ys = []
        weights = []
        modelCtr = 0
        for model in ex.models:
            if model.var['mstar'].sensible(timeIndex=zinds[0])>10**9.7 and model.var['mstar'].sensible(timeIndex=zinds[0])<10.0**11.4:
                modelCtr+=1
                the_filter =  model.var['rx'].sensible(timeIndex=zinds[0])<2.5
                xs = xs + list( np.log10(model.var['colst'].sensible(timeIndex=zinds[0]))[the_filter])
                ys = ys + list( np.log10(model.var['colsfr'].sensible(timeIndex=zinds[0]))[the_filter] )
                weights = weights + list( model.var['colst'].sensible(timeIndex=zinds[0])[the_filter] * model.var['dA'].sensible(timeIndex=zinds[0])[the_filter] )
                #weights = weights + list( model.var['colsfr'].sensible(timeIndex=zinds[0])[the_filter]   *   model.var['dA'].sensible(timeIndex=zinds[0])[the_filter] )
        print ("*****************************************************")
        print ("CALLIBRATION13: number of samples: ", len(xs), "over ", modelCtr, "galaxies")
        print ("*****************************************************")
        ax.hist2d( xs, ys, weights=weights, cmap=colormap_list[i], alpha=0.8, bins=40, range=[[0,3],[-4,-0.5]], norm=pltcolors.LogNorm() )

        #ex.ptMovie(xvar='colst', yvar=['colsfr'], colorby='mstar', prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6, color=colors[i])
    obsDat.datasets['CanoDiaz16_80p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_60p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_40p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_20p'].plot(0, axIn=ax)

    ax.set_xlim(0, 3)
    ax.set_ylim(-4,-0.5)
    ax.set_xlabel(r'$\log_{10} \Sigma_*\ (M_\odot\ \mathrm{pc}^{-2})$')
    ax.set_ylabel(r'$\log_{10} \Sigma_\mathrm{SFR}\ (M_\odot\ \mathrm{kpc}^{-2}\ \mathrm{yr}^{-1})$')
    plt.savefig(basename+'calibration13.pdf')
    plt.close(fig)




    fig,ax = plt.subplots(1,1, figsize=(6,6))
    for i,ex in enumerate(experiment_list):
        xs = []
        ys = []
        weights = []
        for model in ex.models:
            xs = xs + list( np.log10(model.var['colst'].sensible(timeIndex=zinds[0])) )
            ys = ys + list( np.log10(model.var['colsfr'].sensible(timeIndex=zinds[0])) )
            weights = weights + list( model.var['dA'].sensible(timeIndex=zinds[0]) * model.var['colst'].sensible(timeIndex=zinds[0]) )
        ax.hist2d( xs, ys, weights=weights, cmap=colormap_list[i], alpha=1.0, bins=50, range=[[0,3],[-4,-0.5]], norm=pltcolors.LogNorm() )

        #ex.ptMovie(xvar='colst', yvar=['colsfr'], colorby='mstar', prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6, color=colors[i])
    obsDat.datasets['CanoDiaz16_80p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_60p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_40p'].plot(0, axIn=ax)
    obsDat.datasets['CanoDiaz16_20p'].plot(0, axIn=ax)

    ax.set_xlim(0, 3)
    ax.set_ylim(-4,-0.5)
    ax.set_xlabel(r'$\log_{10} \Sigma_*\ (M_\odot\ \mathrm{pc}^{-2})$')
    ax.set_ylabel(r'$\log_{10} \Sigma_\mathrm{SFR}\ (M_\odot\ \mathrm{kpc}^{-2}\ \mathrm{yr}^{-1})$')
    plt.savefig(basename+'calibration17.pdf')
    plt.close(fig)




#    fig,ax = plt.subplots(1,2, figsize=(12,7))
#    for i,ex in enumerate(experiment_list):
#        ex.radialPlot(timeIndex=[zinds[0]], variables=['ageRadial'], colorby='Mh0', percentiles=None, logR=False, scaleR=True, light=True, movie=False, axIn=ax[0], color=colors[i])
#        ex.radialPlot(timeIndex=[zinds[0]], variables=['colst'], colorby='Mh0', percentiles=None, logR=False, scaleR=True, light=True, movie=False, axIn=ax[1], color=colors[i])
#    plt.savefig(basename+'calibration14.pdf')
#    plt.close(fig)


    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    plasma = plt.get_cmap('plasma')
    cNorm = colors.Normalize(vmin=10,vmax=12.5)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=plasma)

    fig,ax = plt.subplots(2,4, figsize=(12,7))
    fig.subplots_adjust(wspace=0.04, hspace=0.05)
    for j in range(4):
        for i,ex in enumerate(experiment_list):
            for model in ex.models:
                ax[0,j].set_title( r'$z=$'+str(j))
                cmod = np.max([np.min( [ (np.log10(model.p['accScaleLength']/0.13) + 0.1 ) * 5, .8 ]), 0])
                colorVal = scalarMap.to_rgba(np.log10(model.var['Mh'].sensible(timeIndex=zinds[j])))
                ax[0,j].plot( model.var['r'].cgs(timeIndex=zinds[j])/(0.015*model.var['Rvir'].cgs(timeIndex=zinds[j])), model.var['colstNormalizedKravtsov'].cgs(timeIndex=zinds[j]), color=colorVal, alpha=.4, lw=1 )
                kravtsovr = np.linspace(0,6, 200)
                kravtsovcolst = 0.6308*np.exp(-kravtsovr/(0.011/0.015))
                kravtsovcol = 0.099944*np.exp(-kravtsovr/(0.029/0.015))
                ls = '-'
                if j>0:
                    ls=':'
                ax[1,j].plot( model.var['r'].cgs(timeIndex=zinds[j])/(0.015*model.var['Rvir'].cgs(timeIndex=zinds[j])), model.var['colNormalizedKravtsov'].cgs(timeIndex=zinds[j]), color=colorVal, alpha=.4, lw=1 )
                if j==0:
                     ax[0,j].plot( kravtsovr, kravtsovcolst, lw=3, c='b', ls=ls)
                     ax[1,j].plot( kravtsovr, kravtsovcol, lw=3, c='b', ls=ls)
                else:
                     ax[0,j].plot( kravtsovr, kravtsovcolst, lw=2, c='gray', ls=ls)
                     ax[1,j].plot( kravtsovr, kravtsovcol, lw=2, c='gray', ls=ls)
        ax[0,j].set_xlim(0,5.97)
        ax[1,j].set_xlim(0,5.97)
        ax[0,j].set_ylim(1.0e-4, 10.0)
        ax[1,j].set_ylim(1.0e-4, 10.0)
        ax[0,j].get_xaxis().set_ticklabels([])
        ax[0,j].set_yscale('log')
        ax[1,j].set_yscale('log')
        ax[1,j].set_xlabel(r'$r/(0.015 R_\mathrm{vir})$')
        if j>0:
            ax[0,j].get_yaxis().set_ticklabels([])
            ax[1,j].get_yaxis().set_ticklabels([])
    ax[0,0].set_ylabel(experiment_list[0].models[0].var['colstNormalizedKravtsov'].texString)
    ax[1,0].set_ylabel(experiment_list[0].models[0].var['colNormalizedKravtsov'].texString)
    plt.savefig(basename+'calibration15.pdf')
    plt.close(fig)





    fig,ax = plt.subplots( figsize=(5,5))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    import behroozi,halo
    b = behroozi.behroozi()
    b.readSmmr()
    mstars = np.power(10.0, np.linspace(6,10.5, 100))
    etas = 3.6 * np.power(mstars/1.0e10, -0.35)
    ax.plot(mstars, etas, c='cyan', label='Muratov15')

    Mhs = np.power(10.0, np.linspace(10,13,100))
    mstars = []
    MhUsed = []
    for Mh in Mhs:
        try:
            ret,flag = b.getSmmr(np.log10(Mh), 0, eff=False)
            mstars.append( np.power(10.0, ret[1]) )
            MhUsed.append( Mh )
        except:
            pass
    
    cos = halo.Cosmology()
    etasZahid = []
    etasVogelsberger = []
    for Mh in MhUsed:
        thisHalo = halo.halo(Mh, 0, cos, 0)
        vcirc = thisHalo.vcirc
        alpha = .57 
        beta = .41
        etasZahid.append( alpha * np.power(100.0/vcirc, beta) )
        alpha = 50.0
        beta = 2.0
        etasVogelsberger.append( alpha * np.power(100.0/vcirc, beta) )
    ax.plot(mstars, etasZahid, c='k', label='Zahid14')
    ax.plot(mstars, etasVogelsberger, c='r', label='Vogelsberger13')
    ax.set_xlabel(r'$M_*\ (M_\odot)$')
    ax.set_ylabel(r'$\mu = \dot{M}_\mathrm{out}/\mathrm{SFR}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(basename+'calibration18b.pdf')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(5,5))
    for j in range(1):
        for i,ex in enumerate(experiment_list):
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mlf,_,_,_ = ex.constructQuantity('metalMassCGM',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( mstar, mlf, c=bigcolors[j], lw=0, s=10 )
    ax.errorbar( [10**10.1], [(0.23+0.46+0.50)*1.0e8], xerr=np.array([10**10.1-10**9.3, 10**10.8-10**10.1]).reshape((1,2)), yerr=np.array([np.sqrt((0.23-0.069)**2 + (0.46-0.28)**2 + (0.5-0.16)**2)*1.0e8,  np.sqrt((0.23-0.89)**2 + (0.46-1.1)**2 + (0.5-1.16)**2)*1.0e8]).reshape((1,2)), lw=2, c='k') # Werk14
    fac = 2.06*10400.0/4600.0
    ax.errorbar( [2.0*10**9], [5.0e5*fac], xerr=np.array([2.0e9-1.5e8, 3.0e9-2.0e9]).reshape((1,2)), yerr=fac*np.array([5.0e5-4.5e5,  1.0e7-5.0e5]).reshape((1,2)), lw=2, c='r') # Bordoloi14
    ax.errorbar( [ 6.0e9], [2.5e6*fac], xerr=np.array([6.0e9-3.0e9, 1.0e10-6.0e9]).reshape((1,2)), yerr=fac*np.array([2.5e6-2.4e6,2.0e7-2.5e6]).reshape((1,2)), lw=2, c='r') # Bordoloi14
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$M_* (M_\odot)$')
    ax.set_ylabel(r'$M_{Z,CGM}$')
    plt.savefig(basename+'calibration18f.pdf')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(5,5))
    for j in range(1):
        for i,ex in enumerate(experiment_list):
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mstarInt,_,_,_ = ex.constructQuantity('mstarIntegrated',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( mstar, mstarInt/mstar, c=bigcolors[j], lw=0, s=10 )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$M_* (M_\odot)$')
    ax.set_ylabel(r'$M_{*,\mathrm{integrated}}/M_{*}$')
    plt.savefig(basename+'calibration18i.pdf')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(5,5))
    for j in range(4):
        for i,ex in enumerate(experiment_list):
            sfr,_,_,_ = ex.constructQuantity('sfr', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mdotBulgeG,_,_,_ = ex.constructQuantity('mdotBulgeG',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( sfr, mdotBulgeG, c=bigcolors[j], lw=0, s=10 )
    ax.set_xlim(1.0e-5, 100)
    ax.set_ylim(1.0e-8, 100)
    ax.plot( [1.0e-5, 100], [1.0e-5, 100], lw=2, c='gray', ls='--' )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'SFR $(M_\odot/\mathrm{yr})$')
    ax.set_ylabel(r'$\dot{M}_{g,\mathrm{bulge}}$')
    plt.savefig(basename+'calibration18j.pdf')
    plt.close(fig)



    fig,ax=plt.subplots(figsize=(5,5))
    for j in range(1):
        for i,ex in enumerate(experiment_list):
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            HI,_,_,_ = ex.constructQuantity('MHI',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            HII,_,_,_ = ex.constructQuantity('MHII',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            H2,_,_,_ = ex.constructQuantity('MH2',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            Mh,_,_,_ = ex.constructQuantity('Mh',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mbar,_,_,_ = ex.constructQuantity('mbar',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            Mout,_,_,_ = ex.constructQuantity('mOut',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( Mh, mstar/Mh/0.17, c=bigcolors[0], lw=0, s=10, label=r'$M_*$', alpha=0.6 )
            ax.scatter( Mh, HI/Mh/0.17, c=bigcolors[1], lw=0, s=10, label=r'$M_\mathrm{HI}$', alpha=0.6 )
            ax.scatter( Mh, HII/Mh/0.17, c=bigcolors[2], lw=0, s=10, label=r'$M_\mathrm{HII}$', alpha=0.6 )
            ax.scatter( Mh, H2/Mh/0.17, c=bigcolors[3], lw=0, s=10, label=r'$M_{\mathrm{H}_2}$', alpha=0.6 )
            ax.scatter( Mh, mbar/Mh/0.17, c=bigcolors[4], lw=0, s=10, label=r'$M_{\mathrm{bar}}$', alpha=0.6 )
            ax.scatter( Mh, Mout/Mh/0.17, c=bigcolors[5], lw=0, s=15, marker='s', label=r'$M_{\mathrm{out}}$', alpha=0.6 )
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$M_h (M_\odot)$')
    ax.set_ylabel(r'$M/M_h/0.17$')
    ax.legend()
    plt.savefig(basename+'calibration18g.pdf')
    plt.close(fig)



    fig,ax = plt.subplots( figsize=(5,5))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    import behroozi,halo
    b = behroozi.behroozi()
    b.readSmmr()
    for j in range(4):
        for i,ex in enumerate( experiment_list ):
            #ex.ptMovie( xvar='mstar', yvar=['integratedMLF'], colorby='z', prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax, textsize=6, color=colors[i])
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mlf,_,_,_ = ex.constructQuantity('ZOutZDisk',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( mstar, mlf, c=bigcolors[j], lw=0, s=10 )
    mstars = np.power(10.0, np.linspace(6,10.5, 100))
    etas = 3.6 * np.power(mstars/1.0e10, -0.35)
    ax.plot(mstars, etas, c='cyan', label='Muratov15')

    Mhs = np.power(10.0, np.linspace(10,13,100))
    mstars = []
    MhUsed = []
    for Mh in Mhs:
        try:
            ret,flag = b.getSmmr(np.log10(Mh), 0, eff=False)
            mstars.append( np.power(10.0, ret[1]) )
            MhUsed.append( Mh )
        except:
            pass
    
    cos = halo.Cosmology()
    etasZahid = []
    etasVogelsberger = []
    for Mh in MhUsed:
        thisHalo = halo.halo(Mh, 0, cos, 0)
        vcirc = thisHalo.vcirc
        alpha = .57 
        beta = .41
        etasZahid.append( alpha * np.power(100.0/vcirc, beta) )
        alpha = 50.0
        beta = 2.0
        etasVogelsberger.append( alpha * np.power(100.0/vcirc, beta) )
    ax.plot(mstars, etasZahid, c='k', label='Zahid14')
    ax.plot(mstars, etasVogelsberger, c='r', label='Vogelsberger13')
    ax.set_xlabel(r'$M_*\ (M_\odot)$')
    ax.set_ylabel(r'$\mu_Z = Z_w \dot{M}_\mathrm{out}/ (Z_\mathrm{disk} \mathrm{SFR})$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(basename+'calibration18d.pdf')
    plt.close(fig)




    fig,ax = plt.subplots( figsize=(5,5))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    import behroozi,halo
    b = behroozi.behroozi()
    b.readSmmr()
    for j in range(4):
        for i,ex in enumerate( experiment_list ):
            #ex.ptMovie( xvar='mstar', yvar=['integratedMLF'], colorby='z', prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax, textsize=6, color=colors[i])
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mlf,_,_,_ = ex.constructQuantity('ZOutZSF',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( mstar, mlf, c=bigcolors[j], lw=0, s=10 )
    mstars = np.power(10.0, np.linspace(6,10.5, 100))
    etas = 3.6 * np.power(mstars/1.0e10, -0.35)
    ax.plot(mstars, etas, c='cyan', label='Muratov15')

    Mhs = np.power(10.0, np.linspace(10,13,100))
    mstars = []
    MhUsed = []
    for Mh in Mhs:
        try:
            ret,flag = b.getSmmr(np.log10(Mh), 0, eff=False)
            mstars.append( np.power(10.0, ret[1]) )
            MhUsed.append( Mh )
        except:
            pass
    
    cos = halo.Cosmology()
    etasZahid = []
    etasVogelsberger = []
    for Mh in MhUsed:
        thisHalo = halo.halo(Mh, 0, cos, 0)
        vcirc = thisHalo.vcirc
        alpha = .57 
        beta = .41
        etasZahid.append( alpha * np.power(100.0/vcirc, beta) )
        alpha = 50.0
        beta = 2.0
        etasVogelsberger.append( alpha * np.power(100.0/vcirc, beta) )
    ax.plot(mstars, etasZahid, c='k', label='Zahid14')
    ax.plot(mstars, etasVogelsberger, c='r', label='Vogelsberger13')
    ax.set_xlabel(r'$M_*\ (M_\odot)$')
    ax.set_ylabel(r'$\mu_Z = Z_w \dot{M}_\mathrm{out}/(y \mathrm{SFR})$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(basename+'calibration18c.pdf')
    plt.close(fig)




    fig,ax = plt.subplots( figsize=(5,5))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    import behroozi,halo
    b = behroozi.behroozi()
    b.readSmmr()
    for j in range(4):
        for i,ex in enumerate( experiment_list ):
            #ex.ptMovie( xvar='mstar', yvar=['integratedMLF'], colorby='z', prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax, textsize=6, color=colors[i])
            mstar,_,_,_ = ex.constructQuantity('mstar', timeIndex=[zinds[j]], locIndex=None, flatten=False)
            mlf,_,_,_ = ex.constructQuantity('integratedMLF',timeIndex=[zinds[j]], locIndex=None, flatten=False)
            ax.scatter( mstar, mlf, c=bigcolors[j], lw=0, s=10 )
    mstars = np.power(10.0, np.linspace(6,10.5, 100))
    etas = 3.6 * np.power(mstars/1.0e10, -0.35)
    ax.plot(mstars, etas, c='cyan', label='Muratov15')

    Mhs = np.power(10.0, np.linspace(10,13,100))
    mstars = []
    MhUsed = []
    for Mh in Mhs:
        try:
            ret,flag = b.getSmmr(np.log10(Mh), 0, eff=False)
            mstars.append( np.power(10.0, ret[1]) )
            MhUsed.append( Mh )
        except:
            pass
    
    cos = halo.Cosmology()
    etasZahid = []
    etasVogelsberger = []
    for Mh in MhUsed:
        thisHalo = halo.halo(Mh, 0, cos, 0)
        vcirc = thisHalo.vcirc
        alpha = .57 
        beta = .41
        etasZahid.append( alpha * np.power(100.0/vcirc, beta) )
        alpha = 50.0
        beta = 2.0
        etasVogelsberger.append( alpha * np.power(100.0/vcirc, beta) )
    ax.plot(mstars, etasZahid, c='k', label='Zahid14')
    ax.plot(mstars, etasVogelsberger, c='r', label='Vogelsberger13')
    ax.set_xlabel(r'$M_*\ (M_\odot)$')
    ax.set_ylabel(r'$\mu = \dot{M}_\mathrm{out}/\mathrm{SFR}$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    plt.savefig(basename+'calibration18.pdf')
    plt.close(fig)



    fig,ax = plt.subplots(2,4, figsize=(10,6))
    fig.subplots_adjust(wspace=0.01, hspace=0.45)
    #cax = fig.add_axes([.9, .1, .05, .12])
    for j in range(4):
        ax[0,j].set_title( r'$z=$'+str(j))
        for i,ex in enumerate(experiment_list):
        #self.ptMovie(xvar='sfr', yvar=['LXProxy'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6)
            ex.ptMovie(xvar='MHI', yvar=['broeilsHI'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[0,j], textsize=6, plotObs=(i==0))#, color=colors[i])
            ex.ptMovie(xvar='sfr', yvar=['sfsig'], colorby=colorby, prev=0, timeIndex=[zinds[j]], movie=False, axIn=ax[1,j], textsize=6, plotObs=(i==0))#, color=colors[i])
        #ax[1,j].text(1.0e7, 2.0, r'$z=$'+str(j))
        markAs(ax[0,j], 1-Theta(j))
        markAs(ax[1,j], 1)
    for j in range(4):
        for i in range(2):
            if j>0:
                ax[i,j].set_ylabel('')
                ax[i,j].get_yaxis().set_ticks([])
    plt.savefig(basename+'calibration16.pdf')
    plt.close(fig)


    fig,ax = plt.subplots(1,1, figsize=(8,8))
    ax.plot([2.0e-2,5.0],[2.0e-2,5.0], c='gray', ls='--', lw=3)
    for i,ex in enumerate(experiment_list):
        ex.ptMovie( xvar='stZ', yvar=['sfZ'], colorby=colorby, prev=0, timeIndex=[zinds[0]], movie=False, axIn=ax, textsize=6, plotObs=(i==0)) 
    gas_z_datasets,gas_z_colors = obsDat.identifyApplicableDatasets('mstar','sfZ',0) 
    star_z_datasets,star_z_colors = obsDat.identifyApplicableDatasets('mstar','stZ',0) 
    mstars = np.power(10, np.linspace(8, 11, 100))
    counter=0
    for k, gds in enumerate(gas_z_datasets):
        # only plot datasets valid at z=0
        if obsDat.datasets[gds].zmin <=0 and obsDat.datasets[gds].zmax>=0:
            y16, y50, y84, logy = obsDat.datasets[gds].returnCachedQuantiles( mstars ) 
            for l, sds in enumerate(star_z_datasets):
                if obsDat.datasets[sds].zmin <=0 and obsDat.datasets[sds].zmax>=0:
                    x16, x50, x84, logx = obsDat.datasets[sds].returnCachedQuantiles( mstars ) 
                    # no extrapolations! Only plot cases where the two datasets have actual overlap in the stellar masses probed.
                    filt = np.logical_and( mstars > np.min(obsDat.datasets[sds].xval), mstars < np.max(obsDat.datasets[sds].xval))
                    filt = np.logical_and( filt, mstars > np.min(obsDat.datasets[gds].xval) )
                    filt = np.logical_and( filt, mstars < np.max(obsDat.datasets[gds].xval) )
                    if np.any(filt):
                        ax.plot(np.power(10,x50[filt]), np.power(10,y50[filt]), c=bigcolors[counter], lw=2, label=obsDat.datasets[sds].label+'+'+obsDat.datasets[gds].label)
                        ax.plot(np.power(10,x84[filt]), np.power(10,y50[filt]), c=bigcolors[counter], lw=1)
                        ax.plot(np.power(10,x50[filt]), np.power(10,y84[filt]), c=bigcolors[counter], lw=1)
                        ax.plot(np.power(10,x16[filt]), np.power(10,y50[filt]), c=bigcolors[counter], lw=1)
                        ax.plot(np.power(10,x50[filt]), np.power(10,y16[filt]), c=bigcolors[counter], lw=1)
                        counter+=1
    ax.set_xlim(2.0e-2,5.0)
    ax.set_ylim(2.0e-2,5.0)
    ax.legend()
    plt.savefig(basename+'calibration19.pdf')
    plt.close(fig)

    ## To add: age-velocity dispersion correlation?


if __name__=='__main__':
    #experiment_list = []
    #for k in range(10):
    #    experiment_list.append( Experiment('rf150_'+str(k)) )
    #experiment_list = [Experiment('rf153'), Experiment('rf154'), Experiment('rf151'), Experiment('rf152') ]
    #experiment_list = [Experiment('rf161'), Experiment('rf162'), Experiment('rf163') ]
#    experiment_list = [Experiment('rf212'), Experiment('rf213')]
    #experiment_list = [ Experiment('rf285')]
    #experiment_list = [ Experiment('rf313')]
    #experiment_list = [ Experiment('rf290'), Experiment('rf291'), Experiment('rf292'), Experiment('rf293') ]
    #experiment_list = [Experiment('rf118'), Experiment('rf119'), Experiment('rf120')]
    #experiment_list = [Experiment('hkm21'), Experiment('hkm22'), Experiment('hkm23')] # hopkins no heating
    #experiment_list = [Experiment('fe22'), Experiment('fe23'), Experiment('fe24')] # hopkins with heating
    #experiment_list = [Experiment('fe20')] # MAP no heating
    experiment_list = [Experiment('fe21')] # MAP with heating
    MONDargs = ['gbar', 'gtot', 'hGas', 'sSFRRadial', 'rxl', 'colstNormalizedKravtsov', 'colNormalizedKravtsov', 'colHI', 'colH2', 'colst', 'fH2', 'vPhi', 'sigstR', 'sigstZ', 'ageRadial', 'colsfr', 'Z', 'sig', 'col', 'vPhiGasRadial', 'vPhiStarsRadial', 'kappaZ', 'kappaZlimit', 'kappaZconservativeLimit']
    for ex in experiment_list:
        ex.read(MONDargs, keepStars=True)
    quickCheck(experiment_list)




