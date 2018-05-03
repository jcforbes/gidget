from readoutput import *
import matplotlib.colors as pltcolors


tW = np.array([ 1.1545454471140513, 1.4562262504912922, 1.6174725466389603, 1.7283910560022944, 1.7965734639435877, 1.8777958523969331, 1.9626902734826572, 2.040115430494593, 2.108906312357657, 2.180016761710267, 2.253524992309725, 2.316671738324658, 2.394787853464255, 2.4892586052022443, 2.5590109183314973, 2.6452984955626317, 2.7648914130211915, 2.8898910798240056, 3.003892903558393, 3.2098854713099376, 3.5066815356965195, 3.830920293793729, 4.374347699826242, 5.337385567099643, 6.128288264102687, 7.075387313922233, 8.034520783900259, 9.225089476610762, 10.592078624278916, 12.571709892264574 ])
sW = np.array([ 10.15947094711839, 10.119366364989414, 11.807187869366663, 9.80417679462689, 11.394270225407617, 12.678456671422186, 11.039371014586864, 12.043018110237995, 10.909152427547689, 12.043018110237995, 13.560276180165056, 15.884716272874831, 13.400321414471529, 15.884716272874831, 13.614017609237196, 14.793111980906472, 12.779149095990569, 13.189979561162446, 15.029019138748827, 15.512179549310623, 14.618614888860135, 16.72285843631423, 17.12447052897269, 20.379765023842037, 20.86920049387745, 24.83635342658302, 24.44650273214034, 27.526432797204478, 30.26749576398466, 30.99439256572452])
sWp = np.array([10.952387051506143, 10.823194512848927,12.678456671422186,10.48608269127704,12.235069260658351,13.667972023735317,11.900960668746963,12.931689153515185,11.714153946296314,12.931689153515185,14.676550674587839,17.12447052897269    ,14.389149754895975,17.056871562617136,14.676550674587839,15.947669810219324,13.776523188728843,14.219417734176508,16.201988811116035,16.72285843631423,15.759553745770539,18.028027163560893,18.460983870288672,21.970344294667754,22.49797872882931,    26.774755997810324,26.458925422884747,29.674787960072617,32.75909940704369,33.545833954643506])

tU = np.array([ 1.1545454471140513, 1.4481996267037256, 1.6174725466389603, 1.7283910560022944, 1.8065309442770299, 1.8777958523969331, 1.9626902734826572, 2.040115430494593, 2.108906312357657, 2.1920994705262253, 2.253524992309725, 2.316671738324658, 2.394787853464255, 2.5030552821005823, 2.5449058176283663, 2.6452984955626317, 2.7496515125055185, 2.8898910798240056, 3.003892903558393, 3.192192792668354, 3.4873529365723437, 3.830920293793729, 4.350236581591234, 5.337385567099643, 6.128288264102687, 7.036388242108341, 8.079051989720647, 9.17424138716694, 10.650785054440707, 12.571709892264574 ])
sU = np.array([ 23.03828470297126,26.25044414355292,24.83635342658302,21.118308554196396,26.04360557809312,26.35447863100062,24.640657091469063,26.774755997810324,27.745047990670184,21.455084220864002,30.38745043564325,31.738746343902054,28.523964850804024,23.12958891586919,31.488663079506047,32.24488713788331,33.41341154210631,33.94625838061393,25.94079826418102 ,32.24488713788331,30.38745043564325,36.740695743681265,34.62428217248994,37.326603098484604,39.45184274341137,41.533482242701545,45.13064926230107,48.46090259220439,53.497788341709885,48.26960284262397] )
sUp = np.array([ 24.83635342658302, 28.299212235621507, 26.774755997810324, 22.676657891822398, 27.965399428100582, 28.411366301434906, 26.458925422884747, 28.864445050736656, 29.910465410861065, 22.947340914087935, 32.75909940704369, 34.080792545905844, 30.50788050584042, 24.83635342658302, 33.94625838061393, 34.761503448799075, 35.87903361647077, 36.4511999637219, 27.855005819617247, 34.761503448799075, 32.75909940704369, 39.45184274341137, 37.326603098484604, 40.239832032642134, 42.53094023006707, 44.775045421823734, 48.65296049167699, 52.243130062368444, 57.673129572180216, 51.831483919616105 ])

tV = np.array([ 1.154781984689458, 1.4489643232788298, 1.6096474622839103, 1.7201919979837987, 1.8080545883896824, 1.879483258284322, 1.9645777937496076, 2.042190040422458, 2.0994976101348772, 2.1824399740999354, 2.2561365673928373, 2.3194478660758224, 2.4110795441388793, 2.5063312063142997, 2.548296747979347, 2.648969287610528, 2.753618977957703, 2.8942661247167507, 3.0253054571024536, 3.1974789370674914, 3.4935682400446213, 3.8382619912089933, 4.3593634327528035, 5.350142932297119,  6.144146161448629, 7.095149298539764, 8.013941241504929, 9.254354972049763, 10.627773930104956, 12.687188271180933 ])
sV = np.array([ 16.855671538516322,14.618614888860135,14.332348486294974,14.618614888860135,15.268688330380503,15.759553745770546,15.884716272874831,16.330665161755658,13.996218708819354,17.192337400274095,16.656844836483952,14.851739325048687,15.697342797897065,16.460363461207066,19.667049629525607,20.299315759857507,17.60522438531374,20.139368700347838,19.823245526441635,19.901808051959897,18.099474966706556,20.70476306688229,21.118308554196396,24.158135846053955,23.96778347962331, 25.736399532087084,27.309540167088848,28.864445050736656,28.97883921983696,38.223029701586235] )
sVp = np.array([ 18.099474966706556,15.822011245056265,15.389952399644592,15.697342797897065,16.395386062114245,16.92247311869177,17.056871562617147,17.53572763977108,15.08858141870079,18.460983870288672,17.956861401114423,15.947669810219333, 16.92247311869177,17.60522438531374,21.20200361071367,21.797230480984602,18.979259116958353,21.71118576226263,21.455084220864002,21.455084220864002,19.512084467414365,22.320707455986533,22.676657891822398,26.04360557809312,25.838396782854662,27.745047990670184,29.557646528486323,31.11722803954335,31.364361363287067,41.20622203183994 ])



def plotHolmberg(ax, sigU=False, sigV=False, sigW=False, ratioU=False, ratioV=False, alpha=1):
    assert sigU or sigV or sigW or ratioU or ratioV
    if sigU:
        ax.errorbar( tU*1.0e9, sU, yerr=sUp-sU, label=r'MW $\sigma_U$ (Holmberg09)', c='k', fmt='o' , alpha=alpha)
    if sigV:
        ax.errorbar( tV*1.0e9, sV, yerr=sVp-sV, label=r'MW $\sigma_V$ (Holmberg09)', c='k', fmt='o', alpha=alpha )
    if sigW:
        ax.errorbar( tW*1.0e9, sW, yerr=sWp-sW, label=r'MW $\sigma_W$ (Holmberg09)', c='k', fmt='o', alpha=alpha )
    if ratioU:
        su = sUp-sU
        sv = sVp-sV
        sw = sWp-sW
        u = sU
        v = sV
        w = sW
        yE = np.sqrt( su*su / (w*w) + sw*sw*u*u/(w*w*w*w) )
        #ax.errorbar( tU*1.0e9, np.sqrt(2.0)*sW/np.sqrt(sU*sU + sV*sV), yerr=yE, label=r'Holmberg09', c='b', ls='--' )
        ax.errorbar( tU*1.0e9, sW/sU, yerr=yE, label=r'Holmberg09', c='k', fmt='o', alpha=alpha )
    if ratioV:
        #su = sUp-sU
        sv = sVp-sV
        sw = sWp-sW
        #u = sU
        v = sV
        w = sW
         # predict u from u = w/(2\Omega/\kappa)
        u = v * np.sqrt(2.0)
        su = sv * np.sqrt(2.0)
        yE = np.sqrt( su*su / (w*w) + sw*sw*u*u/(w*w*w*w) )
        ax.errorbar( tU*1.0e9, sW/sU, yerr=yE, label=r'Holmberg09', c='r', fmt='o', alpha=alpha )

def collectStars(ex, rEff=-1, rFixed=-1):
    # need one of rEff or rFixed to be >0, but not both.
    assert (rEff>0 or rFixed>0) and not (rEff>0 and rFixed>0)

    vars = sorted(ex.models[0].var.keys()) # Get everything
    cols = []
    sigRs = []
    sigZs = []
    Zs = []
    VZs = []
    ages =[]
    startingAges=[]
    endingAges=[]
    N=100 # meh.
    for i in range(N):
        sti = str(i).zfill(2)
        if 'colst'+sti in vars:
            cols.append('colst'+sti)
        if 'sigstR'+sti in vars:
            sigRs.append('sigstR'+sti)
        if 'sigstZ'+sti in vars:
            sigZs.append('sigstZ'+sti)
        if 'Zst'+sti in vars:
            Zs.append('Zst'+sti)
        if 'VarZst'+sti in vars:
            VZs.append('VarZst'+sti)
        if 'ageSt'+sti in vars:
            ages.append('ageSt'+sti)
        if 'stargingAgeSt'+sti in vars:
            startingAges.append('startingAgeSt'+sti)
        if 'endingAgeSt'+sti in vars:
            endingAges.append('endingAgeSt'+sti)
    # Now we'll make some plots.


    allSigRs = []
    allSigZs = []
    allAges = []
    allZs = []
    allZVars = []
    for j,model in enumerate(ex.models):
        r= model.var['r'].sensible(timeIndex=-1)
        theSigRs = []
        theSigZs = []
        theAges = []
        theZs = []
        theZVars = []
        #for j,model in enumerate(self.models):
        rSample = rFixed # kpc
        if rEff>0:
            rSample = rEff*model.var['halfMassStars'].sensible(-1) 
        assert rSample > 0
        for i in range(len(sigRs)):
            theSigRs.append( model.var[ sigRs[i] ].atR([rSample], r, -1)  )
            theSigZs.append( model.var[ sigZs[i] ].atR([rSample], r, -1)  )
            theAges.append( model.var[ ages[i] ].sensible(timeIndex=-1) )
            theZs.append( model.var[ Zs[i] ].atR([rSample], r, -1) )
            #theZVars.append( model.var[ VZs[i] ].atR([rSample], r, -1) )
        theSigRs = np.array(theSigRs)
        theSigZs = np.array(theSigZs)
        theAges = np.array(theAges)
        theZs = np.array(theZs)
        #theZVars = np.array(theZVars)
        #ax.errorbar( [np.mean(theAges)], [np.mean(theSigRs)], xerr=[np.std(theAges)], yerr=[np.std(theSigRs)], ecolor='blue')
        #ax.errorbar( [np.mean(theAges)], [np.mean(theSigZs)], xerr=[np.std(theAges)], yerr=[np.std(theSigZs)], ecolor='red')
        #ax.plot( theAges, theSigRs, c='blue', lw=1, ls='-', alpha=0.5)
        #ax.plot( theAges, theSigZs, c='red', lw=1, ls='-', alpha=0.5)
        allSigRs.append(theSigRs)
        allSigZs.append(theSigZs)
        allAges.append(theAges)
        allZs.append(theZs)
        #allZVars.append(theZVars)
    allSigRs = np.array(allSigRs)
    allSigZs = np.array(allSigZs)
    allAges = np.array(allAges)
    allZs = np.array(allZs)

    return allSigRs, allSigZs, allAges, allZs

def plot(experlist, colors, linestyles, axIn, sigR=False, sigZ=False, ratio=False, vsAge=True, jmax=100, rEff=-1, rFixed=-1, labels=None ):
    assert len(experlist)==len(colors) and len(experlist)==len(linestyles)
    assert sigZ or sigR or ratio # otherwise we're not plotting anything..
    for i,exn in enumerate(experlist):
        try:
            ex = Experiment(exn)
            ex.read(keepStars=True)
            allSigRs, allSigZs, allAges, allZs = collectStars(ex, rEff=rEff, rFixed=rFixed)
            
            x = allAges
            y = allSigRs
            if not vsAge:
                x = allZs
            if sigZ:
                y = allSigZs
            if ratio:
                y = allSigZs/allSigRs
            label = None
            if labels is not None:
                label = labels[i]
            for j in range(len(ex.models)):
                if j==jmax or (j<jmax and jmax==100):
                    axIn.plot( x[j,:],y[j,:], c=colors[i], ls=linestyles[i], alpha=0.6, lw=3, label=label)
        except:
            print "Failed for experiment ",exn

def plotSigR_nicer():
    rEff = 0.5/1.7
    fig,ax = plt.subplots()
    plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, rEff=rEff )
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_sigR_age_mass_history_reff.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, ratio=True, rEff=rEff )
    ax.axhline( 0.8, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.81, r'Lacey84 Eq.')
    ax.axhline( 0.5, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.45, r'Spiral Heating')
    plotHolmberg(ax, ratioU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$ ')
    plt.savefig('rf18x_ratio_age_mass_history_reff.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, rFixed = 8.0 )
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_sigR_age_mass_history_8kpc.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, vsAge=False, rEff=rEff )
    ax.set_xlabel(r'$Z/Z_\odot$ ')
    ax.set_xscale('log')
    ax.set_xlim(1.0e-1, 3.0)
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_sigR_Z_mass_history_reff.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, rFixed=8.0, vsAge=False )
    ax.set_xscale('log')
    ax.set_xlim(1.0e-1, 3.0)
    ax.set_xlabel(r'$Z/Z_\odot$ ')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_sigR_Z_mass_history_8kpc.pdf')
    plt.close(fig)

    fig, ax = plt.subplots()
    plot( ['rf194', 'rf195'], ['b','r'], ['-','-'], ax, sigR=True, rEff=0.5/1.7, vsAge=True)
    ax.set_xlabel(r'Age')
    ax.set_ylabel(r'$\sigma_r$ (km/s)' )
    plt.savefig('rf19x_sigR_age_mass_reff.pdf')
    plt.close(fig)


    exlistHi = ['rf185', 'rf181', 'rf183'] + ['rf187', 'rf189'] + ['rf191', 'rf193']
    exlistLo = ['rf184', 'rf180', 'rf182'] + ['rf186', 'rf188'] + ['rf190', 'rf192']
    colors = ['b', 'k', 'r'] + ['lightblue','green'] + ['purple', 'orange']
    linestyles = ['-']*3 + ['--',':'] + ['-.', '-']
    labels = [ r'$f=0.1$', r'$f=1.0$', r'$f=10$', r'no spiral or GI', r'no SN turb', r'no gas GI', r'$f=0$' ]
    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    plot( exlistHi, colors, linestyles, ax, sigR=True, jmax=1, labels=labels, rEff=rEff )
    #plot( exlistLo, colors, linestyles, ax, sigR=True, jmax=0 )
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    ax.legend()
    plt.savefig('rf18x_sigR_age_physics_reff_1e12.pdf')
    plt.close(fig)

    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    #plot( exlistHi, colors, linestyles, ax, sigR=True, jmax=0 )
    plot( exlistLo, colors, linestyles, ax, sigR=True, jmax=1, labels=labels, rEff=rEff )
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    ax.legend()
    plt.savefig('rf18x_sigR_age_physics_reff_1e11.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    plot( exlistHi, colors, linestyles, ax, ratio=True, jmax=1, labels=labels, rEff=rEff )
    #plot( exlistLo, colors, linestyles, ax, sigR=True, jmax=0 )
    ax.axhline( 0.8, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.81, r'Lacey84 Eq.')
    ax.axhline( 0.5, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.45, r'Spiral Heating')
    plotHolmberg(ax, ratioU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    ax.legend()
    plt.savefig('rf18x_ratio_age_physics_reff_1e12.pdf')
    plt.close(fig)

    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    plot( exlistHi, colors, linestyles, ax, jmax=1, labels=labels, vsAge=False, sigR=True, rEff=rEff )
    #plot( exlistLo, colors, linestyles, ax, sigR=True, jmax=0 )
    ax.set_xscale('log')
    ax.set_xlim(1.0e-1, 3.0)
    ax.set_xlabel(r'$Z/Z_\odot$')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    ax.legend()
    plt.savefig('rf18x_sigR_Z_physics_reff_1e12.pdf')
    plt.close(fig)

    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    plot( exlistHi, colors, linestyles, ax, ratio=True, jmax=1, labels=labels, vsAge=False, rEff=rEff )
    #plot( exlistLo, colors, linestyles, ax, sigR=True, jmax=0 )
    ax.axhline( 0.8, ls='--', lw=2, c='k' )
    ax.text(1, 0.81, r'Lacey84 Eq.')
    ax.axhline( 0.5, ls='--', lw=2, c='k' )
    ax.text(1, 0.45, r'Spiral Heating')
    plotHolmberg(ax, ratioU=True)
    ax.set_xscale('log')
    ax.set_xlim(1.0e-1, 3.0)
    ax.set_xlabel(r'$Z/Z_\odot$')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    ax.legend()
    plt.savefig('rf18x_ratio_Z_physics_reff_1e12.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    #plot( ['rf180', 'rf181'], [ 'b', 'r'], ['-','-'], ax, sigR=True, jmax=0 )
    plot( exlistLo, colors, linestyles, ax, ratio=True, jmax=1, labels=labels, rEff=rEff )
    ax.axhline( 0.8, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.81, r'Lacey84 Eq.')
    ax.axhline( 0.5, ls='--', lw=2, c='k' )
    ax.text(1e10, 0.45, r'Spiral Heating')
    plotHolmberg(ax, ratioU=True)
    ax.set_xlabel(r'Age (yr)')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    ax.legend()
    plt.savefig('rf18x_ratio_age_physics_reff_1e11.pdf')
    plt.close(fig)


def plotSigR_all():
    exlistHi = ['rf185', 'rf181', 'rf183'] + ['rf187', 'rf189'] + ['rf191', 'rf193']
    exlistLo = ['rf184', 'rf180', 'rf182'] + ['rf186', 'rf188'] + ['rf190', 'rf192']
    colors = ['b', 'k', 'r'] + ['lightblue','green'] + ['purple', 'orange']
    linestyles = ['-']*3 + ['--',':'] + ['-.', '-']

    fig,ax = plt.subplots()
    plot( exlistHi, colors, linestyles, ax, sigR=True )
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_1e12_sigR_age.pdf')
    plt.close()

    fig,ax = plt.subplots()
    plot( exlistLo, colors, linestyles, ax, sigR=True )
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$ (km/s)')
    plt.savefig('rf18x_1e11_sigR_age.pdf')
    plt.close()


    fig,ax = plt.subplots()
    plot( exlistHi, colors, linestyles, ax, ratio=True )
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    plotHolmberg(ax,ratioU=True)
    plotHolmberg(ax,ratioV=True)
    plt.savefig('rf18x_1e12_ratio_age.pdf')
    plt.close()

    fig,ax = plt.subplots()
    plot( exlistLo, colors, linestyles, ax, ratio=True )
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    plt.savefig('rf18x_1e11_ratio_age.pdf')
    plt.close()

    fig,ax = plt.subplots()
    plot( exlistHi, colors, linestyles, ax, sigR=True, vsAge=False )
    ax.set_xlabel(r'Z (Gyr)')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf18x_1e12_sigR_Z.pdf')
    plt.close()

    fig,ax = plt.subplots()
    plot( exlistLo, colors, linestyles, ax, sigR=True, vsAge=False )
    ax.set_xlabel(r'Z (Gyr)')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf18x_1e11_sigR_Z.pdf')
    plt.close()


    fig,ax = plt.subplots()
    plot( exlistHi, colors, linestyles, ax, ratio=True, vsAge=False )
    ax.set_xlabel(r'Z (Gyr)')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    plt.savefig('rf18x_1e12_ratio_Z.pdf')
    plt.close()

    fig,ax = plt.subplots()
    plot( exlistLo, colors, linestyles, ax, ratio=True, vsAge=False )
    ax.set_xlabel(r'Z (Gyr)')
    ax.set_xscale('log')
    ax.set_ylabel(r'$\sigma_z/\sigma_r$')
    plt.savefig('rf18x_1e11_ratio_Z.pdf')
    plt.close()

def plotMassScaling(expers, colors):
    ## M* vs central velocity dispersion
    maxsig0=0
    fig,ax = plt.subplots()
    for j, exn in enumerate(expers):
        ex = Experiment(exn)
        ex.read(keepStars=True, keepOnly=['sigstZ', 'sigstR'])
        for m in ex.models:
            #sig0 = m.var['sigstZ'].sensible(timeIndex=-1, locIndex=1)
            hR = m.var['halfMassStars'].sensible(timeIndex=-1)/1.7
            sig0 = m.var['sigstZ'].atR( [hR/8.0], m.var['r'].sensible(timeIndex=-1), -1 )
            # keep track of the maximum value, since it will be used to color lines in the /next/ plot.
            if sig0>maxsig0:
                maxsig0=sig0
            ax.scatter( m.var['mstar'].sensible(-1), sig0, marker='o', lw=0, c=colors[j])
    ax.set_xscale('log')
    ax.set_xlabel(r'$M_* (M_\odot)$')
    ax.set_ylabel(r'Central $\sigma_{z,*}$ (km/s)')
    ax.set_yscale('log')
    plt.savefig('rf19x_massScaling.pdf')
    plt.close(fig)


    fig,ax = plt.subplots()
    vmin = 100000
    vmax = 0
    print expers
    for j, exn in enumerate(expers):
        ex = Experiment(exn)
        ex.read(keepStars=True, keepOnly=['sigstZ', 'sigstR'])
        print "Number of models: ", len(ex.models)
        for i, m in enumerate(ex.models):
            hR = m.var['halfMassStars'].sensible(timeIndex=-1)/1.7
            sig0 = m.var['sigstZ'].atR( [hR/8.0], m.var['r'].sensible(timeIndex=-1), -1 )
            #sig0 = m.var['sigstZ'].sensible(timeIndex=-1, locIndex=1)
            thisv = m.var['vPhi22'].sensible(-1)
            print thisv, sig0, i,j, exn
            # keep track of the maximum value, since it will be used to color lines in the /next/ plot.
            ax.scatter( thisv, sig0, marker='o', lw=0, c=colors[j])
            if thisv>vmax:
                vmax=thisv
            if thisv<vmin:
                vmin=thisv
    vs = np.linspace( vmin, vmax, 100 )
    ax.plot( vs, 0.25 * vs, lw=1, c='k' )
    ax.fill_between( vs, (0.248-0.038)*vs, (0.248+0.038)*vs, facecolor='gray', alpha=0.3 )
    ax.set_xlabel(r'$v_{2.2}$ (km/s)')
    ax.set_ylabel(r'Central $\sigma_{z,*}$ (km/s)')
    plt.savefig('rf19x_vScaling.pdf')
    plt.close(fig)


    #linestyles=['--']*3+['-']
    fig,ax = plt.subplots()
    for j, exn in enumerate(expers):
        ex = Experiment(exn)
        ex.read(keepStars=True, keepOnly=['sigstZ', 'sigstR'])
        for m in ex.models:
            mstar = m.var['mstar'].sensible(timeIndex=-1)
            #sig0 = m.var['sigstZ'].sensible(timeIndex=-1, locIndex=1)
            hR = m.var['halfMassStars'].sensible(timeIndex=-1)/1.7
            sig0 = m.var['sigstZ'].atR( [hR/8.0], m.var['r'].sensible(timeIndex=-1), -1 )
            sigZ = m.var['sigstZ'].sensible(timeIndex=-1)
            r = m.var['r'].sensible(timeIndex=-1)
            ax.plot(r/hR, sigZ/sig0, c=(0.99*sig0/maxsig0,0,1-.99*sig0/maxsig0),  lw=1, alpha=0.5)
            ax.plot(r/hR, np.exp(-r/(2.0*hR)), c='k', lw=1, alpha=0.1 )
    ax.set_xlabel(r'$r/h_R$')
    ax.set_ylabel(r'$\sigma_{z,*}/(\sigma_{z,*}|_{r=0})$')
    ax.set_ylim(0.05, 3.0)
    ax.set_xlim(0,3)
    ax.set_yscale('log')
    plt.savefig('rf19x_radialProfiles.pdf')
    plt.close(fig)

def appendletter(theList, theLetter):
    for i,ele in enumerate(theList):
        theList[i] = ele+theLetter
    return theList

def plotadhoc():


    ## VARY the cloud heating rate
    #rf203 = NewSetOfExperiments(rf202, 'rf203')
    #rf204 = NewSetOfExperiments(rf202, 'rf204')
    #rf205 = NewSetOfExperiments(rf202, 'rf205')
    # Vary spiral heating
    #rf206 = NewSetOfExperiments(rf202, 'rf206') # weak spirals & GI
    #rf207 = NewSetOfExperiments(rf202, 'rf207') # strong spirals & GI
    #rf208 = NewSetOfExperiments(rf202, 'rf208') # default spirals & weak GI

    fig,ax = plt.subplots()
    exps = ['rf202a', 'rf202b', 'rf202c', 'rf202d'] + appendletter(['rf203', 'rf204', 'rf205']+[ 'rf206', 'rf207'], 'b')
    colors = ['gray', 'k', 'gray', 'gray'] + [ (0.8,0,0), (0.9,0,0), (1,0,0) ] + [(0,0,0.8), (0,0,1)]
    linestyles = ['--','-','--','--'] + ['-.','--','-'] + ['-.', '-']
    labels = [None, 'Fiducial', 'Different Posterior Draw', None] + [None,None,'Adjusted Cloud Heating'] + [None, 'Adjusted Spiral Heating']
    plot( exps, colors, linestyles, ax, sigR=True, vsAge=True, rEff=1.0/1.7, labels=labels)
    plotHolmberg(ax, sigU=True, alpha=0.3)
    ax.set_xlabel(r'Stellar Age (yr)')
    ax.set_ylabel(r'$\sigma_r\ (\mathrm{km}\ \mathrm{s}^{-1})$')
    ax.set_yscale('log')
    ax.set_ylim(9.5, 140)
    ax.legend(frameon=False, fontsize=9)
    ax.set_xlim(0,1.28e10)
    ax.text(1.05e10, 20.0, 'No Spirals', color='b')
    ax.text(1.0e10, 119.0, r'$Q_\mathrm{crit}=5.0$', color='b')
    ax.text(1.0e10, 60.0, r'$Q_\mathrm{crit}=2.4$', color='k')
    plt.savefig('rf20x_mass_history_sigR_age.pdf')
    plt.close()

    fig,ax = plt.subplots()
    exps = ['rf202b'] + appendletter(['rf203', 'rf204', 'rf205']+[ 'rf206', 'rf207'], 'b')
    colors = ['k'] + [ (0.8,0,0), (0.9,0,0), (1,0,0) ] + [(0,0,0.8), (0,0,1)]
    linestyles = ['-'] + ['-.','--','-'] + ['-.', '-']
    labels = [ 'Fiducial'] + [None,None,'Adjusted Cloud Heating'] + [None, 'Adjusted Spiral Heating']
    plot( exps, colors, linestyles, ax, sigR=True, vsAge=True, rEff=1.0/1.7, labels=labels)
    plotHolmberg(ax, sigU=True, alpha=0.3)
    ax.set_xlabel(r'Stellar Age (yr)')
    ax.set_ylabel(r'$\sigma_r\ (\mathrm{km}\ \mathrm{s}^{-1})$')
    ax.set_yscale('log')
    ax.set_ylim(9.5, 140)
    ax.legend(frameon=False, fontsize=10)
    ax.set_xlim(0,1.28e10)
    ax.text(1.05e10, 20.0, 'No Spirals', color='b')
    ax.text(1.0e10, 119.0, r'$Q_\mathrm{crit}=5.0$', color='b')
    ax.text(1.0e10, 60.0, r'$Q_\mathrm{crit}=2.4$', color='k')
    plt.savefig('rf20x_mass_history_sigR_age_nodraw.pdf')
    plt.close()
    return

    fig,ax = plt.subplots()
    plot( ['rf198'], ['r'], ['-'], ax, sigR=True, vsAge=True, rEff=1.0/1.7)
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf198_mass_history_sigR_age.pdf')
    plt.close()


    fig,ax = plt.subplots()
    plot( ['rf198'], ['r'], ['-'], ax, ratio=True, vsAge=True, rEff=1.0/1.7)
    plotHolmberg(ax, ratioU=True)
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf198_mass_history_ratio_age.pdf')
    plt.close()


    fig,ax = plt.subplots()
    plot( ['rf200'], ['r'], ['-'], ax, sigR=True, vsAge=True, rEff=1.0/1.7)
    plotHolmberg(ax, sigU=True)
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf200_mass_history_sigR_age.pdf')
    plt.close()


    fig,ax = plt.subplots()
    plot( ['rf200'], ['r'], ['-'], ax, ratio=True, vsAge=True, rEff=1.0/1.7)
    plotHolmberg(ax, ratioU=True)
    ax.set_xlabel(r'Age (Gyr)')
    ax.set_ylabel(r'$\sigma_r$')
    plt.savefig('rf200_mass_history_ratio_age.pdf')
    plt.close()

if __name__=='__main__':
    # everything
    #plotMassScaling(['rf194', 'rf195', 'rf196','rf197', 'rf198', 'rf199'], ['purple','orange','green','k','r', 'b'])
    # fiducial + ad-hoc Q modification for comparison
    #plotMassScaling(['rf194','rf195','rf196','rf198'], ['k','k','k','r'])
    # just the ad-hoc version
    #plotMassScaling(['rf198'], ['r'])
    plotadhoc()
    #plotSigR_nicer()
    #plotSigR_all()

