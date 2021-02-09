from readoutput import *
import pdb
import math
from collections import Counter

def isprime(N):
    ''' Return True if N is prime, False if it's composite.'''
    if N<2:
        return False
    if N==2:
        return True
    if N%2==0:
        return False
    for i in range(3,int(math.sqrt(float(N)))+1,2):
        if(N%i==0):
            return False
    return True

def factors(N):
    ''' Return the set of factors of an integer N'''
    #return set(reduce(list.__add__,([i,N/i] for i in range(1, int(math.sqrt(N))+1) if N%i==0)))
    facs = []
    for i in range(1, int(np.sqrt(N))+1):
        if N%i==0:
            facs.append(int(i))
            facs.append(int(N/i))
    return list(Counter(facs).keys())

def getSubplots(N):
    if(N==1):
        return 1,1
    if isprime(N) and N>3:
        N+=1
    fac = sorted(list(factors(N)))
    if(len(fac) % 2 == 1): # perfect square
        return (fac[1], fac[1])
    col,row = (fac[int(len(fac)/2)], fac[int(len(fac)/2-1)])
    if(float(col)/float(row) > 2.0):
        return getSubplots(N+1)
    else:
        return col,row


def budgetAM(models, timeIndex=None, name=None, sortby=None, ncols=None, nrows=None):
    ts = len(models[0].getData('t'))
    if(timeIndex is None):
        timeIndex = range(1,ts)
        timeIndex.append(ts-1)
    nmodels=len(models)
    if ncols is None and nrows is None:
        ncols,nrows = getSubplots(nmodels)
    #sqr = int(math.ceil(sqrt(float(nmodels))))
    dirname = ''
    if(name is None):
        if(len(models)<5):
            for model in models:
                dirname+=model.name+'_'
        else:
            dirname+=models[0].name+'_'+str(len(models))+'_'
    else:
        dirname+=name+'_'
    dirname+='budgetAM'
    if(not os.path.exists(dirname)):
        os.makedirs(dirname)
    if sortby is None or sortby not in models[0].p.keys():
        sortby = 'Mh0'
    stmodels = sorted(models, key=lambda sb: sb.p[sortby])
    counter=0
    lw = 4.0/(1.0+np.log10(nmodels))
    fs = 20.0/(1.0+np.log10(nmodels))
    fig,ax = plt.subplots(nrows,ncols,sharey=False,figsize=(16,12), squeeze=False)
    #if(nrows==ncols and nrows==1):
    #    ax=np.array([[ax]])
    
    for i,model in enumerate(stmodels):
        t = model.getData('t')
        z = model.getData('z')
        row = int(i/ncols)
        col = i-row*ncols
        JStars = model.getData('JStars')
        JH2 = model.getData('JH2')
        JHI = model.getData('JHI')
        #shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)

        bottomBound = JStars*0
        starsBound = JStars
        H2Bound = JStars+JH2
        HIBound = JStars+JH2+JHI



        if(isinstance(ax,np.ndarray)):
            theAx = ax[row,col]
        else:
            theAx = ax

        theAx.fill_between(t,bottomBound,starsBound,facecolor='red')
        theAx.fill_between(t,starsBound,H2Bound,facecolor='blue')
        #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
        #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=lightblue)#(mdot<=0.0 or stagnation<=0) )
        theAx.fill_between(t,H2Bound,HIBound,facecolor='orange')



        #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)


        theAx.set_xlim(np.min(t),np.max(t))
        #theAx.set_ylim(-1.0,1.0)

    
        dispsb = "%.3e" % model.p[sortby]

        plt.text(.5,.9,sortby+'='+dispsb,transform=ax[row,col].transAxes,fontsize=fs)

        if(row==nrows-1):
            theAx.set_xlabel('t (Gyr)',fontsize=fs)
        if(col==0):
            theAx.set_ylabel(r'$J (M_\odot\ \mathrm{km}/\mathrm{s}\ \mathrm{kpc})$')
        #ax[row,col].set_ylabel(r'Share of $\partial\Sigma/\partial t$')
    
    # end loop over models

    #ax[0,0].set_title(r'Share of $\partial\Sigma/\partial t$',fontsize=fs)
    plt.tight_layout()
    plt.savefig(dirname+'_'+str(counter).zfill(4)+'.png')
    plt.close(fig)
    # end loop over time.





def balanceAM(models, timeIndex=None, name=None, sortby=None, logR=False, ncols=None, nrows=None):
    ts = len(models[0].getData('t'))
    if(timeIndex is None):
        timeIndex = range(1,ts)
        timeIndex.append(ts-1)
    nmodels=len(models)
    if ncols is None and nrows is None:
        ncols,nrows = getSubplots(nmodels)
    #sqr = int(math.ceil(sqrt(float(nmodels))))
    dirname = ''
    if(name is None):
        if(len(models)<5):
            for model in models:
                dirname+=model.name+'_'
        else:
            dirname+=models[0].name+'_'+str(len(models))+'_'
    else:
        dirname+=name+'_'
    dirname+='balanceAM'
    if(not os.path.exists(dirname)):
        os.makedirs(dirname)
    if sortby is None or sortby not in models[0].p.keys():
        sortby = 'Mh0'
    stmodels = sorted(models, key=lambda sb: sb.p[sortby])
    counter=0
    lw = 4.0/(1.0+np.log10(nmodels))
    fs = 20.0/(1.0+np.log10(nmodels))
    fig,ax = plt.subplots(nrows,ncols,sharey=False,figsize=(16,12), squeeze=False)
    #if(nrows==ncols and nrows==1):
    #    ax=np.array([[ax]])
    
    for i,model in enumerate(stmodels):
        t = model.getData('t')
        z = model.getData('z')
        row = int(i/ncols)
        col = i-row*ncols
        JdotAccr = model.getData('JAccr',cgs=True)
        JdotTrPlus = np.sum( np.clip(model.getData('colTr',cgs=True),0,np.inf) * model.getData('r',cgs=True) * model.getData('vPhi', cgs=True) * model.getData('dA', cgs=True), axis=1 )
        JdotTrMinus = np.sum( np.clip(model.getData('colTr',cgs=True),-np.inf,0) * model.getData('r',cgs=True) * model.getData('vPhi', cgs=True) * model.getData('dA', cgs=True), axis=1 )
        JdotSFR = model.getData('JSFR',cgs=True)*model.p['RfREC'] + model.getData('JOut',cgs=True)
        #shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)

        sfrBound = np.clip(JdotTrMinus, -np.inf, 0)
        bottomBound = -1.0*JdotSFR + sfrBound
        trBound = np.clip(JdotTrPlus, 0, np.inf)
        upperBound = trBound + JdotAccr


        if(isinstance(ax,np.ndarray)):
            theAx = ax[row,col]
        else:
            theAx = ax

        theAx.fill_between(t,bottomBound,sfrBound,facecolor='red')
        theAx.fill_between(t,sfrBound,trBound,facecolor='blue')
        #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
        #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=lightblue)#(mdot<=0.0 or stagnation<=0) )
        theAx.fill_between(t,trBound,upperBound,facecolor='orange')

        theAx.plot(t,-upperBound, color='gray', lw=2, ls='--')


        #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)


        theAx.set_xlim(np.min(t),np.max(t))
        #theAx.set_ylim(-1.0,1.0)

    
        dispsb = "%.3e" % model.p[sortby]

        plt.text(.5,.9,sortby+'='+dispsb,transform=ax[row,col].transAxes,fontsize=fs)

        if(row==nrows-1):
            theAx.set_xlabel('t (Gyr)',fontsize=fs)
        if(col==0):
            theAx.set_ylabel(r'$\dot{J}_\mathrm{gas} (M_\odot/\mathrm{yr}\ \mathrm{km}/\mathrm{s}\ \mathrm{kpc})$')
        #ax[row,col].set_ylabel(r'Share of $\partial\Sigma/\partial t$')
    
    # end loop over models

    ax[0,0].set_title(r'Share of $\partial J/\partial t$',fontsize=fs)
    plt.tight_layout()
    plt.savefig(dirname+'_'+str(counter).zfill(4)+'.png')
    plt.close(fig)
    # end loop over time.



def balancesig(models, timeIndex=None, name=None, sortby=None, logR=False, ncols=None, nrows=None):
    ts = len(models[0].getData('t'))
    if(timeIndex is None):
        timeIndex = range(1,ts)
        timeIndex.append(ts-1)
    nmodels=len(models)
    if (ncols is None and nrows is None) or ncols*nrows<nmodels:
        ncols,nrows = getSubplots(nmodels)
    #sqr = int(math.ceil(sqrt(float(nmodels))))
    dirname = 'movie_'
    if(name is None):
        if(len(models)<5):
            for model in models:
                dirname+=model.name+'_'
        else:
            dirname+=models[0].name+'_'+str(len(models))+'_'
    else:
        dirname+=name+'_'
    dirname+='balancesig'
    if(not os.path.exists(dirname)):
        os.makedirs(dirname)
    if sortby is None or sortby not in models[0].p.keys():
        sortby = 'Mh0'
    stmodels = sorted(models, key=lambda sb: sb.p[sortby])
    counter=0
    lw = 4.0/(1.0+np.log10(nmodels))
    fs = 20.0/(1.0+np.log10(nmodels))
    for ti in timeIndex:
        fig,ax = plt.subplots(nrows,ncols,sharey=True,figsize=(15,13), squeeze=False)
        #if(nrows==ncols and nrows==1):
        #    ax=np.array([[ax]])
        
        for i,model in enumerate(stmodels):
            r = model.getData('r',timeIndex=ti)
            rb = model.getData('rb',timeIndex=ti)
            row = int(i/ncols)
            col = i-row*ncols

            sigLoss = model.getData('dsigdtLoss', cgs=True, timeIndex=ti)
            sigGI = model.getData('dsigdtGI', cgs=True, timeIndex=ti)
            sigSN = model.getData('dsigdtSN', cgs=True, timeIndex=ti)
            sigAccr = model.getData('dsigdtAccr', cgs=True, timeIndex=ti)
            sigAdv = model.getData('dsigdtAdv', cgs=True, timeIndex=ti)

            theNorm = np.abs(sigLoss) + np.abs(sigGI) + np.abs(sigSN) + np.abs(sigAccr) + np.abs(sigAdv)




            colAccr = model.getData('colAccr',timeIndex=ti,cgs=True)
            colREC = model.getData('colREC', timeIndex=ti, cgs=True)
            colTr = model.getData('colTr',timeIndex=ti,cgs=True)
            colSFR = model.getData('colsfr',timeIndex=ti,cgs=True)*(model.var['MassLoadingFactor'].sensible(timeIndex=ti)+0.77)
            dcoldt = model.getData('dcoldt',timeIndex=ti,cgs=True)
            shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)+np.abs(colREC)

            colsfrBound = np.clip(colTr/shareNorm, -1, 0)
            bottomBound = -1.0*colSFR/shareNorm + colsfrBound
            colTrBound = np.clip(colTr/shareNorm, 0, 1)
            colRecBound = colTrBound + colREC/shareNorm
            upperBound = colRecBound + colAccr/shareNorm

            fullmdot = model.getData('Mdot',timeIndex=ti)
            mdot =fullmdot[0:-1]
            stagnation = fullmdot[1:]*mdot

            mdSign = np.sign(fullmdot)
            # positive Mdot => inwards
            leftarrows=[]
            rightarrows=[]
            currentDirection = mdSign[0]
            lastChange = 0
            for k,el in enumerate(mdSign):
                newDirection = currentDirection
                if el==currentDirection:
                    pass
                elif currentDirection==1:
                    leftarrows.append([lastChange,k])
                    newDirection = el
                    lastChange = k
                elif currentDirection==-1:
                    rightarrows.append([lastChange,k])
                    newDirection = el
                    lastChange = k
                elif currentDirection==0:
                    newDirection = el
                    lastChange = k
                else:
                    print ("Problem in balanceplots.py!")
                    pdb.set_trace()
                currentDirection = newDirection
            if currentDirection==1:
                leftarrows.append([lastChange,k])
            if currentDirection==-1:
                rightarrows.append([lastChange,k])


            #lightblue = []
            #for i,stg in enumerate(stagnation):
            #    lightblue.append(mdot[i]<=0.0 or stg<=0)
            #lightblue=np.array(lightblue)

            if(isinstance(ax,np.ndarray)):
                theAx = ax[row,col]
            else:
                theAx = ax

            theAx.fill_between(r, sigLoss/theNorm, np.clip(sigAdv/theNorm, None, 0), facecolor='red', label='Dissipation')
            theAx.fill_between(r, np.clip(sigAdv/theNorm, None, 0), np.clip(sigAdv/theNorm, 0, None), facecolor='gray', label='Advection')
            theAx.fill_between(r,  np.clip(sigAdv/theNorm, 0, None) , np.clip(sigAdv/theNorm, 0, None)+ sigGI/theNorm, facecolor='blue', label='GI Heating')
            theAx.fill_between(r, np.clip(sigAdv/theNorm, 0, None)+sigGI/theNorm, np.clip(sigAdv/theNorm, 0, None) + (sigGI+sigSN)/theNorm, facecolor='pink', label='SN')
            theAx.fill_between(r, np.clip(sigAdv/theNorm, 0, None)+(sigGI+sigSN)/theNorm,  np.clip(sigAdv/theNorm, 0, None) + (sigGI+sigSN+sigAccr)/theNorm, facecolor='orange', label='Accretion')


###            theAx.fill_between(r,bottomBound,colsfrBound,facecolor='red')
###            theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue')
###            #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
###            #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=lightblue)#(mdot<=0.0 or stagnation<=0) )
###            theAx.fill_between(r,colTrBound,colRecBound,facecolor='pink')
###            theAx.fill_between(r,colRecBound,upperBound,facecolor='orange')
            for la in leftarrows:
                delt = rb[la[1]] - rb[la[0]]
                try:
                    theAx.arrow(rb[la[1]],.96, -delt*.9,0, head_width=.05, head_length=.1*delt, fc="DarkSlateGray", ec="DarkSlateGray",lw=3)
                except:
                    print ("Failed to draw left arrow")
            for ra in rightarrows:
                delt = rb[ra[1]] - rb[ra[0]]
                try:
                    theAx.arrow(rb[ra[0]],.96, delt*.9,0, head_width=.05, head_length=.1*delt, fc="DarkGray", ec="DarkGray",lw=2)
                except:
                    print ("Failed to draw right arrow")


            for kk in range(len(r)):
                assert np.abs(upperBound[kk] - bottomBound[kk] - 1.0) < 1.0e-5

            #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)

            theAx.plot(r,r*0+0.5,'--',c=(0.5,0.5,0.5),lw=lw)
            theAx.plot(r,r*0-0.5,'--',c=(0.5,0.5,0.5),lw=lw)

            racc = model.getData('accretionRadius',timeIndex=ti)
            theAx.plot([racc,racc],[-1,1],'--',lw=lw)

            theAx.set_xlim(np.min(r),np.max(r))
            theAx.set_ylim(-1.0,1.0)

            if logR:
                theAx.set_xscale('log')
        
            dispsb = "%.3e" % model.p[sortby]

            plt.text(.5,.9,sortby+'='+dispsb,transform=ax[row,col].transAxes,fontsize=fs)

            if(row==nrows-1):
                theAx.set_xlabel('r (kpc)',fontsize=fs)
            #ax[row,col].set_ylabel(r'Share of $\partial\Sigma/\partial t$')
        
        # end loop over models

        ax[0,0].set_title(r'Share of $\partial\sigma/\partial t$',fontsize=fs)
        ax[0,0].legend(loc='lower right')
        z=models[0].getData('z')
        dispz = "%.3f" % z[ti]
        plt.text(.6,0.1,'z='+dispz,transform=ax[nrows-1,ncols-1].transAxes,fontsize=fs)
        plt.tight_layout()
        plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
        counter = counter+1
        plt.close(fig)
    # end loop over time.
    makeMovies(dirname[7:])



def balance(models, timeIndex=None, name=None, sortby=None, logR=False, ncols=None, nrows=None):
    ts = len(models[0].getData('t'))
    if(timeIndex is None):
        timeIndex = range(1,ts)
        timeIndex.append(ts-1)
    nmodels=len(models)
    if (ncols is None and nrows is None) or ncols*nrows<nmodels:
        ncols,nrows = getSubplots(nmodels)
    #sqr = int(math.ceil(sqrt(float(nmodels))))
    dirname = 'movie_'
    if(name is None):
        if(len(models)<5):
            for model in models:
                dirname+=model.name+'_'
        else:
            dirname+=models[0].name+'_'+str(len(models))+'_'
    else:
        dirname+=name+'_'
    dirname+='balance'
    if(not os.path.exists(dirname)):
        os.makedirs(dirname)
    if sortby is None or sortby not in models[0].p.keys():
        sortby = 'Mh0'
    stmodels = sorted(models, key=lambda sb: sb.p[sortby])
    counter=0
    lw = 4.0/(1.0+np.log10(nmodels))
    fs = 20.0/(1.0+np.log10(nmodels))
    for ti in timeIndex:
        fig,ax = plt.subplots(nrows,ncols,sharey=True,figsize=(15,12), squeeze=False)
        if(nrows==ncols and nrows==1):
            ax=np.atleast_2d( np.array([[ax]]) )
        
        for i,model in enumerate(stmodels):
            r = model.getData('r',timeIndex=ti)
            rb = model.getData('rb',timeIndex=ti)
            row = int(i/ncols)
            col = i-row*ncols
            colAccr = model.getData('colAccr',timeIndex=ti,cgs=True)
            colREC = model.getData('colREC', timeIndex=ti, cgs=True)
            colTr = model.getData('colTr',timeIndex=ti,cgs=True)
            colSFR = model.getData('colsfr',timeIndex=ti,cgs=True)*(model.var['MassLoadingFactor'].sensible(timeIndex=ti)+0.77)
            dcoldt = model.getData('dcoldt',timeIndex=ti,cgs=True)
            shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)+np.abs(colREC)

            colsfrBound = np.clip(colTr/shareNorm, -1, 0)
            bottomBound = -1.0*colSFR/shareNorm + colsfrBound
            colTrBound = np.clip(colTr/shareNorm, 0, 1)
            colRecBound = colTrBound + colREC/shareNorm
            upperBound = colRecBound + colAccr/shareNorm

            fullmdot = model.getData('Mdot',timeIndex=ti)
            mdot =fullmdot[0:-1]
            stagnation = fullmdot[1:]*mdot

            mdSign = np.sign(fullmdot)
            # positive Mdot => inwards
            leftarrows=[]
            rightarrows=[]
            currentDirection = mdSign[0]
            lastChange = 0
            for k,el in enumerate(mdSign):
                newDirection = currentDirection
                if el==currentDirection:
                    pass
                elif currentDirection==1:
                    leftarrows.append([lastChange,k])
                    newDirection = el
                    lastChange = k
                elif currentDirection==-1:
                    rightarrows.append([lastChange,k])
                    newDirection = el
                    lastChange = k
                elif currentDirection==0:
                    newDirection = el
                    lastChange = k
                else:
                    print ("Problem in balanceplots.py!")
                    pdb.set_trace()
                currentDirection = newDirection
            if currentDirection==1:
                leftarrows.append([lastChange,k])
            if currentDirection==-1:
                rightarrows.append([lastChange,k])


            #lightblue = []
            #for i,stg in enumerate(stagnation):
            #    lightblue.append(mdot[i]<=0.0 or stg<=0)
            #lightblue=np.array(lightblue)

            if(isinstance(ax,np.ndarray)):
                theAx = ax[row,col]
            else:
                theAx = ax

            theAx.fill_between(r,bottomBound,colsfrBound,facecolor='red')
            theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue')
            #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
            #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=lightblue)#(mdot<=0.0 or stagnation<=0) )
            theAx.fill_between(r,colTrBound,colRecBound,facecolor='pink')
            theAx.fill_between(r,colRecBound,upperBound,facecolor='orange')
            for la in leftarrows:
                delt = rb[la[1]] - rb[la[0]]
                try:
                    theAx.arrow(rb[la[1]],.96, -delt*.9,0, head_width=.05, head_length=.1*delt, fc="DarkSlateGray", ec="DarkSlateGray",lw=3)
                except:
                    print ("Failed to draw left arrow")
            for ra in rightarrows:
                delt = rb[ra[1]] - rb[ra[0]]
                try:
                    theAx.arrow(rb[ra[0]],.96, delt*.9,0, head_width=.05, head_length=.1*delt, fc="DarkGray", ec="DarkGray",lw=2)
                except:
                    print ("Failed to draw right arrow")


            for kk in range(len(r)):
                assert np.abs(upperBound[kk] - bottomBound[kk] - 1.0) < 1.0e-5

            #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)

            theAx.plot(r,r*0+0.5,'--',c=(0.5,0.5,0.5),lw=lw)
            theAx.plot(r,r*0-0.5,'--',c=(0.5,0.5,0.5),lw=lw)

            racc = model.getData('accretionRadius',timeIndex=ti)
            theAx.plot([racc,racc],[-1,1],'--',lw=lw)

            theAx.set_xlim(np.min(r),np.max(r))
            theAx.set_ylim(-1.0,1.0)

            if logR:
                theAx.set_xscale('log')
        
            dispsb = "%.3e" % model.p[sortby]

            plt.text(.5,.9,sortby+'='+dispsb,transform=ax[row,col].transAxes,fontsize=fs)

            if(row==nrows-1):
                theAx.set_xlabel('r (kpc)',fontsize=fs)
            #ax[row,col].set_ylabel(r'Share of $\partial\Sigma/\partial t$')
        
        # end loop over models

        ax[0,0].set_title(r'Share of $\partial\Sigma/\partial t$',fontsize=fs)
        z=models[0].getData('z')
        dispz = "%.3f" % z[ti]
        plt.text(.6,0.1,'z='+dispz,transform=ax[nrows-1,ncols-1].transAxes,fontsize=fs)
        plt.tight_layout()
        plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
        counter = counter+1
        plt.close(fig)
    # end loop over time.
    makeMovies(dirname[7:])
        




