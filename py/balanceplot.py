from readoutput import *
import math

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
    return set(reduce(list.__add__,([i,N/i] for i in range(1, int(math.sqrt(N))+1) if N%i==0)))

def getSubplots(N):
    if(N==1):
        return 1,1
    if isprime(N) and N>3:
        N+=1
    fac = sorted(list(factors(N)))
    if(len(fac) % 2 == 1): # perfect square
        return (fac[1], fac[1])
    col,row = (fac[len(fac)/2], fac[len(fac)/2-1])
    if(float(col)/float(row) > 2.0):
        return getSubplots(N+1)
    else:
        return col,row


def balance(models, timeIndex=None, name=None, sortby=None):
    ts = len(models[0].getData('t'))
    if(timeIndex is None):
        timeIndex = range(1,ts)
        timeIndex.append(ts-1)
    nmodels=len(models)
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
        fig,ax = plt.subplots(nrows,ncols,sharey=True,figsize=(16,12))
        if(nrows==ncols and nrows==1):
            ax=np.array([[ax]])
        
        for i,model in enumerate(stmodels):
            r = model.getData('r',timeIndex=ti)
            row = int(i/ncols)
            col = i-row*ncols
            colAccr = model.getData('colAccr',timeIndex=ti,cgs=True)
            colTr = model.getData('colTr',timeIndex=ti,cgs=True)
            colSFR = model.getData('colsfr',timeIndex=ti,cgs=True)*(model.p['mu']+model.p['RfREC'])
            dcoldt = model.getData('dcoldt',timeIndex=ti,cgs=True)
            shareNorm = np.abs(colAccr)+np.abs(colTr)+np.abs(colSFR)
            colsfrBound = np.clip(colTr,None,0)/shareNorm
            bottomBound = -1.0*colSFR/shareNorm + colsfrBound
            colTrBound = np.clip(colTr,0,None)/shareNorm
            upperBound = colTrBound + colAccr/shareNorm
            fullmdot = model.getData('Mdot',timeIndex=ti)
            mdot =fullmdot[0:-1]
            stagnation = fullmdot[1:]*mdot

            lightblue = []
            for i,stg in enumerate(stagnation):
                lightblue.append(mdot[i]<=0.0 or stg<=0)
            lightblue=np.array(lightblue)

            if(isinstance(ax,np.ndarray)):
                theAx = ax[row,col]
            else:
                theAx = ax
            theAx.fill_between(r,bottomBound,colsfrBound,facecolor='red')
            theAx.fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
            #theAx.fill_between(r,colsfrBound,colTrBound,facecolor='Gray',where=mdot
            theAx.fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=lightblue)#(mdot<=0.0 or stagnation<=0) )
            theAx.fill_between(r,colTrBound,upperBound,facecolor='orange')

            #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)

            theAx.plot(r,r*0+0.5,'--',c=(0.5,0.5,0.5),lw=lw)
            theAx.plot(r,r*0-0.5,'--',c=(0.5,0.5,0.5),lw=lw)

            racc = model.getData('scaleRadius',timeIndex=ti)
            theAx.plot([racc,racc],[-1,1],'--',lw=lw)

            theAx.set_xlim(np.min(r),np.max(r))
            theAx.set_ylim(-1.0,1.0)
        
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
        




