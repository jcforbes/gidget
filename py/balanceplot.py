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
    if sortby is None:
        sortby = 'Mh0'
    stmodels = sorted(models, key=lambda sb: sb.p[sortby])
    counter=0
    for ti in timeIndex:
        fig,ax = plt.subplots(nrows,ncols,sharey=True,figsize=(16,12))
        
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
            mdot = model.getData('Mdot',timeIndex=ti)[0:-1]

            ax[row,col].fill_between(r,bottomBound,colsfrBound,facecolor='red')
            ax[row,col].fill_between(r,colsfrBound,colTrBound,facecolor='blue',where=mdot>0.0)
            ax[row,col].fill_between(r,colsfrBound,colTrBound,facecolor='lightblue',where=mdot<=0.0)
            ax[row,col].fill_between(r,colTrBound,upperBound,facecolor='orange')

            #ax[row,col].plot(r,dcoldt/shareNorm,'yellow',lw=2)

            ax[row,col].plot(r,r*0+0.5,'--',c=(0.5,0.5,0.5))
            ax[row,col].plot(r,r*0-0.5,'--',c=(0.5,0.5,0.5))

            racc = model.getData('scaleRadius',timeIndex=ti)
            ax[row,col].plot([racc,racc],[-1,1],'--')

            ax[row,col].set_xlim(np.min(r),np.max(r))
            ax[row,col].set_ylim(-1.0,1.0)
        
            dispmh0 = "%.3e" % model.p['Mh0']

            plt.text(.5,.9,'Mh0='+dispmh0,transform=ax[row,col].transAxes,fontsize=10)

            if(row==nrows-1):
                ax[row,col].set_xlabel('r')
            #ax[row,col].set_ylabel(r'Share of $\partial\Sigma/\partial t$')

        ax[0,0].set_title(r'Share of $\partial\Sigma/\partial t$')
        z=models[0].getData('z')
        dispz = "%.3f" % z[ti]
        plt.text(.6,0.1,'z='+dispz,transform=ax[nrows-1,ncols-1].transAxes)
        plt.tight_layout()
        plt.savefig(dirname+'/frame_'+str(counter).zfill(4)+'.png')
        counter = counter+1
        plt.clf()
        del fig
        del ax




