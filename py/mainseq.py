from readoutput import *

def generateMockSampleCoords(N,MSparams,z):
    masses = drawFromSchechter(N,-1.3,10.0**10.9,10.0**9.0)
    coords = []
    locs=[]
    for mass in masses:
        rnd = random.gauss(0,1)*MSparams[2]
        coord = (mass, 10.0**((MSparams[0]*np.log10(mass)+MSparams[1])+rnd), z)
        locs.append(rnd/MSparams[2])
        coords.append(coord)
    return coords,locs

def imposeCuts(coords,sigmaMS,zrange,massrange,deltamsRange):
    newCoords = []
    newSigmaMS = []
    for i,coord in enumerate(coords):
        if(massrange[0]<coord[0] and coord[0] < massrange[1] and \
           zrange[0]<coord[2] and coord[2] < zrange[1] and \
           deltamsRange[0] < sigmaMS[i] and sigmaMS[i] < deltamsRange[1]):
            newCoords.append(coord)
            newSigmaMS.append(sigmaMS[i])
    return newCoords,newSigmaMS

def drawFromSchechter(N,alpha,Mcutoff, Mmin):
    vals = []
    attempts=0
    while len(vals)<N:
        draw = np.random.gamma(alpha+2,Mcutoff)
        if draw < Mmin:
            pass
        else:
            if(random.random() > Mmin/draw):
                pass
            else:
                vals.append(draw)
        attempts+=1
    print "To generate N= ",N," draws from a Schechter fn required ",attempts," attempts."
    return vals

def testStackAnalysis():
    test2 = Experiment('rg78')
    test2.read()
    z=1.0
    msCoords = test2.msCoords(z)
    msparams, residuals = test2.fitMS(msCoords)
    sampleCoords, locs = generateMockSampleCoords(1000,msparams,z)
    massbounds = [1.0e9,3.14e9,1.0e10,3.14e10,1.0e11]
    masscuts = [(massbounds[i],massbounds[i+1]) for i in range(len(massbounds)-1)]
    for i,cut in enumerate(masscuts):
        sampleCoordsStack, locsStack = imposeCuts(sampleCoords,locs, (-.1,3.0), cut, (-5,5))
        makeStackFrames(test2, sampleCoordsStack, "massStack_"+str(i))
    mscuts = [(-50,-1),(-1,0),(0,1),(1,50)]
    for i,cut in enumerate(mscuts):
        sampleCoordsStack, locsStack = imposeCuts(sampleCoords,locs, (-.1,3.0), (1.0e10,5.0e10), mscuts)
        makeStackFrames(test2, sampleCoordsStack, "msStack_"+str(i)+"_")


def makeStackFrames(expt,coords,stackname):
#    avgProfR = 10.0 * np.power((200.0 - np.arange(200)[::-1])/200.0) # 10 kpc, logarithmic (untested!)
    # N is the number of cells in the radial direction.
    N=500
    avgProfR = 10.0*np.arange(N)/(N-1.0) # 10 kpc linear
    avgProfColSFR = np.zeros(N)
    # For every given coordinate in M* - SFR - z space, find the closest model galaxy
    for i in range(len(coords)):
        if(i % 50 == 0):
            print "Adding stack frame "+str(i)+" of "+str(len(coords))
        mstar,sfr,z = coords[i]
        model,timeIndex = expt.modelClosestTo(mstar,sfr,z)
        model.TwoDimSFR(expt.name+'_'+stackname+'_frame_'+str(i).zfill(4), random.random()*pi/2.0, random.random()*2.0*pi,timeIndex)
        for k in range(N):
            qa = model.quantityAt(avgProfR[k],'colsfr',timeIndex)
            if(type(qa) != type(np.float64(43))):
                pdb.set_trace()
            avgProfColSFR[k] += model.quantityAt(avgProfR[k],'colsfr',timeIndex) / float(len(coords))
    with open(expt.name+'_'+stackname+'_prof.dat','w') as plotfile:
        for k in range(N):
            plotfile.write(str(avgProfR[k])+' '+str(avgProfColSFR[k])+' \n')




if __name__=='__main__':
    testStackAnalysis()

