import numpy as np
import pdb
class bolshoireader:
    def __init__(self,fn, minmass, maxmass, bolshoidir):
        self.filename = fn
        self.bolshoidir = bolshoidir
        trees=[]
        with open(bolshoidir+fn,'r') as f:
            line = f.readline()
            while line!='':
                sline = line.split()
                treenum = int(sline[0])
                mh0 = float(sline[1])
                weightedmdot = float(sline[2])
                if minmass<mh0 and mh0<maxmass:
                    trees.append([treenum, mh0, weightedmdot])
                line=f.readline()
        self.trees = sorted(trees, key=lambda tree: tree[2])
        self.keys = [tree[2] for tree in self.trees]
        self.masses = np.log10([tree[1] for tree in self.trees]) # log10(M_h)
        self.minimumMass = np.min(self.masses)
        self.maximumMass = np.max(self.masses)
        self.stored = False

    def storeAll(self):
        self.xs = np.zeros( (len(self.trees), 1000) )
        with open(self.bolshoidir+'rf_data4.txt', 'r') as f:
            all_lines = f.readlines()
        treecounter = 0
        bad_rows = []
        for treei, tree in enumerate(self.trees):
            treecounter+=1
            if treecounter % 100 == 0:
                print ("Loading up tree ", treecounter, ' of ', len(self.trees))
            treestring = 'tree '+repr(tree[0])+'\n'
            treeind = all_lines.index(treestring)
            try:
                self.xs[treei, : ] = [float( all_lines[lineind].split()[0] ) for lineind in range(treeind+1,treeind+1001) ]
            except:
                print ("Ran into trouble!")
                bad_rows.append(treei)
        for k, row in enumerate(bad_rows):
            random_ind = int((len(self.trees)-1)*np.random.random())
            if not random_ind in bad_rows:
                self.xs[row, :] = self.xs[random_ind, :]
        self.stored = True
                    

    def returnnearest(self, dirname, filename, value, logMh0):
        ''' copy a file with the closest value of tree[2]'''
        #ind = np.searchsorted(self.keys, value)
        assert value>=0 and value<=1
        if logMh0 < self.minimumMass+0.25:
            acceptable = self.masses<self.minimumMass+0.5
        elif logMh0 > self.maximumMass-0.25:
            acceptable = self.masses>self.maximumMass-0.5
        else:
            acceptable = np.logical_and( self.masses>logMh0-0.25, self.masses<logMh0+0.25 )

        nAcceptable = np.sum(np.ones(len(self.keys))[acceptable])
        indAccept = int(nAcceptable*value) # pick one of the acceptable trees
        ind = np.array(range(len(self.keys)))[acceptable][indAccept] 

        if ind == len(self.keys):
            ind=ind-1
        assert ind<len(self.keys) and ind>-1
        #rffile = bolshoidir+ repr(self.trees[ind][0])+'_rf.txt'
        #shutil.copy(rffile, dirname)
        # This chunk of code searches the file rf_data.txt for 
        #  the appropriate tree number from the register as identified by ind above.
        # Once found, the subsequent lines (containing the values of the random factors
        #  for the main progenitor accretion history of that tree) are copied to a file
        #  in the working directory of the gidget run.

        # If we've already cached the x values, use that!
        if self.stored:
            return self.xs[ind,:]
        else:
            xs = []
            with open(self.bolshoidir+'rf_data4.txt','r') as f:
                line = f.readline()
                while line!='':
                    if 'tree '+repr(self.trees[ind][0]) in line:
                        print ("Copying tree with index ",ind," and Bolshoi ID ",self.trees[ind][0])
                        break # we found the tree we're looking for
                    line = f.readline()
                if line=='':
                    print ("Failed to find line tree "+repr(self.trees[ind][0]))
                line = f.readline()
                collectlines=[]
                while not 'tree' in line and line!='':
                    xs.append( float(line.split()[0]) )
                    collectlines.append(line)
                    line = f.readline()
    #            with open(dirname+'/'+filename,'w') as ff:
    #                assert len(collectlines)>0
    #                ### Major change! Add in a line at the beginning which tells us the z=0 halo mass of the bolshoi tree from which these random factors are drawn
    #                ff.write(str(self.trees[ind][1])+'\n')
    #                for line in collectlines:
    #                    ff.write(line)
            return xs[:-1]
    def copynearest(self, dirname, filename, value, logMh0):
        ''' copy a file with the closest value of tree[2]'''
        #ind = np.searchsorted(self.keys, value)
        assert value>=0 and value<=1
        if logMh0 < self.minimumMass+0.25:
            acceptable = self.masses<self.minimumMass+0.5
        elif logMh0 > self.maximumMass-0.25:
            acceptable = self.masses>self.maximumMass-0.5
        else:
            acceptable = np.logical_and( self.masses>logMh0-0.25, self.masses<logMh0+0.25 )

        nAcceptable = np.sum(np.ones(len(self.keys))[acceptable])
        indAccept = int(nAcceptable*value) # pick one of the acceptable trees
        ind = np.array(range(len(self.keys)))[acceptable][indAccept] 

        if ind == len(self.keys):
            ind=ind-1
        assert ind<len(self.keys) and ind>-1
        #rffile = bolshoidir+ repr(self.trees[ind][0])+'_rf.txt'
        #shutil.copy(rffile, dirname)
        # This chunk of code searches the file rf_data.txt for 
        #  the appropriate tree number from the register as identified by ind above.
        # Once found, the subsequent lines (containing the values of the random factors
        #  for the main progenitor accretion history of that tree) are copied to a file
        #  in the working directory of the gidget run.
        with open(self.bolshoidir+'rf_data4.txt','r') as f:
            line = f.readline()
            while line!='':
                if 'tree '+repr(self.trees[ind][0]) in line:
                    print ("Copying tree with index ",ind," and Bolshoi ID ",self.trees[ind][0])
                    break # we found the tree we're looking for
                line = f.readline()
            if line=='':
                print ("Failed to find line tree "+repr(self.trees[ind][0]))
            line = f.readline()
            collectlines=[]
            while not 'tree' in line and line!='':
                collectlines.append(line)
                line = f.readline()
            with open(dirname+'/'+filename,'w') as ff:
                assert len(collectlines)>0
                ### Major change! Add in a line at the beginning which tells us the z=0 halo mass of the bolshoi tree from which these random factors are drawn
                ff.write(str(self.trees[ind][1])+'\n')
                for line in collectlines:
                    ff.write(line)
    def test(self):
        for i in range(len(self.keys)):
            globalBolshoiReader.copynearest('~/bolshoitest/', 'test_'+str(i).zfill(5)+'_inputRandomFactors.txt', i)


