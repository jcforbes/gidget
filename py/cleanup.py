import shutil,os,glob



def GetDirnames(listOfFiles):
    ''' Given a list of all pdf, eps, mpg, or png in a directory,
    Construct a directory for each, as classified by'''
    dirs=set()
    # For each file
    for f in listOfFiles:
        # Split the filename by _
        parsed = f.rsplit('_')
        dirname=parsed[0]
        dirs.add(dirname)
        ##i=1
        ###i=0
        #### For each part of the filename,
        ###for p in parsed:
        ###    i=i+1
        ###    # iterate until the component is a pure digit - this means we've reached
        ###    # the end of the first part of the filename which tells us which experiments
        ###    # contributed to the data plotted in these files.
        ###    if(p.isdigit()):
        ###        break
        ##if i!=len(parsed):
        ##    dirname = ''
        ##    for j in range(i):
        ##        if(j!=0 or parsed[j] != 'crp'): # exclude crp_ from the dirname
        ##           dirname+=parsed[j]+'_'
        ##    dirname=dirname[0:-1] # remove trailing _
        ##    dirs.add(dirname)
    return list(dirs)


if __name__=="__main__":
    pdfs = glob.glob('*pdf')
    epss = glob.glob('*eps')
    mpgs = glob.glob('*mpg')
    pngs = glob.glob('*png')

    allFiles = pdfs+epss+mpgs+pngs

    ldirs = GetDirnames(allFiles)

    print "Directories: ",ldirs

    for dirname in ldirs:
        if(not os.path.exists(dirname)):
            os.makedirs(dirname)
        for f in allFiles:
            if(dirname in f):
                destname = dirname+'/'+f
                if(os.path.exists(destname)):
                    os.remove(destname)
                shutil.move(f,destname)
                print "The plan is to move ",f," to ",destname









