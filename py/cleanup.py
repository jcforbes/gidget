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
        if 'fake' not in dirname:
            dirs.add(dirname)
    return list(dirs)

def moveFiles(keyname=None):
    if(keyname is None):
        base = '*'
    else:
        base = '*'+keyname+'*'


    pdfs = glob.glob(base+'pdf')
    epss = glob.glob(base+'eps')
    mpgs = glob.glob(base+'mpg')
    movs = glob.glob(base+'mov')
    gifs = glob.glob(base+'gif')
    pngs = glob.glob(base+'png')

    allFiles = pdfs+epss+mpgs+pngs+movs+gifs

    ldirs = GetDirnames(allFiles)

    if keyname is None:
        print ("Directories: ",ldirs)

    for dirname in ldirs:
        if(not os.path.exists(dirname)):
            os.makedirs(dirname)
        for f in allFiles:
            if(dirname in f):
                destname = dirname+'/'+f
                if(os.path.exists(destname)):
                    os.remove(destname)
                shutil.move(f,destname)
                if keyname is None:
                    print ("Moving ",f," to ",destname)


if __name__=="__main__":

    moveFiles()






