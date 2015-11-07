import glob
import pdb
import os
import subprocess
import time

# from exper.py
def HowManyStillRunning(procs):
    ''' Given a list of processes created with subprocess.Popen, (each of which
          in this case will be a single run of the gidget code), count up how
          many of those runs are still going.'''
    nStillRunning=0
    for proc in procs:
        if(proc.poll()==None):
            nStillRunning+=1
    return nStillRunning 


def makeMovies(keyname = None):
    print "running makeMovies in makeGidgetMovies.py with keyname ",keyname
    # What subdirectory are we in?
    subdirname = os.getcwd()[os.getcwd().rfind("/")+1:] # analysis typically
    # Find all the directories in which the IDL code has put frames for movies
    if(keyname is None):
        movieDirs = glob.glob('movie_*')
    else:
        movieDirs = glob.glob('movie*'+keyname+'*')
    

    movieNames = []
    procs = []
    ctr=0

    nulfp = open(os.devnull, "w")

    # For each directory, make a movie of the same name sans the prefix "movie_"
    # Only use 12 processors at a time.
    for movieDir in movieDirs:
        print "WORKING ON MOVIE: ",movieDir
        ctr=ctr+1
        regstring = movieDir+'/frame_????.png'
        movieName = movieDir[6:]+".gif"
        subprocess.call(["rm","-f",movieName])
        if(keyname is None):
            print "Producing movie #",ctr,"of",len(movieDirs)
        #print "processing movieDir ",movieDir
        #print 'running ffmpeg with regstring ',regstring
        #print 'to create movie ',movieName
        #procs.append(subprocess.Popen(["ffmpeg","-loglevel","quiet","-i",regstring,"-vcodec","qtrle",movieName],stderr=nulfp))
        procs.append(subprocess.Popen(["convert","-delay","10","-loop","0",regstring,movieName],stderr=nulfp))
        print "running convert -delay 3 -loop 0 "+regstring+" "+movieName
    #    procs.append(subprocess.Popen(["ffmpeg","-loglevel","quiet","-f","image2","-qscale","1","-i",regstring,movieName],stderr=nulfp))
    #    procs.append(subprocess.Popen(["ffmpeg","-f","image2","-qscale","0","-i",regstring,movieName]))
    #    pdb.set_trace()
        movieNames.append(movieName)
        while True:
            nStillRunning = HowManyStillRunning(procs)
            if(nStillRunning >= 12):
                time.sleep(5) # wait a moment
            else:
                break
    # Do not let the IDL code move on until all movies are made. This way, if this script is called again
    # immediately we don't get a situation where the new script remakes the movies this instance has 
    # already made.
    nPrev=0
    while True:
        nStillRunning=HowManyStillRunning(procs)
        # has anything changed since the last time we checked?
        if(nStillRunning == nPrev and nStillRunning != 0):
            # do nothing except wait a little bit
            time.sleep(5)
        else:
            nPrev=nStillRunning
            if(nPrev == 0):
                break # we're done!
            if(keyname is None):
                print "Still waiting for ",nPrev, " processes to finish; I'll check every few seconds for changes."


    if(keyname is None):
        print 
        print
        print "Movies produced!"
        print "Movies produced: ",movieNames
    
        print "Removing the directories containing movie frames."
    for movieDir in movieDirs:
        pass
        #subprocess.call(["rm","-rf",movieDir])

if __name__=='__main__':
    makeMovies()

