;; This is the same as variability2 (forked August 18, 2012), except that I'll
;; attempt to add support for an arbitrary number of data arrays.

;; Given a 4-d array of the form discussed in the procedures below,
;; generate a 2-d array of the form [[min,max], [min,max], ...] for each variable
;; where min and max should encompass the entire data range (for all times and models)
FUNCTION simpleranges,data,wrtxlog
  ranges=dblarr(2,n_elements(data[0,0,*,0]))
  FOR k=0,n_elements(data[0,0,*,0])-1 DO BEGIN ;;; loop over variables
    ;; do the simplest thing: look at the min and the max in the data.
    ranges[*,k] = [MIN(data[*,*,k,*]),MAX(data[*,*,k,*])]

    ;; If the variable is supposed to be plotted logarithmically but the data 
    ;; isn't plottable that way (i.e. <=0 in some places), only plot the valid data.
    IF(wrtXlog[k] EQ 1) THEN BEGIN
      flat = data[*,*,k,*]
      ind = WHERE(flat GT 0.0, ct)
      IF(ct NE 0) THEN ranges[0,k] = MIN(flat[ind]) ELSE ranges[0,k]=.1
      IF(ct NE 0) THEN med = median(flat[ind]) ELSE med =1
      IF(ct EQ 0) THEN ranges[*,k] = [.1,1]
      fac = 1.0e9
      ;; Set a maximum dynamic range so that a few bad apples don't spoil the readability of the plot
      IF(ranges[1,k] / ranges[0,k] GT fac) THEN BEGIN
	  IF(med GT ranges[1,k]/fac) THEN ranges[0,k] = ranges[1,k]/fac ELSE ranges[1,k]=ranges[0,k]*fac
      ENDIF
    ENDIF

    ;; Give us a little breathing room: add 10% in each direction for each axis.
    IF(wrtXlog[k] EQ 1) THEN BEGIN
      dexRange = alog10(ranges[1,k]/ranges[0,k])
      ranges[1,k] = ranges[1,k]*(10.0 ^ (dexRange/10.0))
      ranges[0,k] = ranges[0,k]/(10.0 ^ (dexRange/10.0))
    ENDIF
    IF(wrtXlog[k] EQ 0) THEN BEGIN
      linRange = ranges[1,k]-ranges[0,k]
      ranges[1,k] = ranges[1,k] + linRange/10.0
      ranges[0,k] = ranges[0,k] - linRange/10.0
    ENDIF

 
  ENDFOR
  return,ranges
END

;; given an array "data" whose first column indexes time, second column indexes position,
;; fourth column indexes models, and a vector of time labels, and axis labels
;; make one or a few movies! The third column indexes which variable we're plotting.
;; 
;; colors is a (# timesteps) by (# models) array.
PRO simpleMovie,data,time,labels,colors,styles,wrtXlog,name,sv,strt,prev=prev,psym=psym,axislabels=axislabels,taillength=taillength,saveFrames=saveFrames,plotContours=plotContours,whichFrames=whichFrames,texLabels=texLabels,horizontal=horizontal,timeText=timeText,NIndVarBins=NIndVarBins,percentiles=percentiles
  lth=1
  chth=1
  IF(sv EQ 3 || sv EQ 4 || sv EQ 5) THEN cg=1 ELSE cg=0
  IF(sv EQ 1) THEN cs=1 ELSE cs=2.8
  IF(sv EQ 4) THEN BEGIN 
    cs=2.5
    chth=2
    lth=5
  ENDIF

  ;; exclude the first time step from setting the range.
  ranges= simpleranges(data[2:n_elements(time)-1,*,*,*],wrtxlog)
  
  IF(n_elements(horizontal) EQ 0) THEN horizontal = 0
  IF(horizontal GE 2 OR horizontal LT 0) THEN horizontal = 1
  IF(n_elements(axislabels) EQ 0) THEN axislabels=labels[*]
  IF(n_elements(prev) EQ 0) THEN prev=0
  IF(n_elements(taillength) EQ 0) THEN taillength=1
  IF(n_elements(saveFrames) EQ 0) THEN saveFrames=0
  IF(n_elements(whichFrames) EQ 0) THEN whichFrames = indgen(n_elements(time))
  IF(n_elements(texLabels) EQ 0) THEN texLabels = axisLabels
  IF(n_elements(timeText) EQ 0) THEN timeText="z = "
  IF(n_elements(NIndVarBins) EQ 0) THEN NIndVarBins=0
  IF(n_elements(percentiles) EQ 0) THEN percentiles=[.025,.5,.975]

  ;; loop over y-axis variables to be plotted (except the first one, which is just the independent var!)
  FOR k=1,n_elements(labels)-1 DO BEGIN    
    fn=name+"_"+labels[k]
    dn="movie_"+name+"_"+labels[k]
    set_plot,'x'
    IF(sv EQ 3 || sv EQ 4) THEN set_plot,'z'
    IF(cg NE 1) THEN WINDOW,1,XSIZE=1024,YSIZE=1024,RETAIN=2
    IF(sv EQ 3 || sv EQ 4) THEN cgDisplay,1024,1024
    IF(sv EQ 1) THEN mpeg_id=MPEG_OPEN([1024,1024],FILENAME=(fn+".mpg"))
    IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN BEGIN
      FILE_MKDIR,dn
    ENDIF
    IF(sv EQ 5) THEN BEGIN
        wf = n_elements(whichFrames)
	IF (wf GT 10) THEN message,"You have requested an unreasonable number of frames for a single plot."
	cs = .3
        figureInit,(fn+"_timeSeries"),2,1+horizontal,1+(1-horizontal)
	multiplot,[0,wf*horizontal+1,wf*(1-horizontal)+1,0,0]
	!p.noerase=0
    ENDIF

    symsize = cs* 2 / (alog10(n_elements(data[0,0,0,*]))+1)

    ;; counter of frames.
    count=0 
    ; length of a tail to put on each point which represents a galaxy when prev=1
    ; tailLength = fix(3.0/(float(n_elements(styles))/928.0)^(.25)) 

    ;; loop over time: each iteration of this loop will save a frame in a movie
    FOR tii=0,n_elements(whichFrames)-1 DO BEGIN  
      IF(sv EQ 5 and tii NE 0) THEN !p.noerase=1
      ti = whichFrames[tii]
      count=count+1
      SETCT,1,n_elements(styles),2

      ;; If this variable k is included in strt, draw a dashed line y=x.
      ;; Either way, this long statement draws the first plot for this frame.
      ;; Everything else will be overplotted on top in a few lines. 
      xt = axisLabels[0]
      IF(sv EQ 5 AND tii NE n_elements(whichFrames)-1) THEN xt=""
      yt = axisLabels[k]
      ;IF(sv EQ 5) THEN yt=yt+" at z="+str(time)
      IF(strt[k] EQ 0) THEN $
        PLOT,[0],[0],COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=xt, $
         YTITLE=axislabels[k], XRANGE=ranges[*,0],YRANGE=ranges[*,k],ylog=wrtXlog[k], $
         xlog=wrtXlog[0],CHARTHICK=chth,CHARSIZE=cs,THICK=lth,XTHICK=lth,YTHICK=lth $
       ELSE $
        PLOT,findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
	 findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
         COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=xt, $
         YTITLE=axislabels[k],XRANGE=ranges[*,0],YRANGE=ranges[*,k],ylog=wrtXlog[k], $
         xlog=wrtXlog[0],linestyle=2,CHARSIZE=cs,CHARTHICK=chth,THICK=lth,XTHICK=lth,YTHICK=lth
      ;; that was all one line! It's now over, and our plotting space is set up.

      ;; Print the time in the upper left of the screen
      
      theTimeText = timeText+string(time[ti],Format='(D0.3)')
      IF(sv NE 5) THEN XYOUTS,.3,.87,theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth $
      ELSE IF(horizontal EQ 0) THEN XYOUTS,.3,.95-float(tii)*(1.0 - 2.0/8.89)/(float(n_elements(whichFrames))),theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth $
      ELSE XYOUTS,.1+float(tii)/float(n_elements(whichFrames)),theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth

      ;; loop over models (4th column of data)
      FOR j=0,n_elements(data[0,0,0,*])-1 DO BEGIN 

        ;; Set the color table.	
;        IF(n_elements(styles) GT 80) THEN IF(sv NE 3 AND sv NE 4) THEN setct,5,n_elements(styles),j ELSE setct,3,n_elements(styles),0
	setct,1,n_elements(styles),j
        
        ;; Overplot the data for this model. If prev is set, draw a tail to include previous values of the data.
        OPLOT, data[ti,*,0,j], data[ti,*,k,j], COLOR=colors[ti,j],linestyle=styles[j],PSYM=psym,SYMSIZE=symsize,THICK=lth
        IF(prev EQ 1) THEN BEGIN
          OPLOT,data[MAX([0,ti-tailLength]):ti,0,0,j],data[MAX([0,ti-tailLength]):ti,0,k,j], COLOR=colors[ti,j],linestyle=0,THICK=1
        ENDIF
      ENDFOR ;; end loop over models

      ;; All of the models have been plotted. Now if the user has requested it, we would like to 
      ;; bin the data by its independent variable, and then its dependent variable!
      IF(NIndVarBins GT 0) THEN BEGIN
          ;;;; ASSUMING that there are some galaxies in each color bin,
          ncolors = MAX(colors)+1
          theCurves = dblarr(NIndVarBins, n_elements(percentiles), ncolors)
          ;; each color gets its own percentile curves
          FOR aColor=0, ncolors-1 DO BEGIN
              subset = WHERE(aColor EQ colors[ti,*], NInThisSubset)
              binCenters = dblarr(NIndVarBins)
              IF(NInThisSubset GT 0) THEN BEGIN
                  FOR indVarIndex=0, NIndVarBins-1 DO BEGIN
                      ;; 
                      thisIndBin = GetIthBin(indVarIndex, data[ti,0,0,*], NIndVarBins, subset)
                      IF(wrtXlog[0] EQ 1 ) THEN binCenters[indVarIndex] = 10.0^average(alog10(data[ti,0,0,subset])) ELSE binCenters[indVarIndex] = average(data[ti,0,0,subset])
                      ;; The following is the list of
                      dependentData = sort((data[ti,0,k,subset])[thisIndBin])
                      FOR percentileIndex=0, n_elements(percentiles)-1 DO BEGIN
                          theCurves[indVarIndex, percentileIndex, aColor] = depndentData[fix(percentiles[percentileIndex]*n_elements(depndentData))]
                      ENDFOR ; end loop over percentiles
                  ENDFOR ; end loop over independent variable bins.
                  FOR percentileIndex=0,n_elements(percentiles)-1 DO OPLOT, binCenters, theCurves[*,percentileIndex,aColor], COLOR=aColor
              ENDIF ; end check that there are models w/ this color
          ENDFOR ; end loop over color
      ENDIF


      ;; We've created a single frame. What do we do with it?
      IF(sv EQ 1) THEN MPEG_PUT,mpeg_id,WINDOW=1,FRAME=count,/ORDER
      IF(sv EQ 2) THEN WRITE_PNG, dn+'/frame_'+string(count,FORMAT="(I4.4)") + '.png',TVRD(TRUE=1)
      IF(sv EQ 3 || sv EQ 4) THEN dummy=cgSnapshot(filename=(dn+'/frame_'+STRING(count-1,FORMAT="(I4.4)")),/png,/true,/nodialog) 
      IF(sv EQ 0) THEN wait,.05
      IF(sv EQ 5 AND tii NE n_elements(whichFrames)-1) THEN multiplot
    ENDFOR ;; end loop over time indices
    IF(sv EQ 1) THEN MPEG_SAVE,mpeg_id
    IF(sv EQ 1) THEN MPEG_CLOSE,mpeg_id
    ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN spawn,("rm " + fn+".mpg") 
    ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN spawn,("ffmpeg -f image2 -qscale 4 -i "+dn+"/frame_%04d.png " + fn + ".mpg")
    ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN IF(saveFrames EQ 0) THEN spawn,("rm -rf "+dn)
    IF(sv EQ 5) THEN BEGIN
	figureClean,(fn+"_timeSeries"),2
        latexify,(fn+"_timeSeries.eps"),[axisLabels[0],axisLabels[k]],[texLabels[0],texLabels[k]],[.6,.6],height=8.89*(2-horizontal),width=8.89*(1+horizontal)
        multiplot,/default
    ENDIF
    set_plot,'x'
  ENDFOR ;; end loop over dependent variables
  IF(sv EQ 2|| sv EQ 3 || sv EQ 4) THEN spawn,"python makeGidgetMovies.py"
END



;; Given an array of mdot's vs. other variables, compute the cross-correlations between mdot and the others
;; and also the autocorrelation of each variable.
PRO ComputeCorrelations,vsmdot,colors,time,labels,names,name,sv=sv,nt0=nt0,thicknesses=thicknesses,logarithms=logarithms
    IF(n_elements(nt0) EQ 0) THEN nt0 = 1
    IF(n_elements(sv) EQ 0 ) THEN sv=1
    IF(n_elements(thicknesses) EQ 0) THEN thicknesses=intarr(n_elements(vsMdot[0,0,0,*]))+1 ; set each model's thickness to 1
    IF(n_elements(logarithms) EQ 0) THEN logarithms=intarr(n_elements(vsMdot[0,0,*,0])) ; set each variable to be not logarithmic

    ;; Construct some cross-correlations
    nt = n_elements(vsMdot[*,0,0,0])
    ;;; (# of delay times to test) x (# of cross-correlations) x (# of models)
    crossCorrelations = dblarr(nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    autoCorrelations = dblarr(nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    ;;; loop over models
    FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
        ;; loop over variables
        FOR k=0, n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
            laggingVar = vsMdot[nt0:(nt-1),0,k,j]
            IF(logarithms[k] EQ 1) THEN laggingVar = alog10(laggingVar)
            mdot = vsMdot[nt0:(nt-1),0,0,j]
            IF(logarithms[0] EQ 1) THEN mdot = alog10(mdot)
            crossCorrelations[*,k,j] = c_correlate(vsmdot[nt0:(nt-1),0,0,j],vsmdot[nt0:(nt-1),0,k,j],indgen(nt-1-nt0))
            autoCorrelations[*,k,j] = c_correlate(vsmdot[nt0:(nt-1),0,k,j],vsmdot[nt0:(nt-1),0,k,j],indgen(nt-1-nt0))
        ENDFOR

    ENDFOR

    ;; loop over variables
    FOR k=0,n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
        setct,1,0,2
        FIGUREINIT,(name+"_cc_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,CHARSIZE=1.4,THICK=1,CHARTHICK=1,XTITLE="Lag Time (Ga)",YTITLE="x-corr: Accr with "+labels[k],XLOG=1
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
            OPLOT, time[1:(nt-1-nt0)],crossCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j]
        ENDFOR
        FIGURECLEAN,(name+"_cc_"+names[k]),sv

        FIGUREINIT,(name+"_ac_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,CHARSIZE=1.4,THICK=1,CHARTHICK=1,XTITLE="Lag Time (Ga)",YTITLE="autocorr: "+labels[k],XLOG=1
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
             OPLOT, time[1:(nt-1-nt0)],autoCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j]
        ENDFOR
        FIGURECLEAN,(name+"_ac_"+names[k]),sv
    ENDFOR

END







;; given an array whose first column indexes time, second column is irrelevant, third column
;;  indexes different data sets fourth column indexes models, plus corresponding time
;; axis labels, colors, and styles, make a series of simple plots
PRO simplePlot,data,time,labels,names,colors,name,sv,styles=styles
  IF(n_elements(styles) EQ 0) THEN styles = intarr(n_elements(data[0,0,0,*]))
  FOR k=0,n_elements(labels)-1 DO BEGIN;; loop over y-axis variables to be plotted against time
    set_plot,'x'
    figureInit,name+"_"+names[k],sv,2,2
    
    PLOT,[MAX(time)],[MAX(data[*,0,k,*])],COLOR=0, BACKGROUND=255,THICK=2,XTITLE="time",YTITLE=labels[k],XSTYLE=1,YSTYLE=1
    FOR j=0,n_elements(data[0,0,0,*])-1 DO BEGIN
      ;PLOT,time,data[*,0,k,j],COLOR=colors[0,j],BACKGROUND=255,linestyle=styles[j],XTITLE="time",YTITLE=labels[k],XSTYLE=1,YSTYLE=1
      OPLOT,time,data[*,0,k,j],COLOR=colors[0,j],linestyle=styles[j]
    ENDFOR

    figureClean,name+"_"+names[k]
  ENDFOR
END

FUNCTION getElements,nvardim,i
  telements=(unpack(i,nVarDim))
  srtd=nVarDim[reverse(sort(nvardim))]
  elements=telements[reverse(sort(nvardim))]
  IF(n_elements(elements) LT 3) THEN BEGIN
    qelements=intarr(3)
    qelements[0:n_elements(elements)-1]=elements[0:n_elements(elements)-1]
    elements=qelements
  ENDIF

  RETURN,elements
END






;;; Meant to do an analysis of many runs at once, particularly when those 
;;; correspond to different accr. histories
;;; expName is an experiment name, like one would specify in exper.py. This 
;;;   code assumes that you are 
;;;     1) in the gidget/analysis directory, and 2) you have a bunch of 
;;;     runs in a folder with names of the form
;;;      expName/expName<identifiers>_<filenames>, 
;;;     where identifiers can be anything as long as it is unique to each run
;;;     and <filenames> encompasses all the files output by the code, e.g. 
;;;     evolution.dat, radial.dat, etc.
;;; keys is an array which lets you specify which variable to plot from the 
;;;   list in wrtXyt. This is useful for very large experiments, where even 
;;;   storing the modest-size arrays here which take a limited amount of
;;;   data for each run can become cumbersome. If you input an empty array,
;;;   by default the code will analyze all the variables
;;; N controls how many runs in the experiment to analyze. For large number 
;;;   runs, sometimes it's useful to just look at a few runs. If N is larger 
;;;   than the number of runs in the experiment, all the runs are analyzed.
;;; sv- behaves like other instances of sv; 0 to just show the plots on the 
;;;   screen, 2 to save them as png's and subsequently movies. 4 saves 
;;;   movies in a format conducive to showing on a low-resolution projector.
PRO variability3,expNames,keys,N,sv
  compile_opt idl2
  set_plot,'x'
  DEVICE,DECOMPOSED=0
  DEVICE,RETAIN=2

  proftimes=dblarr(5)
  proftimes[0]=systime(1)

  nameList2=ListValidModels(expNames[0],N)
  modelcounter = [0,  n_elements(nameList2) ]
  FOR expInd=1,n_elements(expNames)-1 DO BEGIN
      nameList3 = ListValidModels(expNames[expInd],N)
      nameList2 = [nameList2,nameList3]
      modelcounter = [modelcounter, n_elements(namelist2)]
  ENDFOR

  ctr=0

  ;; variables to plot- title, column index, log plot, name, y-range
  wrtXyt=['R (kpc)','Gas Column Density (Msun/pc^2)','Velocity Dispersion (km/s)','Stellar Column Density (Msun/pc^2)','Stellar Velocity Dispersion (km/s)','Gas Fraction','[Z/Zsun]','SFR Column Density (Msun/yr/kpc^2)','Q','Gas Q','Stellar Q','Molecular Fraction','Age At z=0 (Ga)','Age (Ga)','Depletion Time (Ga)','Viscous Time (Ga)','Depletion/Viscous Time','Mass Flux Through Disk (Msun/yr)']
  wrtXyy=[9,10,11,12,13,17,22,14,12,24,23,48,31,32,35,36,37,39] ;; index
  wrtXyp=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1,0] ;; offset? 
  wrtXyn=['r', 'col','sig','colst','sigst','fg','Z','colsfr','Q','Qg','Qst','fH2','ageAtz0','age','tdep','tvisc','tdepOverTvisc','mdotDisk']
;  wrtXyr=[[0,20],[.1,1000],[5,300],[.1,1000],[5,300],[0,1],[-1,1.0],$
;    [1d-6,10],[1.5,5.0],[.5,50],[.5,50],[0,1],[1d9,14d9],[1d5,14d9],[1.0e8,1.0e10],[1.0e8,1.0e10],[.01,100]]
  wrtXyl=[0,1,1,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1] ;; log plot 18

  IF(n_elements(keys) EQ 0) THEN keys=indgen(n_elements(wrtXyt)-1)+1

  dummy = where( keys EQ 19, n1)
  IF(n1 EQ 2) THEN n1=1

  ;; only plot the variables present in 'keys'	
  tk = [[0],keys]
  wrtXyt=wrtXyt[tk]
  wrtXyy=wrtXyy[tk]
  wrtXyp=wrtXyp[tk]
  wrtXyn=wrtXyn[tk]
  ;wrtXyr=wrtXyr[*,tk]
  wrtXyl=wrtXyl[tk]  

  ;; now loop over every model in the list - 
  ;; all we're doing here is reading in each model then freeing 
  ;; the associated memory
  ;; so that we can iterate over a very large number. For each model, 
  ;; record the data specified in the arrays above (e.g. yy)
  modelPtr=ReadModels([nameList2[0]])
  model = (*(modelPtr[0]))
  theNTS = model.NOutputsNominal+2

  time = model.evArray[1,*]  
  z = model.evArray[10-1,*]
  

  ;; (# timesteps) x (nx) x (# of columns) x (# models)
  theData= dblarr(theNTS,n_elements(model.dataCube[0,*,0]),n_elements(wrtXyy),n_elements(nameList2))

  NVS = 28
  ;; (# timesteps) x (# pts/model frame = 1) x (one thing vs another = 20) x (# models)
  vsMdot = dblarr(theNTS,1,NVS,n_elements(nameList2))
  nameList4 = strarr(n_elements(nameList2))


  proftimes[1]=systime(1)

;  FOR i=0,n_elements(nameList2)-1 DO BEGIN ;; loop over every valid model
  ctr=0
  modelcounter2 = modelcounter*0
  FOR expInd=1, n_elements(expNames) DO BEGIN
    low = modelcounter[expInd-1]
    high= modelcounter[expInd]-1
    ctr2=0
    FOR i=low,high DO BEGIN
 
       IF(i MOD 50 EQ 0) THEN PRINT,"Processing model #",i," of ",n_elements(nameList2)

       ;; the number of colors necessary is srtd[0], 
       ;; which color set to use is elements[2]
       IF(i NE 0) THEN modelPtr=ReadModels([nameList2[i]])
       IF(i NE 0) THEN model = (*(modelPtr[0]))
       zi=nearest("redshift",1,model,0.0)
       offsets=[0,model.ncolstep,model.ncolstep+model.npostprocess]
       IF(n_elements(model.evArray[1,*]) EQ theNTS) THEN BEGIN
         time = model.evArray[1,*]  
         z = model.evArray[10-1,*]
       ENDIF
       x = model.dataCube[0,*,0]

       ;; if this model has the appropriate number of time steps..
       ;;PRINT,"Number of time steps",n_elements(model.dataCube[*,0,0]),n_elements(theData[*,0,0,0])
       ;IF(n_elements(model.dataCube[*,0,0]) EQ n_elements(theData[*,0,0,0])) THEN BEGIN
       IF(n_elements(model.dataCube[*,0,0]) EQ theNTS) THEN BEGIN
         FOR k=0,n_elements(wrtXyy)-1 DO BEGIN	
           theData[*,*,k,ctr] = model.dataCube[*,*,wrtXyy[k]+offsets[wrtXyp[k]]-1]
         ENDFOR

;         vsMdot[*,0,0,ctr] =  model.dataCube[*,model.nx-1,39-1]  ;; external accretion rate
         
         vsMdot[*,0,0,ctr] = model.evArray[20-1,*] ;; Cosmological accretion rate.
         vsMdot[*,0,1,ctr] = model.evArray[10,*] ;; SFR [in the disk!]

         FOR zi=0, n_elements(model.evArray[9,*])-1 DO BEGIN ; loop over redshift
           modelInfo = diskStats(model,z=model.evArray[9,zi])
           vsMdot[zi,0,3,ctr] = modelInfo[1] ;; peak radius of column density (kpc)
           vsMdot[zi,0,4,ctr] = modelInfo[2] ;; edge of SF region (kpc)
           vsMdot[zi,0,5,ctr] = modelInfo[3] ;; col @ r=8kpc (solar masses/pc^2)
           vsMdot[zi,0,6,ctr] = modelInfo[7] ;; sSFR (yr^-1)
           vsMdot[zi,0,7,ctr] = modelInfo[8] ;; gas mass in bulge (MSol)
           vsMdot[zi,0,8,ctr] = modelInfo[9] ;; stellar mass in bulge
           vsMdot[zi,0,9,ctr] = modelInfo[10] ;; total f_H2
           vsMdot[zi,0,10,ctr] = modelInfo[11] ;; mdot bulge gas
           vsMdot[zi,0,11,ctr] = modelInfo[12] ;; mdot bulge stars
           vsMdot[zi,0,2,ctr] = modelInfo[19] ;; stellar mass in disk + bulge
           vsMdot[zi,0,12,ctr] = modelInfo[8]/modelInfo[19] ;; B:T
           vsMdot[zi,0,13,ctr] = modelInfo[18] ;; SFR + mdotbulgeGas
           vsMdot[zi,0,14,ctr] = modelInfo[17] ;; efficiency
           vsMdot[zi,0,15,ctr] = alog10(modelInfo[16]) ;; metallicity
           vsMdot[zi,0,16,ctr] = modelInfo[5] ;; fG in SF region.
           vsMdot[zi,0,17,ctr] = model.evArray[20-1,zi] ;; (modelInfo[18] - modelInfo[11])/(modelInfo[18]+1.0d-10) ;; fraction of SF in the bulge
           vsMdot[zi,0,18,ctr] = modelInfo[20] ;; fgL
           vsMdot[zi,0,19,ctr] = alog10(modelInfo[21]) ;; ZL
           vsMdot[zi,0,20,ctr] = modelInfo[15] ;; fg
	   vsMdot[zi,0,21,ctr] = alog10(model.evArray[5-1,zi]/.02) ;; ZBulge
	   vsMdot[zi,0,22,ctr] = modelInfo[23-1];; vrAvg (mass-averaged, km/s)
	   vsMdot[zi,0,23,ctr] = modelInfo[24-1];; GIin (kpc)
	   vsMdot[zi,0,24,ctr] = modelInfo[25-1];; Giout (kpc)
	   vsMdot[zi,0,25,ctr] = modelInfo[26-1];; velocity dispersion (mass-averaged, km/s)
	   vsMdot[zi,0,26,ctr] = modelInfo[27-1];; vrAvgGt0
	   vsMdot[zi,0,27,ctr] = modelInfo[28-1];; tdep = gas mass / SFR
         ENDFOR ;; end loop over redshift

         nameList4[ctr] = nameList2[i]

         ctr2 = ctr2+1
         ctr = ctr+1
       ENDIF ELSE BEGIN
          print,"Warning: discarding model i: ",i," named ",model.name," because it has ",n_elements(model.dataCube[*,0,0])," time steps instead of ",n_elements(theData[*,0,0,0])

       END ;; end check of whether this model has the correct # of time steps.

       IF(i EQ 0) THEN print,model.name,diskStats(model)
       ptr_free,modelPtr[0]
    ENDFOR ;; end loop over specific models (i).
    modelcounter2[expInd] = ctr
  ENDFOR ;; end loop over expNames


;['mdot','DiskSFR','Mst','rPeak','rHI','colsol','sSFR','BulgeGasMass','BulgeStMass','fH2','mdotBulgeGas','mdotBulgeStars',"BT","SFRplusMdot_b_gas","efficiency","Z","f_gInSF","mdot","fgL","ZL","fg","ZBulge","vrgAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg"]

  vsMdotLabels = ["Mdot (Msun/yr)","DiskSFR (Msun/yr)","Stellar Mass (Msun)","Peak Radius of Column Density (kpc)","Radius of HI transition (kpc)","Column Density at r=8 kpc (Msun/pc^2)","Specific SFR (yr^-1)","Bulge Gas Mass (Msun)","Bulge Stellar Mass (Msun)","H_2 Fraction","Mdot Gas into Bulge (Msun/yr)","Mdot Stars into Bulge (Msun/yr)","Bulge to Total Ratio","SFR Including Gas Flux Into Bulge (Msun/yr)","efficiency (Mstar / (f_b M_h))","[Z/Zsun]","Gas Fraction in SF Region","Mdot (Msun/yr)","Gas Fraction in Optical Region","[Z/Zsun] in Optical Region","Gas Fraction","[Z Bulge/Zsun]","Average inward gas radial velocity (km/s)","Inner Edge of GI Region (kpc)","Outer Edge of GI Region (kpc)","Average velocity dispersion (km/s)","Average radial velocity for non-outward velocities only (km/s)","Depletion Time (Ga)"]
  vsMdotTexLabels = ["$\dot{M} (M_\odot/yr)$","Disk SFR ($M_\odot$/yr)","Stellar Mass ($M_\odot$)","Peak Radius (kpc)","Radius of HI trans. (kpc)","$\Sigma$ at $r=8$ kpc ($M_\odot/pc^2$)","sSFR ($yr^{-1}$)","Bulge Gas Mass ($M_\odot$)","Bulge Stellar Mass ($M_\odot$)","$H_2$ Fraction","$\dot{M}_g$ into Bulge ($M_\odot$/yr)","$\dot{M}_*$ into Bulge ($M_\odot$/yr)","Bulge:Total Ratio","SFR + $\dot{M}_{g,\rightarrow\mathrm{bulge}}$ ($M_\odot$/yr)","$\epsilon$ ($M_*$ / ($f_b M_h$))","$[Z/Z_\odot]$","$f_g$ in SF Region","$\dot{M}$ ($M_\odot$/yr)","$f_g$ in Optical Region","$[Z/Z_\odot]$ in Optical Region","$f_g$","$[Z_\mathrm{bulge}/Z_\odot]$","Average $-v_r$ (km/s)","Innermost GI Region (kpc)","Outermost GI Region (kpc)","$\langle\sigma\rangle$ (km/s)","Average $-v_r\ge 0$ (km/s)","$t_{dep} = M_g f_{H_2} / \dot{M}^{SF}$ (Ga)"]

  vsMdotNames  = ["mdot","DiskSFR","Mst","rPeak","rHI","colsol","sSFR","BulgeGasMass","BulgeStMass","fH2","mdotBulgeGas","mdotBulgeStars","BT","SFRplusMdotIn","efficiency","Z","fgInSF","mdot","fgL","ZL","fg","ZBulge","vrAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg"]
  vsMdotStrt =  [1,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  vsMdotToLog = [1,1,1,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1] 

  modelcounter = modelcounter2[*]
  nameList2 = nameList4[*]

  proftimes[2]=systime(1)

  theData = (temporary(theData))[*,*,*,0:ctr-1]
  colors = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))

  ;; having accumulated theData, let's do some analysis!
  sortdData = theData * 0.0
  sortdData[*,*,0,*] = theData[*,*,0,*] ;; don't sort radii..
  ;sortdVsMdot = vsMdot*0.0
  FOR expInd=1, n_elements(expNames) DO BEGIN
    low = modelcounter[expInd-1]
    high= modelcounter[expInd]-1
    FOR j=0,n_elements(theData[*,0,0,0])-1 DO BEGIN ;; loop over time
      FOR k=1,n_elements(wrtXyy)-1 DO BEGIN ;; loop over variables
        FOR xx=0,n_elements(theData[0,*,0,0])-1 DO BEGIN ;; loop over x-coord
	  ;;; Sort the models, within each experiment
	  srtd = sort(theData[j,xx,k,low:high])
          sortdData[j,xx,k,low:high] = theData[j,xx,k,low+srtd]
          IF(xx EQ 4 AND j EQ n_elements(theData[*,0,0,0])-1) THEN BEGIN
;            print, "At n_x=4, the models in order from low to high values of variable ", wrtXyt[k] ," are ",nameList2[srtd]
          ENDIF
	ENDFOR ;; end loop over x-coord 
      ENDFOR ;; end loop over variables
;    ENDFOR ;; end loop over time
      ;st = sort(vsMdot[0,0,0,low:high])
      ;FOR k2=0,n_elements(vsMdot[0,0,*,0])-1 DO BEGIN	;; loop over variables
      ;  sortdVsMdot[j,0,k2,low:high] = vsMdot[j,0,k2,low+st]
      ;ENDFOR ;; end loop over variables
    ENDFOR ;; end loop over time
  ENDFOR ;; end loop over experiments



  proftimes[3]=systime(1)

  intervals = dblarr(n_elements(theData[*,0,0,0]),  n_elements(theData[0,*,0,0]), n_elements(theData[0,0,*,0]), 3*n_elements(expNames))
  FOR expInd=1, n_elements(expNames) DO BEGIN
    low=modelcounter[expInd-1]
    high=modelcounter[expInd]-1
    delta = high-low
    intervals[*,*,*, (expInd-1)*3+0]=sortdData[*,*,*, fix(low+.025*delta)] ;; fix(low+.025*(high-low))
;    intervals[*,*,*, (expInd-1)*5+1]=sortdData[*,*,*, fix(low+.16*delta)]
    intervals[*,*,*, (expInd-1)*3+1]=sortdData[*,*,*, fix(low+.5*delta)]
;    intervals[*,*,*, (expInd-1)*5+3]=sortdData[*,*,*, fix(low+.84*delta)]
    intervals[*,*,*, (expInd-1)*3+2]=sortdData[*,*,*, fix(low+.975*delta)]
    colors[*,low:high] = expInd-1
  ENDFOR


  nmodels = n_elements(vsMdot[0,0,0,*])
  ncolors = 5
  ;; sort according to z=0 sSFR
  IF(n_elements(expNames) EQ 1) THEN BEGIN
    zi0 = n_elements(vsMdot[*,0,0,0])-1

    FOR zi = 0, zi0 DO BEGIN
      st = sort(vsMdot[zi0,0,6,*])

      FOR nc=0, ncolors-1 DO BEGIN
        wh = WHERE (st GE nc*(float(nmodels)/float(ncolors)) and st LT (nc+1)*(float(nmodels)/float(ncolors)))
        ll = fix(float(nc)*float(nmodels)/float(ncolors))
        hh = fix(float(nc+1)*(float(nmodels)/float(ncolors)))-1
        IF(hh LT 0) THEN hh=0
        colors[zi,st[ll:hh]] = nc
        ;STOP
      ENDFOR
    ENDFOR
  ENDIF

  thicknesses=intarr(nmodels)+1
  st= sort(vsMdot[theNTS-1,0,2,*]/vsMdot[theNTS-1,0,14,*]) ;; sort by halo mass at redshift 0.
  FOR nc=0, ncolors-1 DO BEGIN
      wh = WHERE(st GE nc*(float(nmodels)/float(ncolors)) and st LT (nc+1)*(float(nmodels)/float(ncolors)))
      ll = fix(float(nc)*float(nmodels)/float(ncolors))
      hh = fix(float(nc+1)*(float(nmodels)/float(ncolors)))-1
      IF(hh LT 0) THEN hh=0
      thicknesses[st[ll:hh]] = 1.0 + .5*float(nc) ;; thicker lines have higher Mh0.
  ENDFOR






  expName2 = expNames[0]
  IF(n_elements(expNames) LT 5) THEN BEGIN
      FOR i=1, N_elements(expNames)-1 DO BEGIN
          expName2 = expName2 + '_' + expNames[i];; +'_'+strcompress(string(MIN([n_elements(nameList2),N])),/remove) ;; e.g. rk5_113
      ENDFOR
  ENDIF ELSE expName2 +='_' + strcompress(string(n_elements(expNames)),/remove)
  expName2 = expName2 + '_' + strcompress(string(MIN([n_elements(nameList2),N*n_elements(expNames)])),/remove) ;; e.g. rk5_rk6_113



  IF(n_elements(nameList2) GT 50) THEN BEGIN
      setct,1,0,0
      linestyles = intarr(n_elements(expNames)*3)+2
      midpoints = where(indgen(n_elements(expNames)*3) MOD 3 EQ 1)
      linestyles[midpoints] = 0
      intervalsColors0 = indgen(n_elements(expNames)*3)/3
      intervalsColors = intarr(n_elements(time),n_elements(intervalsColors0))
      FOR tt=0,n_elements(time)-1 DO intervalsColors[tt,*] = intervalsColors0 
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervals",sv,intarr(n_elements(wrtXyl)),axislabels=wrtXyt
  ENDIF



  setct,3,n_elements(theData[0,0,0,*]),0
;  IF(n_elements(nameList2) LT 105) THEN BEGIN 	
    simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsorted",sv,intarr(n_elements(wrtXyl)),axislabels=wrtXyt
;  ENDIF
	

  setct,3,n_elements(vsMdot[0,0,0,*]),0
  IF(n1 EQ 1) THEN BEGIN
;    simpleMovie, vsMdot,z,vsMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsmdot",sv,vsMdotStrt,PSYM=1,prev=1,axisLabels=vsMdotLabels
  ENDIF

  ;; average mdot over a number of time steps nt
  nt = 40
  setct,3,n_elements(vsMdot[0,0,0,*]),0
  avgdVsMdot = vsMdot[*,*,*,*]

  FOR ti=0,n_elements(time)-1 DO BEGIN
    FOR mm=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
      avgdVsMdot[ti,0,0,mm] = TOTAL(vsMdot[MAX([0,ti-nt]):ti,0,0,mm]) / (ti - MAX([0,ti-nt]) + 1)
    ENDFOR
  ENDFOR
  IF(n1 EQ 1) THEN BEGIN
     vsAvgMdotLabels = vsMdotLabels[*]
     vsAvgMdotLabels[0] = "Average Mdot Over Past 40 Outputs (Msun/yr)"
     vsAvgMdotTexLabels = vsMdotTexLabels[*]
     vsAvgMdotTexLabels[0] = "$\dot{M}_{40 \mathrm{outputs}} (M_\odot/yr)$"
     vsAvgMdotNames = vsMdotNames[*]
     vsAvgMdotNames = "avgMdot"
     toLog = [1,1,0,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,1,1]
     strt =  [1,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
     simpleMovie, avgdVsMdot,z,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot",sv,vsMdotStrt,PSYM=1,prev=1,axisLabels=vsAvgMdotLabels
  ENDIF


  ;; Plot these quantities vs. the star formation rate.
  setct,3,n_elements(vsMdot[0,0,0,*]),0
  vsSFR = vsMdot[*,*,1:n_elements(vsMdot[0,0,*,0])-1,*] ;; cut out the mdot slice of the array
  IF(n1 EQ 1) THEN BEGIN 

    vsSFRLabels = vsMdotLabels[1:(n_elements(vsMdotLabels)-1)]
    vsSFRTexLabels = vsMdotTexLabels[1:n_elements(vsMdotTexLabels)-1]
    vsSFRNames = vsMdotNames[1:(n_elements(vsMdotNames)-1)]
    vsSFRToLog = vsMdotToLog[1:(n_elements(vsMdotToLog)-1)]
    vsSFRstrt= vsMdotstrt[1:(n_elements(vsMdotStrt)-1)]

;    simpleMovie, vsSFR, z,vsSFRNames, colors, colors*0, vsSFRtoLog,expName2+"_vsSFR",sv,vsSFRstrt,PSYM=1,prev=1,axisLabels=vsSFRLabels
  ENDIF


  ;; Now do a series of plots vs. M*
  vsMstar = vsMdot[*,*,*,*]
  vsMstar[*,*,0,*] = vsMstar[*,*,2,*] ;replace Mdot with M* as the zeroeth column of the data.
  vsMstar[*,*,2:(NVS-2),*] =vsMstar[*,*,3:(NVS-1),*] ;move over the other columns so that we don't produce an M* vs M* plot
  vsMstar=vsMstar[*,*,0:(NVS-2),*] ; eliminate the extraneous last column.

  vsMstarLabels = vsMdotLabels[*]
  vsMstarLabels[0] = vsMstarLabels[2]
  vsMstarLabels[2:(NVS-2)] = vsMstarLabels[3:(NVS-1)]
  vsMstarLabels=vsMstarLabels[0:(NVS-2)]

  vsMstarTexLabels = vsMdotTexLabels[*]
  vsMstarTexLabels[0] = vsMstarTexLabels[2]
  vsMstarTexLabels[2:(NVS-2)] = vsMstarTexLabels[3:(NVS-1)]
  vsMstarTexLabels=vsMstarTexLabels[0:(NVS-2)]

  vsMstarNames = vsMdotNames[*]
  vsMstarNames[0] = vsMstarNames[2]
  vsMstarNames[2:(NVS-2)] = vsMstarNames[3:(NVS-1)]
  vsMstarNames=vsMstarNames[0:(NVS-2)]

  vsMstarStrt = vsMdotStrt[*]
  vsMstarStrt[0] = vsMstarStrt[2]
  vsMstarStrt[2:(NVS-2)] = vsMstarStrt[3:(NVS-1)]
  vsMstarStrt=vsMstarStrt[0:(NVS-2)]

  vsMstarToLog = vsMdotToLog[*]
  vsMstarToLog[0] = vsMstarToLog[2]
  vsMstarToLog[2:(NVS-2)] = vsMstarToLog[3:(NVS-1)]
  vsMstarToLog=vsMstarToLog[0:(NVS-2)]



  setct,3,n_elements(vsMstar[0,0,0,*]),0
  IF(n1 EQ 1) THEN BEGIN
    theLabels =  ["Mst","DiskSFR","rPeak","rHI",'colsol','sSFR','BulgeGasMass','BulgeStMass','fH2','mdotBulgeGas','mdotBulgeStars',"BT","SFRplusMdot_b_gas","efficiency","Z","f_gInSF","mdot","fgL","ZL","fg","ZBulge","vrgAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg"]
    toLog = [1,1,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1]
    strt =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    simpleMovie, vsMstar, z,vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstar",sv,vsMstarStrt,PSYM=1,prev=1,saveFrames=0,axisLabels=vsMstarLabels
  ENDIF



  setct,3,n_elements(vsMdot[0,0,0,*]),0
  vsTime = vsMdot[*,*,*,*]
  vsTimeLabels = vsMdotLabels[*]
  vsTimeLabels[0] = "Time Since z=2 (Ga)"
  vsTimeTexLabels = vsMdotTexLabels[*]
  vsTimeTexLabels[0] = "Time since $z=2$ (Ga)"
  vsTimeNames = vsMdotNames[*]
  vsTimeToLog = vsMdotToLog[*]
  vsTimeToLog[0] = 0
  vsTimeStrt =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]

  ;; for each model, set the independent variable equal to the time: 
  FOR j=0, n_elements(vsTime[0,0,0,*]) -1 DO BEGIN
    vsTime[*,0,0,j] = time[*]
  ENDFOR
  IF(n1 EQ 1) THEN BEGIN
    simpleMovie, vsTime,z,vsTimeNames,colors,colors*0,vsTimeToLog,expName2+"_vstime",sv,vsTimeStrt,PSYM=3,prev=1,taillength=1000,saveFrames=0,axisLabels=vsTimeLabels
  ENDIF


  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2,sv=4,nt0=20,thicknesses=thicknesses

  ;; Now make some standalone plots.
  
  simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstar",5,vsMstarSTrt,PSYM=1,prev=1,saveFrames=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,20,51,104,201]
  

  Mh = vsMdot[*,0,2,*]/(.18*vsMdot[*,0,14,*])
  zn = n_elements(Mh[*,0,0,0])
	
  ;; Plot the accretion rates.
  ;SETCT,3,n_elements(vsMdot[0,0,0,*]),0
  FIGUREInit,(expName2+'_accRates'),1,3,3
  yr=[MIN(vsMdot[zn-1,0,0,*])*.5,MAX(vsMdot[0,0,0,*])*2];[.01,1000];[MIN(sortdVsMdot[*,0,0,*]),MAX(sortdVsmdot[*,0,0,*])]
  Plot,[time[0]],[vsMdot[0,0,0,0]],COLOR=0,BACKGROUND=255,XRANGE=[MIN(time),MAX(time)],YRANGE=yr,YLOG=1,XSTYLE=1,XTITLE="Time",YTITLE="Accretion Rate",CHARSIZE=2,THICK=2,CHARTHICK=1.5
  FOR m=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
    ;setct,5,n_elements(vsMdot[0,0,0,*]),m
    OPLOT,time,vsMdot[*,0,0,m],COLOR=colors[0,m],THICK=thicknesses[m]
	
  ENDFOR
  FIGURECLEAN,(expName2+'_accRates'),1



  ;; Plot the halo mass
  ;SETCT,3,n_elements(vsMdot[0,0,0,*]),0
  FIGUREInit,(expName2+'_Mh'),1,3,3
  ;Plot,time,Mh[*,0,0,0],COLOR=0,BACKGROUND=255,YRANGE=[Min(Mh[zn-1,0,0,*])*.1,Max(Mh[zn-1,0,0,*])*2],YLOG=1,XSTYLE=1,XTITLE="Time",YTITLE="Halo Mass",CHARSIZE=3,THICK=thicknesses[0],CHARTHICK=2
  PLOT,[time[0]],[Mh[0,0,0,0]],COLOR=0,BACKGROUND=255,XRANGE=[MIN(time),MAX(time)],YRANGE=[MIN(Mh[zn-1,0,0,*])*.1,MAX(Mh[zn-1,0,0,*])*2],YLOG=1,XSTYLE=1,YSTYLE=1,XTITLE="Time (Ga)",YTITLE="Halo Mass",CHARSIZE=2,THICK=2,CHARTHICK=1.5
  FOR m=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
    ;setct,5,n_elements(vsMdot[0,0,0,*]),m
    OPLOT,time,Mh[*,0,0,m],COLOR=colors[0,m],THICK=thicknesses[m]
	
  ENDFOR
  FIGURECLEAN,(expName2+'_Mh'),1



;  FIGUREInit,(expName2+'_z0AccRates'),1,1,1
;  cgHistoplot, alog10(vsMdot[n_elements(vsMdot[*,0,0,0])-1,0,0,*]), BinSize=.1
;  FIGURECLEAN,(expName2+'_z0AccRates'),1
	
  proftimes[4]=systime(1)

  PRINT,proftimes[1]-proftimes[0]," seconds to check model validities and initialize"
  PRINT,proftimes[2]-proftimes[1]," seconds to read in data"
  PRINT,proftimes[3]-proftimes[2]," seconds to sort the data"
  PRINT,proftimes[4]-proftimes[3]," seconds to generate movies"

  ;STOP

END
