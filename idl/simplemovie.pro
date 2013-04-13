;; given an array "data" whose first column indexes time, second column indexes position,
;; fourth column indexes models, and a vector of time labels, and axis labels
;; make one or a few movies! The third column indexes which variable we're plotting.
;; 
;; time is a 1D array assumed to contain the "time", or redshift values 
;;    for the corresponding cubes in data, indexed by the 1st column.
;; labels is a 1D array of strings, one for each of the variables indexed by the 3rd column in DATA,
;;    which specify part of the filename produced by simpleMovie.
;; colors is a (# timesteps) by (# models) array.
;; styles is a 1D array with (# models) elements specifying the style for each line.
;; wrtXlog specifies which variables (column 3 of data) should be plotted in log space
;; name is the base name for all files/movies produced in this procedure
;; sv specifies how the movie will be produced:
;;    1: IDL's built-in MPG creator (NOT RECOMMENDED)
;;    2: print to screen, save PNG, compile PNG's w/ ffmpeg (not recommended if you're producing a lot of movies)
;;    3: print to z buffer, save PNG, compile w/ ffmpeg (encouraged)
;;    4: same as 3, but try to adjust parameters so that the video is playable on a low-res projector
;;;;; KEYWORD ARGUMENTS:
;; prev: at a given time step, plot the data from both this time step and previous ones (# of previous timesteps = taillength)
;; psym: just the regular IDL keyword PSYM (see e.g. the documentation for PLOT)
;; axislabels: a vector of strings for labelling the axes appropriate for each variable
;;    indexed by the 3rd column in DATA- defaults to labels if unspecified
;; taillength (see "prev" above)
;; plotContours - not yet implemented: at each time step plot contours instead of lots of points.
;; whichFrames - specify a subset of all time indices to use for movies/plots, e.g. [2,17,30]
;; texLabels - LaTeX-ified version of axislabels. Only used if svSinglePlot == 2 (meaning save plots as .eps)
;; horizontal - for multi-paned plots, arrange them horizontally, as opposed to vertically
;; NIndVarBins - when plotting many galaxies, this variable allows you to divide the galaxies into this many
;;   bins according to their value of the independent variable, such that equal numbers will appear in each bin.
;; percentileList - after dividing galaxies up into bins by their independent value, within each bin find 
;;   the value of each dependent variable at each percentile specified here. For instance the default plots 
;;   2-sigma contours+the median
;; svSinglePlot - the value of sv to use if we're producing a single plot. See figureinit.pro; typically 4 is a good choice
;; strt: for each variable, specify whether we should draw a straight line along x=y
;; yranges: Set the y range for every variable to be the same range, e.g. yranges=[0,1].
;; ranges: manually specify the data ranges. Useful if the default (simpleranges) yields too much whitespace.
;; replacementText: a list the same length as whichFrames with a label for each frame which will
;;          be placed on the plot as an annotation.
;; fill: fill between each line. I don't use this anymore, as the functionality is a bit specialized
;;         (see balanceplots.pro )
;; 

PRO simpleMovie,data,labels,colors,styles,wrtXlog,name,sv, $
    prev=prev,psym=psym,axislabels=axislabels,taillength=taillength, $
    plotContours=plotContours,whichFrames=whichFrames,texLabels=texLabels, $
    horizontal=horizontal,NIndVarBins=NIndVarBins, $
    percentileList=percentileList,svSinglePlot=svSinglePlot,strt=strt, $
    thicknesses=thicknesses,yranges=yranges,ranges=ranges, $
    replacementText=replacementText,fill=fill,includeRandom=includeRandom, $
    additionalLines=additionalLines,additionalLineColors=additionalLineColors, $
    additionalLineLabels=additionalLineLabels,leg=leg, grouping=grouping
  IF(sv EQ 3 || sv EQ 4 || sv EQ 5) THEN cg=1 ELSE cg=0

  PRINT, "Starting to make the files of the form: ", name,"*"

  ;; exclude the first time step from setting the range.
;  IF(n_elements(time) NE n_elements(data[*,0,0,0])) THEN message,"unexpected # of entries for the time vector in simplemovie"
  nTimeStepsInput = n_elements(data[*,0,0,0])
  IF(n_elements(ranges) EQ 0) THEN  $
      IF(n_elements(data[*,0,0,0]) GT 3 ) THEN $
        ranges= simpleranges(data[2:nTimeStepsInput-1,*,*,*],wrtxlog) $
      ELSE $
        ranges=simpleranges(data[0:nTimeStepsInput-1,*,*,*],wrtxlog)
  IF(n_elements(yranges) EQ 2) THEN BEGIN
    ranges[0,1:n_elements(ranges[0,*])-1] = yranges[0]
    ranges[1,1:n_elements(ranges[1,*])-1] = yranges[1]
  ENDIF
  IF(n_elements(horizontal) EQ 0) THEN horizontal = 0
  IF(horizontal GE 2 OR horizontal LT 0) THEN horizontal = 1
  IF(n_elements(axislabels) EQ 0) THEN axislabels=labels[*]
  IF(n_elements(prev) EQ 0) THEN prev=0
  IF(n_elements(taillength) EQ 0) THEN taillength=2
  IF(n_elements(saveFrames) EQ 0) THEN saveFrames=0
  IF(n_elements(whichFrames) EQ 0) THEN whichFrames = indgen(nTimeStepsInput)
  IF(n_elements(texLabels) EQ 0) THEN texLabels = axisLabels[*]
  IF(n_elements(NIndVarBins) EQ 0) THEN NIndVarBins=0
  IF(n_elements(percentileList) EQ 0) THEN percentileList=[.16,.5,.84]
  IF(n_elements(svSinglePlot) EQ 0) THEN svSinglePlot = 4
  IF(n_elements(strt) EQ 0) THEN strt = intarr(n_elements(labels))
  IF(n_elements(thicknesses) EQ 0) THEN thicknesses = intarr(n_elements(data[0,0,0,*]))+2
  IF(n_elements(fill) EQ 0) THEN fill=0
  IF(n_elements(includeRandom) EQ 0) THEN includeRandom = 0
  IF(size(colors,/n_dimensions) EQ 1) THEN BEGIN
      tmp = dblarr(nTimeStepsInput,n_elements(colors))
      FOR ti=0,nTimeStepsInput-1 DO tmp[ti,*]=colors[*]
      colors=tmp[*,*]
  ENDIF
  IF(sv EQ 4) THEN thicknesses= 3+temporary(thicknesses)
  chth = thicknesses[0]
  cs =0.8* thicknesses[0]
  ;; We're making EPS plots in this case. That means that charthick and charsize need to be much smaller
  wf = n_elements(whichFrames)-1  
  IF(sv EQ 5 AND svSinglePlot EQ 2) THEN BEGIN
    chth=.5
    cs = .5
    ;thicknesses = [thicknesses[0]/1.7, thicknesses[1:n_elements(thicknesses)-1]/10.0]
    ;thicknesses=[3,thicknesses[1:n_elements(thicknesses)-1]*0+1]
    thicknesses=thicknesses*0+2.5
    IF(wf EQ 0) THEN BEGIN
        chth=.3
        cs=.3
    ENDIF
  ENDIF


  SETCT,1,n_elements(styles),2

  ngrp = n_elements(grouping)
  IF(ngrp GT 0 and horizontal NE 1) THEN message,"At this point, grouping and !horizontal is not supported"
  IF(ngrp GT 5) THEN message,"Probably too many plots: ",ngrp
  corr =0
  IF(ngrp GT 0) THEN corr=1
  theHeight = 1+(1-horizontal)*wf + ngrp -corr
  theWidth = 1+horizontal*wf
  IF(ngrp GT 3) THEN cs=cs*1.3

  ;; loop over y-axis variables to be plotted (except the first one, which is just the independent var!)
  labelsIT=labels[*]
  IF(ngrp GT 0) THEN labelsIT=[labels[0],labels[grouping]]
  FOR kx=1,n_elements(labelsIT)-1 DO BEGIN    
    k=kx
    IF(ngrp GT 0) THEN k=grouping[kx-1]
    IF(ngrp EQ 0) THEN fn=name+"_"+labels[k] ELSE fn=name+"_"+labelsIT[1]
    dn="movie_"+fn

    IF(ngrp GT 0 and kx EQ 1 or ngrp EQ 0) THEN BEGIN
        set_plot,'x'
        IF(sv EQ 3 || sv EQ 4) THEN set_plot,'z'
        IF(cg NE 1) THEN WINDOW,1,XSIZE=1024,YSIZE=1024,RETAIN=2
        IF(sv EQ 3 || sv EQ 4 ) THEN cgDisplay,1024,1024,color=255
        IF(sv EQ 1) THEN mpeg_id=MPEG_OPEN([1024,1024],FILENAME=(fn+".mpg"))
        IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN BEGIN
          FILE_MKDIR,dn
        ENDIF
        IF(sv EQ 5 ) THEN BEGIN
        	IF (wf GT 6) THEN message,"You have requested an unreasonable number of frames for a single plot."
    	    ;cs = .3
            corr = 0
            figureInit,(fn+"_timeSeries"),svSinglePlot,theWidth,theHeight,sm=1
    ;	    IF(wf NE 0) THEN 
            modmultiplot,[0,theWidth,theHeight,0,0],/square,additionalMargin=0.2;,mxtitle=axisLabels[0],mytitle=axislabels[k],gap=.003,mxtitsize=cs,mytitsize=cs
    ;	    IF(wf NE 0) THEN 
            !p.noerase=0
        ENDIF
    ENDIF

    symsize = cs* 2 / (alog10(n_elements(data[0,0,0,*]))+1)
    IF(NIndVarBins GT 0) THEN symsize = symsize/3.0

    ;; counter of frames.
    count=0 
    ; length of a tail to put on each point which represents a galaxy when prev=1
    ; tailLength = fix(3.0/(float(n_elements(styles))/928.0)^(.25)) 

    ;; loop over time: each iteration of this loop will save a frame in a movie 
    ;; or a panel in a multipanel plot
    FOR tii=0,n_elements(whichFrames)-1 DO BEGIN  
      IF(sv EQ 5 and tii NE 0 and wf NE 0) THEN !p.noerase=1
      ti = whichFrames[tii]
      count=count+1

      ;; If this variable k is included in strt, draw a dashed line y=x.
      ;; Either way, this long statement draws the first plot for this frame.
      ;; Everything else will be overplotted on top in a few lines. 
      xt = axisLabels[0]
      IF(sv EQ 5 and horizontal EQ 1 and (ngrp GT 0 and kx NE ngrp)) THEN xt=""
      yt = axisLabels[k]
      IF(sv EQ 5 AND tii NE 0 AND horizontal EQ 1 and wf NE 0) THEN yt=""
      ;IF(sv EQ 5) THEN yt=yt+" at z="+str(time)


      ;;; SET UP THE PLOTTING SPACE
      IF(strt[k] EQ 0) THEN $
        PLOT,[0],[0],COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=xt, $
         YTITLE=yt, XRANGE=ranges[*,0],YRANGE=ranges[*,k],ylog=wrtXlog[k], $
         xlog=wrtXlog[0],CHARTHICK=chth,CHARSIZE=cs,THICK=thicknesses[0], $
         XTHICK=thicknesses[0],YTHICK=thicknesses[0],/nodata $
       ELSE $
        PLOT,findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
	 findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
         COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=xt, $
         YTITLE=yt,XRANGE=ranges[*,0],YRANGE=ranges[*,k],ylog=wrtXlog[k], $
         xlog=wrtXlog[0],linestyle=2,CHARSIZE=cs,CHARTHICK=chth, $
         THICK=thicknesses[0],XTHICK=thicknesses[0],YTHICK=thicknesses[0]
      ;; that was all one line! It's now over, and our plotting space is set up.

      ;; Print the time in the upper left of the screen
      
      IF(N_ELEMENTS(replacementText) EQ 0) THEN theTimeText = "";timeText+string(abs(time[ti]),Format='(D0.2)')
      IF(N_ELEMENTS(replacementText) NE 0) THEN theTimeText = replacementText[tii]

      ; Having set up the plotting space, loop over models and plot each one.



      ;; FILL IN REGIONS BETWEEN LINES - most plots don't do this, and typically the input
      ;; lines need to be set up very carefully to show the appropriate regions in the same colors
      ;; In practice I tend to do this separately using balanceplots and makethebalanceplots. 
      IF(fill EQ 1 and NIndVarBins EQ 0) THEN BEGIN
          aa = data[ti,*,k,*]
          FOR xx=0,n_elements(data[ti,*,k,0])-1 DO aa[0,xx,0,*] = aa[0,xx,0,SORT(aa[0,xx,0,*])]
          FOR jj=1,n_elements(data[0,0,0,*])-1 DO BEGIN
              POLYFILL, [[data[ti,*,0,jj-1]], [REVERSE(data[ti,*,0,jj],2)]], [[aa[0,*,0,jj-1]], [REVERSE(aa[0,*,0,jj],2)]], COLOR=colors[ti,jj];,line_fill=1
	      
          ENDFOR
      ENDIF


      ; if not ( binning and x-dependence)
      ; i.e. if not binning or if no x-dependence 
      ;;; PLOT INDIVIDUAL MODELS, either as functions of the x-coordinate or a single (x,y) point.
      nmodels = n_elements(data[0,0,0,*])
      nr=1
      IF(n_elements(data[0,*,0,0]) EQ 1) $
          THEN nr = MIN([nmodels,100]) $
          ELSE nr = MIN([nmodels,5])
      IF(nr NE nmodels) THEN randomSet = cgRandomIndices(nmodels,nr) ELSE randomSet=indgen(nmodels)
      IF(includeRandom EQ 1 or NIndVarBins EQ 0 or n_elements(data[0,*,0,0]) EQ 1) THEN BEGIN
          ;; loop over models (4th column of data)
          nplot = nmodels*(1-includeRandom) + n_elements(randomSet)*includeRandom
          mindices = intarr(nplot)
          FOR jj=0,nplot-1 DO BEGIN 
            IF(includeRandom NE 1) THEN j=nmodels-1-jj $ ;; reverse order!
                ELSE j=randomSet[jj]
            mindices[jj]=j
            ;; Set the color table.	
    ;        IF(n_elements(styles) GT 80) THEN IF(sv NE 3 AND sv NE 4) THEN setct,5,n_elements(styles),j ELSE setct,3,n_elements(styles),0
            
            ;; Overplot the data for this model. If prev is set, draw a tail to include previous values of the data.
            OPLOT, data[ti,*,0,j], data[ti,*,k,j], COLOR=colors[ti,j]+50*includeRandom,linestyle=styles[j]*(1-includerandom)+includeRandom*jj,PSYM=psym,SYMSIZE=symsize,THICK=thicknesses[j]
            IF(prev EQ 1) THEN BEGIN
              OPLOT,data[MAX([0,ti-tailLength]):ti,0,0,j],data[MAX([0,ti-tailLength]):ti,0,k,j], COLOR=colors[ti,j],linestyle=styles[j],THICK=thicknesses[j]
            ENDIF


          ENDFOR ;; end loop over models
          IF(n_elements(leg) NE 0 and tii EQ wf and (labels[k] NE 'tdep' and labels[k] NE 'fH2' and labels[k] NE 'tdepH2')) THEN BEGIN
              legend,leg[mindices],COLORS=colors[ti,mindices],THICK=thicknesses[mindices],/right,charsize=cs,charthick=chth,linestyle=intarr(n_elements(mindices)),box=0,horizontal=(labels[k] EQ 'colPerCrit')
          ENDIF

      ENDIF




      ;; All of the models have been plotted. Now if the user has requested it, we would like to 
      ;; bin the data by its independent variable, and then its dependent variable!
      ;; FOR THE CASE THAT THERE IS NO X-DEPENDENCE
      IF(NIndVarBins GT 0 and n_elements(data[0,*,0,0]) EQ 1) THEN BEGIN
          ;;;; ASSUMING that there are some galaxies in each color bin,
          ncolors = MAX(colors)+1
          theCurves = dblarr(NIndVarBins, n_elements(percentileList), ncolors)
          ;; each color gets its own percentile curves
          FOR aColor=0, ncolors-1 DO BEGIN
              subset = WHERE(aColor EQ colors[ti,*], NInThisSubset)
              binCenters = dblarr(NIndVarBins)
              IF(NInThisSubset GT 0) THEN BEGIN
                  FOR indVarIndex=0, NIndVarBins-1 DO BEGIN
                      ;; 
                      thisIndBin = GetIthBin(indVarIndex, data[ti,0,0,*], NIndVarBins, subset=subset)
                      IF(wrtXlog[0] EQ 1 ) THEN binCenters[indVarIndex] = 10.0^avg(alog10((data[ti,0,0,subset])[thisIndBin]))$
                        ELSE binCenters[indVarIndex] = avg((data[ti,0,0,subset])[thisIndBin])
                      
                      ; for this bin of the independent variable, for this color (i.e. subset of all models)
                      ; record the values of the dependent variable corresponding to the percentiles in
                      ; percentileList.
		              theCurves[indVarIndex, *, aColor] = percentiles((data[ti,0,k,subset])[thisIndBin],percentileList,wrtXlog[k])

                  ENDFOR ; end loop over independent variable bins.
                  FOR percentileIndex=0,n_elements(percentileList)-1 DO OPLOT, binCenters, theCurves[*,percentileIndex,aColor], COLOR=aColor,THICK=thicknesses[0],PSYM=-2,SYMSIZE=symsize*3
              ENDIF ; end check that there are models w/ this color
          ENDFOR ; end loop over color
          
      ENDIF

      ; FOR THE CASE THAT THERE /IS/ X-DEPENDENCE WHICH IS THE SAME FROM MODEL TO MODEL
      ; For now, we ignore the actual value of NIndVarBins and just find the percentiles for every single "x" value.
      nper = n_elements(percentileList)
      orientations = [0,0,30,20,0]
      linefills = [0,0,1,1,0]
      IF(NIndVarBins GT 0 and n_elements(data[0,*,0,0]) NE 1) THEN BEGIN
          ;;;; ASSUMING that there are some galaxies in each color bin,
          ncolors = MAX(colors)+1
          theCurves = dblarr(n_elements(data[0,*,0,0]), nper, ncolors)
          ;; each color gets its own percentile curves
          FOR aColorRev=0, ncolors-1 DO BEGIN
              aColor = ncolors-1-aColorRev ; reverse order! This is so the default model is plotted last.
              subset = WHERE(aColor EQ colors[ti,*], NInThisSubset)
;              binCenters = dblarr(NIndVarBins)
              IF(NInThisSubset GT 0) THEN BEGIN
                  binCenters = data[0,*,0,MIN(subset)]
                  FOR indVarIndex=0, n_elements(binCenters)-1 DO BEGIN
                      ;; the differences with the above IF block begin here. 
                      ;thisIndBin = GetIthBin(indVarIndex, data[ti,*,0,0], NIndVarBins, subset=subset)
                      ;IF(wrtXlog[0] EQ 1 ) THEN binCenters[indVarIndex] = 10.0^avg(alog10((data[ti,thisIndBin,0,subset]))) $
                      ;  ELSE binCenters[indVarIndex] = avg((data[ti,thisIndBin,0,subset]))
                      
                      ; for this bin of the independent variable, for this color (i.e. subset of all models)
                      ; record the values of the dependent variable corresponding to the percentiles in
                      ; percentileList.
		              theCurves[indVarIndex, *, aColor] = percentiles((data[ti,IndVarIndex,k,subset]),percentileList,wrtXlog[k])

                  ENDFOR ; end loop over independent variable bins.
		  ;stop
                  FOR percentileIndex=0,nper-1 DO BEGIN
                      theDist = ABS(percentileIndex-nper/2)
                      OPLOT, binCenters, theCurves[*,percentileIndex,aColor], COLOR=aColor, $
                        THICK=max(thicknesses)*1.3/(theDist+1.0),linestyle=0
                  ENDFOR
                  IF(fill EQ 1) THEN BEGIN
                      aa = theCurves[*,*,aColor]  ;data[ti,*,k,*]
                      wh = WHERE(aa LT ranges[0,k], ct)
                      IF(ct GT 0) THEN aa[wh] = ranges[0,k]
                      wh = WHERE(aa GT ranges[1,k], ct)
                      IF(ct GT 0) THEN aa[wh] = ranges[1,k]
                      FOR xx=0,n_elements(binCenters)-1 DO aa[xx,*] = aa[xx,SORT(aa[xx,*])]
                      FOR jj=1,n_elements(aa[0,*])-1 DO BEGIN
                          xp = [[binCenters],[REVERSE(binCenters,2)]]
                          yp = [aa[*,jj-1], REVERSE(aa[*,jj])]
                          IF(linefills[jj] EQ 1 ) THEN POLYFILL, xp[*], yp[*], COLOR=aColor,line_fill=linefills[jj], orientation=orientations[jj], thick=0.3,spacing=0.02
                          ;IF(linefills[jj] EQ 1 and acolor EQ 1) THEN stop
                      ENDFOR
                  ENDIF
              ENDIF ; end check that there are models w/ this color
          ENDFOR ; end loop over color
      ENDIF


      IF(n_elements(additionalLines) NE 0) THEN BEGIN
            FOR pcaI=0,n_elements(additionalLines[0,0,*])-1 DO BEGIN
                arrow,additionalLines[tii,k-1,0,pcaI],additionalLines[tii,k-1,1,pcaI], $
                    additionalLines[tii,k-1,2,pcaI],additionalLines[tii,k-1,3,pcaI], $
                    /data,color=additionalLineColors[pcaI],THICK=max(thicknesses)/2.0, $
                    hsize=-0.1
                XYOUTS,additionalLines[tii,k-1,2,pcaI], additionalLines[tii,k-1,3,pcaI], $
                    additionalLineLabels[pcaI],COLOR=additionalLineColors[pcaI], $
                    CHARSIZE=cs*.45
            ENDFOR    

      ENDIF

      yfiddle = 0.0
      IF(labels[k] EQ 'equilibrium') THEN yfiddle = -0.26
      IF(labels[k] EQ 'colPerCrit') THEN yfiddle = -0.26
      ;;; PRINT TEXT ON THE PLOTS
      IF(ngrp EQ 0 or (ngrp GT 0 and kx EQ 1)) THEN BEGIN
          IF(wf NE 0) THEN BEGIN ; if there's more than one frame, label them:
            IF(sv NE 5) THEN XYOUTS,.3,.87,theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth  $
              ELSE IF(horizontal EQ 0) THEN XYOUTS,.3,.95-float(tii)*(1.0 - 2.0/8.89)/(float(n_elements(whichFrames))),theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth $
              ELSE IF(wf EQ 1 and labels[k] NE 'tdep' and labels[k] NE 'tdepH2' and labels[k] NE 'fH2') THEN XYOUTS,.36+(.36+.01*ngrp)*float(tii)/float(n_elements(whichFrames)),.69-ngrp*.01,theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth $
              ELSE IF(wf EQ 2) THEN XYOUTS,.31+.36*float(tii)/float(n_elements(whichFrames)),.66+yfiddle,theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth $
              ELSE IF(wf EQ 3) THEN xyouts,.10+.67*float(tii)/float(n_elements(whichFrames)),.7,theTimeText,/NORMAL,COLOR=0,CHARSIZE=cs,CHARTHICK=chth
          ENDIF
      ENDIF

      IF(ngrp EQ 0 or (ngrp GT 0 and ngrp EQ kx)) THEN BEGIN
          ;; We've created a single frame. What do we do with it?
          IF(sv EQ 1) THEN MPEG_PUT,mpeg_id,WINDOW=1,FRAME=count,/ORDER
          IF(sv EQ 2) THEN WRITE_PNG, dn+'/frame_'+string(count,FORMAT="(I4.4)") + '.png',TVRD(TRUE=1)
          IF(sv EQ 3 || sv EQ 4) THEN dummy=cgSnapshot(filename=(dn+'/frame_'+STRING(count-1,FORMAT="(I4.4)")),/png,/true,/nodialog) 
          IF(sv EQ 0) THEN wait,.05
    ;      IF(sv EQ 5 AND tii NE wf AND wf NE 0) THEN 
      ENDIF
      IF(sv EQ 5 and (not (tii EQ wf and kx EQ ngrp) and ngrp GT 0) or (ngrp EQ 0 and tii NE wf)) THEN modmultiplot

    ENDFOR ;; end loop over time indices
    IF(ngrp EQ 0 or (ngrp GT 0 and ngrp EQ kx)) THEN BEGIN
        IF(sv EQ 1) THEN MPEG_SAVE,mpeg_id
        IF(sv EQ 1) THEN MPEG_CLOSE,mpeg_id
        ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN spawn,("rm " + fn+".mpg") 
        ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN spawn,("ffmpeg -f image2 -qscale 4 -i "+dn+"/frame_%04d.png " + fn + ".mpg")
        ;IF(sv EQ 2 || sv EQ 3 || sv EQ 4) THEN IF(saveFrames EQ 0) THEN spawn,("rm -rf "+dn)
        IF(sv EQ 5) THEN BEGIN
        	figureClean,(fn+"_timeSeries"),svSinglePlot
            scale = 1.5
            ;IF(ngrp GT 0) THEN scale = 1.0
            IF(svSinglePlot EQ 2) THEN latexify,(fn+"_timeSeries.eps"),axisLabels,texLabels,replicate(scale*cs,n_elements(texLabels)),height=8.89*theHeight,width=8.89*theWidth,tempname=("TMP_"+fn+"_timeSeries")
            modmultiplot,/default
            IF(sv EQ 5 and svSinglePlot EQ 2) THEN BEGIN

                spawn,"cat "+fn+"_timeSeries.eps | ps2eps --loose --preserveorientation > crp_"+fn+"_timeSeries.eps"
                spawn,"ps2pdf -dEPScrop crp_"+fn+"_timeSeries.eps "+fn+"_timeSeries.pdf"

            ENDIF
            IF(sv EQ 5 and svSinglePlot EQ 2) THEN spawn,"rm TMP_"+fn+"_timeSeries.eps"
    ;        IF(sv EQ 5 and svSinglePlot EQ 2) THEN spawn,"epstopdf "+fn+"_timeSeries.eps"
        ENDIF
        set_plot,'x'
    ENDIF
  ENDFOR ;; end loop over dependent variables
  IF(sv EQ 2|| sv EQ 3 || sv EQ 4) THEN spawn,"python makeGidgetMovies.py"
END

