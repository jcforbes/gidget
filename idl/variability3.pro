

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


FUNCTION vsMstarTransform,original,NVS
    vsMstarLabels = original[*]
    vsMstarLabels[0] = vsMstarLabels[2]
    vsMstarLabels[2:(NVS-2)] = vsMstarLabels[3:(NVS-1)]
    vsMstarLabels=vsMstarLabels[0:(NVS-2)]
    RETURN,vsMstarLabels
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
;;;   data for each run can become cumbersome. By default the code will
;;;   analyze all the variables
;;; N controls how many runs in the experiment to analyze. For large number 
;;;   runs, sometimes it's useful to just look at a few runs. If N is larger 
;;;   than the number of runs in the experiment, all the runs are analyzed.
;;; sv- behaves like other instances of sv; 0 to just show the plots on the 
;;;   screen, 2 to save them as png's and subsequently movies. 4 saves 
;;;   movies in a format conducive to showing on a low-resolution projector.
PRO variability3,expNames,N=N,sv=sv,keys=keys,annotations=annotations,integratedQuantities=integratedQuantities,useColors=useColors
  IF(n_elements(N) EQ 0) THEN N=1000000
  IF(n_elements(sv) EQ 0) THEN sv=3
  IF(n_elements(useColors) EQ 0) THEN useColors=[6] ; colorByExperiment. (each expt. gets its own color)
  IF(n_elements(annotations) EQ 0) THEN annotations=expNames[*]
  IF(n_elements(annotations) NE n_elements(expNames)) THEN $
      message,"expNames and annotations should have the same number of elements"

  notations=annotations
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

;  IF(n_elements(keys) EQ 0) THEN keys=[[indgen(38)+1],47,48,49,50]
;  tk=[0,keys[*]]
  IF(n_elements(keys) EQ 0) THEN BEGIN
      keys = [ $ ; compare to the data in getlabels.pro
           0, 1, 2, 3, 4, $  ; r, col, sig, colst, sigst
           5, 6, 7, 8, 9, $  ; f_g, Z, colsfr, Q, Qg
          10,11,12,   14, $  ; Qst, fH2, ageAtz0,   tdep
          15,   17,18,19, $  ; tvisc, mdotDisk, vrg, vst
             21,   23,    $  ; tdepH2, sfrPerAccr
          25,   27,28,29, $  ; accr, BB4, beta, u
                      34, $  ; mdotGtr0
             36,          $  ; sfPerAll
                47,   49, $  ; equilibrium, colPerCrit
          50,51           $  ; sigPerAn, MJeans
          ]

  ENDIF

  ;; variables to plot- title, column index, log plot, name, y-range
  dummy = GetLabels(keys=keys,tex=wrtXytex, labels=wrtXyt, names=wrtXyn,$
      log=wrtXyl, base=wrtXyy, offset=wrtXyp)

  IF(n_elements(integratedQuantities) EQ 0) THEN integratedQuantities=1
  n1 = IntegratedQuantities



  ;; now loop over every model in the list - 
  ;; all we're doing here is reading in each model then freeing 
  ;; the associated memory
  ;; so that we can iterate over a very large number. For each model, 
  ;; record the data specified in the arrays above (e.g. yy)
  modelPtr=ReadModels([nameList2[0]])
  model = (*(modelPtr[0]))
  theNTS = model.NOutputsNominal+2

  whichFrames = indgen(theNTS) ;; plot every output
  whichFrames = indgen((theNTS+2)/3)*3 ; plot every 3rd output.
  whichFrames3 = [1,51,201] ;; roughly z=2,1,0

  time = model.evArray[1,*]  
  z = model.evArray[10-1,*]

  ;; (# timesteps) x (nx) x (# of columns) x (# models)
  theData= dblarr(theNTS,n_elements(model.dataCube[0,*,0]),n_elements(wrtXyy),n_elements(nameList2))


  nameList4 = strarr(n_elements(nameList2))

  ; so far we have set up two data structures, namely one for quantities which vary in radius & time, 
  ; and one for quantities which vary only in time. Now in analogy to the first structure, we will 
  ; construct an array for stellar properties at redshift zero which vary with radius & age.
  byAge = dblarr(model.npassive+1,n_elements(model.dataCube[0,*,0]),5+1,n_elements(nameList2))
  vsAge = dblarr(3,model.npassive+1,5+1,n_elements(nameList2))

  proftimes[1]=systime(1)

;  FOR i=0,n_elements(nameList2)-1 DO BEGIN ;; loop over every valid model
  ctr=0
  modelcounter2 = modelcounter*0

  dummy=GetTimeLabels(names=vsMdotNames,tex=vsMdotTexLabels,labels=vsMdotLabels,$
        log=vsMdotToLog, offset=offset, base=base, ranges=vsMdotRanges, strt=vsMdotStrt)

  NVS = n_elements(vsMdotLabels)
  ;; (# timesteps) x (# pts/model frame = 1) x (one thing vs another = 20) x (# models)
  vsMdot = dblarr(theNTS,1,NVS,n_elements(nameList2))
  mh0s = dblarr(n_elements(nameList2))
  lambdas = dblarr(n_elements(nameList2))

  xBig = dblarr(theNTS,n_elements(theData[0,*,0,0]),1,n_elements(namelist2))

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
       r = x[*]*model.Radius

       ;; if this model has the appropriate number of time steps..
       ;;PRINT,"Number of time steps",n_elements(model.dataCube[*,0,0]),n_elements(theData[*,0,0,0])
       IF(n_elements(model.dataCube[*,0,0]) EQ theNTS) THEN BEGIN
         Mh0s[ctr] = model.mh0
         lambdas[ctr] = model.accScaleLength/(model.mh0 ^ (1.0/3.0))
         FOR k=0,n_elements(wrtXyy)-1 DO BEGIN	
           theData[*,*,k,ctr] = model.dataCube[*,*,wrtXyy[k]+offsets[wrtXyp[k]]-1]
         ENDFOR
         xBig[*,*,0,ctr] = model.dataCube[*,*,0]*model.Radius/model.accScaleLength
         rText=['r=1 kpc','r=8 kpc','r=15 kpc']
         positionIndices=[nearest("position",1,model,1.0/model.Radius), $
            nearest("position",1,model,8.0/model.Radius), $
            nearest("position",1,model,15.0/model.Radius)]
;         IF(model.npassive GT n_elements(byAge[*,0,0,0])) THEN BEGIN
;            byAgeT=temporary(byAge)[*,*,*,*]
;            byAge=dblarr(model.npassive+1,n_elements(model.dataCube[0,*,0]),5+1,n_elements(namelist2))
;            byAge[0:n_elements(byAgeT[*,0,0,0])-1,*,*,*] = temporary(byAgeT)[*,*,*,*]
;         ENDIF
         FOR popI=0, model.npassive DO BEGIN
            byAge[popI, *, 0, ctr] = model.dataCube[0,*,0]*model.Radius
            byAge[popI, *, 1:5, ctr] = model.stellarData[theNTS-1, popI, *, 0:4]
            FOR m=0,n_elements(positionIndices)-1 DO BEGIN
                vsAge[m,popI,1:5, ctr] = model.stellarData[theNTS-1,popI, positionIndices[m], 0:4]
            ENDFOR
         ENDFOR
         ages = model.stellarAges[theNTS-1,*]
         FOR m=0,n_elements(positionIndices)-1 DO vsAge[m,*,0,ctr] = ages[*]




         FOR zi=0, n_elements(model.evArray[9,*])-1 DO BEGIN ; loop over redshift
           modelInfo = diskStats(model,z=model.evArray[9,zi])
           
           wh = where(offset EQ 0, ct)
           IF(ct GT 0) THEN vsMdot[zi,0,wh,ctr] = model.evArray[base[wh]-1,zi]
           wh = where(offset EQ 1, ct)
           IF(ct GT 1) THEN vsMdot[zi,0,wh,ctr] = modelInfo[base[wh]-1]

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




  byAgeLabels=["Radius (kpc)","Column Density","r Vel Dis","z Vel Dis","Z","var Z"]
  byAgeTexLabels=["Radius (kpc)","$\Sigma_{*,i}$","$\sigma_{*,r}$","$\sigma_{*,z}$","Z","var Z"]
  byAgeNames=["r","stCol","sigr","sigz","stZ","stVZ"]
  byAgeToLog=[1,1,1,1,0,0]

  vsAgeLabels=byAgeLabels[*]
  vsAgeTexLabels=byAgeTexLabels[*]
  vsAgeNames=byAgeNames[*]
  vsAgeToLog=byAgeToLog[*]
  vsAgeLabels[0]="Ages (Ga)"
  vsAgeTexLabels[0]="Ages (Ga)"
  vsAgeNames[0]="ages"

  modelcounter = modelcounter2[*]
  nameList2 = nameList4[0:ctr-1]
  lambdas=lambdas[0:ctr-1]
  Mh0s=Mh0s[0:ctr-1]

  proftimes[2]=systime(1)

  theData = (temporary(theData))[*,*,*,0:ctr-1]
  xBig = (temporary(xBig))[*,*,*,0:ctr-1]
  theDataScaleR = theData[*,*,*,*]
  theDataScaleR[*,*,0,*] = xBig[*,*,*,*]
  resample,theData
  resample,theDataScaleR

  colors = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))

  colorsByMh0 = intarr(n_elements(nameList2))
  colorsByLambda = intarr(n_elements(nameList2))
  colorsByExperiment = intarr(n_elements(nameList2))

  colorsByMS = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))
  colorsByMSandMass = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))
  colorsByDeltaAccr = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))
  colorsByDeltaAccrAndMass = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))

  vsMdot = (temporary(vsMdot))[*,*,*,0:ctr-1]
  byAge = (temporary(byAge))[*,*,*,0:ctr-1]
  vsAge = (temporary(vsAge))[*,*,*,0:ctr-1]

  thePercentiles = [0.025,.16,.5,.84,.975]
  nper = n_elements(thePercentiles)



  proftimes[3]=systime(1)

  FOR expInd=1, n_elements(expNames) DO BEGIN
    low=modelcounter[expInd-1]
    high=modelcounter[expInd]-1
    colors[*,low:high] = expInd-1
    colorsByExperiment[low:high] = expInd-1
  ENDFOR


  nmodels = n_elements(vsMdot[0,0,0,*])
  ncolors = 4
  zi0 = n_elements(vsMdot[*,0,0,0])-1
  FOR zi = 0, zi0 DO BEGIN
    FOR nc=0, ncolors-1 DO BEGIN
      ; sort according to z=0 sSFR - feel free to change this!
      ;; examples include structural parameters, e.g. CentralDensity, rPeak, rHalfSFR,..
      ;; also metallicity, initial mass,...
      ;colors[zi,GetIthBin(nc, vsMdot[zi0,0,6,*], ncolors)] = nc
      
    ENDFOR
  ENDFOR

  FOR nc=0, ncolors-1 DO BEGIN
    ; sort according to z=0 sSFR - feel free to change this!
    ;; examples include structural parameters, e.g. CentralDensity, rPeak, rHalfSFR,..
    ;; also metallicity, initial mass,...
    ;colors[zi,GetIthBin(nc, vsMdot[zi0,0,6,*], ncolors)] = nc
    
    colorsByMh0[GetIthBin(nc, Mh0s,ncolors)] = nc
    colorsByLambda[GetIthBin(nc, lambdas,ncolors)] = nc
  ENDFOR




  thicknesses=intarr(nmodels)+1
  Mh0 = vsMdot[theNTS-1,0,2,*]/vsMdot[theNTS-1,0,14,*]
  FOR nc=0, ncolors-1 DO BEGIN
      thicknesses[GetIthBin(nc,Mh0,ncolors)] = 12.0 + .5*float(nc) ;; thicker lines have higher Mh0.
  ENDFOR


  expName2 = expNames[0]
;  IF(n_elements(expNames) LT 5) THEN BEGIN
      FOR i=1, N_elements(expNames)-1 DO BEGIN
          expName2 = expName2 + '_' + expNames[i];; +'_'+strcompress(string(MIN([n_elements(nameList2),N])),/remove) ;; e.g. rk5_113
      ENDFOR
;  ENDIF ELSE expName2 +='_' + strcompress(string(n_elements(expNames)),/remove)
  expName2 = expName2 + '_' + strcompress(string(MIN([n_elements(nameList2),N*n_elements(expNames)])),/remove) ;; e.g. rk5_rk6_113




  unsRanges = simpleRanges(thedata[1:n_elements(time)-1,*,*,*],wrtXyl)
  wh = WHERE(keys EQ 1 or keys EQ 3, ct)
  IF(ct GT 0) THEN unsRanges[0,wh] = 0.3 ;; manually set minimum col and colst to 1 Msun/pc^2
  wh = WHERE(keys EQ 7, ct)
  IF(ct GT 0) THEN unsRanges[0,7] =  1.0d-5 ;; manually set min SFR col density to 10^-5 Msun/yr/kpc^2
  IF(ct GT 0) THEN unsRanges[1,7] =  1.0d2 ;; manually set max SFR col density to 10^2 Msun/yr/kpc^2
  wh = WHERE(keys EQ 8 or keys EQ 9 or keys EQ 10, ct)
  IF(ct GT 0) THEN unsRanges[1,wh] = 50.0 ;; manually set max Q,Qgas,Qst
  wh = WHERE(keys EQ 17 , ct)
  IF(ct GT 0) THEN unsRanges[0,wh] = -2.8 ;; mdotDisk
  IF(ct GT 0) THEN unsRanges[1,wh] = 5.8
  wh = WHERE(keys EQ 3-1 , ct)
  IF(ct GT 0) THEN  unsRanges[0,wh] = 0.0 ;; velocity dispersion
  IF(ct GT 0) THEN  unsRanges[1,wh] = 35.0
  wh = WHERE(keys EQ 50-1, ct)
  IF(ct GT 0) THEN  unsRanges[0,wh] = .09 ;; colPerCrit
  IF(ct GT 0) THEN  unsRanges[1,wh] = 3.1
  wh = WHERE(keys EQ 28-1 , ct)
  IF(ct GT 0) THEN  unsRanges[0,wh] = 0.099 ; col/BB profile
  IF(ct GT 0) THEN  unsRanges[1,wh] = 21.0
  wh = WHERE(keys EQ 52-1 , ct)
  IF(ct GT 0) THEN  unsRanges[0,wh] = 3.0e6
  IF(ct GT 0) THEN  unsRanges[1,wh] = 3.0e9 ; 2D Jeans Mass
  
  unsRangesLin = unsRanges[*,*]
  unsRangesLin[0,0] = - 0.5
  unsRangesLin[1,0] = 39.999

  IF(unsRanges[1,0] LT unsRangesLin[1,0]) THEN unsRangesLin[1,0] = unsRanges[1,0]
;  unsRanges[1,0] = unsRanges[1,0]-.3

  unsRangesScaled = unsRangesLin[*,*]
  unsRangesScaled[0,0] = -.01
  unsRangesScaled[1,0] = 1.9999


  linestyles = intarr(n_elements(expNames)*nper)+2;nper-1
  midpoints = where(indgen(n_elements(expNames)*nper) MOD nper EQ nper/2)
  linestyles[midpoints] = 0
  intervalsColors0 = indgen(n_elements(expNames)*nper)/nper
  intervalsColors = intarr(n_elements(time),n_elements(intervalsColors0))
  FOR tt=0,n_elements(time)-1 DO intervalsColors[tt,*] = intervalsColors0 


  unsThicknesses=colors*0+2
  unsThicknesses[0]=12
	


  ;; average mdot over a number of time steps nt
;  nt = 40
;  avgdVsMdot = vsMdot[*,*,*,*]
;
;  FOR ti=0,n_elements(time)-1 DO BEGIN
;    FOR mm=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
;      avgdVsMdot[ti,0,0,mm] = TOTAL(vsMdot[MAX([0,ti-nt]):ti,0,0,mm]) / (ti - MAX([0,ti-nt]) + 1)
;    ENDFOR
;  ENDFOR



;     vsAvgMdotLabels = vsMdotLabels[*]
;     vsAvgMdotLabels[0] = "Average Mdot Over Past 40 Outputs (Msun/yr)"
;     vsAvgMdotTexLabels = vsMdotTexLabels[*]
;     vsAvgMdotTexLabels[0] = "$\dot{M}_{40 \mathrm{outputs}} (M_\odot/yr)$"
;     vsAvgMdotNames = vsMdotNames[*]
;     vsAvgMdotNames[0] = "avgMdot"



;    vsSpecificSFR = vsMdot[*,*,*,*]
;    vsSpecificSFR[*,*,0,*] = vsMdot[*,*,6,*]
;    vssSFRLabels = vsMdotLabels[*]
;    vssSFRLabels[0] = vsMdotLabels[6]
;    vssSFRTexLabels = vsMdotTexLabels[*]
;    vssSFRTexLabels[0] =vsMdotTexLabels[6]
;    vssSFRNames = vsMdotNames[*]
;    vssSFRNames[0] = vsMdotNames[6]
;    vssSFRToLog = vsMdotToLog[*]
;    vssSFRToLog[0] = vsMdotToLog[6]


 ; ;; Plot these quantities vs. the star formation rate.
;  vsSFR = vsMdot[*,*,1:n_elements(vsMdot[0,0,*,0])-1,*] ;; cut out the mdot slice of the array

;    vsSFRLabels = vsMdotLabels[1:(n_elements(vsMdotLabels)-1)]
;    vsSFRTexLabels = vsMdotTexLabels[1:n_elements(vsMdotTexLabels)-1]
;    vsSFRNames = vsMdotNames[1:(n_elements(vsMdotNames)-1)]
;    vsSFRToLog = vsMdotToLog[1:(n_elements(vsMdotToLog)-1)]
;    vsSFRstrt= vsMdotstrt[1:(n_elements(vsMdotStrt)-1)]



  ;; Now do a series of plots vs. M*
  vsMstar = vsMdot[*,*,*,*]
  vsMstar[*,*,0,*] = vsMstar[*,*,2,*] ;replace Mdot with M* as the zeroeth column of the data.
  vsMstar[*,*,2:(NVS-2),*] =vsMstar[*,*,3:(NVS-1),*] ;move over the other columns so that we don't produce an M* vs M* plot
  vsMstar=vsMstar[*,*,0:(NVS-2),*] ; eliminate the extraneous last column.

  vsMstarLabels = vsMstarTransform(vsMdotLabels,NVS)
  vsMstarTexLabels = vsMstarTransform(vsMdotTexLabels,NVS)
  vsMstarNames = vsMstarTransform(vsMdotNames,NVS)
  vsMstarStrt = vsMstarTransform(vsMdotStrt,NVS)
  vsMstarToLog = vsMstarTransform(vsMdotToLog,NVS)

  vsMstarRanges = simpleRanges(vsMstar[2:n_elements(vsMstar[*,0,0,0])-1,*,*,*],vsmstartolog )
  vsMstarRanges[0,0] = 1.0e9 ; just for the rg70ish series.


  deltams,vsMstar,binnedMS,binnedAvgAccr,theDeltas,deltaAccr,deltaAvgAccr,widths,widthAccr,widthsAvgAccr,msText,NMassBins=1,NDeltaBins=3
  colorsByMS = binnedMS[*,*]
  colorsByDeltaAccr = binnedAvgAccr[*,*]
  deltams,vsMstar,binnedMS,binnedAvgAccr,theDeltas,deltaAccr,deltaAvgAccr,widths,widthAccr,widthsAvgAccr,msText,NMassBins=1,NDeltaBins=5
  colorsByZ0MS = dblarr(n_elements(binnedMS[0,*]))
  colorsByZ0MS[*] = binnedMS[zi0,*]
  deltams,vsMstar,binnedMS,binnedAvgAccr,theDeltas,deltaAccr,deltaAvgAccr,widths,widthAccr,widthsAvgAccr,msText,NMassBins=2,NDeltaBins=2
  colorsByMSAndMass = binnedMS[*,*]
  colorsByDeltaAccrAndMass = binnedAvgAccr[*,*]
;  deltams,vsMstar,binnedMS,binnedAvgAccr,theDeltas,deltaAccr,deltaAvgAccr,widths,widthAccr,widthsAvgAccr,msText,NMassBins=3,NDeltaBins=1
  

;  deltaMS,vsMstar, per , msColors, theDeltas, deltaAccr, widths, widthAccr, msText, NMassBins = 2
;  IF(useMsColors EQ 1) THEN colors=msColors
  vsMstar[*,0,60-2,*] = theDeltas
  vsMstar[*,0,61-2,*] = deltaAccr

  vsTimeLabels = vsMdotLabels[*]
  vsTimeLabels[0] = "Time Since z=2 (Ga)"
  vsTimeTexLabels = vsMdotTexLabels[*]
  vsTimeTexLabels[0] = "Time since $z=2$ (Ga)"
  vsTimeNames = vsMdotNames[*]
  vsTimeToLog = vsMdotToLog[*]
  vsTimeToLog[0] = 0

  nmodels = n_elements(vsMdot[0,0,0,*])
  vsTime = dblarr(1,n_elements(vsMdot[*,0,0,0])-1,n_elements(vsMdot[0,0,*,0]),nmodels)
  FOR tt=1, n_elements(vsMdot[*,0,0,0])-1 DO BEGIN
    vsTime[0,tt-1,*,*] = vsMdot[tt,0,*,*]
    vsTime[0,tt-1,0,*] = time[tt]
  ENDFOR

  vsTimeStyles=colors*0
  vsTimeColors=colors[*,*]
  IF(n_elements(vsTime[0,0,0,*]) GT 50) THEN vsTimeStyles=temporary(vsTimeStyles)+1
  vsTimeRanges = simpleRanges(vstime,vstimetolog)
  vsTimeRanges[0,53-1] = 0.95e11 ; Lower Mh
  vsTimeRanges[1,53-1] = 1.05e12 ; upper Mh
  vsTimeRanges[0,[56,57]-1] = 0.0
  vsTimeRanges[1,[56,57]-1] = 44.9
  vsTimeRanges[0,58-1] = 0.0
  vsTimeRanges[1,58-1] =  2.99
  vsTimeRanges[0,13-1] = 0.0
  vsTimeRanges[1,13-1] = 0.5
  vsTimeRanges[0,59-1] = 0.0
  vsTimeRanges[1,59-1] = 1.0 ; bt ratio

  svSinglePlot=2 ;; 2=eps, 4=png

  setct,3,n_elements(vsMdot[0,0,0,*]),0
;  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2,sv=4,nt0=2,thicknesses=thicknesses
  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2+"_log",sv=svSinglePlot,nt0=2,thicknesses=thicknesses,logarithms=vsMdotToLog,normalize=(n_elements(expNames) NE n_elements(vsMdot[0,0,0,*])),texLabelsIn=vsMdotTexLabels



  ;;;; MOVIES ;;;;;;
;;  simpleMovie,intervals,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervals",sv,axislabels=wrtXyt,whichFrames=whichFrames,ranges=unsRanges
;;  IF(n_elements(vsSpecificSFR[0,0,0,*]) LE 50) THEN simpleMovie, vsSpecificSFR,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFR",sv,PSYM=1,prev=1,axisLabels=vssSFRLabels,taillength=2,whichFrames=whichFrames
;;  IF(n_elements(vsSpecificSFR[0,0,0,*]) GT 50) THEN simpleMovie, vsSpecificSFR,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFRDist",sv,PSYM=1,prev=0,axisLabels=vssSFRLabels,NIndVarBins=5,whichFrames=whichFrames
;;  IF(n_elements(vsMstar[0,0,0,*]) LE 50) THEN simpleMovie, vsMstar, vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstar",sv,PSYM=1,prev=0,axisLabels=vsMstarLabels,whichFrames=whichFrames
;;  simpleMovie, avgdVsMdot,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot",sv,PSYM=1,prev=1,axisLabels=vsAvgMdotLabels,whichFrames=whichFrames,thicknesses=unsThicknesses
;;  simpleMovie,theData,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsorted",sv,axislabels=wrtXyt,whichFrames=whichFrames,ranges=unsRanges
;;;    simpleMovie, vsSFR, vsSFRNames, colors, colors*0, vsSFRtoLog,expName2+"_vsSFR",sv,vsSFRstrt,PSYM=1,prev=1,axisLabels=vsSFRLabels
;;;    simpleMovie, vsMdot,vsMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsmdot",sv,vsMdotStrt,PSYM=1,prev=1,axisLabels=vsMdotLabels




   ;;;;; STANDALONE PNG's (or eps as the case may be) ;;;;;
   ;vs radius. First two are log radius, second two are identical with linear radius
  TheNumberOfModels = n_elements(vsMstar[0,0,0,*])

  whichFramesz0 = [n_elements(z)-1]
  whichFramesz2 = [1]
  whichFrames5 = [1,20,51,105,201]; a guess for z=2,1.5,1,.5,0
  whichFramesThird = [indgen(n_elements(z)/3) * 3 +1,n_elements(z)-1] ; every third element, starting with 2nd frame and ending w/ last frame
  whichFramesAll = indgen(n_elements(z)-1)+1
  z5 = strarr(5)  ;z[whichFrames5]
  z3 = strarr(3)
  zAll= strarr(n_elements(z)-1)
  zThird = strarr(n_elements(whichFramesThird))
  FOR zi=0,3-1 DO z3[zi] = 'z='+string(abs(z[whichFrames3[zi]]),Format='(D0.2)')
  FOR zi=0,5-1 DO z5[zi] = 'z='+string(abs(z[whichFrames5[zi]]),Format='(D0.2)')
  FOR zi=0,n_elements(zAll)-1 DO zAll[zi] = 'z='+string(abs(z[whichFramesAll[zi]]),Format='(D0.2)')
  FOR zi=0,n_elements(zThird)-1 DO zThird[zi] = 'z='+string(abs(z[whichFramesThird[zi]]),Format='(D0.2)')

  setct,1,n_elements(linestyles),2



    ;; Try this out.
  FOR ccc=0,n_elements(useColors)-1 DO BEGIN
    cc=useColors[ccc]
    IF(cc EQ 0) THEN BEGIN
        colors= colorsByMh0
        app = "_cbmh0"
    ENDIF
    IF(cc EQ 1) THEN BEGIN
        colors = colorsByLambda
        app = "_cbl"
    ENDIF
    IF(cc EQ 2) THEN BEGIN
        colors=colorsByMS
        app = "_cbms"
    ENDIF
    IF(cc EQ 3) THEN BEGIN
        colors=colorsByMsAndMass
        app="_cbmsm"
    ENDIF
    IF(cc EQ 4) THEN BEGIN
        colors=colorsByDeltaAccr
        app="_cba"
    ENDIF
    IF(cc EQ 5) THEN BEGIN
        colors=colorsByDeltaAccrAndMass
        app="_cbam"
    ENDIF
    IF(cc EQ 6) THEN BEGIN
        colors=colorsByExperiment
        app="_cbe"
    ENDIF
    IF(cc EQ 7) THEN BEGIN
        colors=colorsByZ0MS 
        app="_cbz0ms"
    ENDIF

    IF(thenumberOfModels LT 20) THEN BEGIN
      simpleMovie,theData,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsortedLogR"+app, $
          5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRanges,horizontal=1,thicknesses=unsThicknesses, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
      simpleMovie,theData,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsorted"+app, $
          5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRangesLin,horizontal=1,thicknesses=unsThicknesses, $
           svSinglePlot=svSinglePlot,texLabels=wrtXytex
  
      simpleMovie,theData,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedz2"+app, $
          5,axislabels=wrtXyt,whichFrames=whichFramesz2,ranges=unsRangesLin,horizontal=1,thicknesses=unsThicknesses, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex


    ENDIF
    IF(thenumberofmodels GT 100) THEN BEGIN
    simpleMovie,theData,wrtXyn,colors,colors*0+2,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedDistFill"+app, $
      5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRangesLin,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex,NIndVarBins=20,percentileList=[.025,.16,.5,.84,.975], $
      replacementText=z3,fill=1
    simpleMovie,theDataScaleR,wrtXyn,colors,colors*0+2,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_scaleDistFill"+app, $
      5,axislabels=['x',wrtXyt[1:n_elements(wrtXyt)-1]],whichFrames=whichFrames3,ranges=unsRangesScaled,horizontal=1, $
      thicknesses=unsThicknesses, svSinglePlot=svSinglePlot,texLabels=["$x = r/R$",wrtXytex[1:n_elements(wrtXytex)-1]], $
      NIndVarBins=20,percentileList=[.025,.16,.5,.84,.975], replacementText=z3,fill=1
    ENDIF


;  IF(TheNumberOfModels GT 10) THEN $
;;      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDist"+app, $
;;          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
;;          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
;;          ranges=vsTimeRanges
;  IF(TheNumberOfModels GT 10) THEN $
;      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroup"+app, $
;          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
;          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
;          ranges=vsTimeRanges, grouping = [53,18]-1,horizontal=1 ; Mh
;  IF(TheNumberOfModels GT 10) THEN $
;      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroup"+app, $
;          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
;          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
;          ranges=vsTimeRanges, grouping = [57,56,58]-1,horizontal=1;sigmax
;  IF(TheNumberOfModels GT 10) THEN $
;      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroup"+app, $
;          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
;          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
;          ranges=vsTimeRanges, grouping = [54,52]-1,horizontal=1 ;bbparams

  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistFill"+app, $
          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
          ranges=vsTimeRanges,fill=1
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroupFill"+app, $
          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
          ranges=vsTimeRanges, grouping = [53,18]-1,horizontal=1,fill=1 ; Mh
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroupFill"+app, $
          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
          ranges=vsTimeRanges, grouping = [57,56,58]-1,horizontal=1,fill=1;sigmax
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimecolors*0+2,vsTimeToLog,expName2+"_vstimeDistGroupFill"+app, $
          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975], $
          ranges=vsTimeRanges, grouping = [54,52]-1,horizontal=1,fill=1 ;bbparams



      IF(TheNumberOfModels LT 20) THEN $
  simpleMovie, vsTime,vsTimeNames,vsTimecolors,vsTimeStyles,vsTimeToLog,expName2+"_vstime"+app, $
      5,axisLabels=vsTimeLabels,whichFrames=[0],thicknesses=unsThicknesses,svSinglePlot=svSinglePlot, $
      texLabels=vsTimeTexLabels, ranges=vsTimeRanges
;  simpleMovie, avgdVsMdot,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot", $
;      5,PSYM=1,axisLabels=vsAvgMdotLabels,horizontal=1,whichFrames=whichFrames3,thicknesses=unsThicknesses
  ;simpleMovie, byAge,byAgeNames,colors,colors*0,byAgeToLog,expName2+"_byAgeLogR", $
  ;    5,axisLabels=byAgeLabels,whichFrames=[0,1,5,10],horizontal=1,thicknesses=unsThicknesses, $
  ;    svSinglePlot=svSinglePlot,texLabels=byAgeTexLabels,timeText="age = "
  IF(n_elements(ages) GT 3) THEN whichAges=[0,n_elements(ages)/2,n_elements(ages)-1]
  IF(n_elements(ages) LE 3) THEN whichAges=indgen(n_elements(ages))
  ageText=strarr(n_elements(whichAges))
  FOR st=0,n_elements(ageText)-1 DO ageText[st] = "Age ="+ string(ages[st])
  simpleMovie, byAge,byAgeNames,colors,colors*0,[0,byAgeToLog[1:n_elements(byAgeToLog)-1]],expName2+"_byAge"+app, $
      5,axisLabels=byAgeLabels,whichFrames=whichAges,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=byAgeTexLabels,replacementText=ageText
  simpleMovie, vsAge,vsAgeNames,colors,colors*0,vsAgeToLog,expName2+"_vsAge", $
      5,axisLabels=vsAgeLabels,whichFrames=[0,1,2],horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=vsAgeTexLabels,replacementText=rText


  IF(TheNumberOfModels GT 30) THEN $
  vsMstarPCA,vsMstar,whichFrames=[1,51,201],log=vsMStarToLog,names=vsMstarNames,colors=colors, $
      expName2=expName2+app,labels=vsMstarLabels,texLabels=vsMstarTexLabels,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,replacementText=replacementText


    
  IF(TheNumberOfModels GT 10) THEN $
       simpleMovie, vsMstar, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarH"+app, $
           5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201], $
           horizontal=1,thicknesses=unsThicknesses,svSinglePlot=svSinglePlot



  ;; We've made some plots where each panel is a specific time and each line is a specific model. Maybe we can do the converse.
  ;; Right now unsorted ~ intervals ~ (time) (pos) (variable) (model).
  ;; The thing simpleMovie really does is ~ (pane) (pos) (plot) (line)
  IF(TheNumberOfModels LE 3) THEN BEGIN
    singlePane = dblarr(TheNumberOfModels,n_elements(theData[0,*,0,0]),n_elements(theData[0,0,*,0]),5)
    FOR j=0,TheNumberOfModels-1 DO BEGIN
      FOR ii=0,5-1 DO BEGIN
	    singlePane[j,*,*,ii] = theData[whichFrames5[ii],*,*,j]
      ENDFOR
    ENDFOR
    spColors= indgen(5)
    FOR j=1,TheNumberOfModels-1 DO spColors=[[spColors],[indgen(5)]]
    simpleMovie,singlePane,wrtXyn,transpose(spcolors), $
        transpose(spcolors)*0,wrtXyl,expName2+"_unsortedSPlogR"+app, $
    	5,axisLabels=wrtXyt,whichFrames=indgen(TheNumberOfModels),ranges=unsRanges, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations, $
        leg = z5
    simpleMovie,singlePane,wrtXyn,transpose(spcolors), $
        transpose(spcolors)*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedSP"+app, $
     	5,axisLabels=wrtXyt,whichFrames=indgen(TheNumberOfModels),ranges=unsRangesLin, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations, $
        leg = z5
    wh = WHERE(keys EQ 8-1 or keys EQ 15-1 or keys EQ 12-1 or keys EQ 22-1, ct)
    simpleMovie,singlePane,wrtXyn,transpose(spcolors), $
        transpose(spcolors)*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedSPGroup"+app, $
     	5,axisLabels=wrtXyt,whichFrames=indgen(TheNumberOfModels),ranges=unsRangesLin, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations, $
        leg = z5, grouping = wh

  ENDIF



 IF(TheNumberOfModels GT 10) THEN BEGIN
    singlePaneAvg = dblarr(n_elements(expNames),n_elements(theData[0,*,0,0]),n_elements(theData[0,0,*,0]),5)
    FOR j=0,n_elements(expNames)-1 DO BEGIN
      FOR ii=0,5-1 DO BEGIN
        IF(modelCounter[j] NE modelCounter[j+1]-1) THEN $
            singlePaneAvg[j,*,*,ii] = avg(theData[whichFrames5[ii], *,*, modelCounter[j]:modelCounter[j+1]-1],3) $
            ELSE singlePaneAvg[j,*,*,ii] = theData[whichFrames5[ii], *, *, modelCounter[j]]
      ENDFOR
    ENDFOR
    spColors= indgen(5)
    FOR j=1,n_elements(expNames)-1 DO spColors=[[spColors],[indgen(5)]]
    simpleMovie,singlePaneAvg,wrtXyn,transpose(spColors), $
        transpose(spColors)*0,wrtXyl,expName2+"_intervalsSPAvglogR"+app, $
    	5,axisLabels=wrtXyt,whichFrames=indgen(n_elements(expNames)),ranges=unsRanges, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText='Average of '+notations
    simpleMovie,singlePaneAvg,wrtXyn,transpose(spColors), $
        transpose(spColors)*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_intervalsSPAvg"+app, $
	    5,axisLabels=wrtXyt,whichFrames=indgen(n_elements(expNames)),ranges=unsRangesLin,$
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText='Average of '+notations

  ENDIF




;  ;;; These 6 are movies
;    simpleMovie,theData,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedmv", $
;      4,axislabels=wrtXyt,ranges=unsRangesLin,horizontal=1,thicknesses=unsThicknesses, $
;      texLabels=wrtXytex,whichFrames=whichFramesAll,replacementText=zAll

IF(TheNumberOfModels GT 10) THEN BEGIN

    simpleMovie,theData,wrtXyn,colors,colors*0+2,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedDistFillmv"+app, $
      4,axislabels=wrtXyt,ranges=unsRangesLin,horizontal=1,thicknesses=unsThicknesses*.1, $
      texLabels=wrtXytex,NIndVarBins=20,percentileList=[.025,.16,.5,.84,.975], $
      replacementText=zAll,fill=1,whichFrames=whichFramesAll

  simpleMovie,theDataScaleR,wrtXyn,colors,colors*0+2,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_scaleDistFillmv"+app, $
      4,axislabels=['x (r/racc)',wrtXyt[1:n_elements(wrtXyt)-1]],ranges=unsRangesScaled,horizontal=1,thicknesses=unsThicknesses*.1, $
      texLabels=["$x = r/R$",wrtXytex[1:n_elements(wrtXytex)-1]],NIndVarBins=20,percentileList=[.025,.16,.5,.84,.975], $
      replacementText=zAll,fill=1,whichFrames=whichFramesAll
  ENDIF

 ; IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, vsMstarNames, colors, colors*0, $
;      vsMstarToLog,expName2+"_vsMstarDist",4,PSYM=1,prev=0,axisLabels=vsMstarLabels,NIndVarBins=5, $
;      percentileList=[.17,.5,.83],whichFrames=whichFramesAll,replacementText=zAll

      ;; This one is a movie.
  IF(TheNumberOfModels GT 10) THEN $
       simpleMovie, vsMstar, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstar"+app, $
           4,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels, $
           horizontal=1,thicknesses=unsThicknesses*.1,whichFrames=whichFramesAll, replacementText=zAll, $
           ranges=vsMstarRanges
   ;; as is this one
  IF(TheNumberOfModels GT 10) THEN $
       simpleMovie, vsMstar, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarPrev"+app, $
           4,PSYM=1,prev=1,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels, $
           horizontal=1,thicknesses=unsThicknesses,tailLength=2,whichFrames=whichFramesAll, $
           replacementText=zAll,ranges=vsMstarRanges


  ENDFOR

	
  proftimes[4]=systime(1)

  PRINT,proftimes[1]-proftimes[0]," seconds to check model validities and initialize"
  PRINT,proftimes[2]-proftimes[1]," seconds to read in data"
  PRINT,proftimes[3]-proftimes[2]," seconds to sort the data"
  PRINT,proftimes[4]-proftimes[3]," seconds to generate movies"

  STOP

END
