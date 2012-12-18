
;; Given an array of mdot's vs. other variables, compute the cross-correlations between mdot and the others
;; and also the autocorrelation of each variable.
PRO ComputeCorrelations,vsmdot,colors,time,labelsIn,names,name,sv=sv,nt0=nt0,thicknesses=thicknesses,logarithms=logarithms,normalize=normalize,percentileList=percentileList
    IF(n_elements(nt0) EQ 0) THEN nt0 = 1
    IF(n_elements(sv) EQ 0 ) THEN sv=4
    IF(n_elements(normalize) EQ 0) THEN normalize=1
    IF(n_elements(thicknesses) EQ 0) THEN thicknesses=intarr(n_elements(vsMdot[0,0,0,*]))+1 ; set each model's thickness to 1
    IF(n_elements(logarithms) EQ 0) THEN logarithms=intarr(n_elements(vsMdot[0,0,*,0])) ; set each variable to be not logarithmic
    IF(n_elements(percentileList) EQ 0) THEN percentileList=[.16,.5,.84]
    labels=labelsIn[*]
    FOR j=0, n_elements(logarithms)-1 DO IF(logarithms[j] EQ 1) THEN labels[j]="log "+labelsIn[j]
    cs = 2
    ct = 1
    IF(sv EQ 4) THEN cs = 3
    IF(sv EQ 4) THEN ct = 3
    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)*3
    ls=0
    nmodels = n_elements(vsMdot[0,0,0,*])
    nvars = n_elements(vsMdot[0,0,*,0])
    IF(nmodels GT 50) THEN ls =1 


    normalizations = dblarr(n_elements(vsMdot[*,0,0,0]), n_elements(vsMdot[0,0,*,0]))
    FOR k=0, n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
        FOR ti=0, n_elements(vsMdot[*,0,0,0])-1 DO BEGIN
            IF(normalize EQ 1) THEN normalizations[ti,k] = median(vsMdot[ti,0,k,*]) ELSE medians[ti,k]=1.0
        ENDFOR
    ENDFOR


    ;; Construct some cross-correlations
    nt = n_elements(vsMdot[*,0,0,0])
    ;;; (# of delay times to test) x (# of cross-correlations) x (# of models)
    crossCorrelations = dblarr(1,nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    autoCorrelations = dblarr(1,nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    ccPercentiles = dblarr(nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(percentileList))
    acPercentiles = dblarr(nt-1,nt0,n_elements(vsMdot[0,0,*,0]),n_elements(percentileList))
    ;;; loop over models
    FOR j=0, nmodels-1 DO BEGIN
        ;; loop over variables
        FOR k=0, n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
            laggingVar = vsMdot[nt0:(nt-1),0,k,j]/normalizations[nt0:(nt-1),k]
            IF(logarithms[k] EQ 1) THEN laggingVar = alog10(laggingVar)
            mdot = vsMdot[nt0:(nt-1),0,0,j]/normalizations[nt0:(nt-1),0]
            IF(logarithms[0] EQ 1) THEN mdot = alog10(mdot)
            crossCorrelations[0,*,k,j] = c_correlate(mdot,laggingVar,indgen(nt-1-nt0))
            autoCorrelations[0,*,k,j] = c_correlate(laggingVar,laggingVar,indgen(nt-1-nt0))
        ENDFOR

    ENDFOR

    wh = where(crossCorrelations NE crossCorrelations, ct)
    crossCorrelations[wh] = -1.1
    wh = where(autoCorrelations NE autoCorrelations, ct)
    autoCorrelations[wh] = -1.1

    ccVsTime = dblarr(1,nt-1-nt0,nvars+1,nmodels)
    acVsTime = dblarr(1,nt-1-nt0,nvars+1,nmodels)
    ccVsTime[0,*,1:nvars,*] = crossCorrelations[*,*,*,*]
    acVsTime[0,*,1:nvars,*] = autoCorrelations[*,*,*,*]
    FOR j=0, nmodels-1 DO BEGIN
        ccVsTime[0,*,0,j] = time[1:(nt-1-nt0)]
        acVsTime[0,*,0,j] = time[1:(nt-1-nt0)]
    ENDFOR

    qualifier = ""
    IF(normalize EQ 1) THEN qualifier="normalized "
    simpleMovie, crossCorrelations[0,1:nt-2-nt0,*,*], [0.0], names, colors, intarr(nmodels)+ls, intarr(nvars), name+"_xcc", 5, axisLabels="x-corr "+qualifier+labelsIn, whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses
    simpleMovie, autoCorrelations[0,1:nt-2-nt0,*,*],  [0.0], names, colors, intarr(nmodels)+ls, intarr(nvars), name+"_xac", 5, axisLabels="autocorr "+qualifier+labelsIn, whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses


    simpleMovie, ccVsTime[0,1:nt-2-nt0,*,*], [0.0], ['time',names], colors, intarr(nmodels)+ls, [1,intarr(nvars)], name+"_cc", 5, axisLabels=['lag time',"x-corr "+qualifier+labelsIn], whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses
    simpleMovie, acVsTime[0,1:nt-2-nt0,*,*], [0.0], ['time',names], colors, intarr(nmodels)+ls, [1,intarr(nvars)], name+"_ac", 5, axisLabels=['lag time',"autocorr "+qualifier+labelsIn], whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses

    ;; loop over variables to make plots
;   FOR k=0,n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
;       ;; for each lag time, compute the percentiles
;       FOR tt=0,n_elements(crossCorrelations[*,0,0])-1 DO BEGIN
;           ccPercentiles[tt,k,*] = percentiles(crossCorrelations[tt,k,*], percentileList, 0)
;           acPercentiles[tt,k,*] = percentiles( autoCorrelations[tt,k,*], percentileList, 0)
;       ENDFOR
;       setct,1,0,2
;       FIGUREINIT,(name+"_cc_"+names[k]),sv,2,2
;       PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="Lag Time (Ga)",YTITLE="x-corr: "+labels[0]+" with "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
;       FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO $
;           OPLOT, time[1:(nt-1-nt0)],crossCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
;       FOR j=0, n_elements(percentileList)-1 DO $
;           OPLOT, time[1:(nt-1-nt0)],ccPercentiles[*,k,j],COLOR=colors[0,
;       
;       FIGURECLEAN,(name+"_cc_"+names[k]),sv
;
;       FIGUREINIT,(name+"_ac_"+names[k]),sv,2,2
;       PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="Lag Time (Ga)",YTITLE="autocorr: "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
;       FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
;            OPLOT, time[1:(nt-1-nt0)],autoCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
;       ENDFOR
;       FIGURECLEAN,(name+"_ac_"+names[k]),sv
;
;       FIGUREINIT,(name+"_xcc_"+names[k]),sv,2,2
;       PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[-1.0,1.0],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="mdot autocorrelation",YTITLE="x-corr: "+labels[0]+" with "+labels[k],CHARSIZE=cs,CHARTHICK=ct
;       FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
;           OPLOT, crossCorrelations[*,0,j] , crossCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
;       ENDFOR
;       FIGURECLEAN,(name+"_xcc_"+names[k]),sv
;
;       FIGUREINIT,(name+"_xac_"+names[k]),sv,2,2
;       PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[-1.0,1.0],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="mdot autocorrelation",YTITLE="autocorr: "+labels[k],CHARSIZE=cs,CHARTHICK=ct
;       FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
;           OPLOT, crossCorrelations[*,0,j] , autoCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
;       ENDFOR
;       FIGURECLEAN,(name+"_xac_"+names[k]),sv
;
;
;   ENDFOR
    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)/3.0
END


