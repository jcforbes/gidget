
;; Given an array of mdot's vs. other variables, compute the cross-correlations between mdot and the others
;; and also the autocorrelation of each variable.
PRO ComputeCorrelations,vsmdot,colors,time,labelsIn,names,name,sv=sv,nt0=nt0,thicknesses=thicknesses,logarithms=logarithms,normalize=normalize,percentileList=percentileList,texLabelsIn=texLabelsIn
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
    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)
    ls=0
    nmodels = n_elements(vsMdot[0,0,0,*])
    nvars = n_elements(vsMdot[0,0,*,0])
    IF(nmodels GT 50) THEN ls =1 


    lowestColor = MIN(colors)
    highestColor = MAX(colors)
    ncolors = highestcolor-lowestcolor+1
    normalizations = dblarr(n_elements(vsMdot[*,0,0,0]), n_elements(vsMdot[0,0,*,0]),ncolors)
    ; An array to store the # of times each color appears in the list.
    listOfColors=intarr(highestColor-lowestColor+1)+lowestColor
    FOR aColor=lowestColor, highestColor DO BEGIN
        FOR ti=0, n_elements(vsMdot[*,0,0,0])-1 DO BEGIN
            wh = WHERE(colors[ti,*] EQ aColor, count)
            listOfColors[aColor-lowestColor] = count
            IF(count NE 0) THEN BEGIN
                FOR k=0, n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
                    IF(normalize EQ 1) THEN normalizations[ti,k,acolor-lowestColor] = median(vsMdot[ti,0,k,wh]) ELSE normalizations[ti,k,wh]=1.0
                ENDFOR
            ENDIF
        ENDFOR
    ENDFOR

    ;stop,"Normalizations produced"

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
            FOR aColor=lowestColor, highestColor DO BEGIN
;                wh = WHERE(colors[0,*] EQ aColor, count)
;                IF(count GT 0) THEN BEGIN
                IF(colors[0,j] EQ aColor) THEN BEGIN
                    laggingVar = vsMdot[nt0:(nt-1),0,k,j];/normalizations[nt0:(nt-1),k,aColor]
		    ;;; Note that if normalize is false (0), then all values of normalizations[*,*,*] will be 1, in which case we're either subtracting zero or one from every variable, which will make no difference in the correlations.
                    IF(logarithms[k] EQ 1) THEN laggingVar = alog10(laggingVar/normalizations[nt0:(nt-1),k,aColor-lowestColor]) ELSE laggingVar = laggingVar - normalizations[nt0:(nt-1),k,aColor-lowestColor]
                    mdot = vsMdot[nt0:(nt-1),0,0,j];/normalizations[nt0:(nt-1),0,aColor]
                    IF(logarithms[0] EQ 1) THEN mdot = alog10(mdot/normalizations[nt0:(nt-1),0,aColor-lowestColor]) ELSE mdot=mdot - normalizations[nt0:(nt-1),0,aColor-lowestColor]
                    crossCorrelations[0,*,k,j] = c_correlate(mdot,laggingVar,indgen(nt-1-nt0))
                    autoCorrelations[0,*,k,j] = c_correlate(laggingVar,laggingVar,indgen(nt-1-nt0))
                    ;IF(j EQ 1) THEN STOP, "in the loop"
                ENDIF
            ENDFOR
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
;    simpleMovie, crossCorrelations[0,1:nt-2-nt0,*,*], names, colors, intarr(nmodels)+ls, intarr(nvars), name+"_xcc", 5, axisLabels="x-corr "+qualifier+labelsIn, whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses
;    simpleMovie, autoCorrelations[0,1:nt-2-nt0,*,*],  names, colors, intarr(nmodels)+ls, intarr(nvars), name+"_xac", 5, axisLabels="autocorr "+qualifier+labelsIn, whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses


    simpleMovie, ccVsTime[0,1:nt-2-nt0,*,*], ['time',names], colors, $
        intarr(nmodels)+ls, [1,intarr(nvars)], name+"_cc", 5, $
        axisLabels=['lag time',"x-corr "+qualifier+labelsIn], $
        whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses, $
        svSinglePlot=2,percentileList=[.025,.16,.5,.84,.975],fill=1,$
        texLabels=['$t_{lag}$ (Gyr)','$P_{t_{lag}}(\dot{M}_{ext},$ '+texLabelsIn+'$)$']
    simpleMovie, acVsTime[0,1:nt-2-nt0,*,*], ['time',names], colors, $
        intarr(nmodels)+ls, [1,intarr(nvars)], name+"_ac", 5, $
        axisLabels=['lag time',"autocorr "+qualifier+labelsIn], $
        whichFrames=[0], NIndVarBins=20,thicknesses=thicknesses, $
        svSinglePlot=2,percentileList=[.025,.16,.5,.84,.975],fill=1,$
        texLabels=['$t_{lag}$ (Gyr)','$P_{t_{lag}}(\dot{M}_{ext},$ '+texLabelsIn+'$)$']

    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)/3.0
END


