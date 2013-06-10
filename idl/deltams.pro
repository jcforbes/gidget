; Take vsMstar as an input.

PRO deltams,vsMstar,binnedMS,binnedAvgAccr,deltas,deltaAccr,deltaAvgAccr,widths,widthsAccr,widthsAvgAccr,msText,NMassBins=NMassBins,NDeltaBins=NDeltaBins
    if(n_elements(nmassbins) EQ 0) THEN nmassbins=1
    if(n_elements(ndeltabins) EQ 0) THEN ndeltabins=1
    ; Loop over time.
    ; At each time, 1) find main sequence. 2) for each model, compute deltaMS, 3) bin it up
    nt = n_elements(vsMstar[*,0,0,0])
    nmodels = n_elements(vsMstar[0,0,0,*])
    deltas = dblarr(nt,nmodels)
    deltaAccr = dblarr(nt,nmodels)
    deltaAvgAccr = dblarr(nt,nmodels)
    binnedMS = intarr(nt,nmodels)
    binnedAvgAccr = intarr(nt,nmodels)
    widths = dblarr(nt)
    widthsAccr=dblarr(nt)
    widthsAvgAccr=dblarr(nt)
    msText = strarr(nt)
    FOR ti=0, nt-1 DO BEGIN
        theIndex = 2-1 ; 14-2 ; incl bulge   2-1 ; disk SFR only
        ms = dblarr(nmodels)
        logMstar = alog10(vsMstar[ti,0,0,*])
        logSFR = alog10(vsMstar[ti,0,theIndex,*])
        fitparams = linfit(logMstar,logSFR, yfit = ms)
        msText[ti] = 'logSFR = '+string(fitparams[1],Format='(D0.2)')+' logM* + ' $
            +string(fitparams[0],Format='(D0.2)')
        deltas[ti,*] = logSFR - ms
        width = SQRT(VARIANCE(deltas[ti,*]))
        logAccr = alog10(vsMstar[ti,0,17-1,*])
        nAvg = ti-MAX([0,ti-20])
        avgAccr = vsMstar[ti,0,17-1,*]/MAX([float(nAvg),1.0])
        FOR tii=MAX([0,ti-20]), ti-1 DO BEGIN
            avgAccr += vsMstar[tii,0,17-1,*]/float(nAvg)
        ENDFOR

        fitparamsAcc = linfit(logMstar, logAccr, yfit = msAccr)
        deltaAccr[ti,*] = logAccr-msAccr
        widthAccr = SQRT(VARIANCE( deltaAccr[ti,*] ))

        widthsAccr[ti] = widthAccr
        widths[ti]=width

        deltaAvgAccr[ti,*] = alog10(avgAccr) - msAccr
        widthsAvgAccr[ti] = SQRT(Variance( deltaAvgAccr[ti,*]))

        percenMstar = PERCENTILES(logMstar[*], [findgen(NMassBins+1)/float(NMassBins)], 0)
        percenMstar[0] = MIN(logMstar)-.1
        percenMstar[NMassBins] = MAX(logMstar)+.1

        percenAvgAccr = PERCENTILES(deltaAvgAccr[ti,*], [findgen(NDeltaBins+1)/float(NDeltaBins)], 0)
        percenAvgAccr[0] = MIN(deltaAvgAccr[ti,*])-.1
        percenAvgAccr[NDeltaBins] = MAX(deltaAvgAccr[ti,*])+.1

        percenMS= PERCENTILES(deltas[ti,*], [findgen(NDeltaBins+1)/float(NDeltaBins)], 0)
        percenMS[NDeltaBins] = MAX(deltas[ti,*])+.1
        percenMS[0] = MIN(deltas[ti,*])-.1
        
        
        counts = [0]
        FOR k=1, n_elements(percenMstar)-1 DO BEGIN
            FOR i=1, n_elements(percenMS)-1 DO BEGIN
                wh = WHERE(deltas[ti,*] LT percenms[i] and deltas[ti,*] GE percenms[i-1] $
                    and logMstar LT percenMstar[k] and logMstar GE percenMstar[k-1], ct)
                IF(ct GT 0) THEN binnedMS[ti,wh] = i-1 + (k-1)*(n_elements(percenms)-1)  ; ELSE message,"Problem in deltams"
                counts=[counts,ct]
            ENDFOR
        ENDFOR


        countsAccr = [0]
        FOR k=1, n_elements(percenMstar)-1 DO BEGIN
            FOR i=1, n_elements(percenAvgAccr)-1 DO BEGIN
                wh = WHERE(deltaAvgAccr[ti,*] LT percenAvgAccr[i] and deltaAvgAccr[ti,*] GE percenAvgAccr[i-1] $
                    and logMstar LT percenMstar[k] and logMstar GE percenMstar[k-1], ct)
                IF(ct GT 0) THEN binnedAvgAccr[ti,wh] = i-1 + (k-1)*(n_elements(percenAvgAccr)-1)  ; ELSE message,"Problem in deltams"
                countsAccr=[countsAccr,ct]
            ENDFOR
        ENDFOR


    ENDFOR


END
