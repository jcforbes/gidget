; Take vsMstar as an input.

PRO deltams,vsMstar,percentileList,binned,deltas
    ; Loop over time.
    ; At each time, 1) find main sequence. 2) for each model, compute deltaMS, 3) bin it up
    nt = n_elements(vsMstar[*,0,0,0])
    nmodels = n_elements(vsMstar[0,0,0,*])
    deltas = dblarr(nt,nmodels)
    binned = intarr(nt,nmodels)
    FOR ti=0, nt-1 DO BEGIN
        theIndex = 14-2 ; incl bulge   2-1 ; disk SFR only
        ms = dblarr(nmodels)
        logMstar = alog10(vsMstar[ti,0,0,*])
        logSFR = alog10(vsMstar[ti,0,theIndex,*])
        fitparams = linfit(logMstar,logSFR, yfit = ms)
        deltas[ti,*] = logSFR - ms
        width = SQRT(VARIANCE(deltas[ti,*]))

        percens = percentileList[*]*0.0
        percens = PERCENTILES(deltas[ti,*], percentileList, 0)
        percens=[min(deltas[ti,*])-1,percens, max(deltas[ti,*])+1]
        counts = [0]
        FOR i=1, n_elements(percens)-1 DO BEGIN
            wh = WHERE(deltas[ti,*] LT percens[i] and deltas[ti,*] GE percens[i-1], ct)
            IF(ct GT 0) THEN binned[ti,wh] = i-1 ; ELSE message,"Problem in deltams"
            counts=[counts,ct]
        ENDFOR

    ENDFOR


END
