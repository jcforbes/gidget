;; Very simple.
;; Return the values of theArray at the given percentiles.
;; In particular, theArray is a 1d vector of data.
;; percentileArray is the list of percentiles we'd like to compute,
;;    e.g. .025, .975 for "2-sigma"
;; and finally log tells us whether this data naturally "prefers" a 
;;    logarithmic spacing- this will make our sorting more accurate.
;; WARNING: using w/ the log option will significantly alter the given matrix!
;;    If you need it back unaltered, use a copy when calling the function.
FUNCTION percentiles, theArray, percentileArray, log
    theResults = dblarr(n_elements(percentileArray))

    theTotal = n_elements(theArray)

    wh = where(theArray GT 1.0e30, count)
    IF(count GT 0) THEN theArray[wh] = 1.0e30


    IF(log EQ 1) THEN BEGIN
        wh = where(theArray GT 0.0, count, complement=comp)
        IF(count GT 0) THEN theValidMin = MIN(theArray[wh]) ELSE theValidMin=1.
        IF(count NE theTotal) THEN theArray[comp] = theValidMin
        theArray = alog10(temporary(theArray))

    ENDIF

    thenbins = 300
    themin = min(theArray,/nan)
    spacing = (max(theArray)-theMin)/(thenbins-1.)
    theLinearHist = histogram(theArray,min=theMin,binsize=spacing,nbins=theNbins)
    cumulativeLinear = TOTAL(theLinearHist, /CUMULATIVE)

    FOR i=0, n_elements(theResults)-1 DO BEGIN
    	wh = WHERE(cumulativeLinear GE theTotal*percentileArray[i],count)

        IF(count GT 0) THEN theResults[i] = wh[0]*spacing+theMin+spacing/2.0 ;$
            ;ELSE stop,"Problem in percentiles.pro; this is usually caused by a NaN in the data."
        
    ENDFOR
    
    IF(log EQ 1) THEN theResults=10.0^temporary(theResults)
    
    RETURN,theResults

END

