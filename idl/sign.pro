FUNCTION sign,array
    absArray = abs(array)
    wh = where(absArray GT 0, ct) 
    signArray = absArray*0+1
    IF(ct GT 0) THEN BEGIN
        signArray[wh] = array[wh]/absArray[wh]
    ENDIF

    RETURN,signArray

END




