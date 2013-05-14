FUNCTION posonly,arr
    wh = where(arr LT 0, ct)
    IF(ct GT 0) THEN arr[wh]=0.0
    ; The following line is a quick diagnostic. In the context of bulgeExcess,
    ; we hope that ct is smallish, i.e. the exp fit to the column density is
    ; never an over-estimate of the line, except at large radii. So, we hope that
    ; ct is small and wh is only at large radii.
;    print, "PosOnly count, min(wh): ",ct," ",MIN(wh)
    return,arr

END

FUNCTION bulgeExcess,r,col,scaleRadius
    ; first step: locate the index closest to f*scale radius
    doneFlag = 0
    excess=-1.0
    whichTry=0
    nx=n_elements(r)
    maxTries = 200
    FOR try=0,maxTries DO BEGIN
        IF(doneFlag EQ 0) THEN BEGIN
            f = 1.5 ; - float(try)*.01
            theMin = MIN(abs(r-scaleRadius*f), ind)
            failflag = 0
            failLoc=-1
            logCol = alog10(col)
            slope1 =  (logcol[ind-1]-logcol[ind])/(r[ind-1]-r[ind]) * float(maxTries/2-try)/(maxTries/2)
;            slope2 =  (logcol[ind-2]-logcol[ind])/(r[ind-2]-r[ind])
;            slope3 =  (logcol[ind-3]-logcol[ind])/(r[ind-3]-r[ind])
;            slope4 =  (logcol[ind-4]-logcol[ind])/(r[ind-4]-r[ind])
;            slope5 =  (logcol[ind-5]-logcol[ind])/(r[ind-5]-r[ind])

            yGuess1 = logcol[ind] + (r[*] - r[ind])*slope1
;            yGuess2 = logcol[ind] + (r[*] - r[ind])*slope2
;            yGuess3 = logcol[ind] + (r[*] - r[ind])*slope3
;            yGuess4 = logcol[ind] + (r[*] - r[ind])*slope4
;            yGuess5 = logcol[ind] + (r[*] - r[ind])*slope5

            fail1 = PosOnly( yGuess1[0:ind] - logcol[0:ind] )
;            fail2 = PosOnly( yGuess2[0:ind] - logcol[0:ind] )
;            fail3 = PosOnly( yGuess3[0:ind] - logcol[0:ind] )
;            fail4 = PosOnly( yGuess4[0:ind] - logcol[0:ind] )
;            fail5 = PosOnly( yGuess5[0:ind] - logcol[0:ind] )

;            PRINT, "failure stats: current try = ",try, " index = ",ind
;            PRINT, "fail1: ",MAX(fail1,failloc1), fail1[0]
;            PRINT, "fail2: ",MAX(fail2,failloc1), fail2[0]
;            PRINT, "fail3: ",MAX(fail3,failloc1), fail3[0]
;            PRINT, "fail4: ",MAX(fail4,failloc1), fail4[0]
;            PRINT, "fail5: ",MAX(fail5,failloc1), fail5[0]

            failFlag = (MAX(fail1) GT 0.0); 1.0e-2) 

;            print,"Try, fail: ",try,failflag,failloc1
            IF(failFlag EQ 0) THEN BEGIN
                theFit = 10.0^(logCol[ind] + (r[*]-r[ind])*slope1)
                excess1 = TOTAL(r[*]*r[*]*PosOnly(col[*] - 10.0^(logcol[ind] + (r[*] - r[ind])*slope1)) )
                excR = TOTAL(r[ind:nx-1]*r[ind:nx-1]*PosOnly(col[ind:nx-1] - theFit[ind:nx-1]))
                excL = TOTAL(r[0:ind-1]*r[0:ind-1]*PosOnly(col[0:ind-1] - theFit[0:ind-1]))
                stellar = TOTAL(r[*]*r[*]*col[*])
                fitTotal = TOTAL(r[*]*r[*]* theFit *(1.0+sign(col[*] - theFit+1.0e-8))/2.0  $
                    + r[*]*r[*]*col[*] * (1.0+sign(theFit-col[*]-1.0e-8))/2.0)
;                excess2 = TOTAL(r[*]*r[*]*PosOnly(col[*] - 10.0^(logcol[ind] + (r[*] - r[ind])*slope2)) )
;                excess3 = TOTAL(r[*]*r[*]*PosOnly(col[*] - 10.0^(logcol[ind] + (r[*] - r[ind])*slope3)) )
;                excess4 = TOTAL(r[*]*r[*]*PosOnly(col[*] - 10.0^(logcol[ind] + (r[*] - r[ind])*slope4)) )
;                excess5 = TOTAL(r[*]*r[*]*PosOnly(col[*] - 10.0^(logcol[ind] + (r[*] - r[ind])*slope5)) )
                doneFlag = 1
                whichTry=try
                ;IF(excL LT excR and excess1 GT .1*stellar) THEN stop
            ENDIF
            ;IF(doneFlag EQ 1) THEN print,"Try: ",try
            
        ENDIF
    ENDFOR
;    PRINT, "bulgeExcess: try = ",whichtry
;    PRINT, "bulgeExcess: excesses: ",excess1,excess2,excess3,excess4,excess5
    ;; Assertion 1:  fitTotal + excess1 = stellar
;    PRINT, "Assertion 1: ", fitTotal, excess1, fitTotal+excess1, stellar
;    PRINT, "BT guess in bulgeExcess: ",excess1/stellar
;    PRINT, "right-hand excess vs total: ",excL,excR,excL+excR,excess1,(excess1+excess2+excess3+excess4+excess5)/5.0
;    PRINT, "right-hand excess vs total: ",excL/(excL+excR),excess1/(excL+excR),excess1/((excess1+excess2+excess3+excess4+excess5)/5.0)
;    RETURN,excess1
   ; stop
   IF(doneFlag EQ 0) THEN BEGIN
      PRINT,"bulgeExcess did not converge"
      stop
      return,0.0
   ENDIF
    RETURN,excL
;    RETURN,-1 ; hope you never get here, but no guarantee.


END

