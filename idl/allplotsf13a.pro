; Make all of the plots which will be in F13a.
; numbers allows you to make each set of plots one-at-a-time
; since the poorly-written code may require a fresh IDL session each time.

PRO allplotsf13a,numbers
    wh = WHERE(numbers EQ 0,ct)
    IF(ct GT 0) THEN variability3,['rg01','rg20b']
    wh = WHERE(numbers EQ 1,ct)
    IF(ct GT 0) THEN variability3,['rg01','rg36d','rg36e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 2,ct)
    IF(ct GT 0) THEN variability3,['rg01','rg36e'],annotations=['Fiducial','No GI']
    wh = WHERE(numbers EQ 3,ct)
    IF(ct GT 0) THEN postagestamps
    wh = WHERE(numbers EQ 4,ct)
    IF(ct GT 0) THEN balanceplots,['rg01/rg01','rg36d/rg36d','rg36e/rg36e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 5,ct)
    IF(ct GT 0) THEN variability3,['rg01','rg36d'],annotations=['Fiducial','No SF']

END
