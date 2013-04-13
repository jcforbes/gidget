; Make all of the plots which will be in F13a.
; numbers allows you to make each set of plots one-at-a-time
; since the poorly-written code may require a fresh IDL session each time.

PRO allplotsf13a,numbers
    wh = WHERE(numbers EQ 0,ct)
    IF(ct GT 0) THEN variability3,['rf01','rf20b']
    wh = WHERE(numbers EQ 1,ct)
    IF(ct GT 0) THEN variability3,['rf01','rf36d','rf36e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 2,ct)
    IF(ct GT 0) THEN variability3,['rf01','rf36e'],annotations=['Fiducial','No GI']
    wh = WHERE(numbers EQ 3,ct)
    IF(ct GT 0) THEN postagestamps
    wh = WHERE(numbers EQ 4,ct)
    IF(ct GT 0) THEN balanceplots,['rf01/rf01','rf36d/rf36d','rf36e/rf36e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 5,ct)
    IF(ct GT 0) THEN variability3,['rf01','rf36d'],annotations=['Fiducial','No SF']

END
