; Make all of the plots which will be in F13a.
; numbers allows you to make each set of plots one-at-a-time
; since the poorly-written code may require a fresh IDL session each time.

PRO allplotsf13a,numbers
    fid = 'rg01'
    stoch='rg69b';'rg20b'
    alt='rg36'


    wh = WHERE(numbers EQ 0,ct)
    IF(ct GT 0) THEN variability3,[fid,stoch]
    wh = WHERE(numbers EQ 1,ct)
    IF(ct GT 0) THEN variability3,[fid,alt+'d',alt+'e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 2,ct)
    IF(ct GT 0) THEN variability3,[fid,alt+'e'],annotations=['Fiducial','No GI']
    wh = WHERE(numbers EQ 3,ct)
    IF(ct GT 0) THEN postagestamps
    wh = WHERE(numbers EQ 4,ct)
    IF(ct GT 0) THEN balanceplots,[fid+'/'+fid,alt+'d/'+alt+'d',alt+'e/'+alt+'e'],annotations=['Fiducial','No SF','No GI']
    wh = WHERE(numbers EQ 5,ct)
    IF(ct GT 0) THEN variability3,[fid,alt+'d'],annotations=['Fiducial','No SF']

END
