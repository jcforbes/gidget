; Take as input the unsorted data hyper-cube. The first variable will be the radius of each object in the cube.
; Some functions can handle this, i.e. those which simply plot radius vs. whatever
; But others (or rather, the often the same function with a different set of options), specifically
; those which attempt to examine the distrubtion, e.g. the median and central 68%, will do something completely
; and totally wrong, namely treat each y-value at a fixed index as if it comes from the same x-value.
; In order for that to be the correct behavior, the data must be interpolated
PRO resample, data
    nx = n_elements(data[0,*,0,0]) ; number of radial data points per model
    nt = n_elements(data[*,0,0,0]) ; number of timesteps
    nm = n_elements(data[0,0,0,*]) ; number of models
    nv = n_elements(data[0,0,*,0]) ; number of variables

    ; the range over which all models actually have data
    restrictive = [MAX(data[*,0,0,*]),MIN(data[*,nx-1,0,*])]
    IF(restrictive[0] GT restrictive[1]) THEN message,"No points qualify for resampling!"
    ; the full range over which we have data:
    full = [MIN(data[*,0,0,*]),MAX(data[*,nx-1,0,*])]
    IF(full[0] EQ restrictive[0] and full[1] EQ restrictive[1]) THEN RETURN ; no need to resample

    newGrid = 10.0^(alog10(restrictive[0]) + findgen(nx)/(nx-1) * alog10(restrictive[1]/restrictive[0]))

    print, "Resampling to a universal grid; reducing the x-range from ",full," to ",restrictive

    for k=1, nv-1 DO BEGIN
        for ti=0, nt-1 DO BEGIN
            for j=0, nm-1 DO BEGIN
                data[ti,*,k,j] = INTERPOL(data[ti,*,k,j],data[ti,*,0,j], newGrid )
            ENDFOR
        ENDFOR
    ENDFOR

    FOR ti=0, nt-1 DO BEGIN
        FOR j=0, nm-1 DO BEGIN
            data[ti,*,0,j] = newGrid
        ENDFOR
    ENDFOR

END

