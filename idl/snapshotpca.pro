;; Find the principal components of a particular snapshot (in redshift) of a population of galaxies.
;; Designed to be called by vsMstarPCA.
PRO snapshotPCA, vsMstar , whichVars, log, frame,frameI, means, variances, eigenvect, proj_obj, proj_atr, theArrows
    nvars = n_elements(whichVars)
    nmodels = n_elements(vsMstar[0,0,0,*])

    dataMatrix = dblarr(nvars + 1, nmodels)
    FOR j=0, nmodels-1 DO BEGIN
        dataMatrix[0:nvars-1,j] = vsMstar[frame,0,whichVars-1,j]
        wh = WHERE(log[whichVars-1] EQ 1, ct) 
        IF(ct GT 0) THEN dataMatrix[wh,j] = alog10(dataMatrix[wh,j])
        dataMatrix[nvars,j] = alog10(vsMstar[0,0,51,j])

    ENDFOR

    theMeans = dblarr(nvars+1)
    theVariances = dblarr(nvars+1)
    FOR v=0, nvars DO theMeans[v] = MEAN(dataMatrix[v,*])
    FOR v=0, nvars DO theVariances[v] = VARIANCE(dataMatrix[v,*])


    ;; If the variance of one quantity is zero, bad things will happen
    ;; In particular, when we normalize the dataMatrix by the standard deviation, we'll get NaN.
    ;; So, here we artificially put in some variance, hopefully such that it will not affect the results.
    ;; The idea is that if a quantity with truly zero variance has some small non-zero variance, then
    ;; it will allow the procedure to go forward and this quantity won't appear in the first few
    ;; principal components.
    wh = WHERE(theVariances LE 0, ct, complement=comp)
    IF(ct GT 0) THEN BEGIN
        FOR k=0, ct-1 DO BEGIN
            minvar = MIN(theVariances[comp])
            dataMatrix[wh[k],*] = theMeans[wh[k]] + randomn(seed,nmodels)*minvar*1.0e-5
            theVariances[wh[k]] = VARIANCE(dataMatrix[wh[k],*])

        ENDFOR
    ENDIF

    FOR v=0, nvars DO dataMatrix[v,*] = (dataMatrix[v,*] - theMeans[v])/(sqrt(theVariances[v]*double(nmodels-1.0)))

    pca, transpose(dataMatrix), eigenval, eigenvect, percentages, proj_obj, proj_atr,matrix=matrix,/covariance, /silent
	
    FOR v=0, nvars DO BEGIN
        dataMatrix[v,*] = dataMatrix[v,*]*sqrt(theVariances[v]*double(nmodels-1.0))+theMeans[v]
        eigenvect[*,v] = eigenvect[*,v] * sqrt(theVariances[v])
        
    ENDFOR


    FOR pcaI=0,5 DO eigenvect[pcaI,*] =10.0* eigenvect[pcaI,*] / SQRT(TOTAL(eigenvect[pcaI,*]*eigenvect[pcaI,*]))
    FOR varI=0,nvars-2 DO BEGIN
        theArrows[frameI,varI,0,*] = 10.0 ^ (theMeans[0]-.5*sqrt(theVariances[0]))
        theArrows[frameI,varI,1,*] = theMeans[varI+1]+sqrt(theVariances[varI+1])

        theArrows[frameI,varI,2,*] = theArrows[frameI,varI,0,*]*10.0 ^ ( eigenvect[0:5,0])
        theArrows[frameI,varI,3,*] = theArrows[frameI,varI,1,*] + eigenvect[0:5,varI+1] 

        IF(log[whichVars[varI+1]-1] EQ 1) THEN theArrows[frameI,varI,3,*] = 10.0^theArrows[frameI,varI,3,*]
        IF(log[whichVars[varI+1]-1] EQ 1) THEN theArrows[frameI,varI,1,*] = 10.0^theArrows[frameI,varI,1,*]
    ENDFOR





END

