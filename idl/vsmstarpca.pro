; Assume vsMdot is the same 4-d matrix where as constructed in variability3:
; 1- which time
; 2- which x-position
; 3- which variable
; 4- which model
PRO vsMstarPCA,vsMstar,whichFrames=whichFrames,log=log,names=names,colors=colors,expName2=expName2,labels=labels,texLabels=texLabels,thicknesses=thicknesses,svSinglePlot=svSinglePlot,replacementText=replacementText
    IF(n_elements(whichFrames) EQ 0) THEN whichFrames = indgen(n_elements(vsMstar[*,0,0,0]))
    IF(n_elements(log) EQ 0) THEN log=intarr(n_elements(vsMstar[0,0,*,0]))
    nmodels = n_elements(vsMstar[0,0,0,*])
    ; these are indices of GetTimeLabels (indexed from 1!)
    ;;; M_*, SFR (disk only), mdot, Z, fg, 
    ;;; fgmol, colTrans, Mh, r25d, sigMax, 
    ;;; (Mhzs) 
;    whichVars = [3,2,18,16,21,44,52,53,54,56]
    whichVars = [1,2,17,15,20,43,51,52,53,55]
    nvars = n_elements(whichVars)
    theArrows = dblarr(n_elements(whichFrames),nvars-1,4,6) ;; 3ish frames, nvars plots, 4 pts per arrow, 6 principal components

    ;; First let's construct the matrix upon which we would like to do PCA.
    ;; We want an m by n matrix. m is a variable, n is a model.
    FOR frameI = 0, n_elements(whichFrames)-1 DO BEGIN
        frame = whichFrames[frameI]
        snapshotPCA, vsMstar , whichVars, log, frame,frameI, means, variances, eigenvect, proj_obj, proj_atr, theArrows
    ENDFOR

    ;; STEP ONE: plot the projections (theArrows, stored from above) of each eigenvector
    ;; onto the 2d-subspaces of M_* vs <other quantities>.
    ;; The projections are of the eigenvectors normalized to unit length after scaling (this may be incorrect)
    simpleMovie, vsMstar[*,*,whichVars-1,*],names[whichVars-1], $
        colors,colors*0,log[whichVars-1],expName2+"_vsMstarPCA", $
        5,PSYM=1,prev=0,axisLabels=labels[whichvars-1], $
        texLabels=texLabels[whichvars-1],whichFrames=whichFrames, $
        NIndVarBins=5,horizontal=1,thicknesses=thicknesses,svSinglePlot=svSinglePlot, $
        replacementText=replacementText, additionalLineColors=[0,1,2,3,4,5,6], $
        additionalLines=theArrows,additionalLineLabels=["PC1","PC2","PC3","PC4","PC5","PC6"]



    ;; STEP TWO: plot each model in PCA-space. This requires setting up the following matrix:
    ;; time/redshift, -, PC, model  

END





