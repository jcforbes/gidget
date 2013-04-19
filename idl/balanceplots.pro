;; We'd like to produce some postage stamps.
;; Use the results of variability3, and stitch together
;; N x M set of plots at a particular time.
;; I think this will just involve constructing the 
;; appropriate matrix to send to simpleMovie


PRO balanceplotsCol,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations

    compile_opt idl2
    set_plot,'x'
    device,decomposed=0
    device,retain=2

    IF(n_elements(annotations) EQ 0) THEN BEGIN
	annotations=ListOfModels[*]
	FOR i=0,n_elements(annotations)-1 DO annotations[i]=ExtractSlashName(annotations[i],1)
    ENDIF
    IF(n_elements(texAnnotations) EQ 0) THEN texAnnotations=annotations[*]

    whichRedshifts=[1,51,201]
    zs=['z=2.00','z=1.01','z=0.00']
    columns = n_elements(ListOfModels)
    rows = n_elements(zs)
    setct,1,0,0
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2


    dummy = GetLabels(keys=[0,31,35,36,17], names=names, $
        labels=labels, tex=texLabels, base=vars, $
        log=logs, offset=offsets,ranges=ranges)



    ;; collect data from disk.
    ptrListOfModels = ptrarr(n_elements(ListOfModels), /allocate_heap)
    FOR k=0, n_elements(ListOfModels)-1 DO BEGIN

        model = readOutput(ListOfModels[k])
        theData = dblarr(n_elements(whichRedshifts),model.nx,n_elements(vars))
        theData[*,*,*] = model.dataCube[whichRedshifts,*,vars+offsets*model.ncolstep-1]
        ;; at this point, per the keys requested of GetLabels, we have
        ;; radius, SigTr/All, accr/All, SF/all, mdotDisk
        theProcessedData=dblarr(n_elements(whichRedshifts),model.nx,6)
        theProcessedData[*,*,0]=theData[*,*,0] ;; radius
        FOR xx=0,n_elements(theData[0,*,0])-1 DO BEGIN
            FOR wh=0,n_elements(theData[*,0,0])-1 DO BEGIN;WHERE(abs(theData[*,xx,4]) GT 0, ct)
                IF(abs(theData[wh,xx,4]) GT 0) THEN BEGIN
                    theProcessedData[wh,xx,1]=MIN([theData[wh,xx,1], 0])
                    theProcessedData[wh,xx,2]=MAX([theData[wh,xx,1], 0])
                    theProcessedData[wh,xx,3]=(theData[wh,xx,4]/abs(theData[wh,xx,4]) + 1)/2 *theProcessedData[wh,xx,2] + (theData[wh,xx,4]/abs(theData[wh,xx,4]) - 1)/(-2.0) * theProcessedData[wh,xx,1]

                ENDIF
            ENDFOR
        ENDFOR
        theProcessedData[*,*,4]=theData[*,*,2] + (theData[*,*,1] +abs(theData[*,*,1]))/2.0
        theProcessedData[*,*,5]=-theData[*,*,3] + (theData[*,*,1] -abs(theData[*,*,1]))/2.0

        stamp = {name:ListOfModels[k], theData:theProcessedData}
        *(ptrListOfModels[k]) = stamp

    ENDFOR

    colors=indgen(n_elements(theProcessedData[0,0,*]))
    colors=[0,1,2,5,3]
    labels=["r (kpc)","Fraction of dcoldt"]
    texLabels=["$r$ (kpc)","Share of $\partial \Sigma/\partial t$"]
    baseName="balance_col"

    textArr = intarr(n_elements(whichRedshifts), n_elements(ListOfModels),3)
    textArr[2,1,[0,1,2]]=1
    textArr[2,2,0]=1

    theText = strarr(n_elements(whichRedshifts), n_elements(ListOfModels),3)
    theText[2,1,0] = "Inward Transport"
    theText[2,1,1] = "Outward Transport"
    theText[2,1,2] = "Cosmological Accretion"
    theText[2,2,0] = "SF + galactic winds"

    textPos = dblarr(n_elements(whichRedshifts), n_elements(ListOfModels),3,2)
    textPos[2,1,1,0] = 19.0
    textPos[2,1,2,0] = 6.0
    textPos[2,1,0,0]= 2.0
    textPos[2,1,0,1] = -.7
    textPos[2,1,1,1] = -.4
    textPos[2,1,2,1] = 0.28
    textPos[2,2,0,0] = 4.0
    textPos[2,2,0,1] = -0.26

    textCol = intarr(n_elements(whichRedshifts), n_elements(ListOfModels),3)
    textCol[2,1,0] = 1




;    FOR a=0, 100 DO BEGIN
        makeBalancePlots, ptrListOfModels, whichRedshifts,zs,th,annotations,texAnnotations,svSinglePlot,cs,chth,logs,ranges,colors,labels,texLabels,baseName,fill=1,textArr=textArr,theText=theText,textPos=textPos
;	    stop
;    ENDFOR

END






PRO balanceplotsCol2,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations

    compile_opt idl2
    set_plot,'x'
    device,decomposed=0
    device,retain=2

    IF(n_elements(annotations) EQ 0) THEN BEGIN
	annotations=ListOfModels[*]
	FOR i=0,n_elements(annotations)-1 DO annotations[i]=ExtractSlashName(annotations[i],1)
    ENDIF
    IF(n_elements(texAnnotations) EQ 0) THEN texAnnotations=annotations[*]

    whichRedshifts=[1,51,201]
    zs=['z=2.00','z=1.01','z=0.00']
    columns = n_elements(ListOfModels)
    rows = n_elements(zs)
    setct,1,0,0
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2

    ; incoming, outgoing, accr, sfPerAccr
    dummy = GetLabels(keys=[0,45,46,25,23], names=names, $
        labels=labels, tex=texLabels, base=vars, $
        log=logs, offset=offsets,ranges=ranges)



    ;; collect data from disk.
    ptrListOfModels = ptrarr(n_elements(ListOfModels), /allocate_heap)
    FOR k=0, n_elements(ListOfModels)-1 DO BEGIN

        model = readOutput(ListOfModels[k])
        theData = dblarr(n_elements(whichRedshifts),model.nx,n_elements(vars))
        theData[*,*,*] = model.dataCube[whichRedshifts,*,vars+offsets*model.ncolstep-1]
        theData[*,*,3] = model.dataCube[whichRedshifts,*,30-1]
        theData[*,*,4] = model.dataCube[whichRedshifts,*,29-1]

        theProcessedData=dblarr(n_elements(whichRedshifts),model.nx,6)
        theProcessedData[*,*,0]=theData[*,*,0] ;; radius
        theNorm = TOTAL(abs(theData[*,*,1:4]),3)
        theProcessedData[*,*,1]=theData[*,*,1]/theNorm
        theProcessedData[*,*,2]=-theData[*,*,2]/theNorm
        theProcessedData[*,*,3]=(theData[*,*,1]+theData[*,*,3])/theNorm
        theProcessedData[*,*,4]=-(theData[*,*,2]+theData[*,*,4])/theNorm
        theProcessedData[*,*,5]=theData[*,*,1]*0

        stamp = {name:ListOfModels[k], theData:theProcessedData}
        *(ptrListOfModels[k]) = stamp

    ENDFOR

    colors=indgen(n_elements(theProcessedData[0,0,*]))
    colors=[0,1,2,5,3]
    labels=["r (kpc)","Fraction of dcoldt"]
    texLabels=["$r$ (kpc)","Share of $\partial \Sigma/\partial t$"]
    baseName="balance_col2"

;    FOR a=0, 100 DO BEGIN
        makeBalancePlots, ptrListOfModels, whichRedshifts,zs,th,annotations,texAnnotations,svSinglePlot,cs,chth,logs,ranges,colors,labels,texLabels,baseName,fill=1
;	    stop
;    ENDFOR

END








PRO balanceplotsSig,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations

    compile_opt idl2
    set_plot,'x'
    device,decomposed=0
    device,retain=2

    IF(n_elements(annotations) EQ 0) THEN BEGIN
	annotations=ListOfModels[*]
	FOR i=0,n_elements(annotations)-1 DO annotations[i]=ExtractSlashName(annotations[i],1)
    ENDIF
    IF(n_elements(texAnnotations) EQ 0) THEN texAnnotations=annotations[*]

    whichRedshifts=[1,51,201]
    zs=['z=2.00','z=1.01','z=0.00']
    columns = n_elements(ListOfModels)
    rows = n_elements(zs)
    setct,1,0,0
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2




    dummy = GetLabels(keys=[0,37,38,39,40], names=names, $
        labels=labels, tex=texLabels, base=vars, $
        log=logs, offset=offsets,ranges=ranges)



    ;; collect data from disk.
    ptrListOfModels = ptrarr(n_elements(ListOfModels), /allocate_heap)
    FOR k=0, n_elements(ListOfModels)-1 DO BEGIN

        model = readOutput(ListOfModels[k])
        theData = dblarr(n_elements(whichRedshifts),model.nx,n_elements(vars))
        theData[*,*,*] = model.dataCube[whichRedshifts,*,vars+offsets*model.ncolstep-1]
        ;; at this point, per the keys requested of GetLabels, we have
        ;; radius, dsigTr, dsigDdx, dsigHeat, dsigCool
        ;; Now we normalize things (besides the radius!)
;        thenorm = abs(theData[*,*,1]) + abs(theData[*,*,2]) + abs(theData[*,*,3])+abs(theData[*,*,4])
        thenorm = TOTAL(abs(theData[*,*,1:n_elements(vars)-1]),3)
        
        FOR i=1,n_elements(vars)-1 DO theData[*,*,i] = theData[*,*,i]/theNorm
        theProcessedData=dblarr(n_elements(whichRedshifts),model.nx,7)
        theProcessedData[*,*,0]=theData[*,*,0] ;; radius
        FOR xx=0,n_elements(theData[0,*,0]) -1 DO BEGIN
            FOR ti=0,n_elements(theData[*,0,0]) -1 DO BEGIN
                theProcessedData[ti,xx,1] = MIN([0,theData[ti,xx,2]]) ;; closest to zero- bound from below (2)
                theProcessedData[ti,xx,2] = MAX([0,theData[ti,xx,2]]) ;; closest to zero- bound from above (3)
                theProcessedData[ti,xx,3] = MAX([0,theProcessedData[ti,xx,2],theData[ti,xx,2]+theData[ti,xx,1]]) ; (4)
                theProcessedData[ti,xx,4] = MIN([0,theProcessedData[ti,xx,1],theData[ti,xx,2]+theData[ti,xx,1]]) ; (1)
                theProcessedData[ti,xx,5] = theProcessedData[ti,xx,4] + theData[ti,xx,4] ; 0
                theProcessedData[ti,xx,6] = theProcessedData[ti,xx,3] + theData[ti,xx,3] ; 5

            ENDFOR
        ENDFOR



        stamp = {name:ListOfModels[k], theData:theProcessedData}
        *(ptrListOfModels[k]) = stamp

    ENDFOR

    colors=indgen(n_elements(theProcessedData[0,0,*]))
    colors=[0,1,2,3,2,4]
    labels=["r (kpc)","Fraction of dsigdt"]
    texLabels=["$r$ (kpc)","Share of $\partial \sigma/\partial t$"]
    baseName="balance_sig"
    leg=["Viscous heating","Turbulent dissipation","Mass flux gradient","Velocity dispersion gradient"]
    legColors=[4,1,2,3]

;    FOR a=0, 100 DO BEGIN
        makeBalancePlots, ptrListOfModels, whichRedshifts,zs,th,annotations,texAnnotations,svSinglePlot,cs,chth,logs,ranges,colors,labels,texLabels,baseName,fill=1,legName=leg,legColors=legColors
;	    stop
;    ENDFOR

END


PRO balanceplotsZ,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations

    compile_opt idl2
    set_plot,'x'
    device,decomposed=0
    device,retain=2

    IF(n_elements(annotations) EQ 0) THEN BEGIN
	annotations=ListOfModels[*]
	FOR i=0,n_elements(annotations)-1 DO annotations[i]=ExtractSlashName(annotations[i],1)
    ENDIF
    IF(n_elements(texAnnotations) EQ 0) THEN texAnnotations=annotations[*]

    whichRedshifts=[1,51,201]
    zs=['z=2.00','z=1.01','z=0.00']
    columns = n_elements(ListOfModels)
    rows = n_elements(zs)
    setct,1,0,0
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2


    dummy = GetLabels(keys=[0,41,42,43,44], names=names, $
        labels=labels, tex=texLabels, base=vars, $
        log=logs, offset=offsets,ranges=ranges)



    ;; collect data from disk.
    ptrListOfModels = ptrarr(n_elements(ListOfModels), /allocate_heap)
    FOR k=0, n_elements(ListOfModels)-1 DO BEGIN

        model = readOutput(ListOfModels[k])
        theData = dblarr(n_elements(whichRedshifts),model.nx,n_elements(vars))
        theData[*,*,*] = model.dataCube[whichRedshifts,*,vars+offsets*model.ncolstep-1]
        ;; at this point, per the keys requested of GetLabels, we have
        ;; radius,"dZdtAdv","dZdtDiff","dZdtDil","dZdtSF" 
        ;; Now we normalize things (besides the radius!)
;        thenorm = abs(theData[*,*,1]) + abs(theData[*,*,2]) + abs(theData[*,*,3])+abs(theData[*,*,4])
        thenorm = TOTAL(abs(theData[*,*,1:n_elements(vars)-1]),3)
        
        FOR i=1,n_elements(vars)-1 DO theData[*,*,i] = theData[*,*,i]/theNorm
        theProcessedData=dblarr(n_elements(whichRedshifts),model.nx,7)
        theProcessedData[*,*,0]=theData[*,*,0] ;; radius
        FOR xx=0,n_elements(theData[0,*,0]) -1 DO BEGIN
            FOR ti=0,n_elements(theData[*,0,0]) -1 DO BEGIN
                theProcessedData[ti,xx,1] = MIN([0,theData[ti,xx,2]]) ;; closest to zero- bound from below (2)
                theProcessedData[ti,xx,2] = MAX([0,theData[ti,xx,2]]) ;; closest to zero- bound from above (3)
                theProcessedData[ti,xx,3] = MAX([0,theProcessedData[ti,xx,2],theData[ti,xx,2]+theData[ti,xx,1]]) ; (4)
                theProcessedData[ti,xx,4] = MIN([0,theProcessedData[ti,xx,1],theData[ti,xx,2]+theData[ti,xx,1]]) ; (1)
                theProcessedData[ti,xx,5] = theProcessedData[ti,xx,4] + theData[ti,xx,3] ; 0
                theProcessedData[ti,xx,6] = theProcessedData[ti,xx,3] + theData[ti,xx,4] ; 5

            ENDFOR
        ENDFOR



        stamp = {name:ListOfModels[k], theData:theProcessedData}
        *(ptrListOfModels[k]) = stamp

    ENDFOR

    colors=indgen(n_elements(theProcessedData[0,0,*]))
    colors=[0,3,2,4,2,1]
    labels=["r (kpc)","Fraction of dsigdt"]
    texLabels=["$r$ (kpc)","Share of $\partial Z/\partial t$"]
    baseName="balance_Z"

;    FOR a=0, 100 DO BEGIN
        makeBalancePlots, ptrListOfModels, whichRedshifts,zs,th,annotations,texAnnotations,svSinglePlot,cs,chth,logs,ranges,colors,labels,texLabels,baseName,fill=1
;	    stop
;    ENDFOR

END






PRO balanceplots,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations
    balanceplotsCol,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations
    balanceplotsCol2,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations
    balanceplotsSig,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations
    balanceplotsZ,ListOfModels,Annotations=Annotations,TexAnnotations=TexAnnotations
END
