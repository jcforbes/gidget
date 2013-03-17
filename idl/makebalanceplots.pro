PRO makeBalancePlots, ptrListOfModels, whichRedshifts,zs,th,annotations,texAnnotations,svSinglePlot,cs,chth,logs,ranges,colors,labels,texLabels,baseName,fill=fill

    cs=0.53

    filename = ""
    FOR k=0,n_elements(ptrListOfModels)-1 DO filename=filename+ExtractSlashName((*(ptrListOfModels[k])).name,1)+"_"
    filename=filename+basename+"_all"
    FigureInit,filename,svSinglePlot,n_elements(whichRedshifts),n_elements(ptrListOfModels)
    myTitle = "Fraction of";labels[1]
    modmultiplot,[n_elements(whichRedshifts),n_elements(ptrListOfModels)],mxtitle=labels[0],mxtitSize=cs, myTitle=myTitle,myTitSize=cs,additionalMargin=.11
    
    !p.noerase=0
    FOR k=0,n_elements(ptrListOfModels)-1 DO BEGIN
        theData = (*(ptrListOfModels[k])).theData
        theName = ExtractSlashName((*(ptrListOfModels[k])).name,1)
        gap =0.15; 0.2
        yoff = 1.0 - 1.3*gap - (1.0-2.35*gap)* double(k)/double(n_elements(ptrListOfModels))
        FOR i=0,n_elements(whichRedshifts)-1 DO BEGIN
            IF(k NE 0 OR i NE 0) THEN !p.noerase=1
            yt=""
            IF(i EQ 0 and k EQ n_elements(ptrListOfModels)/2 ) THEN yt=labels[1]
            PLOT,[0],[0],/nodata,COLOR=0,BACKGROUND=255,XRANGE=ranges[*,0],YRANGE=[-1.1,1.1],$
                XSTYLE=1,YSTYLE=1,ylog=logs[1],CHARSIZE=cs,CHARTHICK=chth,THICK=th,XTHICK=th, $
                YTHICK=th,LINESTYLE=0,ticklen=.05,YTITLE=yt;,XTITLE=labels[0],YTITLE=labels[1]
            IF(fill EQ 1) THEN BEGIN
                aa = theData[i,*,1:n_elements(theData[0,0,*])-1]
                FOR xx=0, n_elements(aa[0,*,0])-1 DO aa[0,xx,*] = aa[0,xx,SORT(aa[0,xx,*])]
                FOR j=1,n_elements(aa[0,0,*])-1 DO BEGIN
                    POLYFILL, [[theData[i,*,0]], [REVERSE(theData[i,*,0],2)]], [[aa[0,*,j-1]],[REVERSE(aa[0,*,j],2)]], COLOR= colors[j];line_fill=1
                ENDFOR
            ENDIF

            xoff = .075+gap + (1.0-2.8*gap) * (.9/double(n_elements(zs)) + double(i)/double(n_elements(zs)))
            IF(k EQ 0) THEN XYOUTS,xoff,yoff,zs[i],/normal,color=0,charsize=cs*1.2,charthick=chth
            IF(i EQ 0) THEN XYOUTS,.87,yoff-.07+.07*k/double(n_elements(ptrListOfModels)),annotations[k],color=0,charsize=cs*1.5,charthick=chth,/normal

            modmultiplot

        ENDFOR
;        XYOUTS,.9,yoff,theName,color=0,charsize=cs,charthick=chth,/normal
    ENDFOR
    FigureClean,filename,svSinglePlot
    latexify,(filename+".eps"),[myTitle,labels],[myTitle,texLabels],[replicate(1.0,n_elements(labels)+1)],tempname=("TMP_"+filename);height=8.89*rows,width=8.89*columns,tempname=("TMP_"+filename)
    modmultiplot,/default
    spawn,"ps2pdf -dEPSCrop "+filename+".eps"
    spawn,"rm "+filename+".eps"







    FOR k=0,n_elements(ptrListOfModels)-1 DO BEGIN
        theData = (*(ptrListOfModels[k])).theData
        theName = ExtractSlashName( (*(ptrListOfModels[k])).name, 1)
        filename = theName+"_"+baseName+"_only"
        FigureInit,filename,svSinglePlot,n_elements(whichRedshifts),1
        modmultiplot,[n_elements(whichRedshifts),1],mxtitle=labels[0],mxtitSize=cs, myTitle=myTitle,myTitSize=cs,/square,mxTitOffset=-.35,myTitOffset=-.35
; myTitle=myTitle,myTitSize=cs;,mxtitle=labels[0],mxtitSize=cs,/square,additionalMargin=0.2       
        !p.noerase=0
        
        gap =0.15; 0.2
        yoff = 1.0 - 1.7*gap ;- (1.0-2.05*gap)* double(k)/double(rows)
        FOR i=0,n_elements(whichRedshifts)-1 DO BEGIN
            IF( i NE 0) THEN !p.noerase=1
            yt=""
            IF(i EQ 0) THEN yt=labels[1] 
            PLOT,[0],[0],/nodata,COLOR=0,BACKGROUND=255,XRANGE=ranges[*,0],YRANGE=[-1.1,1.1],$
                XSTYLE=1,YSTYLE=1,ylog=logs[1],CHARSIZE=cs,CHARTHICK=chth,THICK=th,XTHICK=th, $
                YTHICK=th,LINESTYLE=0,ticklen=.05,YTITLE=yt
            IF(fill EQ 1) THEN BEGIN
                aa = theData[i,*,1:n_elements(theData[0,0,*])-1]
                FOR xx=0, n_elements(aa[0,*,0])-1 DO aa[0,xx,*] = aa[0,xx,SORT(aa[0,xx,*])]
                FOR j=1,n_elements(aa[0,0,*])-1 DO BEGIN
                    POLYFILL, [[theData[i,*,0]], [REVERSE(theData[i,*,0],2)]], [[aa[0,*,j-1]],[REVERSE(aa[0,*,j],2)]], COLOR= colors[j];line_fill=1
                ENDFOR
            ENDIF

            xoff = .828*gap + (1.0-3.5*gap) * (.9/double(n_elements(zs)) + double(i)/double(n_elements(zs)))
            XYOUTS,xoff,yoff,zs[i],/normal,color=0,charsize=cs*1.2,charthick=chth

            modmultiplot

        ENDFOR
;        XYOUTS,.9,yoff,theName,color=0,charsize=cs,charthick=chth,/normal
        FigureClean,filename,svSinglePlot
        latexify,(filename+".eps"),labels,texLabels,[replicate(1.3,n_elements(labels))],tempname=("TMP_"+filename);height=8.89*rows,width=8.89*columns,tempname=("TMP_"+filename)
        modmultiplot,/default
        spawn,"ps2pdf -dEPSCrop "+filename+".eps"
        spawn,"rm "+filename+".eps"
    ;    spawn,"rm TMP_"+filename+".eps"
    ;    spawn,"rm TMP_"+filename+".pdf"
    ENDFOR


    set_plot,'x'
END
