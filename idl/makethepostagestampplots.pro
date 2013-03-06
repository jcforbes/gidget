PRO makeThePostageStampPlots, stampList, whichRedshifts,vars,zs,th, $
    labels,stampLabels,texLabels,texStampLabels, $
    lowRange, hiRange, $
    names,columns,rows,svSinglePlot,cs,chth,logs,ranges

    cs=0.53

    ;; now that we have the data, let's make the plots.
    ; loop over variables
    FOR j=0,n_elements(whichRedshifts)-1 DO BEGIN
        FOR i=1,n_elements(vars)-1 DO BEGIN
            filename = "stamps_z"+zs[j]+"_"+names[i]
            FigureInit,filename,svSinglePlot,columns,rows
            modmultiplot,[columns,rows],mxtitle=labels[0],myTitle=labels[i],mxTitSize=cs*2,myTitSize=cs*2,mxTitOffset=-0.6,myTitOffset=-1.9;,xgap=.001,ygap=.001
            !p.noerase=0
            FOR k=0, n_elements(stampList)-1 DO BEGIN
                IF(k NE 0) THEN !p.noerase = 1
                theData = (*(stampList[k])).theData
                numberOfModels = (*(stampList[k])).numberOfModels

                PLOT,[0],[0],/nodata,COLOR=0,BACKGROUND=255,XRANGE=ranges[*,0],YRANGE=ranges[*,i], $
                    XSTYLE=1,YSTYLE=1,ylog=logs[i],CHARSIZE=cs,CHARTHICK=chth,THICK=th,XTHICK=th, $
                    YTHICK=th,LINESTYLE=0,ticklen=.05
                FOR exper=0,n_elements(numberOfModels)-1 DO BEGIN
                    offsetModel = 0
                    tthh = th
                    IF(exper NE 0) THEN tthh=1.0
                    IF(exper NE 0) THEN offsetModel = TOTAL(numberOfModels[0:exper-1])
                    FOR m=0,numberOfModels[exper]-1 DO OPLOT, theData[j,*,0,m+offsetModel],theData[j,*,i,m+offsetModel],COLOR=exper,THICK=tthh

                ENDFOR
                gap =0.15; 0.2
                xoff = gap + (1.0-2.3*gap) * (.8/double(columns) + double((k MOD columns))/double(columns))
                yoff = 1.0 - 1.05*gap - (1.0-2.05*gap)* double(k/columns)/double(rows)
                XYOUTS,xoff,yoff,stampLabels[k],/normal,color=0,charsize=cs*.9,charthick=chth
                XYOUTS,xoff,yoff-.01,lowRange[k],/normal,color=1,charsize=cs*.8,charthick=chth*.8
                XYOUTS,xoff,yoff-.02,hiRange[k],/normal,color=2,charsize=cs*.8,charthick=chth*.8
                modmultiplot
            ENDFOR
            FigureClean,filename,svSinglePlot
            latexify,(filename+".eps"),[labels,stampLabels],[texLabels,texStampLabels],[replicate(1.2,n_elements(labels)),replicate(.6,n_elements(stampLabels))],tempname=("TMP_"+filename);height=8.89*rows,width=8.89*columns,tempname=("TMP_"+filename)
            modmultiplot,/default
            spawn,"ps2pdf -dEPSCrop "+filename+".eps"
        ENDFOR
    ENDFOR

    set_plot,'x'
END
