PRO makeThePostageStampPlots, stampList, whichRedshifts,vars,zs,th,labels,stampLabels,texLabels,texStampLabels,names,columns,rows,svSinglePlot,cs,chth,logs,ranges

    cs=0.53

    ;; now that we have the data, let's make the plots.
    ; loop over variables
    ranges=dblarr(2,n_elements(vars))
    ranges[*,0]=[-5,45]
    ranges[*,1]=[.5,140]
    ranges[*,2]=[-1.1,.7]
    ranges[*,3]=[4.0,50.0]
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
                    gap =0.15; 0.2
                    xoff = .075+gap + (1.0-2.8*gap) * (.9/double(columns) + double((k MOD columns))/double(columns))
                    yoff = 1.0 - 1.05*gap - (1.0-2.05*gap)* double(k/columns)/double(rows)
                    XYOUTS,xoff,yoff,stampLabels[k],/normal,color=0,charsize=cs,charthick=chth
                ENDFOR

                modmultiplot
            ENDFOR
            FigureClean,filename,svSinglePlot
            latexify,(filename+".eps"),[labels,stampLabels],[texLabels,texStampLabels],[replicate(1.4,n_elements(labels)),replicate(.9,n_elements(stampLabels))],tempname=("TMP_"+filename);height=8.89*rows,width=8.89*columns,tempname=("TMP_"+filename)
            modmultiplot,/default
            spawn,"ps2pdf -dEPSCrop "+filename+".eps"
        ENDFOR
    ENDFOR

    set_plot,'x'
END
