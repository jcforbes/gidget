
PRO plotvsx,model,timeInd,ind,overPlot,_EXTRA=EXTRA_KEYWORDS
        PlotVsXind,model,timeInd,Xind,ind,overPlot,1,_EXTRA=EXTRA_KEYWORDS
END

PRO plotvsxind,model,timeInd,Xind,ind,overPlot,getLabels,_EXTRA=EXTRA_KEYWORDS;range,chars,col
        IF(getLabels EQ 1) THEN BEGIN
                label=GetLabel(ind,model.ncolstep,$
			model.npostprocess,$
			model.npassive,model.stvars)
                xlabel=GetLabel(xind,model.ncolstep,$
			model.npostprocess,$
			model.npassive,model.stvars)
        ENDIF
        log=logvar(ind,model)

        IF overPlot EQ 1 THEN BEGIN
                oplot,model.dataCube[timeInd,*,Xind-1],$
			model.dataCube[timeInd,*,ind-1],$
                        _extra=extra_keywords
        ENDIF ELSE BEGIN
                plot,model.dataCube[timeInd,*,Xind-1],$
			model.dataCube[timeInd,*,ind-1],$
                        XTITLE=xlabel,YTITLE=label,$
                        BACKGROUND=255,ylog=log,xlog=1,$
			_extra=extra_keywords
        END
;	PRINT,"PlotVsXind- min max of x,y: ",min(model.dataCube[timeInd,*,Xind-1]),max(model.dataCube[timeInd,*,Xind-1]),min(model.dataCube[timeInd,*,ind-1]),max(model.dataCube[timeInd,*,ind-1]),overPlot
END
