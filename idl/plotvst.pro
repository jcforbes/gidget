;; here index refers to the column number in the evolution file if position < 0
;; or to the column number in the step files
PRO plotvst,model,index,position,overPlot,_EXTRA=EXTRA_KEYWORDS
        IF (position LT 0) THEN BEGIN
                IF (overPlot EQ 1) THEN $
                        oplot,model.evArray[1,*],$
			model.evArray[index-1,*],$
			_EXTRA=EXTRA_KEYWORDS $ 
                ELSE plot,model.evArray[1,*],$
			model.evArray[index-1,*],$
			_EXTRA=EXTRA_KEYWORDS
                
        ENDIF ELSE BEGIN
                IF(overPlot EQ 1) THEN $
                        oplot,model.evArray[1,*],$
			model.dataCube[*,nearest("position",$
				1,model,position),index-1],$
			_EXTRA=EXTRA_KEYWORDS $
                ELSE plot,model.evArray[1,*],$
			model.dataCube[*,nearest("position",1,$
				model,position),index-1],$
			_EXTRA=EXTRA_KEYWORDS
        END
        
END