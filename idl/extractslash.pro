;; Simple function for getting the base name of a run
;; even if it is in a subdirectory
FUNCTION ExtractSlash,model,after
;	RETURN,strmid((*(model)).name,(strlen((*(model)).name)-1)/2 + 1)
	IF(after EQ 1) THEN RETURN,strmid((*(model)).name,$
		strpos((*(model)).name,'/')+1) $
	ELSE RETURN,strmid((*(model)).name,0,strpos((*(model)).name,'/'))
END