;; Simple function for getting the base name of a run
;; even if it is in a subdirectory
FUNCTION ExtractSlashName,theName,after
;	RETURN,strmid(thename,(strlen(thename)-1)/2 + 1)
	IF(after EQ 1) THEN RETURN,strmid(theName,$
		strpos(theName,'/')+1) $
	ELSE RETURN,strmid((*(model)).name,0,strpos((*(model)).name,'/'))
END