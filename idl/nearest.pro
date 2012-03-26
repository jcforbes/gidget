;;;; A series of functions
;; given a time, find the actual recorded time nearest to the requested time.
FUNCTION nearest,typ,ind,model,var
	IF(typ EQ "redshift") THEN BEGIN
		near=min(abs(model.evArray[9,*]-var),index)
		val=model.evArray[9,index]
	ENDIF
	IF(typ EQ "time") THEN BEGIN
		near=min(abs(model.evArray[1,*]-var),index)
		val=model.evArray[1,index]
	ENDIF
	IF(typ EQ "position") THEN BEGIN
		near=min(abs(model.dataCube[0,*,0]-var),index)
		val=model.dataCube[0,index,0]
	ENDIF
	IF(ind EQ 1) THEN RETURN,index ELSE RETURN,val
	
END
