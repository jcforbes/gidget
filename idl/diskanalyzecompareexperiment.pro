
;; wrapper for diskAnalyzeCompare, where the set of runs to be analyzed is an 'experiment' from the xgrid
PRO diskAnalyzeCompareExperiment,exp,sv,cmt
	IF(n_elements(exp) EQ 0 || n_elements(sv) EQ 0 || n_elements(cmt) EQ 0) THEN message,"Usage: diskAnalyzeCompareExperiment,expName,sv,cmt; expName is the experiment name, i.e. a directory containing the output of a set of related runs produced by the python script. sv=0 -> print plots to screen only, sv=1 -> use IDL's built-in MPEG algorithm (not recommended), sv=2 -> save plots as pngs and stack them together with ffmpeg, sv=3 -> same as sv=2, except use zbuffer rather than printing to the screen first (recommended if you are running this code over a remote connection). cmt is a comment which will be written in every frame of the movies created"
	nameList=file_search(exp+'/*_comment.txt')
	;; in each elements of the list remove '_comment.txt':
	FOR i=0, n_elements(nameList)-1 DO BEGIN
		nameList[i]=strmid(nameList[i],0,strlen(nameList[i])-12)
	ENDFOR

	;; now pick out only those runs which were successful, i.e. with no standard error output
	nameList2=strarr(n_elements(nameList))
	counter=0
	FOR i=0, n_elements(nameList2)-1 DO BEGIN
		IF(success(nameList[i])) THEN BEGIN
			nameList2[counter]=nameList[i]
			counter=counter+1
		ENDIF
	ENDFOR
	nameList2=nameList2[0:counter-1]
	leg=strarr(n_elements(nameList2))
	diskAnalyzeCompare,nameList2,sv,cmt,leg
END
