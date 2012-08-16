FUNCTION listvalidmodels,expName,N
  nm=expName+'/'+expName+'_'
  nameList=file_search(expName+'/*_comment.txt')
  FOR i=0, n_elements(nameList)-1 DO nameList[i] = strmid(nameList[i], 0, strlen(nameList[i])-12)
  nameList2=nameList
  ctr=0
  FOR i=0, n_elements(nameList2)-1 DO BEGIN
    IF(success(nameList[i])) THEN BEGIN
      nameList2[ctr]=nameList[i]
      ctr=ctr+1
    ENDIF
  ENDFOR
  IF(ctr EQ 0) THEN message,"No valid models found"
  nameList2=nameList2[0:ctr-1]

  IF(N LT ctr) THEN BEGIN 
;      nameList2=nameList2[0:N-1]
      FOR k=0,N-1 DO BEGIN
          nameList2[k] = nameList2[k*FIX(float(ctr)/float(N-1))]
      ENDFOR
      nameList2=nameList2[0:N-1]
  ENDIF

  RETURN,nameList2  
  
END

