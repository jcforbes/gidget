;; The goal here is to take an experiment (i.e. a set of gidget runs) where it is assumed that one variable is being changed, in particular the seed of the accretion history. Thus we expect every parameter to have some variability in time and radius.

PRO extremize,extL,extR,model
	newMins = WHERE(model.evArray LT extL.evArray,count)
	IF(count NE 0) THEN extL.evArray[newMins]=model.evArray[newMins]
	newMaxs = WHERE(model.evArray GT extR.evArray,count)
	IF(count NE 0) THEN extR.evArray[newMaxs]=model.evArray[newMaxs]

	newMins = WHERE(model.dataCube LT extL.dataCube,count)
	IF(count NE 0) THEN extL.dataCube[newMins]=model.dataCube[newMins]
	newMaxs = WHERE(model.dataCube GT extR.dataCube,count)
	IF(count NE 0) THEN extR.dataCube[newMaxs]=model.dataCube[newMaxs]


	newMins = WHERE(model.stellarData LT extL.stellarData,count)
	IF(count NE 0) THEN extL.evArray[newMins]=model.evArray[newMins]
	newMaxs = WHERE(model.stellarData GT extR.stellarData,count)
	IF(count NE 0) THEN extR.evArray[newMaxs]=model.evArray[newMaxs]

	newMins = WHERE(model.stellarAges LT extL.stellarAges,count)
	IF(count NE 0) THEN extL.stellarAges[newMins]=model.stellarAges[newMins]
	newMaxs = WHERE(model.stellarAges GT extR.stellarAges,count)
	IF(count NE 0) THEN extR.stellarAges[newMaxs]=model.stellarAges[newMaxs]
	
END

PRO addTo,sum,model
	sum.evArray = temporary(sum.evArray) + model.evArray
	sum.dataCube= temporary(sum.dataCube)+ model.dataCube
	sum.stellarData= temporary(sum.stellarData) + model.stellarData
	sum.stellarAges= temporary(sum.stellarAges) + model.stellarAges
	
END

PRO add2To,sum2,model
	sum2.evArray=temporary(sum2.evArray) + (model.evArray)*(model.evArray)
	sum2.dataCube=temporary(sum2.dataCube) +(model.dataCube)*(model.dataCube)
	sum2.stellarData=temporary(sum2.stellarData) + (model.stellarData)*(model.stellarData)
	sum2.stellarAges=temporary(sum2.stellarAges) + (model.stellarAges)*(model.stellarAges)
END
PRO addconst,model,const
	model.evArray = temporary(model.evArray)+const
	model.dataCube=temporary(model.dataCube)+const
	model.stellarData=temporary(model.stellarData)+const
	model.stellarAges=temporary(model.stellarAges)+const
END
PRO zero,model
	model.evArray=model.evArray*0.0
	model.dataCube=model.dataCube*0.0
	model.stellarData=model.stellarData*0.0
	model.stellarAges=model.stellarAges*0.0
END

PRO linsum,out,in1,in2,c1,c2
	out.evArray = c1*in1.evArray + c2*in2.evArray
	out.dataCube= c1*in1.dataCube+ c2*in2.dataCube
	out.stellarData=c1*in1.stellarData+c2*in2.stellarData
	out.stellarAges=c1*in1.stellarAges+c2*in2.stellarAGes
END

PRO multiply,out,in1,in2
	out.evArray = in1.evArray*in2.evArray
	out.dataCube  = in1.dataCube* in2.dataCube
	out.stellarData= in1.stellarData * in2.stellarData
	out.stellarAges= in1.stellarAges * in2.stellarAges
END

PRO power,out,in,exponent
	out.evArray  = in.evArray ^ exponent
	out.dataCube = in.dataCube ^ exponent
	out.stellarData=in.stellarData ^ exponent
	out.stellarAges=in.stellarAges ^ exponent
END

PRO countbetween, counter, lower, upper,model
	isBetween = WHERE(lower.evArray LE model.evArray AND model.evArray LE upper.evArray)
	counter.evArray[isBetween] = counter.evArray[isBetween] + (lower.evArray[isBetween]*0+1)

	isBetween = WHERE(lower.dataCube LE model.dataCube AND model.dataCube LE upper.dataCube)
	counter.dataCube[isBetween] = counter.dataCube[isBetween] + (lower.dataCube[isBetween]*0+1)

	isBetween = WHERE(lower.stellarData LE model.stellarData AND model.stellarData LE upper.stellarData)
	counter.stellarData[isBetween] = counter.stellarData[isBetween] + (lower.stellarData[isBetween]*0+1)

	isBetween = WHERE(lower.stellarAges LE model.stellarAges AND model.stellarAges LE upper.stellarAges)
	counter.stellarAges[isBetween] = counter.stellarAges[isBetween] + (lower.stellarAges[isBetween]*0+1)
	
END

PRO variability,expName,sv
	nm = expName+'/'+expName+'_'
	nameList = file_search(expName+'/*_comment.txt')
	FOR i=0, n_elements(nameList)-1 DO nameList[i] = strmid(nameList[i], 0, strlen(nameList[i])-12)
	nameList2 = nameList
	ctr=0
	FOR i=0, n_elements(nameList2)-1 DO BEGIN
		IF(success(nameList[i])) THEN BEGIN
			nameList2[ctr]=nameList[i]
			ctr=ctr+1
		ENDIF
	ENDFOR
	nameList2=nameList2[0:ctr-1]

	print,"Of the ",n_elements(nameList), " models, ",ctr," are valid. Beginning analysis of valid models.."

	;; First pass - loop over all models.
	storeModelPtr = ReadModels([nameList2[0],nameList2[0],nameList2[0],nameList2[0],nameList2[0],nameList2[0], nameList2[0],nameList2[0],nameList2[0]])
	avg = (*(storemodelPtr[0]))
	sigL = (*(storemodelPtr[1]))
	sigR = (*(storemodelPtr[2]))
	extL = (*(storemodelPtr[3]))
	extR = (*(storemodelPtr[4]))
	temp=(*(storemodelPtr[5]))
	temp2=(*(storemodelPtr[6]))
	sum=(*(storemodelPtr[7]))
	zero,sum
	sum2=(*(storemodelPtr[8]))
	zero,sum2


	valid = 0
	FOR j=0, n_elements(nameList2)-1 DO BEGIN
		modelPtr = ReadModels([nameList2[j]])
		model = (*(modelPtr[0]))
		;; 
		IF(j EQ 0) THEN BEGIN
			fid = model
;			avg = ReadModels([nameList2
			; try it..
		ENDIF
		
		IF(j GE 3 AND n_elements(avg.evArray) EQ n_elements(model.evArray)) THEN BEGIN
			extremize,extL,extR,model
			addTo,sum,model ;; sum+=model
			add2To,sum2,model;;sum2+=model^2
			valid=valid+1
		ENDIF
		
		IF(j MOD 100 EQ 0) THEN PRINT, "Currently reading in model ",j

		ptr_free,modelPtr[0]
	ENDFOR
	linsum, avg, sum,sum,double(1.0/valid),0.0 ;; compute the average 
	power,temp,avg,2.0 ;; temp = mu^2
	linsum,temp2, sum2,temp, double(1.0/(valid-1.0)), -double(valid)/(valid-1.0) ;; temp = variance
	power,sigL,temp2,.5 ;; sigL = std dev
	linsum,sigR,sigL,sigL,1.0,0.0 ;; sigL=sigR
	
	;; first pass complete: we have extL = min, extR=max, avg=average, sigL=sigR= std dev.
	PRINT, "First pass complete! Valid = ",valid


	bL=avg
	bR=avg
	linsum,bL,avg,sigL,1.0,-1.0  ;; bL = avg - sigL
	linsum,bR,avg,sigR,1.0,1.0   ;; bR = avg + sigR

	;; construct an array of model pointers:
	ptrlist=ptrarr(3, /allocate_heap)
	extL.name="min"
	bL.name="bL"
	avg.name="avg"
	bR.name="bR"
	extR.name="max"
	*(ptrlist[0])=bL
	*(ptrlist[1])=avg
	*(ptrlist[2])=bR
	diskanalyzecompare,ptrlist,sv,"stats after 1st pass",['m-s','m','m+s']

	NL = avg
	NR = avg
	zero,NL
	zero,NR
	addconst,NL,1 ;; this will avoid potential divide-by-zero issues
	addconst,NR,1
	;; next pass; count up number of models between extL and bL=avg-sigL, and between bR=avg+sigR and extR:
	
	FOR j=3, n_elements(nameList2)-1 DO BEGIN
		modelPtr = ReadModels([nameList2[j]])
		model = (*(modelPtr[0]))

		IF(n_elements(avg.evArray) EQ n_elements(model.evArray)) THEN BEGIN
			;;
			CountBetween,NL,extL,bL,model
			CountBetween,NR,bR,extR,model
		ENDIF
		ptr_free,modelptr[0]
	ENDFOR

	alpha=1.0
	delta=avg
	linsum,delta,bL,extL,(.16*valid)^alpha,-(.16*valid)^alpha
	power,temp,NL,-alpha
	multiply,temp2,temp,delta
	linsum,bL,temp2,extL,1.0,1.0

	linsum,delta,extR,bR,-(.16*valid)^alpha,(.16*valid)^alpha
	power,temp,NR,-alpha
	multiply,temp2,temp,delta
	linsum,bR,extR,temp2,1.0,1.0 
	;; second pass complete

	ptrlist=ptrarr(3, /allocate_heap)
	extL.name="min"
	bL.name="bL2"
	avg.name="avg"
	bR.name="bR2"
	extR.name="max"
	*(ptrlist[0])=bL
	*(ptrlist[1])=avg
	*(ptrlist[2])=bR
	diskanalyzecompare,ptrlist,sv,"stats after 2nd pass",['m-s','m','m+s']



	ptr_free,ptrlist[0]
	ptr_free,ptrlist[1]
	ptr_free,ptrlist[2]
	FOR i=0,n_elements(storemodelptr)-1 DO ptr_free,storemodelptr[i]
END

