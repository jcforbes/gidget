
;; an n-dimensional array can be indexed by n indices or just 1. Given the single index, this function
;; will give you the n indices to which it corresponds.
;; index is the single index from 1 to N elements in the multi-dimensional array
;; nVarDim is an n-element array which contains the number of variables in each dimension of the n-dimensional array
;; note that N=n[0]*n[1]*...
FUNCTION unpack,index,nVarDim
	elements=intarr(n_elements(nVarDim))
	;;; index = elements[last] + elements[last-1]*nvardim[last] + elements[last-2]*nvardim[last]*nvardim[last-1] + ... + elements[0]*nvardim[last]*nvardim[last-1]*...*nvardim[1]
	;; also,  elements[i]<nvardim[i]
	;; 
	;; maybe not the most elegent solution, but this is an
	;; implementation of the general formula
	;; elements_i = (index - (index % N_n*...*N_{i+1})
	;;		 - elements_0 * N_n *...*N_1
	;;		 - elements_1 * N_n *...*N_2
	;;		 - ... - elements_{i-1} * N_n*...*N_i)
	FOR i=0,n_elements(nVarDim)-1 DO BEGIN
		den=1
		FOR j=i+1,n_elements(nVarDim)-1 DO den*=nVarDim[j]
		subterms=0
		FOR j=0,i-1 DO BEGIN
			m2=1
			FOR k=j+1,n_elements(nVarDim)-1 DO m2*=nVarDim[k]
			subterms += elements[j]*m2
		ENDFOR
		elements[i] = (index - (index MOD den) - subterms)/den
	ENDFOR

	RETURN,elements
END

FUNCTION autoNVarDim,nameList
	flag=0
	dims=0
	FOR j=0,strlen(nameList[0])-1 DO BEGIN
		;; Go through the first name on this list from the end until the letter is not a.
		jthChar= strmid(nameList[0],j,1,/REVERSE_OFFSET)
		IF((jthChar NE 'a') AND flag EQ 0) THEN BEGIN
			flag=1
			dims=j
		ENDIF
	ENDFOR
	nvardim=intarr(dims)
	alphabet='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
	FOR k=0, dims-1 DO BEGIN
		nvardim[dims-1-k]=strpos(alphabet, strmid((nameList[n_elements(nameList)-1]), k, 1, /REVERSE_OFFSET))+1
	ENDFOR

	RETURN, nvardim
END

;; For a given variable, e.g. column density, sequentially plot that variable from eachmodel
PRO seqDiskPlot,nameList,sv,yind,ytitle,yr,yp,prnt,k
	compile_opt idl2
	
;;;	nVarDim = autoNVarDim,nameList
	nVarDim = [n_elements(nameList)]

;;;	IF(sv EQ 2) THEN figureInit, "seqDiskPlot",sv,1,1

	FOR i=0,n_elements(nameList)-1 DO BEGIN
		telements=(unpack(i,nVarDim))
		srtd=nVarDim[reverse(sort(nvardim))]
		elements=telements[reverse(sort(nvardim))]
		IF(n_elements(elements) LT 3) THEN BEGIN
			qelements=intarr(3)
			qelements[0:n_elements(elements)-1]=elements[0:n_elements(elements)-1]
			elements=qelements
		ENDIF
		;; the number of colors necessary is srtd[0], 
		;; which color set to use is elements[2]
		setct,3,srtd[0],elements[2]
		modelPtr=ReadModels([nameList[i]])
		model = (*(modelPtr[0]))
		zi=nearest("redshift",1,model,0.0)
		offsets=[0,model.ncolstep,model.ncolstep+model.npostprocess]
		
;;;		FOR k=0,n_elements(yind)-1 DO BEGIN
;;;			IF(sv EQ 0) THEN IF(i EQ 0) THEN window,k ELSE wset,k
			PlotVsXInd,model,zi,model.ncolstep+9,yind[k]+offsets[yp[k]],OP(i),0,COLOR=elements[0],LINESTYLE=(elements[1]+fix(elements[1] GE 1)),THICK=1,XTITLE="r [kpc]",YTITLE=ytitle[k],YRANGE=yr[*,k],YSTYLE=1,XSTYLE=1
;;;			STOP
		

;;;		ENDFOR
		IF(prnt EQ 0) THEN print,model.name,diskStats(model)
		ptr_free,modelPtr[0]
	ENDFOR

	IF(sv EQ 2) THEN figureClean, "seqDiskPlot",sv
END


PRO seqDiskAnalyze,expName,sv
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
	nameList2=nameList2[0:ctr-1]

	;; variables to plot- title, column index, log plot, name, y-range
	yt=['Col','Sig','Col_*','Sig_*','fg','Z']
	yy=[10,11,12,13,17,22]
	yp=[ 1, 1, 1, 1, 0, 0]
	yn=['col','sig','colst','sigst','fg','Z']
	yr=[[.1,1000],[5,300],[.1,1000],[5,300],[0,1],[-1,1.0]]
	;; loop over variables

;;;;	seqDiskPlot,nameList2,sv,yy,yt,yr,yp,1

	FOR k=0, n_elements(yy)-1 DO BEGIN
		figureInit,nm+yn[k],sv,1,1
		seqDiskPlot,nameList2,sv,yy,yt,yr,yp,k,k
		figureClean,nm+yn[k],sv
	ENDFOR
END
