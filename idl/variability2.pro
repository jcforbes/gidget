;; given an array whose first column indexes time, second column indexes position,
;; fourth column indexes models, and a vector of time labels, and axis labels
;; make one or a few movies! The third column indexes which variable we're plotting.
PRO simpleMovie,data,time,labels,colors,styles,ranges,wrtXlog,name,sv,strt,prev=prev,psym=psym
	IF(sv EQ 3) THEN cg=1 ELSE cg=0
	IF(sv EQ 1) THEN cs=1 ELSE cs=2.8

	FOR k=1,n_elements(labels)-1 DO BEGIN ;; loop over y-axis variables to be plotted (except the first one, which is just radius!)
		fn=name+"_"+labels[k]
		dn="movie_"+name+"_"+labels[k]
		set_plot,'x'
		IF(sv EQ 3) THEN set_plot,'z'
		IF(cg NE 1) THEN WINDOW,1,XSIZE=1024,YSIZE=1024,RETAIN=2
		IF(sv EQ 3) THEN cgDisplay,1024,1024
		IF(sv EQ 1) THEN mpeg_id=MPEG_OPEN([1024,1024],FILENAME=(fn+".mpg"))
		IF(sv EQ 2 || sv EQ 3) THEN BEGIN
			FILE_MKDIR,dn
		ENDIF
		count=0
		tailLength = fix(3.0/(float(n_elements(styles))/928.0)^(.25)) ; length of a tail to put on each point which represents a galaxy when prev=1
		FOR ti=0,n_elements(time)-1 DO BEGIN  ;; loop over time
			
			count=count+1
			SETCT,1,n_elements(styles)
			IF(strt[k] EQ 0) THEN $
			  PLOT,[0],[0],COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=labels[0],YTITLE=labels[k], $ 
			    XRANGE=ranges[*,0],YRANGE=ranges[*,k],ylog=wrtXlog[k],xlog=wrtXlog[0] ELSE $
			  PLOT,findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
			    findgen(1001)*(ranges[1,0]-ranges[0,0])/1000 + ranges[0,0], $
			    COLOR=0,BACKGROUND=255,XSTYLE=1,YSTYLE=1,XTITLE=labels[0],YTITLE=labels[k],XRANGE=ranges[*,0], $
			    YRANGE=ranges[*,k],ylog=wrtXlog[k],xlog=wrtXlog[0],linestyle=2
			XYOUTS,.1,.9,string(time[ti]),/NORMAL,COLOR=0,CHARSIZE=2	
			FOR j=0,n_elements(data[0,0,0,*])-1 DO BEGIN ;; loop over models (4th column of data)	
				IF(n_elements(styles) GT 20) THEN setct,5,n_elements(styles),j
				OPLOT, data[ti,*,0,j], data[ti,*,k,j], COLOR=colors[j],linestyle=styles[j],PSYM=psym
				IF(n_elements(prev) NE 0) THEN BEGIN
					FOR ti2=MAX([0,ti-tailLength]),ti DO BEGIN ;; loop over previous time steps
						OPLOT, data[ti2,*,0,j],data[ti2,*,k,j], COLOR=colors[j],PSYM=3,symsize=2
					ENDFOR
					OPLOT,data[MAX([0,ti-tailLength]):ti,0,0,j],data[MAX([0,ti-tailLength]):ti,0,k,j], COLOR=colors[j],linestyle=0
				ENDIF
			ENDFOR
			

			IF(sv EQ 1) THEN MPEG_PUT,mpeg_id,WINDOW=1,FRAME=count,/ORDER
			IF(sv EQ 2) THEN WRITE_PNG, dn+'/frame_'+string(count,FORMAT="(I4.4)") + '.png',TVRD(TRUE=1)
			IF(sv EQ 3) THEN dummy=cgSnapshot(filename=(dn+'/frame_'+STRING(count,FORMAT="(I4.4)")),/png,/true,/nodialog)
			IF(sv EQ 0) THEN wait,.05
		ENDFOR
		IF(sv EQ 1) THEN MPEG_SAVE,mpeg_id
		IF(sv EQ 1) THEN MPEG_CLOSE,mpeg_id
		IF(sv EQ 2 || sv EQ 3) THEN spawn,("rm " + fn+".mpg") 
		IF(sv EQ 2 || sv EQ 3) THEN spawn,("ffmpeg -f image2 -qscale 4 -i "+dn+"/frame_%04d.png " + fn + ".mpg")
		IF(sv EQ 2 || sv EQ 3) THEN spawn,("rm -rf "+dn)
	
		set_plot,'x'

	ENDFOR
END

;; given an array whose first column indexes time, second column indexes different data sets
;;  third column indexes models, plus corresponding time
;; axis labels, colors, and styles, make a series of simple plots
PRO simplePlot,data,time,labels,colors,styles,name,sv
	FOR k=0,n_elements(labels)-1 DO BEGIN;; loop over y-axis variables to be plotted against time
		set_plot,'x'
		figureInit,name+"_"+labels[k],1,1
		
		FOR j=0,n_elements(data[0,0,*])-1 DO BEGIN
			IF(j EQ 0) THEN PLOT,time,data[*,k,j],COLOR=colors[j],BACKGROUND=255,linestyle=styles[j],XTITLE="time",YTITLE=labels[k],XSTYLE=1,YSTYLE=1
			IF(j NE 0) THEN OPLOT,time,data[*,k,j],COLOR=colors[j],linestyle=styles[j]
		ENDFOR
		
		figureClean,name+"_"+labels[k]
	ENDFOR
END

FUNCTION getElements,nvardim,i
	telements=(unpack(i,nVarDim))
	srtd=nVarDim[reverse(sort(nvardim))]
	elements=telements[reverse(sort(nvardim))]
	IF(n_elements(elements) LT 3) THEN BEGIN
		qelements=intarr(3)
		qelements[0:n_elements(elements)-1]=elements[0:n_elements(elements)-1]
		elements=qelements
	ENDIF

	RETURN,elements
END

PRO variability2,expName,keys,N,sv
	compile_opt idl2

	;; take the experiment name and make a list of the names of all 
	;; runs in that experiment which were not terminated due to some
	;; fatal error!
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

	IF(N LT ctr) THEN nameList2=nameList2[0:N-1]
	;; Now 
	ctr=0

	;; variables to plot- title, column index, log plot, name, y-range
	wrtXyt=['R','Col','Sig','Col_*','Sig_*','fg','Z','ColSFR','Q','Qg','Qst','fH2']
	wrtXyy=[9,10,11,12,13,17,22,14,12,24,23,48]
	wrtXyp=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0,0]
	wrtXyn=['r', 'col','sig','colst','sigst','fg','Z','colsfr','Q','Qg','Qst','fH2']
	wrtXyr=[[0,20],[.1,1000],[5,300],[.1,1000],[5,300],[0,1],[-1,1.0],[1d-6,10],[1.5,5.0],[.5,50],[.5,50],[0,1]]
	wrtXyl=[0,1,1,1,1,0,0,1,0,1,1,0]


	tk = [[0],keys]
	wrtXyt=wrtXyt[tk]
	wrtXyy=wrtXyy[tk]
	wrtXyp=wrtXyp[tk]
	wrtXyn=wrtXyn[tk]
	wrtXyr=wrtXyr[*,tk]
	wrtXyl=wrtXyl[tk]

	;; now loop over every model in the list - 
	;; all we're doing here is reading in each model then freeing 
	;; the associated memory
	;; so that we can iterate over a very large number. For each model, 
	;; record the data specified in the arrays above (e.g. yy)
	modelPtr=ReadModels([nameList2[0]])
	model = (*(modelPtr[0]))

	time = model.evArray[1,*]

	;; (# timesteps) x (nx) x (# of columns) x (# models)
	theData= dblarr(n_elements(model.dataCube[*,0,0]),n_elements(model.dataCube[0,*,0]),n_elements(wrtXyy),n_elements(nameList2))

	;; (# timesteps) x (# pts/model frame = 1) x (one thing vs another = 4) x (# models)
	vsMdot = dblarr(n_elements(model.dataCube[*,0,0]),1,10,n_elements(nameList2))

;	theTData = dblarr(n_elements(model.dataCube[*,0,0]), n_elements(wrtTyy), n_elements(nameList2))

	FOR i=0,n_elements(nameList2)-1 DO BEGIN ;; loop over every valid model
		IF(i MOD 50 EQ 0) THEN PRINT,"Processing model #",i," of ",n_elements(nameList2)

		;; the number of colors necessary is srtd[0], 
		;; which color set to use is elements[2]
		IF(i NE 0) THEN modelPtr=ReadModels([nameList2[i]])
		IF(i NE 0) THEN model = (*(modelPtr[0]))
		zi=nearest("redshift",1,model,0.0)
		offsets=[0,model.ncolstep,model.ncolstep+model.npostprocess]
		
		x = model.dataCube[0,*,0]

		;; if this model has the appropriate number of time steps..
		IF(n_elements(model.dataCube[*,0,0]) EQ n_elements(theData[*,0,0,0])) THEN BEGIN
			FOR k=0,n_elements(wrtXyy)-1 DO BEGIN	
				theData[*,*,k,i] = model.dataCube[*,*,wrtXyy[k]+offsets[wrtXyp[k]]-1]
				
			ENDFOR
			
			vsMdot[*,0,0,i] = - model.dataCube[*,model.nx-1,3-1] * model.mdotext0 ;; external accretion rate
;			FOR theT=0, n_elements(model.dataCube[*,0,0])-1 DO $
;				vsMdot[theT,0,1,i] = TOTAL(2*!pi*x[*]*x[*]*model.dlnx*model.dataCube[theT,*,model.ncolstep+14-1]*model.Radius*model.Radius)

			vsMdot[*,0,1,i] = model.evArray[10,*] ;; SFR
			vsMdot[*,0,2,i] = model.evArray[5,*]/model.evArray[6,*] * (1.0 - model.evArray[6,*])

			FOR zi=0, n_elements(model.evArray[9,*])-1 DO BEGIN ; loop over redshift
				modelInfo = diskStats(model,z=model.evArray[9,zi])
				vsMdot[zi,0,3,i] = modelInfo[1] ;; peak radius of column density (kpc)
				vsMdot[zi,0,4,i] = modelInfo[2] ;; edge of SF region (kpc)
				vsMdot[zi,0,5,i] = modelInfo[3] ;; col @ r=8kpc (solar masses/pc^2)
				vsMdot[zi,0,6,i] = modelInfo[7] ;; sSFR (yr^-1)
				vsMdot[zi,0,7,i] = modelInfo[8] ;; mass in bulge (MSol)
				vsMdot[zi,0,8,i] = modelInfo[9] ;; fSt in bulge
				vsMdot[zi,0,9,i] = modelInfo[10];; total f_H2
			ENDFOR


;			FOR kk=0,n_elements(wrtTyind)-1 DO BEGIN
;				theTData[*,kk,i] =;; model.evArray[wrtTyind[kk],*]
;			ENDFOR
			ctr = ctr+1
		ENDIF
		

;			PlotVsXInd,model,zi,model.ncolstep+9,yind[k]+offsets[yp[k]],OP(i),0,COLOR=elements[0],LINESTYLE=(elements[1]+fix(elements[1] GE 1)),THICK=1,XTITLE="r [kpc]",YTITLE=ytitle[k],YRANGE=yr[*,k],YSTYLE=1,XSTYLE=1


		IF(i EQ 0) THEN print,model.name,diskStats(model)
		ptr_free,modelPtr[0]
	ENDFOR

	theData = (temporary(theData))[*,*,*,0:ctr-1]
;	theTData = temporary(theTData)[*,*,0:ctr-1]


	;; having accumulated theData, let's do some analysis!
	sortdData = theData * 0.0
	sortdVsMdot = vsMdot*0.0
	FOR j=0,n_elements(theData[*,0,0,0])-1 DO BEGIN ;; loop over time
		FOR k=0,n_elements(wrtXyy)-1 DO BEGIN ;; loop over variables
			FOR xx=0,n_elements(theData[0,*,0,0])-1 DO BEGIN ;; loop over x-coord
				sortdData[j,xx,k,*] = theData[j,xx,k,sort(theData[j,xx,k,*])]
			ENDFOR
			
		ENDFOR
		st = sort(vsMdot[0,0,0,*])
		FOR k2=0,n_elements(vsMdot[0,0,*,0])-1 DO BEGIN	;; loop over variables
			sortdVsMdot[j,0,k2,*] = vsMdot[j,0,k2,st]
		ENDFOR
	ENDFOR


	intervals = dblarr(n_elements(theData[*,0,0,0]),  n_elements(theData[0,*,0,0]), n_elements(theData[0,0,*,0]), 5)
	intervals[*,*,*, 0]=sortdData[*,*,*, fix(.025*ctr)]
	intervals[*,*,*, 1]=sortdData[*,*,*, fix(.16*ctr)]
	intervals[*,*,*, 2]=sortdData[*,*,*, fix(.5*ctr)]
	intervals[*,*,*, 3]=sortdData[*,*,*, fix(.84*ctr)]
	intervals[*,*,*, 4]=sortdData[*,*,*, fix(.975*ctr)]

	expName2 = expName+'_'+strcompress(string(MIN([n_elements(nameList2),N])),/remove)

	setct,1,0,0
	simpleMovie,intervals,time,wrtXyn,[2,1,0,1,2],[2,2,0,2,2],wrtXyr,wrtXyl,expName2+"_intervals",sv,intarr(n_elements(wrtXyl))

	

	setct,3,n_elements(sortdData[0,0,0,*]),0
	simpleMovie,sortdData,time,wrtXyn,indgen(n_elements(sortdData[0,0,0,*])),intarr(n_elements(sortdData[0,0,0,*])),wrtXyr,wrtXyl,expName2+"_sorted",sv,intarr(n_elements(sortdData[0,0,0,*]))

	setct,3,n_elements(theData[0,0,0,*]),0
	simpleMovie,theData,time,wrtXyn,indgen(n_elements(theData[0,0,0,*])),intarr(n_elements(theData[0,0,0,*])),wrtXyr,wrtXyl,expName2+"_unsorted",sv,intarr(n_elements(theData[0,0,0,*]))

	

	setct,3,n_elements(vsMdot[0,0,0,*]),0
	simpleMovie, sortdVsMdot,time,['mdot','SFR','stMass','rPeak','rHI','colsol','sSFR','BulgeMass','BulgeFSt','fH2'],$
		indgen(n_elements(vsMdot[0,0,0,*])), $
		intarr(n_elements(vsMdot[0,0,0,*])), $
		[[.03,500],[1,100],[3d9,3d10],[0.2,20],[.2,22],[6,100],[3d-11,3d-8],[1d7,3d11],[-0.1,1.1],[-0.1,1.1]], $
		[1,1,1,1,0,1,1,1,0,0],expName2+"_vsmdot",sv,[1,1,0,0,0,0,0,0,0,0],PSYM=1,prev=1

	;; average mdot over a number of time steps nt
	nt = 50
	setct,3,n_elements(vsMdot[0,0,0,*]),0
	avgdVsMdot = sortdVsMdot[*,*,*,*]
	FOR ti=0,n_elements(time)-1 DO BEGIN
		FOR mm=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
			avgdVsMdot[ti,0,0,mm] = TOTAL(sortdVsMdot[MAX([0,ti-nt]):ti,0,0,mm]) / (ti - MAX([0,ti-nt]) + 1)
		ENDFOR
	ENDFOR
	simpleMovie, avgdVsMdot,time,['avgMdot','SFR','stMass','rPeak','rHI','colsol','sSFR','BulgeMass','BulgeFSt','fH2'],$
		indgen(n_elements(vsMdot[0,0,0,*])), $
		intarr(n_elements(vsMdot[0,0,0,*])), $
		[[.03,500],[1,100],[3d9,3d10],[0.2,20],[.2,22],[6,100],[3d-11,1d-8],[1d7,3d11],[-0.1,1.1],[-0.1,1.1]], $
		[1,1,1,1,0,1,1,1,0,0],expName2+"_vsavgmdot",sv,[1,1,0,0,0,0,0,0,0,0],PSYM=1,prev=1


	setct,3,n_elements(vsMdot[0,0,0,*]),0
	vsSFR = sortdVsMdot[*,*,1:n_elements(sortdVsMdot[0,0,*,0])-1,*] ;; cut out the mdot slice of the array
	simpleMovie, vsSFR, time, ['SFR','stMass','rPeak','rHI','colsol','sSFR','BulgeMass','BulgeFSt','fH2'], $
		indgen(n_elements(vsMdot[0,0,0,*])), $
		intarr(n_elements(vsMdot[0,0,0,*])), $
		[[1,100],[3d9,3d10],[0.2,20],[.2,22],[6,100],[3d-11,1d-8],[1d7,3d11],[-0.1,1.1],[-0.1,1.1]], $
		[1,1,1,0,1,1,1,0,0],expName2+"_vsSFR",sv,[1,0,0,0,0,0,0,0,0],PSYM=1,prev=1

	
	SETCT,3,n_elements(vsMdot[0,0,0,*]),0
	FIGUREInit,(expName2+'_accRates'),1,1,1
	Plot,time,sortdVsMdot[*,0,0,0],COLOR=0,BACKGROUND=255,YRANGE=[.03,200],YLOG=1,XSTYLE=1,XTITLE="Time",YTITLE="Accretion Rate"
	FOR m=1,n_elements(sortdVsMdot[0,0,0,*])-1 DO BEGIN
		setct,5,n_elements(sortdVsMdot[0,0,0,*]),m
		OPLOT,time,sortdVsMdot[*,0,0,m],COLOR=m
		
	ENDFOR
	FIGURECLEAN,(expName2+'_accRates'),1

	FIGUREInit,(expName+'_z0AccRates'),1,1,1
	cgHistoplot, alog10(sortdVsMdot[n_elements(sortdVsMdot[*,0,0,0])-1,0,0,*]), BinSize=.1
	FIGURECLEAN,(expName2+'_z0AccRates'),1
	


END
