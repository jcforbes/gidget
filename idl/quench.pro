PRO Quench,modelList,sv
	compile_opt idl2
	model=(*(modelList[0]))
	SETCT,1,0,0
	figureInit,(model.name+'_quench1'),sv,1,2
	multi=1
	IF(multi EQ 0) THEN !p.multi=[0,3,1]
	IF(multi EQ 1) THEN BEGIN
		multiplot,[0,1,3,0,0],ygap=.003
		!p.noerase=0
	ENDIF

	COMMON modelBlock,theModel
	theModel=(*(modelList[0]))

	print,"Model 0 parameters: ",model.diskScaleLength,model.MLF,model.Qlim,model.eta

	xt="Time since z=2 [Gyr]"
	cs=.7
	xtvi=[nearest("redshift",1,model,2.0),nearest("redshift",1,model,1.5),nearest("redshift",1,model,1.0),nearest("redshift",1,model,0.5),nearest("redshift",1,model,.001)]
	xtv=model.evArray[1,xtvi]
	PlotVsT,model,11,-1,0,YRANGE=[1,40],ylog=1,YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE="Star Formation Rate [MSol/yr]",CHARSIZE=cs,XSTYLE=9,ymargin=[4,6],THICK=3
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		PlotVsT,*(modelList[i]),11,-1,1,COLOR=i
	ENDFOR

	axis,XSTYLE=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4

	IF(multi EQ 1) THEN BEGIN
		multiplot,ygap=0.0
		!p.noerase=1
	ENDIF

	totalMass = model.evArray[5,*]/model.evArray[6,*] ;; gas mass * 1/f_g	
	Plot,model.evArray[1,*],totalMass,YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE="Total Disk Mass [MSol]",CHARSIZE=cs,XSTYLE=1,THICK=3,YRANGE=[5.0d9,3.0d10],ylog=1
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		oplot,(*(modelList[i])).evArray[1,*],  (*(modelList[i])).evArray[5,*]/(*(modelList[i])).evArray[6,*], COLOR=i
	ENDFOR

	IF(multi EQ 1) THEN multiplot

	sSFR= model.evArray[11-1,*]/totalMass
	Plot,model.evArray[1,*],sSFR,YSTYLE=1,COLOR=0,BACKGROUND=255,XTITLE=xt,YTITLE="SFR/Total Mass [yr^-1]",CHARSIZE=cs,XSTYLE=1,THICK=3,YRANGE=[1.0d-11,1.0d-8],ylog=1
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		oplot,(*(modelList[i])).evArray[1,*], (*(modelList[i])).evArray[11-1,*] * (*(modelList[i])).evArray[6,*] / (*(modelList[i])).evArray[5,*], COLOR=i
	ENDFOR

	figureClean,(model.name+'_quench1'),sv
	IF(sv EQ 2) THEN latexify,(model.name+'_quench1.eps'),["Total Disk Mass [MSol]","Star Formation Rate[MSol/yr]","SFR/Total Mass [yr^-1]","Stellar Disk Mass [MSol]","SFR/Stellar Mass [yr^-1]","Redshift",xt],["Total Disk Mass [$M_\odot$]","SFR [$M_\odot/yr$]","SFR/Total Mass [yr$^{-1}$]","Stellar Disk Mass [$M_\odot$]","SFR/Stellar Mass [yr$^{-1}$]","Redshift",xt],replicate(1.3,7),/full
	IF(multi EQ 0) THEN !p.multi=[0]
	IF(multi EQ 1) THEN multiplot,/default	


	figureInit,(model.name+'_quench2'),sv,1,2
	multi=1
	IF(multi EQ 0) THEN !p.multi=[0,3,1]
	IF(multi EQ 1) THEN BEGIN
		multiplot,[0,1,3,0,0],ygap=.003
		!p.noerase=0
	ENDIF
	PlotVsT,model,11,-1,0,YRANGE=[1,40],ylog=1,YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE="Star Formation Rate [MSol/yr]",CHARSIZE=cs,XSTYLE=9,ymargin=[4,6],THICK=3
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		PlotVsT,*(modelList[i]),11,-1,1,COLOR=i
	ENDFOR

	axis,XSTYLE=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4

	IF(multi EQ 1) THEN BEGIN
		multiplot,ygap=0.0
		!p.noerase=1
	ENDIF


	stellarMass = model.evArray[5,*]/model.evArray[6,*] * (1.0-model.evArray[6,*]) ;; gas mass * (1/f_g)* (1-f_g)
	Plot,model.evArray[1,*],stellarMass,YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE="Stellar Disk Mass [MSol]",CHARSIZE=cs,XSTYLE=1,THICK=3,ylog=1,YRANGE=[5.0d9,3.0d10]
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		fs=- (*(modelList[i])).evArray[6,*] + 1.0
		a= ((*(modelList[i])).evArray[5,*])/((*(modelList[i])).evArray[6,*]) * fs
		oplot,(*(modelList[i])).evArray[1,*], a, COLOR=i
	ENDFOR

	IF(multi EQ 1) THEN multiplot

	SFRoSTM= model.evArray[11-1,*]/stellarMass
	Plot,model.evArray[1,*],SFRoSTM,YSTYLE=1,COLOR=0,BACKGROUND=255,XTITLE=xt,YTITLE="SFR/Stellar Mass [yr^-1]",CHARSIZE=cs,XSTYLE=1,THICK=3,ylog=1,YRANGE=[1.0d-11,1.0d-8]
	FOR i=1,n_elements(modelList)-1 DO BEGIN
		oplot,(*(modelList[i])).evArray[1,*], (*(modelList[i])).evArray[11-1,*] * (*(modelList[i])).evArray[6,*] / ((*(modelList[i])).evArray[5,*] * (1.0 - (*(modelList[i])).evArray[6,*])), COLOR=i
	ENDFOR




	figureClean,(model.name+'_quench2'),sv
	IF(sv EQ 2) THEN latexify,(model.name+'_quench2.eps'),["Total Disk Mass [MSol]","Star Formation Rate[MSol/yr]","SFR/Total Mass [yr^-1]","Stellar Disk Mass [MSol]","SFR/Stellar Mass [yr^-1]","Redshift",xt],["Total Disk Mass [$M_\odot$]","SFR [$M_\odot/yr$]","SFR/Total Mass [yr$^{-1}$]","Stellar Disk Mass [$M_\odot$]","SFR/Stellar Mass [yr$^{-1}$]","Redshift",xt],replicate(1.3,7),/full
	IF(multi EQ 0) THEN !p.multi=[0]
	IF(multi EQ 1) THEN multiplot,/default
END

