;; Some functions related to plotting a quantity of some stellar population,
;; e.g. its column density, vs its age.


;; model - a model object
;; index - refers to S, s, Z, s_Z between 0-3 
;; position - number between 0 and 1 referring to x-position in the disk
;; time - the time at which we should plot this quantity (in orbital periods since the beginning of the simulation
;;  overPlot - 1 to overplot, 0 to make a new plot
;; mult - a constant by which to multiply the dependent variable
PRO PlotVsAge,model,index,position,time,overPlot,mult,_EXTRA=EXTRA_KEYWORDS
        tne = n_elements(model.stellarages[nearest("time",1,model,time),*])
        
        IF(overPlot EQ 1) THEN $
                oplot,model.stellarages[nearest("time",1,model,time),$
			1:(tne-1)]*1e-9,$
			model.stellarData[nearest("time",1,model,time),$
			1:(tne-1),nearest("position",1,model,position),$
			index]*mult,_EXTRA=EXTRA_KEYWORDS $
        ELSE BEGIN
		plot,model.stellarages[nearest("time",1,model,time),$
		  1:(tne-1)]*1e-9,model.stellarData[nearest("time",1,$
		  model,time),1:(tne-1),nearest("position",1,model,$
		  position),index]*mult,_EXTRA=EXTRA_KEYWORDS
	END


;       IF(overPlot EQ 1) THEN $
;               oplot,model.stellarages[nearest("time",1,model,time),0:(tne-1)]*1e-9,model.stellarData[nearest("time",1,model,time),0:(tne-1),nearest("position",1,model,position),index]*mult,_EXTRA=EXTRA_KEYWORDS $
;       ELSE plot,model.stellarages[nearest("time",1,model,time),0:(tne-1)]*1e-9,model.stellarData[nearest("time",1,model,time),0:(tne-1),nearest("position",1,model,position),index]*mult,_EXTRA=EXTRA_KEYWORDS

END

FUNCTION GetStellarRange,modelList,position,index
        minimum=1.e50
        maximum=0.
        FOR j=0,n_elements(modelList)-1 DO BEGIN
                tne=n_elements((*(modelList[j])).stellarages[nearest("time",1,*(modelList[j]),40),*])
                RR = (*(modelList[j])).Radius
                minimum=MIN([MIN((*(modelList[j])).stellarData[*,1:(tne-1),nearest("position",1,*(modelList[j]),position/RR),index]),minimum])
                maximum=MAX([MAX((*(modelList[j])).stellarData[*,1:(tne-1),nearest("position",1,*(modelList[j]),position/RR),index]),maximum])
        ENDFOR        
	RETURN,[minimum,maximum]
END



;; modelList is an array of pointers to models, generally obtained form readModels
;; position is a radius (in kpc) in the disk at which to plot the requested quantities
;; name - a string to include in the filename
;; comment - a string to write on the plot
;; sv - a number 0-2 resp. print to the screen, print to png, or print to eps
;; legend - a list of strings to label the different models
PRO PlotStarsVsAge,modelList,position,name,comment,sv,legend=legend
        fn=""
	IF(n_elements(modelList) LT 6) THEN $
	        FOR j=0,n_elements(modelList)-1 DO $
			fn+=(ExtractSlash((modelList[j]),1)+"_") $
	ELSE fn+= ExtractSlash((modelList[0]),0)+"_"
        fn= fn+"_"+name
        figureInit,fn,sv,1,2
        multi=1

        IF(multi EQ 1) THEN multiplot,[0,1,3,0,0]
        cs=.7
        !p.noerase=0

        ;; Convert from dimensionless to dimensional column density
        conv=dblarr(n_elements(modelList))
        FOR j=0,n_elements(modelList)-1 DO BEGIN
                conv[j] = ((*(modelList[j])).mdotext0)/$
			(((*(modelList[j])).Radius)*$
			((*(modelList[j])).vphiR)) * 977.813952 
        ENDFOR
        poss = dblarr(n_elements(modelList))
        FOR i=0, n_elements(modelList)-1 DO BEGIN
		poss[i]=nearest("position",0,*(modelList[i]),$ 
			(position/((*(modelList[i])).Radius)))
	ENDFOR

        tt = "x="+strcompress(string(position),/remove)+" "+comment+ " "
        FOR j=0,n_elements(legend)-1 DO tt += (", " + legend[j] )

	lsj=intarr(n_elements(modelList))
	FOR j=0,n_elements(lsj)-1 DO BEGIN
		IF((*(modelList[j])).whichAccHistory GE 1) THEN lsj[j] = 2 ELSE lsj[j]=0
	ENDFOR
		

	;;; Make a plot of all the models, namely stellar age vs. column density.
	;;; The structure is: plot the first model, then overplot the rest.
	;;; For each model, plot a line and a collection of *'s, i.e. PSYM=2
	;;; This shows explicitly that the plot is produced from a limited number of populations.

	; model 0: line then symbols
        PlotVsAge,*(modelList[0]),0,poss[0],$
		max((*(modelList[0])).evArray[1,*]),0,conv[0],$
		BACKGROUND=255,COLOR=0,YTITLE="Column Density [MSol/pc^2]",$
		CHARSIZE=cs,YSTYLE=1,XSTYLE=1,YRANGE=[0,4],LINESTYLE=lsj[0]
        PlotVsAge,*(modelList[0]),0,poss[0],$
		max((*(modelList[0])).evArray[1,*]),1,conv[0],$
		COLOR=0,PSYM=2
	;; all other models: symbols then line
        FOR j=1, n_elements(modelList)-1 DO $
		PlotVsAge,(*(modelList[j])),0,poss[j],$
			max((*(modelList[j])).evArray[1,*]),$
			1,conv[j],COLOR=j,PSYM=2
        FOR j=1, n_elements(modelList)-1 DO $
		PlotVsAge,(*(modelList[j])),0,poss[j],$
			max((*(modelList[j])).evArray[1,*]),$
			1,conv[j],COLOR=j,LINESTYLE=lsj[j]
        IF(multi EQ 1) THEN BEGIN
                !p.noerase=1
                multiplot
        ENDIF


	;;; Now make a plot of stellar age vs. velocity dispersion for all models
        PlotVsAge,*(modelList[0]),1,poss[0], $
		max((*(modelList[0])).evArray[1,*]),0,(*(modelList[0])).vphiR,$
		BACKGROUND=255,COLOR=0,YTITLE="Velocity Dispersion [km/s]",$
		CHARSIZE=cs,YRANGE=GetStellarRange(modelList,position,1)* ((*(modelList[0])).vphiR),$
		XSTYLE=1,YSTYLE=1,LINESTYLE=lsj[0]
        PlotVsAge,*(modelList[0]),1,poss[0],max((*(modelList[0])).evArray[1,*]),$
		1,(*(modelList[0])).vphiR,COLOR=0,PSYM=2
        FOR j=1, n_elements(modelList)-1 DO $
		PlotVsAge,*(modelList[j]),1,poss[j],max((*(modelList[j])).evArray[1,*]),$
			1,(*(modelList[0])).vphiR,COLOR=j,PSYM=2
        FOR j=1, n_elements(modelList)-1 DO $
		PlotVsAge,*(modelList[j]),1,poss[j],max((*(modelList[j])).evArray[1,*]),$
			1,(*(modelList[0])).vphiR,COLOR=j,LINESTYLE=lsj[j]

;       oploterror,ageData,vdData,errVdData,errcolor=0,color=0,psym=3,errthick=1,nohat=0

        IF(multi EQ 1) THEN multiplot

	;;; Make a plot of stellar age vs mean metallicity
        PlotVsAge,*(modelList[0]),2,poss[0],$
		max((*(modelList[0])).evArray[1,*]),0,1.0,$
		BACKGROUND=255,COLOR=0,XTITLE="Stellar Age at z=0 [Gyr]",$
		YTITLE="[<Z>]",CHARSIZE=cs,YRANGE=[-1.1,.1],XSTYLE=1,$
		YSTYLE=1,LINESTYLE=lsj[0];GetStellarRange(modelList,position,2)
        PlotVsAge,*(modelList[0]),2,poss[0],$
		max((*(modelList[0])).evArray[1,*]),1,1.0,COLOR=0,PSYM=2


        ;;; overplot other models
        FOR j=1, n_elements(modelList)-1 DO $
		PlotVsAge,*(modelList[j]),2,poss[j],$
			max((*(modelList[j])).evArray[1,*]),1,$
			1.0,COLOR=j,PSYM=2
        FOR j=1, n_elements(modelList)-1 DO $
		 PlotVsAge,*(modelList[j]),2,poss[j],$
			max((*(modelList[j])).evArray[1,*]),1,$
			1.0,COLOR=j,LINESTYLE=lsj[j]
        ns=0.
        ns2=.1

        figureClean,fn,sv
        IF(multi EQ 1) THEN BEGIN
                multiplot,/reset
                multiplot,/default
        ENDIF
        IF(sv EQ 2) THEN latexify,fn+'.eps',['Stellar Age at z=0 [Gyr]',$
		'Column Density [MSol/pc^2]','Velocity Dispersion [km/s]',$
		'[Z]'],['Stellar Age at $z=0$ [Gyr]',$
		'$\Sigma_i$ [$M_\odot/pc^2$]','$\sigma_i$ [km/s]',$
		'[$\ <Z>\ $]'],replicate(1.3,4),/full
END
