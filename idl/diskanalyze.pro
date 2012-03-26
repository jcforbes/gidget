
FUNCTION OP,j
	IF (j EQ 0) THEN v=0 ELSE v=1
	RETURN,v
END

;; functions to fit to stellar column density profiles
PRO gfunct, X, p, F, pder
        A=p[0]
        IF(A LT 0) THEN A=0
        B=p[1]
        IF(B LT 0) THEN B=0
        R0=p[2]
        k=p[3]
        n=p[4]
        F = A*exp(-x/R0) + B*exp(-k*x^(1.0/n))
        IF(n_params() GE 4) THEN BEGIN
                pder=dblarr(n_elements(x),n_elements(p))
                pder[*,0]= exp(-x/R0)
;;              pder[*,0]= 1.0 / (A + B * exp(-k*x^(1.0/n) + x/R0))
                pder[*,1]=exp(-k*x^(1.0/n))
;;              pder[*,1]= 1.0 / (B + A* exp(k * x^(1.0/n)-x/R0))
                pder[*,2]=A*exp(-x/R0)*x/(R0^2)
;;              pder[*,2]=A*exp(k*x^(1.0/n))*x/(R0*R0*(A*exp(k*x^(1.0/n))+B*exp(x/R0)))
                pder[*,3]=-B*exp(-k*x^(1.0/n))*x^(1.0/n)
;;              pder[*,3]=-B*exp(x/R0)*x^(1.0/n)/(A*exp(k*x^(1.0/n))+B*exp(x/R0))
                pder[*,4]=B*exp(-k*x^(1.0/n))*k*x^(1.0/n)*alog(x)/(n*n)
;;              pder[*,4]=B*exp(x/R0)*k*x^(1.0/n)*alog(x)/(n*n*(A*exp(k*x^(1.0/n)) + B*exp(x/R0)))
                
        ENDIF
END

PRO gfunct2, X, p, F, pder        
	A=p[0]
	R0=p[1]
        F=A*exp(-x/R0)
        IF(n_params() GE 4) THEN BEGIN
                pder=dblarr(n_elements(x),n_elements(p))
                pder[*,0] = exp(-x/R0)
                pder[*,1] = A*exp(-x/R0)*x/(R0*R0)
        ENDIF
END

PRO globalMassBudget,model,sv
        figureInit,model.name+'_gasMassBudget',sv,1,1
        time = model.evArray[1,*]
        deltaGasMass = model.evArray[11,*]
        deltaStellarMass = model.evArray[12,*]
        gasThroughIB = -1*model.evArray[13,*]
        starsThroughIB = -1*model.evArray[14,*]
        SFmass = model.evArray[15,*]
        massAcc = model.evArray[16,*]
        mu = model.MLF ;; mass loading factor
        Rf = model.Rf  ;; remnant fraction
;       mu = 0
        ;; massAcc - (1+mu)SFmass - gasThroughIB == deltaGasMass
        ;; SFmass - starsThroughIB == deltaStellarMass
        all=[deltaGasMass,gasThroughIB,(Rf+mu)*SFMass,deltaGasMass+gasThroughIB+(1+mu)*SFMass,massAcc]
        plot,time,massAcc,BACKGROUND=255,COLOR=0,THICK=3,LINESTYLE=2,XTITLE="Time (Ga)",YTITLE="Mass (MSol)",YRANGE=[MIN(all),MAX(all)]
        oplot,time,deltaGasMass,COLOR=1
        oplot,time,gasThroughIB,COLOR=2
        oplot,time,(Rf+mu)*SFmass,COLOR=3
        oplot,time,deltaGasMass+gasThroughIB+(Rf+mu)*SFmass,COLOR=0,THICK=1
        legend,["ExtMassAcc","Curr-Init Gas","Gas Thru IB","(1+mu)SF","Sum"],COLORS=[0,1,2,3,0],LINESTYLE=[2,0,0,0,0],THICK=[3,1,1,1,1],TEXTCOLORS=[0,0,0,0,0]
        figureClean,model.name+'_gasMassBudget',sv

        figureInit,model.name+'_stMassBudget',sv,1,1
        all=[SFmass,deltaStellarMass,starsThroughIB,deltaStellarMass+starsThroughIB]
        plot,time,SFmass,BACKGROUND=255,COLOR=0,THICK=4,LINESTYLE=2,XTITLE="Time (Ga)",YTITLE="Mass (MSol)",YRANGE=[MIN(all),MAX(all)]
        oplot,time,deltaStellarMass,COLOR=1
        oplot,time,starsThroughIB,COLOR=2
        oplot,time,deltaStellarMass+starsThroughIB,COLOR=0,THICK=2
        legend,["SF","Curr-Init Stars","stars thru IB","Sum"],COLORS=[0,1,2,0],THICK=[4,1,1,2],LINESTYLE=[2,0,0,0],TEXTCOLORS=[0,0,0,0]        
	figureClean,model.name+'_stMassBudget',sv

;       newWindow
;       plot,time,(SFmass-starsThroughIB-deltaStellarMass),COLOR=0,BACKGROUND=255,XTITLE="Time (Ga)",YTITLE="Stellar Mass (MSol)"
;       oplot,time,starsThroughIB,COLOR=2       
;       oplot,time,.01*deltaGasMass,COLOR=1
;;       newWindow;       plot,time,(deltaGasMass+gasThroughIB+(1+mu)*SFmass - massAcc),COLOR=0,BACKGROUND=255,XTITLE="Time (Ga)",YTITLE="Gas Mass (MSol)"
;       oplot,time,gasThroughIB,COLOR=2
;       oplot,time,-.01*deltaGasMass,COLOR=1
        
END

PRO stellarDisk,model,sv
        figureInit,(model.name+'_pub_stdisk'),sv,1,1
        zi=nearest("redshift",1,model,0.0)
        fitThese = where(model.dataCube[zi,*,0] LE 2,ctr)
        IF(ctr EQ 0) THEN RETURN
        r=dblarr(n_elements(fitThese))
        r[*]=model.dataCube[zi,fitThese,0]*model.Radius 
        col_st=model.dataCube[zi,fitThese,5]
        weights=r
        ;; initial guess
        p=[col_st[0],3]
        colfit = CURVEFIT(r,col_st,weights,p,FUNCTION_NAME='gfunct2',ITMAX=500)
        PRINT,'Function parameters: ',p
        PLOT,r,col_st,/ylog,background=255,color=0,xtitle="Radius (kpc)",YTITLE="Column Density"
        Oplot,r,p[0]*exp(-r/p[1]),linestyle=3,color=0
        figureClean,(model.name+'_pub_stdisk'),sv
END


PRO bulgeDisk,model,sv
        figureInit,(model.name+'_pub_bd'),sv,1,1
        zi=nearest("redshift",1,model,0.0)
        fitThese = where(model.dataCube[zi,*,22] LE 2.0)
        IF(n_elements(fitThese) EQ 1) THEN IF(fitThese EQ -1L) THEN RETURN
        r=dblarr(n_elements(fitThese))
        r[*]=model.dataCube[zi,fitThese,0]*model.Radius
        col_st=model.dataCube[zi,fitThese,5]
        weights=r
        ;; initial guess
        p=[.5*col_st[0],.5*col_st[0],3.5,1.0,3.0]       
        help,fitThese
        help,r
        help,col_st
        colfit = CURVEFIT(r,col_st,weights,p,FUNCTION_NAME='gfunct',ITMAX=500)
        Print,'Function parameters: ',p
        bu = dblarr(n_elements(r))
        bu= p[1]*exp(-p[3]*r^(1.0/p[4]))
        di = dblarr(n_elements(r))
        di =p[0]*exp(-r/p[2])
        bd = TOTAL(2*!pi*r*r*bu)/Total(2*!pi*r*r*di)
        Print,'Bulge-to-Disk ratio: ',bd
;;      newWindow
        Plot,r,col_st,/ylog,background=255,color=0,xtitle="Radius (kpc)",YTITLE="Column Density"
        Oplot,r,bu+di,linestyle=2,THICK=2,color=0
        Oplot,r,bu,linestyle=3,color=0
        Oplot,r,di,linestyle=3,color=0
        legend,["Stellar Column Density","Fit","Bulge","Disk"],COLORS=[0,0,0,0],LINESTYLE=[0,2,3,3],THICK=[0,2,0,0],TEXTCOLORS=[0,0,0,0]

        figureClean,(model.name+'_pub_bd'),sv
END

PRO GlobalTimeSeries,model,sv
        figureInit,(model.name+'_pub_ts'),sv,1,2
        multi=1
        IF(multi EQ 0) THEN !p.multi=[0,3,1]
        IF(multi EQ 1) THEN BEGIN
                multiplot,[0,1,3,0,0],ygap=.003
                !p.noerase=0
        ENDIF
        xt="Time since z=2 [Gyr]"
        cs=.7
        xtvi=[nearest("redshift",1,model,2.0),nearest("redshift",1,model,1.5),nearest("redshift",1,model,1.0),nearest("redshift",1,model,0.5),nearest("redshift",1,model,.001)]
        xtv=model.evArray[1,xtvi]
        PlotVsT,model,7,-1,0,YRANGE=[0,.6],YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE="Gas Fraction",CHARSIZE=cs,XSTYLE=9,YMARGIN=[4,6]
        axis,XSTYLE=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4

        IF(multi EQ 1) THEN BEGIN
                 multiplot,ygap=0.0
                !p.noerase=1
        ENDIF

        PlotVsT,model,model.ncolstep+8,8/model.Radius,0,YRANGE=[8d6,2d9],YSTYLE=1,COLOR=0,BACKGROUND=255,YTITLE='2D Jeans Mass at r=8 kpc [MSun]',CHARSIZE=cs,XSTYLE=9,ylog=1,ymargin=[4,6]
;       axis,XSTYLE=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4

        IF(multi EQ 1) THEN multiplot

        PlotVsT,model,11,-1,0,YRANGE=[1,25],ylog=1,YSTYLE=1,COLOR=0,BACKGROUND=255,XTITLE="Time since z=2 [Gyr]",YTITLE="Star Formation Rate [MSol/yr]",CHARSIZE=cs,XSTYLE=9,ymargin=[4,6];       
;	axis,XSTYLE=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4

        figureClean,(model.name+'_pub_ts'),sv
        IF(sv EQ 2) THEN latexify,(model.name+'_pub_ts.eps'),['Time since z=2 [Gyr]','Gas Fraction','2D Jeans Mass at r=8 kpc [MSun]','Star Formation Rate [MSol/yr]'],['Time since $z=2$ [Gyr]','Gas Fraction','2D Jeans Mass at $r=8$ kpc [$M_\odot$]','Star Formation Rate [$M_\odot/yr$]'],replicate(1.3,4),/full
        IF(multi EQ 0) THEN !p.multi=[0]        
	IF(multi EQ 1) THEN multiplot,/default
END

PRO MdotSFR,model,sv
	figureInit,(model.name+'_mdotsfr'),sv,1,1
	xt="Time since z=2 [Gyr]"
	cs=1
	xtvi=[nearest("redshift",1,model,2.0),nearest("redshift",1,model,1.5),nearest("redshift",1,model,1.0),nearest("redshift",1,model,0.5),nearest("redshift",1,model,.001)]
	xtv=model.evArray[1,xtvi]

	t=model.evArray[1,*]
	x=model.dataCube[0,*,0]
	dlnx=model.dlnx
        x2=nearest("position",1,model,2.0/model.Radius)
        x4=nearest("position",1,model,4.0/model.Radius)
        x8=nearest("position",1,model,8.0/model.Radius)
        x12=nearest("position",1,model,12.0/model.Radius)
        x20=nearest("position",1,model,20.0/model.Radius)
        xo=nearest("position",1,model,1)

;        PRINT,"Total SFR within 2kpc",TOTAL(2*!pi*x[0:x2-1]*x[0:x2-1]*dlnx*model.dataCube[z2,0:x2-1,model.ncolstep+14-1]*model.Radius*model.Radius)," MSol/yr"
;        PRINT,"Total SFR outside 2kpc",TOTAL(2*!pi*x[x2-1:xo-1]*x[x2-1:xo-1]*dlnx*model.dataCube[z2,x2-1:xo-1,model.ncolstep+14-1]*model.Radius*model.Radius),"MSol/yr"
;        PRINT,"mass flux across 2,4,8,16,20 kpc",-model.dataCube[z2,x2-1,3-1]*model.mdotext0,-model.dataCube[z2,x4-1,3-1]*model.mdotext0,-model.dataCube[z2,x8-1,3-1]*model.mdotext0,-model.dataCube[z2,x16-1,3-1]*model.mdotext0,-model.dataCube[z2,x20-1,3-1]*model.mdotext0," MSol/yr"

	SFR2=dblarr(n_elements(model.dataCube[*,0,0]))
	SFR4=SFR2[*]
	SFR8=SFR2[*]
	SFR12=SFR2[*]
	SFR20=SFR2[*]
	Mdot2=-model.dataCube[*,x2-1,3-1]*model.mdotext0
	Mdot4=-model.dataCube[*,x4-1,3-1]*model.mdotext0
	Mdot8=-model.dataCube[*,x8-1,3-1]*model.mdotext0
	Mdot12=-model.dataCube[*,x12-1,3-1]*model.mdotext0
	Mdot20=-model.dataCube[*,x20-1,3-1]*model.mdotext0
	
	FOR zz=0,n_elements(SFR2)-1 DO BEGIN
		SFR2[zz] = TOTAL(2*!pi*x[0:x2-1]*x[0:x2-1]*dlnx*model.dataCube[zz,0:x2-1,model.ncolstep+14-1]*model.Radius*model.Radius)
		SFR4[zz] = TOTAL(2*!pi*x[0:x4-1]*x[0:x4-1]*dlnx*model.dataCube[zz,0:x4-1,model.ncolstep+14-1]*model.Radius*model.Radius)
		SFR8[zz] = TOTAL(2*!pi*x[0:x8-1]*x[0:x8-1]*dlnx*model.dataCube[zz,0:x8-1,model.ncolstep+14-1]*model.Radius*model.Radius)
		SFR12[zz] = TOTAL(2*!pi*x[0:x12-1]*x[0:x12-1]*dlnx*model.dataCube[zz,0:x12-1,model.ncolstep+14-1]*model.Radius*model.Radius)
		SFR20[zz] = TOTAL(2*!pi*x[0:x20-1]*x[0:x20-1]*dlnx*model.dataCube[zz,0:x20-1,model.ncolstep+14-1]*model.Radius*model.Radius)
;		sSFR2[zz] = TOTAL(2*!pi*x[0:x2-1]*x[0:x2-1]*dlnx*model.dataCube[zz,0:x2-1,model.ncolstep+14-1]*model.Radius*model.Radius)
;		sSFR4[zz] = TOTAL(2*!pi*x[0:x4-1]*x[0:x4-1]*dlnx*model.dataCube[zz,0:x4-1,model.ncolstep+14-1]*model.Radius*model.Radius)
;		sSFR8[zz] = TOTAL(2*!pi*x[0:x8-1]*x[0:x8-1]*dlnx*model.dataCube[zz,0:x8-1,model.ncolstep+14-1]*model.Radius*model.Radius)
;		sSFR12[zz] = TOTAL(2*!pi*x[0:x12-1]*x[0:x12-1]*dlnx*model.dataCube[zz,0:x12-1,model.ncolstep+14-1]*model.Radius*model.Radius)
;		sSFR20[zz] = TOTAL(2*!pi*x[0:x20-1]*x[0:x20-1]*dlnx*model.dataCube[zz,0:x20-1,model.ncolstep+14-1]*model.Radius*model.Radius)

	ENDFOR
	
	PLOT,t,SFR20,THICK=1,XTITLE=xt,YTITLE="SFR, Mass flux [MSol/yr]",COLOR=0,BACKGROUND=255,/ylog,XSTYLE=1,YSTYLE=1,XRANGE=[MIN(t),MAX(t)],YRANGE=[.1,16]
	OPLOT,t,SFR12,THICK=1,COLOR=0,linestyle=2
	OPLOT,t,SFR8,THICK=1,COLOR=0,linestyle=3
	OPLOT,t,SFR4,THICK=1,COLOR=0
	OPLOT,t,SFR2,THICK=1,COLOR=0
	OPLOT,t,Mdot20,COLOR=1,THICK=1
	OPLOT,t,Mdot12,COLOR=1,THICK=1,linestyle=2
	OPLOT,t,Mdot8,COLOR=1,THICK=1,linestyle=3
	OPLOT,t,Mdot4,COLOR=1,THICK=1
	OPLOT,t,Mdot2,COLOR=1,THICK=1

	axis,XSTyle=1,XAXIS=1,XTICKFORMAT='conv_axis',CHARSIZE=cs,XTITLE="Redshift",xtickv=xtv,COLOR=0,xticks=4
	figureClean,(model.name+'_mdotsfr'),sv

	figureInit,(model.name+'_mdotsfr2'),sv,1,1
	PLOT,SFR20,Mdot20,XTITLE="SFR [MSol/yr]",YTITLE="Mdot [MSol/yr]",COLOR=0,BACKGROUND=255,XRANGE=[.2,15],YRANGE=[.2,15],/ylog,/xlog,XSTYLE=1,YSTYLE=1
	OPLOT,SFR12,Mdot12,COLOR=1
	OPLOT,SFR8,Mdot8,COLOR=2
	OPLOT,SFR4,Mdot4,COLOR=3
	OPLOT,SFR2,Mdot2,COLOR=4
	OPLOT,findgen(1000)*.016+.001,findgen(1000)*.016+.001,COLOR=0,linestyle=3
	figureClean,(model.name+'_mdotsfr2'),sv
END

PRO angularMomentum,model,sv
	figureInit,(model.name+'_angMomRad'),sv,1,2
	multi=1
	IF(multi EQ 0) THEN !p.multi=[0,1,3]
	IF(multi EQ 1) THEN BEGIN
		multiplot,[0,1,3,0,0]
		!p.noerase=0
	ENDIF
	xr=fltarr(2,3)
	xr[0,*]=0.0
	xr[1,*]=.99999*model.Radius
	yr=fltarr(2,3)
	yr[0,*]=.1
	yr[1,*]=20
	zj=intarr(4)
	zj[0]=nearest("redshift",1,model,0.0)
	zj[1]=nearest("redshift",1,model,1.0)
	zj[2]=nearest("redshift",1,model,1.5)
	zj[3]=nearest("redshift",1,model,2.0)
	lsj=[0,2,3,1]
	ncs=model.ncolstep
	yind=[ncs+28,ncs+29,ncs+30]
	xtitles=["","","Radius [kpc]"]
	ytitles=["Gas J per unit area","Stellar J per unit area","Cumulative Outflow J per unit area"]
	xind=[0,0,0]+ncs+9
	FOR i=0,2 DO BEGIN
                IF(i NE 0 AND multi EQ 1) THEN !p.noerase=1
		FOR j=0,n_elements(zj)-1 DO BEGIN
			PlotVsXind,model,zj[j],xind[i],yind[i],OP(j),0,xrange=xr[*,i],CHARSIZE=0.85,COLOR=j,XTITLE=xtitles[i],YTITLE=ytitles[i],XSTYLE=1,YSTYLE=1,LINESTYLE=lsj[j],YRANGE=yr[*,i]
		ENDFOR
		IF(multi EQ 1) THEN multiplot
	ENDFOR
	figureClean,(model.name+'_angMomRad'),sv
        IF(multi EQ 1) THEN BEGIN
                multiplot,/reset
                multiplot,/default
        ENDIF

	figureInit,(model.name+'_angMomTime'),sv,1,1
	gasJ=model.dataCube[*,*,ncs+28-1] ;;;; angular momenta per unit area, J/area = S x u
	starJ=model.dataCube[*,*,ncs+29-1]
	outJ=model.dataCube[*,*,ncs+30-1]
	gasS=model.dataCube[*,*,3]
	starS=model.dataCube[*,*,5]
	outS=model.dataCube[*,*,51]*model.mlf
	x=model.dataCube[*,*,0]

	extJ = model.evArray[18-1,*] * model.mdotext0 * 2*!pi*model.Radius*model.Radius ;; MSol/yr * kpc^2	

	macc=model.evArray[17-1,*] ;; cumulative mass accreted in solar masses
	speryear=31556926d
	kmperkpc=3.08568025d16
	totInJ=macc*model.Radius*model.vphiR*speryear/kmperkpc ;; solar masses * radius (kpc) * vphiR (km/s) * speryear / kmperkpc


	nstp=n_elements(model.dataCube[*,0,0])
	t=model.evArray[1,*]
	spGasJ=dblarr(nstp)
	spStarJ=dblarr(nstp)
	spOutJ=dblarr(nstp)
	
	totGasJ=dblarr(nstp)
	totStarJ=dblarr(nstp)
	totOutJ=dblarr(nstp)
	BulgeM=dblarr(nstp)
	BulgeJ=dblarr(nstp)
	FOR i=0,nstp-1 DO BEGIN ;; for every time step..
		spGasJ[i]=TOTAL(x[i,*]*x[i,*]*gasJ[i,*])/TOTAL(x[i,*]*x[i,*]*gasS[i,*])
		spStarJ[i]=TOTAL(x[i,*]*x[i,*]*starJ[i,*])/TOTAL(x[i,*]*x[i,*]*starS[i,*])
		spOutJ[i]=TOTAL(x[i,*]*x[i,*]*outJ[i,*])/TOTAL(x[i,*]*x[i,*]*outS[i,*])

		;; These values are dimensional: solar masses * kpc^2 / year
		totGasJ[i]=TOTAL(2.0*!pi*x[i,*]*x[i,*]*x[i,*]*model.dlnx*gasS[i,*]*model.mdotext0*model.radius*model.radius);; needs modification if non-flat rotation curve, i.e. a factor of u, which is currently not output by the gidget executable.
		totStarJ[i]=TOTAL(2.0*!pi*x[i,*]*x[i,*]*x[i,*]*model.dlnx*starS[i,*]*model.mdotext0*model.radius*model.radius)
		totOutJ[i]=TOTAL(2.0*!pi*x[i,*]*x[i,*]*x[i,*]*model.dlnx*outS[i,*]*model.mdotext0*model.radius*model.radius)

		;; assume that the mass in the bulge is equal to the mass which used to be in the disk but no longer is, taking star formation, and ejection by winds into account
		BulgeM[i]=macc[i] + TOTAL((2*!pi*model.Radius*model.mdotext0*model.dlnx*(kmperkpc/speryear)/model.vphiR) * x[i,*]*x[i,*]*(-outS[i,*]-starS[i,*]-gasS[i,*]+starS[0,*]+gasS[0,*])) ;; kpc * MSol/yr *s/km * (1/speryear) * (kmperkpc) = MSol
		BulgeJ[i]=BulgeM[i] * x[i,0]*model.Radius * model.vphiR *(speryear/kmperkpc) ;; solar masses * kpc * km/s * (speryear/kmperkpc) = solar masses kpc^2 / year
	ENDFOR
	consvd = totStarJ+totGasJ+totOutJ+BulgeJ - (extJ+totInJ)
	initJ = replicate(totStarJ[0]+totGasJ[0],n_elements(totStarJ))

	plot,t,spGasJ,COLOR=0,BACKGROUND=255,XSTYLE=1,THICK=2,XTITLE="time (Ga)",YTITLE="specific angular momentum"
	oplot,t,spStarJ,COLOR=1,linestyle=2
	oplot,t,spOutJ,COLOR=2,linestyle=3
        legend,["Gas","Stars","Outflows"],COLORS=[0,1,2],LINESTYLE=[0,2,3],THICK=[1,1,1],TEXTCOLORS=[0,0,0],pos=[.6,.35],/norm,charsize=.7
	figureClean,(model.name+'_angMomTime'),sv

	figureInit,(model.name+'_totAngMomTime'),sv,1.5,1.5
	plot,t,totGasJ,COLOR=0,BACKGROUND=255,XSTYLE=1,THICK=3,XTITLE="time (Ga)",YTITLE="Ang Mom [MSol kpc^2/yr]",ylog=1,charsize=2,YRANGE=[1d3,1d6]
	oplot,t,totStarJ,COLOR=1,linestyle=2,THICK=2
	oplot,t,totOutJ,COLOR=2,linestyle=3,THICK=2
	oplot,t,totInJ,COLOR=3,linestyle=3,THICK=2
	oplot,t,-extJ,COLOR=4,linestyle=3,THICK=2
	oplot,t,BulgeJ,COLOR=6,linestyle=6,THICK=2
	oplot,t,consvd ,COLOR=5,linestyle=4,THICK=2  ;; should be const.
	oplot,t,initJ,COLOR=7,linestyle=5,THICK=2
	legend,["Gas","Stars","Outflows","Cosmic Accretion","-extJ","Bulge","gas+stars+out+bulge+extTorque+CosAcc","initial"],COLORS=[0,1,2,3,4,6,5,7],LINESTYLE=[0,2,3,3,3,6,4,5],THICK=[2,1,1,1,1,1,1,1],TEXTCOLORS=[0,0,0,0,0,0,0,0],pos=[.5,.4],/norm,charsize=1.1
	figureClean,(model.name+'_totAngMomTime'),sv

	nt=n_elements(totStarJ)
        massErr=(initJ[nt-1]-consvd[nt-1])/(model.Radius*model.vphiR*speryear/kmperkpc)
	fracErr=massErr/TOTAL(2.0*!pi*model.Radius*model.dlnx*model.mdotext0/model.vphiR * (kmperkpc/speryear) * x[nt-1,*]*x[nt-1,*]*starS[nt-1,*])
	PRINT,"Percentage error in J: ", (consvd[nt-1]-initJ[nt-1])/(-initJ[nt-1]), " corresponding to a mass of ",massErr,"solar masses at the outer boundary. This is a fraction ",fracErr,"of the stellar mass in the disk."

END


PRO RadialStellarPops,model,sv
        figureInit,(model.name+'_pub_rst'),sv,1,2
        multi=1
        IF(multi EQ 0) THEN !p.multi=[0,1,3]
        IF(multi EQ 1) THEN BEGIN
                multiplot,[0,1,3,0,0]
                !p.noerase=0
        ENDIF
        xr=fltarr(2,3)
        yr=fltarr(2,3)
        xr[0,*]=0.0 ;; Radius (kpc)
        xr[1,*]=.99999*model.Radius ;; kpc
        yr[*,0]=[.1,1000] ; column density
        yr[*,1]=[7,120]   ; velocity dispersion
        yr[*,2]=[-1.5,1.0]; metallicity
        zj=nearest("redshift",1,model,0.0)
        ncs=model.ncolstep
        npp=model.npostprocess
        stv=model.stvars
        np=model.NPassive
        cs=.7
        yind=[indgen(np+1)*stv+1,indgen(np+1)*stv+2,indgen(np+2)*stv+3]+ncs+npp
        xtitles=["","","Radius [kpc]"]
        ytitles=["Column Density [MSol/pc^2]","Velocity Dispersion [km/s]","[Z]"]
        xind=[0,0,0]+ncs+9
        FOR i=0,2 DO BEGIN
                IF(i NE 0 AND multi EQ 1) THEN !p.noerase=1
                FOR j=0,model.NPassive DO BEGIN
                        PlotVsXInd,model,zj,xind[i],yind[i*(np+1)+j],OP(j),0,xrange=xr[*,i],yrange=yr[*,i],CHARSIZE=cs,COLOR=j,XTITLE=xtitles[i],YTITLE=ytitles[i],XSTYLE=1,YSTYLE=1,LINESTYLE=(1-OP(j))
                ENDFOR
                IF(multi EQ 1) THEN multiplot
        ENDFOR
        figureClean,(model.name+'_pub_rst'),sv
        IF(sv EQ 2) THEN latexify,(model.name+'_pub_rst.eps'),["Radius [kpc]",ytitles],['Radius [kpc]','Column Density [$M_\odot/pc^2$]','Velocity Dispersion [km/s]','[\ $<Z>$\ ]'],replicate(1.3,4),/full
        IF(multi EQ 1) THEN BEGIN
                multiplot,/reset
                multiplot,/default
        ENDIF

END


PRO RadialTimeSeriesGlobal,model,sv
	figureInit,(model.name+'_pub_rtsg1'),sv,1,2
        multi = 1
        IF(multi EQ 0) THEN !p.multi=[0,3,1]
        IF(multi EQ 1) THEN BEGIN
                multiplot,[0,1,3,0,0]
                !p.noerase=0
        ENDIF
        
        xr=fltarr(2,6)
        yr=fltarr(2,6)
        yr[*,0]=[0.0,1.0] ;; gas fraction
        yr[*,1]=[0.0,1.0] ;; H2 fraction
        yr[*,2]=[1.0d-5,5.0 ];; SFR
        yr[*,3]=[1.0d6,2.0d9] ;; Two-Dimensional Jeans Mass
        yr[*,4]=[1.8,4.2] ;; Q
        yr[*,5]=[.1,2000] ;; Cumulative Star Formation
        xr[0,*]=0.0 ;; Radius (kpc)
        xr[1,*]=0.9999*model.Radius ;;(kpc)
        zj=fltarr(4)
        zj[3]=nearest("redshift",1,model,2.0)
        zj[2]=nearest("redshift",1,model,1.5)
        zj[1]=nearest("redshift",1,model,1.0)
        zj[0]=nearest("redshift",1,model,0.0)
        lsj=intarr(4)
        lsj=[0,2,3,1]
        thj=[2,2,2,2]
        yind = [17,48,model.ncolstep+14,model.ncolstep+8,12,model.ncolstep+15]
        xind = yind*0 + model.ncolstep+9
        cs=.7
        xtitles=strarr(6)
        xtitles[[2,5]]="Radius (kpc)"        
	ytitles=["Gas Fraction","H2 Fraction","Star Formation Rate [MSol/yr/kpc^2]",'2D Jeans Mass [MSol]',"Q","Cumulative Star Formation [MSol/pc^2]"]

        FOR i=0,2 DO BEGIN
                IF(i NE 0 AND multi EQ 1) THEN !p.noerase=1
                FOR j=0,n_elements(lsj)-1 DO BEGIN
                        PlotVsXInd,model,zj[j],xind[i],yind[i],OP(j),0,xrange=xr[*,i],yrange=yr[*,i],CHARSIZE=cs,COLOR=j,Linestyle=lsj[j],XTITLE=xtitles[i],YTITLE=ytitles[i],XSTYLE=1,YSTYLE=1,THICK=thj[j]
                ENDFOR
                IF(multi EQ 1) THEN multiplot
        ENDFOR
        figureClean,(model.name+'_pub_rtsg1'),sv
        IF(sv EQ 2) THEN latexify,(model.name+'_pub_rtsg1.eps'),["Radius (kpc)",ytitles[0],ytitles[1],ytitles[2]],['Radius [kpc]','Gas Fraction','$H_2$ Fraction','Star Formation Rate [$M_\odot/yr/kpc^2$]'],replicate(1.3,4),/full

        IF(multi EQ 1) THEN BEGIN
                multiplot,/reset
                multiplot,/default
	ENDIF

	figureInit,(model.name+'_pub_rtsg2'),sv,1,2

	IF(multi EQ 1) THEN BEGIN
                multiplot,[0,1,3,0,0]
                !p.noerase=0
        ENDIF
        FOR i=3,5 DO BEGIN
                IF(i NE 3 AND multi EQ 1) THEN !p.noerase=1
                FOR j=0,n_elements(lsj)-1 DO BEGIN
                        PlotVsXInd,model,zj[j],xind[i],yind[i],OP(j),0,xrange=xr[*,i],yrange=yr[*,i],Charsize=cs,COLOR=j,LINESTYLE=lsj[j],XTITLE=xtitles[i],YTITLE=ytitles[i],XSTYLE=1,YSTYLE=1,THICK=thj[j]
                ENDFOR
                IF(multi EQ 1) THEN multiplot
        ENDFOR

        IF(multi EQ 1) THEN BEGIN
                multiplot,/reset
                multiplot,/default
        ENDIF
        figureClean,(model.name+'_pub_rtsg2'),sv
        IF(sv EQ 2) THEN latexify,(model.name+'_pub_rtsg2.eps'),["Radius (kpc)",ytitles[3],ytitles[4],ytitles[5]],['Radius [kpc]','2d Jeans Mass [$M_\odot$]','$Q$','Cumulative Star Formation [$M_\odot/pc^2$]'],replicate(1.3,4),/full
        IF(multi EQ 0) THEN !p.multi=[0,1,1]
END

PRO RadialTimeSeriesComponents,model,sv        
	figureInit,(model.name+'_pub_rtsc'),sv,1,2
        multiplot,[0,2,6,0,0]
        !p.noerase=0
        xr=fltarr(2,12)
        yr=fltarr(2,12)
        yr[0,0:1] = 0.
                yr[1,0:1] = 3   ;;; log column densities in MSol/pc^2
        yr[0,2:3] = 5                
		yr[1,2:3] = 250 ;; velocity dispersion (km/s)        
	yr[0,4:5]= -1.2
                yr[1,4:5]=1 ;; Metallicity
        yr[0,6:7] = 0
                yr[1,6:7]=5.9 ;; scale height (kpc)

        yr[0,8:9]=0.8
                yr[1,8:9]=30 ;; Q
        xr[0,*]=0
                xr[1,*]=.9999*model.Radius ;; radius (kpc)
        yr[0,10:11]=-3
                yr[1,10:11]=1.0 ;; log v_r (km/s)
        zj=fltarr(4)
	compile_opt strictarr
;	print,Nearest("redshift",1,model,2.0) ;; needed to auto-compile NRI()
        zj[3]=nearest("redshift",1,model,2.0)
        zj[2]=nearest("redshift",1,model,1.5)
        zj[1]=nearest("redshift",1,model,1.0)
        zj[0]=nearest("redshift",1,model,0.0)
        lsj=findgen(4)
        lsj=[0,2,3,1]
        thj=[1,1,1,2]
        colj=[0,1,2,3]
        ;;;; Column Density
        ;;;; Velocity Dispersion
        ;;;; Metallicity
        ;;;; Scale Height
        ;;;; Q_i
        yind = [model.ncolstep+24, model.ncolstep+25,$
                model.ncolstep+11, model.ncolstep+13,$
                22,model.ncolstep+18,$
                model.ncolstep+7,model.ncolstep+6,$
                24,23,$
                model.ncolstep+26,model.ncolstep+27]
        xind=yind*0+model.ncolstep+9
        cs=.7
        titles=strarr(12)
        titles[0:1]=["Gas","Stars"];,"Old Stars","Young Stars"]
        xtitles=strarr(12)
        xtitles[10:11]=Replicate("Radius (kpc)",2)
        ytitles=strarr(12)
        ytitles[0]="log Column Density [MSol/pc^2]"
        ytitles[2]="Velocity Dispersion [km/s]"
        ytitles[4]="[Z]"
        ytitles[6]="Scale Height [kpc]"
        ytitles[8]="Qi"
        ytitles[10]="log Inward Velocity [km/s]"
                FOR i=0,11 DO BEGIN                IF(i NE 0) THEN !p.noerase=1
                FOR j=0,n_elements(lsj)-1 DO BEGIN
                        PlotVsXInd,model,zj[j],xind[i],yind[i],OP(j),0,xrange=xr[*,i],yrange=yr[*,i],CHARSIZE=cs,COLOR=colj[j],Linestyle=lsj[j],XTITLE=xtitles[i],YTITLE=ytitles[i],TITLE=titles[i],XSTYLE=1,YSTYLE=1,THICK=thj[j]
                ENDFOR
                multiplot
        ENDFOR        
	multiplot,/reset        
	multiplot,/default
        figureClean,(model.name+'_pub_rtsc'),sv
        IF(sv EQ 2) THEN latexify,(model.name+'_pub_rtsc.eps'),[ytitles[0],ytitles[2],ytitles[4],ytitles[6],ytitles[8],ytitles[10]],['$\log \Sigma_j$ [$M_\odot/pc^2$]','$\sigma_j$ [km/s]','[\ $<Z>$\ ]','$H_j$ [kpc]','$Q_j$','$\log v_r$ [km/s]'],replicate(.83,6),/full

END



;; sv=0 : just show in a window
;; sv=1 : save as a png
;; sv=2 : save as an eps
PRO diskAnalyze,name,sv
        compile_opt idl2
        SETCT,1,0,0

        modelList=readmodels([name])
        model=(*(modelList[0]))
	modelList2=readmodels([ 'ri32/ri32babe',$
				'ri32/ri32baae',$
				'ri32/ri32bbbe',$
				'ri32/ri32bcbe',$
				'ri32/ri32bcae'])
;	modelList2=readmodels([ 'ri16/ri16'
        modelList2[0]=modelList[0]

        COMMON modelBlock,theModel
        theModel=(*(modelList[0]))

        RadialTimeSeriesComponents,model,sv
        RadialTimeSeriesGlobal,model,sv
	angularMomentum,model,sv
        quench,modelList,sv
	stellarDisk,model,sv
        setct,2,0,0
        RadialStellarPops,model,sv
        setct,1,0,0
;       globalMassBudget,model,sv
        PlotStarsVsAge,modelList2,8,"stpop_sol","",sv
        GlobalTimeSeries,model,sv
	MdotSFR,model,sv        

        ;; Compute gas column density at the solar radius:
        z0 = nearest("redshift",1,model,0.0)
        xs = nearest("position",1,model,8/model.Radius)
        colg=model.dataCube[*,*,model.ncolStep+9]
        sigg=model.dataCube[*,*,model.ncolStep+10]
        PRINT,"Gas Column Density at the solar radius: ", colg[z0,xs],z0,xs

        ;; Compute surface density ratio between inner edge of disk and 10 kpc (gas)
        x10 = nearest("position",1,model,10/model.Radius)
        PRINT,"Gas Column Density Contrast:", MAX(colg[z0,*],maxind)/colg[z0,x10]

        PRINT,"Gas sigma contrast (inner edge of star-forming disk/r=10 kpc): ", sigg[z0,x10]/sigg[z0,maxind]


        ;; Compute stellar mass within 10 kpc
        x=model.dataCube[0,*,0]
        dlnx=-alog(x[0])/(n_elements(x)-1)
        PRINT,"Total Stellar Mass within 10 kpc", TOTAL(2*!pi*x[0:x10-1]*x[0:x10-1]*dlnx*model.dataCube[z0,0:x10-1,model.ncolStep+11]*model.Radius*model.Radius*1d6)," MSol"

	z2=nearest("redshift",1,model,2.0)
	x2=nearest("position",1,model,2.0/model.Radius)
	x4=nearest("position",1,model,4.0/model.Radius)
	x8=nearest("position",1,model,8.0/model.Radius)
	x16=nearest("position",1,model,16.0/model.Radius)
	x20=nearest("position",1,model,20.0/model.Radius)
	xo=nearest("position",1,model,1)
	z1p9=nearest("redshift",1,model,1.9)
	

	PRINT,"At z=2:"
	PRINT,"Total SFR within 2kpc",TOTAL(2*!pi*x[0:x2-1]*x[0:x2-1]*dlnx*model.dataCube[z2,0:x2-1,model.ncolstep+14-1]*model.Radius*model.Radius)," MSol/yr"
	PRINT,"Total SFR outside 2kpc",TOTAL(2*!pi*x[x2-1:xo-1]*x[x2-1:xo-1]*dlnx*model.dataCube[z2,x2-1:xo-1,model.ncolstep+14-1]*model.Radius*model.Radius),"MSol/yr"
	PRINT,"mass flux across 2,4,8,16,20 kpc",-model.dataCube[z2,x2-1,3-1]*model.mdotext0,-model.dataCube[z2,x4-1,3-1]*model.mdotext0,-model.dataCube[z2,x8-1,3-1]*model.mdotext0,-model.dataCube[z2,x16-1,3-1]*model.mdotext0,-model.dataCube[z2,x20-1,3-1]*model.mdotext0," MSol/yr"

	PRINT,"At z=1.9:"
	PRINT,"Total SFR within 2kpc",TOTAL(2*!pi*x[0:x2-1]*x[0:x2-1]*dlnx*model.dataCube[z1p9,0:x2-1,model.ncolstep+14-1]*model.Radius*model.Radius)," MSol/yr"
	PRINT,"Total SFR outside 2kpc",TOTAL(2*!pi*x[x2-1:xo-1]*x[x2-1:xo-1]*dlnx*model.dataCube[z1p9,x2-1:xo-1,model.ncolstep+14-1]*model.Radius*model.Radius),"MSol/yr"
	PRINT,"mass flux across 2,4,8,16,20 kpc",-model.dataCube[z1p9,x2-1,3-1]*model.mdotext0,-model.dataCube[z1p9,x4-1,3-1]*model.mdotext0,-model.dataCube[z1p9,x8-1,3-1]*model.mdotext0,-model.dataCube[z1p9,x16-1,3-1]*model.mdotext0,-model.dataCube[z1p9,x20-1,3-1]*model.mdotext0," MSol/yr"
	
END
