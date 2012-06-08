
;; Kind of an ad-hoc function to take a full simulation and reduce it to a set of simple numbers.
;; In essence, instead of having a bunch of quantities vs. radius, have a bunch of standalone quantities.
FUNCTION diskStats,model,z=z
        info=dblarr(7)
	IF(n_elements(z) EQ 0) THEN z=0.0
        ct=1
        ;; all quantities computed at z=0
        zj=nearest("redshift",1,model,z)
        ncs=model.ncolstep
        ;; Find x_in and x_out between which f_{H_2}>.03
        x=model.dataCube[zj,*,0]
        fH2= model.dataCube[zj,*,48-1]
        Q=model.dataCube[zj,*,12-1]
        col= model.dataCube[zj,*,ncs+10-1]
        colst=model.dataCube[zj,*,ncs+12-1]
;       xin=min(x[where(fH2 GT .1)])
        inr=where(Q GT model.fixedQ*2.5,ct)
        IF(ct GT 0) THEN xin=max(x[inr]) ELSE xin=0.0
	outr=where(fH2 GT .1,ct)
        IF(ct GT 0) THEN xout=max(x[outr]) ELSE xout=0.0
        tmp=max(col,xpeakInd)
        xpeak=x[xpeakInd]
        xinInd = nearest("position",1,model,xin)
        xoutInd= nearest("position",1,model,xout)
        vrg=model.dataCube[zj,*,ncs+17-1]
        vst=model.dataCube[zj,*,ncs+16-1]
        fg=col/(col+colst)
	stMass = model.evArray[5,zj]/model.evArray[6,zj] * (1.0-model.evArray[6,zj])
	sSFR = model.evArray[11-1,zj]/stMass
        metallicity = model.dataCube[zj,*,22-1]


	gasS=model.dataCube[*,*,3] ;; dimensionless gas column density
	starS=model.dataCube[*,*,5] ;; dimensionless stellar column density
	outS=model.dataCube[*,*,51]*model.mlf ;; dimensionless outflow column density = stars formed * mass loading factor
	macc = model.evArray[17-1,*] ;; cumulative mass accreted in solar masses

	kmperkpc=3.08568025d16
	speryear=31556926d


	;; Compute bulge mass by conservation arguments: cosmological input - total (outflowing mass + change in stellar mass + change in gas mass in the disk)
	BulgeM = macc[zj] + TOTAL((2*!pi*model.Radius*model.mdotext0*model.dlnx*(kmperkpc/speryear)/model.vphiR) * x[*]*x[*]*(-outS[zj,*]-starS[zj,*]-gasS[zj,*]+starS[0,*]+gasS[0,*])) ;; kpc* MSol/yr *s/km * (1/speryear) * (kmperkpc) = MSol

	BulgeMSt = TOTAL((2*!pi*model.Radius*model.mdotext0*model.dlnx*(kmperkpc/speryear)/model.vphiR) * x[*]*x[*]*(-starS[zj,*] + starS[0,*] + model.Rf * model.dataCube[zj,*,51]))

	;; hopefully these two numbers are the same.
	sfr = model.evArray[11-1,zj] ;; SFR integrated over the whole galaxy in MSol/yr (I think)
	sfr2= TOTAL((2*!pi*model.Radius*model.Radius*model.dlnx * x[*]*x[*] * model.dataCube[zj,*,ncs+14-1]))

	mdotacc = - model.dataCube[zj,model.nx-1,3-1] * model.mdotext0 ;; infalling mass in MSol/yr
	
	mdotbulgeSt = -model.dataCube[zj,0,40-1] * model.mdotext0
	mdotbulgeG = -model.dataCube[zj,0,3-1] * model.mdotext0

;	stop

	IF(xoutInd LT xinInd) THEN BEGIN 
		xoutInd = xinInd
;		PRINT,"diskstats: Warning: changing xoutInd to xinInd "
	ENDIF
 


        xnuc=x[0:xinInd]
        xsf=x[xinInd:xoutInd]
        xHI=x[xoutInd:model.nx-1]
        fgnuc = TOTAL(xnuc*xnuc*fg[0:xinInd])/TOTAL(xnuc*xnuc)
        fgsf  = TOTAL(xsf*xsf*fg[xinInd:xoutInd])/TOTAL(xsf*xsf)
        fgHI  = TOTAL(xHI*xHI*fg[xoutInd:model.nx-1])/TOTAL(xHI*xHI)
        vrgNuc= TOTAL(xnuc*xnuc*vrg[0:xinInd])/TOTAL(xnuc*xnuc)
        vrgSf = TOTAL(xsf*xsf*vrg[xinInd:xoutInd])/TOTAL(xsf*xsf)
        vrgHI = TOTAL(xHI*xHI*vrg[xoutInd:model.nx-1])/TOTAL(xHI*xHI)
        vstNuc= TOTAL(xnuc*xnuc*vst[0:xinInd])/TOTAL(xnuc*xnuc)
        vstSf = TOTAL(xsf*xsf*vst[xinInd:xoutInd])/TOTAL(xsf*xsf)
        vstHI = TOTAL(xHI*xHI*vst[xoutInd:model.nx-1])/TOTAL(xHI*xHI)   
        colsol= col[nearest("position",1,model,8.0/model.Radius)]
	totfH2 = TOTAL(x[*] * x[*] * fH2[*])/TOTAL(x[*]*x[*])
        totFg = TOTAL(x[*] * x[*] * fg[*]) / TOTAL(x[*] * x[*])
        totZ = TOTAL(x[*] * x[*] * (10.0^metallicity)  )/TOTAL(x[*]*x[*])
        eff = stMass/(model.mh0*.18)

        info=[xin*model.Radius,$    ;; inner edge of SF region  - 1
                xpeak*model.Radius,$ ;; peak of the column density - 2
                xout*model.Radius,$ ;; outer edge of SF region- 3
                colsol,$            ;; column density at r=8 kpc- 4
                fgnuc,fgsf,fgHI,sSFR, $   ;; gas fraction, specific SFR- 5,6,7,8
		BulgeM,BulgeMSt, totfH2, $  ; - 9,10,11
		mdotBulgeG,mdotbulgeSt, $ ; 12,13
		sfr,stMass,totFg,totZ,eff]; - 14,15,16,17,18
;               vrgNuc,vrgSf,vrgHI,$;; radial gas velocity [km/s]
;               vstNuc,vstSf,vstHI] ;; radial stellar velocity [km/s]

        return,info
END
