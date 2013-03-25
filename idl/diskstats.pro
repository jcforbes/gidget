; return the radius at which half of the quantity (assumed to be in units of integral quantity/area) is contained
FUNCTION halfRadius,r,dlnx,quantity,inner=inner
    IF(n_elements(inner) EQ 0) THEN inner = 0.0

    ; simplest option:
    overall = TOTAL((2*!pi*2.0*sinh(dlnx/2.0) * r*r * quantity)) + inner
    FOR i=0, n_elements(r)-1 DO BEGIN
        IF(inner + TOTAL((2*!pi*2.0*sinh(dlnx/2.0)*r[0:i]*r[0:i]*quantity[0:i])) GT .5*overall) THEN RETURN,r[i]
    ENDFOR

    ; Never get here, unless quanitity or inner aren't positive-definite.
    RETURN,r[n_elements(r)-1]

END

FUNCTION generatelinearish,x,upperCutoff
    ; first pass:
    ret = lindgen(n_elements(x))
    near = MIN( abs(x-upperCutoff), cutoffindex)
    IF(cutoffIndex LT 10) THEN cutoffIndex = 10
    dynRange = x[cutoffIndex]/x[0]
    maxDeltaX = x[cutoffIndex] - x[cutoffIndex-1]
    linearRange = x[cutoffIndex] - x[0]
    nrescaled = FIX(linearRange / maxDeltaX )-1
    ret = lonarr(nrescaled)
    counter = 0
    FOR i=1, nrescaled DO BEGIN
        thisX = maxDeltaX*i + MIN(x)
        near = MIN( abs(x - thisX), index)
        IF(thisX LE upperCutoff) THEN BEGIN
            ret[i-1] = index
            counter = counter+1
        ENDIF
    ENDFOR
    ret = ret[0:counter-1]
    return,ret
END

;; Kind of an ad-hoc function to take a full simulation and reduce it to a set of simple numbers.
;; In essence, instead of having a bunch of quantities vs. radius, have a bunch of standalone quantities.
FUNCTION diskStats,model,z=z
    info=dblarr(7)
    kmperkpc=3.08568025d16
    speryear=31556926d
    cmperkpc = kmperkpc*1.0d5
    Msol=1.9891d33 ; grams 

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
    outr=where(fH2 GT .5,ct)
    IF(ct GT 0) THEN xout=max(x[outr]) ELSE xout=model.xmin
    tmp=max(col,xpeakInd)
    xpeak=x[xpeakInd]
    xinInd = nearest("position",1,model,xin)
    xoutInd= nearest("position",1,model,xout)
    vrg=model.dataCube[zj,*,ncs+17-1]
    vst=model.dataCube[zj,*,ncs+16-1]
    fg=col/(col+colst)
    sig = model.dataCube[zj,*,ncs+11-1] ;; velocity dispersion of gas (km/s)
    tdep = model.dataCube[zj,*,ncs+35-1] ;; depletion time := gas mass / SFR

    
    sfrHalfRadius = halfRadius(x*model.Radius,model.dlnx,model.dataCube[zj,*,ncs+14-1],inner = mdotbulgeG)
    gasHalfRadius = halfRadius(x*model.Radius,model.dlnx,col)
    stHalfRadius = halfRadius(x*model.Radius,model.dlnx,colst)

    optInd = nearest("position",1,model,stHalfRadius*2/model.Radius)

    fake = 1

    ;; fit the stellar & gas column density profiles to sersic+exp profile.
    ; for the gas, let's only fit inside a disk defined by where vrg becomes outwards
    vrglt0 = WHERE(vrg LT 0,gCtr)
    IF(gCtr GT 0) THEN fitGas = lindgen(MAX([MIN(vrglt0),10])) ELSE fitGas = lindgen(n_elements(x))
;    fitGas = lindgen(n_elements(x))
    fitGas = GenerateLinearish(x,2.0*gasHalfRadius/model.Radius)
    pg=mpfitprofile(x[fitGas]*model.Radius,col[fitGas],gchisq,single=1,fake=fake)
    gasScaleLength = pg[1] ; in kpc
    ;gchisq=2
    ;stchisq=2


    ; for the stars, only fit where Qst has fallen below Qlim
    fitStars= where(model.dataCube[zj,*,23-1] LE model.qlim, stCtr)
    IF(stCtr LT 10) THEN fitStars = lindgen(n_elements(x))
    colstfit = dblarr(n_elements(fitStars))
;    fitStars = lindgen(n_elements(x))
    fitStars = GenerateLinearish(x,2.0*stHalfRadius/model.Radius)
    pst=mpfitprofile(x[fitStars]*model.Radius,colst[fitStars],stchisq,colfit=colstfit,fake=fake)
    starScaleLengthMeas = pst[2]
    starSersic = pst[4]
    starBTMeas = TOTAL(x[fitStars]*(pst[1]*exp(-(x[fitStars]*model.Radius/pst[3])^(1.0/pst[4]))))/TOTAL(x[fitStars]*colstFit)

    indexGI = where(Q LE model.fixedQ,ct)
    IF(ct GT 0) THEN BEGIN
          GIin = min(x[indexGI])
          GIout= max(x[indexGI])
    ENDIF ELSE BEGIN
          GIin = model.xmin*model.xmin
          GIout= model.xmin*model.xmin
    END

    gasMass = model.dataCube[zj,*,3] * (2.0*!pi*model.Radius*model.mdotext0*(2.0*sinh(model.dlnx/2.0))*(kmperkpc/speryear)/model.vphiR) * x[*] * x[*]
    stellarMass =  model.dataCube[zj,*,5] * (2.0*!pi*model.Radius*model.mdotext0*(2.0*sinh(model.dlnx/2.0))*(kmperkpc/speryear)/model.vphiR) * x[*] * x[*]


    stMass = total(stellarMass) ;;;;; model.evArray[5,zj]/model.evArray[6,zj] * (1.0-model.evArray[6,zj])
    metallicity = model.dataCube[zj,*,22-1] ; = log10(  (mass in metals / mass in gas) / .02)
    stMetallicity = model.dataCube[zj,*,ncs+18-1] ; = log10( (mass in metals / mass in gas) / .02)

    gasS=model.dataCube[*,*,3] ;; dimensionless gas column density
    starS=model.dataCube[*,*,5] ;; dimensionless stellar column density
    outS=model.dataCube[*,*,51]*model.mlf ;; dimensionless outflow column density = stars formed * mass loading factor
    macc = model.evArray[17-1,*] ;; cumulative mass accreted in solar masses



    ;; Compute bulge mass by conservation arguments: cosmological input - total (outflowing mass + change in stellar mass + change in gas mass in the disk)
    BulgeM = macc[zj] + TOTAL((2*!pi*model.Radius*model.mdotext0*(2.0*sinh(model.dlnx/2.0))*(kmperkpc/speryear)/model.vphiR) * x[*]*x[*]*(-outS[zj,*]-starS[zj,*]-gasS[zj,*]+starS[0,*]+gasS[0,*])) ;; kpc* MSol/yr *s/km * (1/speryear) * (kmperkpc) = MSol

    BulgeMSt = TOTAL((2*!pi*model.Radius*model.mdotext0*(2.0*sinh(model.dlnx/2.0))*(kmperkpc/speryear)/model.vphiR) * x[*]*x[*]*(-starS[zj,*] + starS[0,*] + model.Rf * model.dataCube[zj,*,51]))

    ;; hopefully these two numbers are the same.
    sfr = model.evArray[11-1,zj] ;; SFR integrated over the whole galaxy in MSol/yr (I think)
    sfr2= TOTAL((2*!pi*model.Radius*model.Radius*model.dlnx * x[*]*x[*] * model.dataCube[zj,*,ncs+14-1]))


    ;;	mdotacc = - model.dataCube[zj,model.nx-1,3-1] * model.mdotext0 ;; infalling mass in MSol/yr

    mdotbulgeSt = -model.dataCube[zj,0,40-1] * model.mdotext0
    ;;;	mdotbulgeG = -model.dataCube[zj,0,3-1] * model.mdotext0

    mdotbulgeG = model.evArray[9-1,zj] * model.mdotext0  

    sSFR = model.evArray[11-1,zj]/(stMass+BulgeM) ;; yr^-1


    IF(xoutInd LT xinInd) THEN BEGIN 
              xoutInd = xinInd
              ;		PRINT,"diskstats: Warning: changing xoutInd to xinInd "
    ENDIF


    centralDensity = colst[0]
    dimensionlessMbulge = BulgeM * model.vphiR / (model.radius * model.mdotext0) * 1.0e5*speryear/(cmperkpc * Msol) ; = pi S_bulge,eff x_0^2
    sStJ = model.Radius*model.vphiR*(dimensionlessMbulge*model.dataCube[zj,0,16-1]*x[0]+TOTAL(2.0*!pi*(2.0*sinh(model.dlnx/2.0))*x*x*x*model.dataCube[zj,*,16-1]*model.dataCube[zj,*,6-1]))/(dimensionlessMbulge+TOTAL(2*!pi*(2.0*sinh(model.dlnx/2.0))*x*x*model.dataCube[zj,*,5]))
    sGasJ= model.Radius*model.vphiR*TOTAL(x*x*x*model.dataCube[zj,*,16-1]*model.dataCube[zj,*,4-1])/TOTAL(x*x*model.dataCube[zj,*,3])
    sOutJ= model.Radius*model.vphiR*TOTAL(x*x*x*model.dataCube[zj,*,16-1]*model.dataCube[zj,*,52-1])/TOTAL(x*x*model.dataCube[zj,*,51])



    gasZ = TOTAL( gasMass *  10.0^metallicity ) / TOTAL( gasMass ) ; metallicity = log10 ( M_Z/M_g / .02)
    totFg = TOTAL ( gasMass ) / (TOTAL (gasMass + stellarMass) + BulgeM)

    fgNuc = TOTAL( gasMass[0:xinInd] ) / TOTAL( gasMass[0:xinInd] + stellarMass[0:xinInd])
    fgsf = TOTAL (gasMass[xinInd:xoutInd] ) / TOTAL( gasMass[xinInd:xoutInd] + stellarMass[xinInd:xoutInd])
    fgHI = TOTAL (gasMass[xoutInd:model.nx-1] ) / TOTAL(gasMass[xoutInd:model.nx-1] + stellarMass[xoutInd:model.nx-1])

    totfH2 = TOTAL( gasMass * fH2 ) / TOTAL( gasMass )

    fgL = TOTAL( gasMass[0:optInd] ) / (TOTAL(gasMass[0:xoutInd] + stellarMass[0:xoutInd]) + BulgeM)
    ZL = TOTAL( gasMass[0:optInd] * 10.0^metallicity[0:xoutInd] ) / (TOTAL(gasMass[0:xoutInd]))

    vrAvg = TOTAL( gasMass * vrg ) / TOTAL(gasMass)
    sigAvg = TOTAL(sig*gasMass)/TOTAL(gasMass)

    indVrg = where(vrg GE 0.0)
    vrGE0 = TOTAL(gasMass[indVrg] * vrg[indVrg]) / TOTAL(gasMass[indVrg])
    den = model.evArray[11-1,zj]
    IF(den EQ 0) THEN den = model.evArray[11-1,1]
    tdepAvg = TOTAL( gasMass )*totfH2/den ; = molecular gas mass / SFR

    ; other quantities to compute: variance in Z_gas, Z_*, variance in Z_*, stellar age, mol gas frac, tdep for all gas

    bulgeMetallicity = alog10( model.evArray[5-1,zj] / .02)
    stZdelta = MAX([MAX(stmetallicity),alog10(model.evArray[5-1,zj]/.02)])-MIN([MIN(stmetallicity),alog10(model.evArray[5-1,zj]/.02)])
    molfg = TOTAL( gasMass* fH2 )/(TOTAL(gasMass + stellarMass)+BulgeM)
    tdepAllGas = tdepAvg / totfH2
    stZ = .02*((TOTAL(stellarMass * 10.0^stMetallicity) + BulgeM*10.0^bulgeMetallicity)/(stMass + BulgeM))
    gasZdelta = MAX(metallicity)-MIN(metallicity)
    stAge = TOTAL(model.evArray[11-1,0:zj]*(model.evArray[1,zj] - model.evArray[1,0:zj]))/TOTAL(model.evArray[11-1,0:zj]) ;; age(z) = Cumulative (SFR(z') * (lbt(z') - lbt(z))) / Cumulative(SFR) ---- on second thought, this is only an OK approximation because it assumes SFR is const. between outputs. Even worse, does not account for initial stars!
    deltaCol = alog10(col[xpeakind]/col[0])
    ;; Rather than comparing col @ peak to col@0, compare col@peak to col averaged over 10 innermost cells.
    ;; As usual, the two factors of x are the result of using a logarithmic grid.
    deltaColAvg = alog10(col[xpeakind] / (TOTAL(col[0:9] * x[0:9] * x[0:9])/TOTAL(x[0:9]*x[0:9])) )

    xnuc=x[0:xinInd]
    xsf=x[xinInd:xoutInd]
    xHI=x[xoutInd:model.nx-1]
    ;        fgnuc = TOTAL(xnuc*xnuc*fg[0:xinInd])/TOTAL(xnuc*xnuc)
    ;        fgsf  = TOTAL(xsf*xsf*fg[xinInd:xoutInd])/TOTAL(xsf*xsf)
    ;        fgHI  = TOTAL(xHI*xHI*fg[xoutInd:model.nx-1])/TOTAL(xHI*xHI)
    vrgNuc= TOTAL(xnuc*xnuc*vrg[0:xinInd])/TOTAL(xnuc*xnuc)
    vrgSf = TOTAL(xsf*xsf*vrg[xinInd:xoutInd])/TOTAL(xsf*xsf)
    vrgHI = TOTAL(xHI*xHI*vrg[xoutInd:model.nx-1])/TOTAL(xHI*xHI)
    vstNuc= TOTAL(xnuc*xnuc*vst[0:xinInd])/TOTAL(xnuc*xnuc)
    vstSf = TOTAL(xsf*xsf*vst[xinInd:xoutInd])/TOTAL(xsf*xsf)
    vstHI = TOTAL(xHI*xHI*vst[xoutInd:model.nx-1])/TOTAL(xHI*xHI)
    colsol= col[nearest("position",1,model,8.0/model.Radius)]
    ;	totfH2 = TOTAL(x[*] * x[*] * fH2[*])/TOTAL(x[*]*x[*])
    ;        totFg = TOTAL(x[*] * x[*] * fg[*]) / TOTAL(x[*] * x[*])
    ;        totZ = TOTAL(x[*] * x[*] * (10.0^metallicity)  )/TOTAL(x[*]*x[*])
    eff = (stMass+bulgeM)/(model.evArray[19-1,zj]*.18)

    info=[xin*model.Radius,$    ;; inner edge of SF region  - 1
        xpeak*model.Radius,$ ;; peak of the column density - 2
        xout*model.Radius,$ ;; outer edge of SF region- 3
        colsol,$            ;; column density at r=8 kpc- 4
        fgnuc,fgsf,fgHI,sSFR, $   ;; gas fraction, specific SFR- 5,6,7,8
        BulgeM,BulgeMSt, totfH2, $  ; - 9,10,11
        mdotBulgeG,mdotbulgeSt, $ ; 12,13
        sfr,stMass,totFg,gasZ,eff, $ ; - 14,15,16,17,18
        sfr+mdotBulgeG,stMass+BulgeM, $ ;; 19,20
        fgL, ZL, vrAvg, GIin*model.Radius, $ ; 21, 22, 23, 24,
        GIout*model.Radius, sigAvg, vrGE0,$   ; 25, 26, 27
        tdepAvg, gasScaleLength, $ ; 28,29
        starScaleLengthMeas,starSersic,starBTMeas, $ ; 30,31,32
        gChiSq,stChiSq, sfrHalfRadius, $ ; 33,34,35
        gasHalfRadius, stHalfRadius, $ ; 36,37
        centralDensity,sStJ,sGasJ,sOutJ, $ ; 38,39,40,41
        gasZdelta,stZdelta,molfg,tdepAllGas, $ ; 42,43,44,45
        stZ,stAge,deltaCol,deltaColAvg, $ ; 46,47,48,49
        BulgeM/(stMass+BulgeM), gasZ*.02, $ ; 50,51
        ZL*.02, model.x25s[zj]*model.Radius, $ ; 52, 53
        model.x25s_2[zj]*model.Radius, $ ; 54
        model.x25s_3[zj]*model.Radius, $ ; 55
        model.colTranses[zj], $ ;56
        (stMass+BulgeM)/(0.18 * eff), $ ;57
        model.x25s_4[zj]*model.Radius, $ ; 58
        model.sigmath*sqrt(1.0+(model.NN[zj]*model.NN[zj]))^(1.0/3.0), $ ;59
        max(model.dataCube[zj,*,model.ncolstep+11-1]), $ ;60
        0] ;61
    info[61-1]= info[60-1]/info[59-1]
    ;               vrgNuc,vrgSf,vrgHI,$;; radial gas velocity [km/s]
        ;               vstNuc,vstSf,vstHI] ;; radial stellar velocity [km/s]

    ;;; IF(BulgeM LT 0.0) THEN STOP,"Negative BulgeM: ",model.name," ",BulgeM


    wh = where(info NE info, count)
;    IF(count NE 0) THEN STOP, "NaN found in modelInfo!"

    return,info
END
