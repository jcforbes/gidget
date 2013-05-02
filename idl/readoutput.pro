compile_opt strictarr
;; This file contains the functions necessary to read in the data from one run
;; It also contains a number of auxiliary functions which will be useful later.


;; taken from the internet; counts number of lines in a file
;; from http://physics.nyu.edu/grierlab/idl_html_help/F12.html#wp894303
FUNCTION file_lines, filename 
   genOPEN, unit, filename, /GET_LUN 
   str = '' 
   count = 0ll 
   WHILE ~ EOF(unit) DO BEGIN 
      READF, unit, str 
      count = count + 1 
   ENDWHILE 
   FREE_LUN, unit 
   RETURN, count 
END

;; Super-simple function for later convenience.
FUNCTION OP,modelIndex
	IF (modelIndex EQ 0) THEN re=0 ELSE re=1
	RETURN,re
END

;; labels for each column in the dataCube.
;; Not guaranteed to be updated- readOutput post-processes the output from 
;; the C++ program, so the correctness of this function depends on the C++ function
;; WriteOutput() and the IDL function (below) readOutput().
FUNCTION GetLabel,ind,ncolstep,npostprocess,npassive,stvars ;; ind to be indexed by 1
  Labels=["x","tau","tau'","Gas Column Density","Gas Velocity Dispersion",$  ;; 1..5
    "Stellar Column Density","Stellar Velocity Dispersion",$           ;; 6,7
    "d/dt Gas Column Density","d/dT Gas Velocity Dispersion",$         ;; 8,9
    "d/dT Stellar Column Density","d/dT Stellar Velocity Dispersion",$ ;; 10,11
    "Q","h0","h1","h2","H","Gas Fraction","q","Toomre Length",$        ;; 12..19
    "Toomre Mass","dZ/dT","Z","Stellar Q","Gas Q","Rafikov Q",$	   ;; 20..25
    "Wang-Silk Q","Romeo-Wiegert Q","Q(q)","SF Column Density",$       ;; 26..29
    "tau''","dQdS","dQds","dQdS Error","dQds Error","yy",$		   ;; 30..35
    "Error in Torque Eq.","vr gas","Cumulative Stars Out (MSol)",$     ;; 36..38
    "Cumulative Gas Out (MSol)","Cumulative SF (MSol)",$               ;; 39,40
    "Change in Stellar Mass (MSol)","Change in Gas Mass (MSol)",$      ;; 41,42
    "d/dx tau'","d/dx Gas Velocity Dispersion",$                               ;; 43,44
    "Cumulative Stars In (MSol)","Cumulative Gas In (MSol)",$          ;; 45,46 
    "alpha viscosity","fH2",$ ;; 47..48
    "Cumulative Torque Err","Cumulative Torque Err 2","d2taudx2",$     ;; 49..51
     "Cumulative SF"]						   ;; 52 
  tLabels=strarr(ncolstep+npostprocess+stvars*(npassive+1))
  tLabels[0:(ncolstep-1)] = Labels[0:(ncolstep-1)]
  tLabels[ncolstep:(ncolstep+npostprocess-1)]= [$					   ;; ncolev+...
    "Col Timescale","Sig Timescale","Col_st Timescale", $              ;; 1..3
    "Sig_st Timescale","Gas Vel. Disp. / Stellar Vel. Disp.",$         ;; 4..5
    "Stellar Scale Height (kpc)","Gas Scale Height (kpc)",$            ;; 6..7
    "Two-Dim Jeans Mass (MSol)","Radius (kpc)",$			   ;; 8..9
    "Gas Col Density [MSol/pc^2]","Gas Vel Disp. [km/s]",$ 		   ;; 10..11
    "Stellar Col Density [MSol/pc^2]","Stellar Vel Disp. [km/s]",$     ;; 12..13
    "SFR Col Density [MSol/yr/kpc^2]","Cumulative SF [MSol/pc^2]",$    ;; 14..15
    "Inward Stellar Velocity [km/s]","Inward Gas Velocity [km/s]",$    ;; 16..17
    "Z_*","mdot","Gas(In-Out-(Rf+mu)SF)","St(In-Out+Rf*SF)",$          ;; 18..21
    "Err in Mass Cons (Gas)","Err in Mass Cons (Stars)",$              ;; 22..23
    "log col","log colst","log vrg","log vrst","gas J/area",$	   ;; 24..27
    "star J/area","cumulative outflow J/area",$			   ;; 28..30
    "average age at z=0","average present age",$		   ;; 31..32
    "-tau","-tau'","depletion time","viscous time","tDep/tVisc", $ ;; 33..37
    "tvisc^-1 (yr^-1)","Ratio w univ prof","Ratio2 w univ prof", $ ;; 38..40
    "Ratio3 w univ prof","","","","","","","","","" ] ;; 41-50

  FOR i=0, npassive DO BEGIN
    pop=strcompress(STRING(i),/remove)
    tLabels[ncolstep+npostprocess+i*stvars : ncolstep+npostprocess+(i+1)*stvars -1]$
    =["Column Density "+pop,"Velocity Dispersion R "+pop,$
    "Velocity Dispersion Z "+pop,$
    "Metallicity "+pop,"Metallicity disp "+pop,"Scale Height "+pop,$
    "Q_i "+pop,"[Z-sig_Z] "+pop,"[Z+sig_Z] "+pop]
  ENDFOR

  RETURN,tLabels[ind-1]
END


;; Is this index of the dataCube best plotted with
;; the y-axis logarithmic?
FUNCTION logvar, ind, model
  st0=model.ncolstep + model.NPostProcess
  ncs=model.ncolstep
  nsv=model.STVars
  log=0
  IF ( (ind GE 4 && ind LE 7) || ind EQ 16 || ind EQ 19 || ind EQ 23 || ind EQ 20 || ind EQ 24 $
        || ind EQ 38 || ind EQ 39 || ind EQ 40 || ind EQ 45 || ind EQ 46 || ind EQ 52) $
     THEN log=1
  IF (ind EQ ncs+8 || ind EQ ncs+22 || ind EQ ncs+23 || (ind GE ncs+10 AND ind LE ncs+17) $
      || ind EQ 19 || (ind GE ncs+1 AND ind LE ncs+5) || (ind GE ncs+28 AND ind LE ncs+30) $
      || ind EQ ncs+33 || ind EQ ncs+34 || ind EQ ncs+6 || ind EQ ncs+7) THEN log=1
  FOR i=0,model.NPassive DO BEGIN
    IF(ind EQ i*nsv+st0+1 || ind EQ i*nsv+st0+2 || ind EQ i*nsv+st0+4 ) THEN log=1
  ENDFOR
                
  RETURN,log
END

;; Given an already-open file's unit number, read a line from that file
;; and discard everything before the colon. Return a long int constructed
;; from the number following the colon.
FUNCTION ExtractCmtL,lun
  dummy=""
  paramStr=""
  READF,lun,dummy
  paramStr = strsplit(dummy,':',/extract)
  param = long(paramStr[1])
  RETURN,param
END
;; Same as above, but for floats
FUNCTION ExtractCmtFlt,lun
  dummy=""
  paramStr=""
  READF,lun,dummy
  paramStr = strsplit(dummy,':',/extract)
  param = float(paramStr[1])
  RETURN,param
END

;; Read in the output from one run, and do some post-processing.
FUNCTION readOutput,name
  ncolconv=11 ;; number of columns in _convergence.dat file
    speryear=31556926.0
    kmperkpc=3.08568025d16
    pcperkpc=1d3



  ;;; Read the comment file for this run. This contains all
  ;;; of the parameters that the user has the ability to provide.
  comment=""
  dummy=""
  paramStr=""
  genOPEN,lunCom,(name+"_comment.txt"),/get_lun
  READF,lunCom,dummy ;; first line is just dashes - ignore
  READF,lunCom,comment ;; next line is the comment
  nx=ExtractCmtL(lunCom)
  eta=ExtractCmtFlt(lunCom)
  EPS_ff=ExtractCmtFlt(lunCom)
  tauHeat=ExtractCmtFlt(lunCom)
  analyticQ=ExtractCmtL(lunCom)
  cosmologyOn=ExtractCmtL(lunCom)
  xmin=ExtractCmtFlt(lunCom)
  NActive=ExtractCmtL(lunCom)
  NPassive=ExtractCmtL(lunCom)
  vphiR=ExtractCmtFlt(lunCom)
  Radius=ExtractCmtFlt(lunCom)
  GasTemp=ExtractCmtFlt(lunCom)
  Qlim=ExtractCmtFlt(lunCom)
  fg0=ExtractCmtFlt(lunCom)
  tempRatio=ExtractCmtFlt(lunCom)
  zstart=ExtractCmtFlt(lunCom)
  tmax=ExtractCmtFlt(lunCom)
  nstepmax=ExtractCmtL(lunCom)
  TOL=ExtractCmtFlt(lunCom)
  MLF=ExtractCmtFlt(lunCom)
  b=ExtractCmtFlt(lunCom)
  beta0=ExtractCmtFlt(lunCom)
  nRotCurve=ExtractCmtFlt(lunCom)
;  md0=ExtractCmtFlt(lunCom)
  scaleLength=ExtractCmtFlt(lunCom)
  whichAccHistory = ExtractCmtL(lunCom)
  alphaMRI=ExtractCmtFlt(lunCom)
  thick=ExtractCmtFlt(lunCom)
  migPassivePop=ExtractCmtL(lunCom)
  fixedQ=ExtractCmtFlt(lunCom)
  kappaMetals=ExtractCmtFlt(lunCom)
  mh0=ExtractCmtFlt(lunCom)
  minStsig = ExtractCmtFlt(lunCom)
  ndecay = ExtractCmtFlt(lunCom)
  debugParam = ExtractCmtL(lunCom)
  accScaleLength = ExtractCmtFlt(lunCom)
  zquench = ExtractCmtFlt(lunCom)
  zrelax = ExtractCmtFlt(lunCom)
  zetaREC = ExtractCmtFlt(lunCom)
  Rf = ExtractCmtFlt(lunCom)
  deltaOmega = ExtractCmtFlt(lunCom)
  NoutputsNominal = ExtractCmtL(lunCom)
  accNormalization=ExtractCmtFlt(lunCom)
  accAlphaZ = ExtractCmtFlt(lunCom)
  accAlphaMh = ExtractCmtFlt(lunCom)
  accCeiling = ExtractCmtFlt(lunCom)
  fscatter = ExtractCmtFlt(lunCom)
  invMassRatio = ExtractCmtFlt(lunCom)
  fcool = ExtractCmtFlt(lunCom)
  whichAccretionProfile = ExtractCmtFlt(lunCom)
  alphaAccProf = ExtractCmtFlt(lunCom)
  accWidth = ExtractCmtFlt(lunCom)
  fH2Min = ExtractCmtFlt(lunCom)
  tDepSC = ExtractCmtFlt(lunCom)
  ZIGM = ExtractCmtFlt(lunCom)

  ;; Some useful constants
  MSol = 1.98892d33 ;; solar masses in grams
  G=6.673d-8 ;; newton's constant in cgs
  speryear = 31556926d
  cmperkpc = 3.08568025d21	
  cmperpc =  3.08568025d18
  mH = 1.6733d-24 ;; grams
  kB = 1.3806503d-16 ;; erg/K
  cm2perpc2perMsol = cmperpc*cmperpc / MSol
  cmperkpcpers2peryear2 = cmperkpc/(speryear*speryear)

  ;; these next two are special: The user doesn't specify them; they are computed
  ;; at the beginning of a run based on whichAccretionHistory.
  genOPEN,lunComAux,(name+"_aux.txt"),/get_lun
  READF,lunComAux,dummy
  READF,lunComAux,dummy
  md0 = ExtractCmtFlt(lunComAux)
  ND08attempts = ExtractCmtL(lunComAux)
  FREE_LUN,luncomaux 

  ;; Now we need to rescale some values for mh0!=1e12
 ; MLF = MLF * (Mh0/1.0e12) ^ (-1.0/3.0);
 ; vphiR = vphiR * (Mh0/1.0e12)^  (1.0/3.0)
	
  sigmath = sqrt(kB * GasTemp/mH) * 1.0e-5 ;; km/s
  b=b/Radius ;; b is read in in kpc (as is Radius). Now set b to its dimensionless value (incidentally, b is the turnover radius of the rotation curve)
  FREE_LUN,luncom

  ;; The maximum number of outputs from the program we'll read in.
  maxrec = 200000L;
  ;; Open up the file which directly gives time-dependent information
  ;; only, and no radial information. I.e. one line per timestep,
  ;; as opposed to nx lines per time step.
  genOPEN,lunEv,(name+"_evolution.dat"),/get_lun;,/swap_if_big_endian
  ;; Read in all the data from this file, keeping track of the number
  ;; of lines it has, since the larger arrays of data will correspond
  ;; to the same number of outputs.
  ncolev = 0L
  READU,lunEv,ncolev
  evArray = dblarr(ncolev,maxrec)
  evArrayStep = dblarr(ncolev)
  currentStep = 0L
  WHILE(EOF(lunEv) NE 1) DO BEGIN
    READU,lunEv,evArrayStep
    evArray[*,currentStep] = evArrayStep
    currentStep += 1
  ENDWHILE
  FREE_LUN,lunEv

  evArray=(temporary(evArray))[*,0:(currentStep-1)] ;; reduce the array to its natural size

  ;; This is the number of variables output for each stellar population		
  ;; at each timestep at each radius.
  STVars=5
  ;; starsHyperCube stores the passive stars..
  starsHyperCube = dblarr(currentStep,NPassive+1,nx,STVars) 
      ;; (reduced nsteps) x (NAgeBins+1) x (nx) x (4) 
  ;; starsHyperCubeA stores the active stars. The usual scenario
  ;; is that NActive will be 1.
  starsHyperCubeA= dblarr(currentStep,NActive+1,nx,STVars)
	;; (reduced nsteps) x (NActive+1) x (nx) x (4)
  stellarAges = dblarr(currentStep,NPassive+1) ;; these are the stellar ages at z=0
	;; (reduced nsteps) x (NAgeBins+1)
  stellarAgesIS = dblarr(currentStep,NPassive+1) ;; these are the stellar ages In Situ
	;; (reduced nsteps) x (NAgeBins+1)
  dataCube=dblarr(currentStep,nx,1) 
	;; this should be the exact size necessary except for the 3rd dimension


  mdotext0=md0* MSol/speryear ;; in grams/s
  chi = G*mdotext0/(vphiR*vphiR*vphiR*1d15)
  timeOutputs=0

  ;; Read in all the data which depends on radius.
  fname=name+"_radial.dat"
  genOPEN,lunRadial,fname,/GET_LUN
  genOPEN,lunStars,(name+"_stars.dat"),/GET_LUN
  genOPEN,lunStarsA,(name+"_act_stars.dat"),/GET_LUN
  ;; There will be a dataset at every timestep
  FOR stp=0, currentStep-1 DO BEGIN
    step=evArray[0,stp]
    t  = evArray[1,stp]
    dt = evArray[2,stp]

    ;; The files will tell us how many rows and columns
    ;; to read in for each timestep		
    ncolStep = 0L
    nrowStep = 0L
    READU,lunRadial,ncolStep
    READU,lunRadial,nrowStep
		
    ;; Create the corresponding array..
    dataSquare=dblarr(ncolStep,nrowStep)
    IF(stp EQ 0) THEN dataCube=dblarr(currentStep,nx,ncolStep)
    READU,lunRadial,dataSquare

    ;; Read in that 2d array to its appropriate place
    ;; in the 3d array.
    dataCube[stp,0:nrowStep-1,0:ncolStep-1] = TRANSPOSE(temporary(dataSquare))
    timeOut = (1 EQ 1)
    ;; read in the stellar populations at this time.
    IF timeOut THEN BEGIN
;;;;		OPENR,lunSt,(name+"_stars_to20-"+strcompress(string(floor(20.*(t+dt)))$
;;;;			,/remove)+".dat"),/GET_LUN
;;;;		OPENR,lunStA,(name+"_act_stars_to20-"+strcompress(string(floor(20.*(t+dt)))$
;;;;			,/remove)+".dat"),/GET_LUN
      NABp1=0L
      NABp1A=0L
      stsSize=0L
      stsSizeA=0L
      startAge=double(0.)
      endAge=double(0.)
      age=double(0.)
      nx2=0L
      nx2A=0L
      x=dblarr(nx)
      contOfAgeBin=dblarr(nx,STVars)
      contOfAgeBinA = dblarr(nx,STVars)
      READU,lunStars, NABp1,stsSize,nx2
      READU,lunStars,x
      READU,lunStarsA,NABp1A,stsSizeA,nx2A
      READU,lunStarsA,x
      FOR i=0, stsSize-1 DO BEGIN
;;;	IF(EOF(lunSt) NE 1) THEN BEGIN
		READU,lunStars,age,startAge,endAge
		stellarAges[timeOutputs,i]=age
		READU,lunStars,contOfAgeBin
		starsHyperCube[timeOutputs,i,0:nx-1,0:(STVars-1)]=contOfAgeBin
;;;	ENDIF
      ENDFOR
      FOR ii=0, stsSizeA-1 DO BEGIN
;;;	IF(EOF(lunStA) NE 1) THEN BEGIN
		READU,lunStarsA,age,startAge,endAge
		READU,lunStarsA,contOfAgeBinA
		starsHyperCubeA[timeOutputs,ii,0:nx-1,0:(STVars-1)]=contOfAgeBinA
;;;	ENDIF
      ENDFOR
					
;;    FREE_LUN,lunSt
;;    FREE_LUN,lunStA
			timeOutputs+=1
    ENDIF

    ENDFOR ;; end loop over steps in the evolution file
    FREE_LUN,lunRadial
    FREE_LUN,lunStars
    FREE_LUN,lunStarsA


    starsHyperCube=(temporary(starsHyperCube))[0:(timeOutputs-1),0:(NABp1-1),0:(nx-1),0:STVars-1] ;; trim down the hypercube
    starsHyperCubeA=(temporary(starsHyperCubeA))[0:(timeOutputs-1),0:(NABp1A-1),0:(nx-1),0:STVars-1] ;; trim down the active hypercube
    

    stvpp=4 ;; number of stellar properties calculated in this routine from the raw output
    STVars=STVars+stvpp ;; add room for scale height & Q_i & Z-s_Z & Z+s_Z

    shc=dblarr(timeOutputs,NABp1,nx,STVars)
    shcA=dblarr(timeOutputs,NABp1A,nx,STVars)
    shc[*,*,*,0:(STVars-stvpp-1)] = (temporary(starsHyperCube))[0:(timeOutputs-1),0:(NABp1-1),0:(nx-1),0:(STVars-stvpp-1)]
    shcA[*,*,*,0:(STVars-stvpp-1)] = (temporary(starsHyperCubeA))[0:(timeOutputs-1),0:(NABp1A-1),0:(nx-1),0:(STVars-stvpp-1)]

    x=dataCube[*,*,0]
    uu = dataCube[*,*,16-1]
    bbeta =  beta0 - beta0/(1.0+(b/x)^(beta0*nRotCurve));

    FOR i=0,NPassive DO BEGIN
	    shc[*,i,*,STVars-4] = shc[*,i,*,2]*dataCube[*,*,0]*Radius/(uu*vphiR);; scale height = sigma_r / Omega
    	shc[*,i,*,STVars-3] = shc[*,i,*,1]*sqrt(2*(bbeta+1))*uu/(!pi*chi*dataCube[*,*,0]*shc[*,i,*,0]) ;; Q_i
	    shc[*,i,*,STVars-2] = shc[*,i,*,3] - shc[*,i,*,4] ;; Z - sigma_Z
    	shc[*,i,*,STVars-1] = shc[*,i,*,3] + shc[*,i,*,4] ;; Z + sigma_Z
        shc[*,i,*,0] = shc[*,i,*,0]*md0*kmperkpc/(vphiR*Radius*speryear*1.0e6) ; Column density
        shc[*,i,*,1] = shc[*,i,*,1]*vphiR ;; r-velocity dispersion
        shc[*,i,*,2] = shc[*,i,*,2]*vphiR ;; z velocity dispersion
        ;shc[*,i,*,3] is metallicity in absolute units
        ;shc[*,i,*,4] is the square root of the variance of Z in absolute units
    ENDFOR
    FOR ii=0,NActive DO BEGIN
    	shcA[*,ii,*,STVars-4] = shcA[*,ii,*,2]*dataCube[*,*,0]*Radius/(uu*vphiR);; scale height
	    shcA[*,ii,*,STVars-3] = shcA[*,ii,*,1]*sqrt(2*(bbeta+1))*uu/(!pi*chi*dataCube[*,*,0]*shcA[*,ii,*,0]) ;; Q_i
    	shcA[*,ii,*,STVars-2] = shcA[*,ii,*,3] - shcA[*,ii,*,4] ;; Z - sigma_Z
	    shcA[*,ii,*,STVars-1] = shcA[*,ii,*,3] + shcA[*,ii,*,4] ;; Z + sigma_Z
        shcA[*,ii,*,0] = shc[*,ii,*,0]*md0*kmperkpc/(vphiR*Radius*speryear*1.0e6)
        shcA[*,ii,*,1] = shc[*,ii,*,1]*vphiR
        shcA[*,ii,*,2] = shc[*,ii,*,2]*vphiR
    ENDFOR
    starsHyperCube=(temporary(shc))[*,*,*,*]
    starsHyperCubeA=(temporary(shcA))[*,*,*,*]
    stellarAges=(temporary(stellarAges))[0:timeOutputs-1,0:NABp1-1] ;; trim down the stellarAges array


    conv=dblarr(ncolconv)


    NPostProcess= 56
    tdc=dblarr(n_elements(dataCube[*,0,0]),n_elements(dataCube[0,*,0]),n_elements(dataCube[0,0,*])+NPostProcess)
    tdc[*,*,0:(ncolStep-1)] = (temporary(dataCube))[*,*,*]
    ; ZstNum and ZstDen:  ( time steps ) x 1 x ( nx ) x 1
    ZstNum = dblarr(n_elements(starsHyperCubeA[*,0,0,0]),1,n_elements(starsHyperCubeA[0,0,*,0]),1)
    ZstDen = dblarr(n_elements(starsHyperCubeA[*,0,0,0]),1,n_elements(starsHyperCubeA[0,0,*,0]),1)
    FOR i=0,NPassive DO BEGIN
;    	tdc[*,*,ncolstep+npostprocess+i*STVars] = starsHyperCube[*,i,*,0] ;;; * mdotext0/(vphiR*1d5*Radius*cmperkpc) * cmperpc*cmperpc/MSol ;; g/cm^2 * (cm/pc)^2 *msol/g
;    	tdc[*,*,ncolstep+npostprocess+i*STVars+1] = starsHyperCube[*,i,*,1] ;;; * vphiR
;        tdc[*,*,ncolstep+npostprocess+i*STVars+2:ncolstep+npostprocess+(i+1)*STVars-1] = starsHyperCube[*,i,*,2:STVars-1]
    ENDFOR
    FOR ii=0,NActive DO BEGIN
        ZstNum[*,0,*,0] += starsHyperCubeA[*,ii,*,3]*starsHyperCubeA[*,ii,*,0]
       	ZstDen[*,0,*,0] += starsHyperCubeA[*,ii,*,0]
    ENDFOR

    tdc[0,*,39-1] = tdc[1,*,39-1] ;; currently the initial value is zero.
    evArray[9-1,0] = evArray[9-1,1]

    tdc[*,*,ncolStep]=  abs(tdc[*,*,3])/(abs(tdc[*,*,7]) ) * 2*!pi*Radius*kmperkpc/(vphiR*speryear)
    tdc[*,*,ncolStep+1]=abs(tdc[*,*,4])/(abs(tdc[*,*,8]) ) * 2*!pi*Radius*kmperkpc/(vphiR*speryear)
    tdc[*,*,ncolStep+2]=abs(tdc[*,*,5])/(abs(tdc[*,*,9]) ) * 2*!pi*Radius*kmperkpc/(vphiR*speryear)
    tdc[*,*,ncolStep+3]=abs(tdc[*,*,6])/(abs(tdc[*,*,10]) )* 2*!pi*Radius*kmperkpc/(vphiR*speryear)
    tdc[*,*,ncolStep+4]=abs(tdc[*,*,4])/(abs(tdc[*,*,6])) ;; gas / stellar temp
    ;; scale heights: 
    tdc[*,*,ncolStep+5]=tdc[*,*,6]*(tdc[*,*,6]/(tdc[*,*,5]+tdc[*,*,3]))*Radius/(chi*!pi) ;; stellar scale height
    tdc[*,*,ncolStep+6]=tdc[*,*,4]*(tdc[*,*,4]/(tdc[*,*,5]*tdc[*,*,4]/tdc[*,*,6]+tdc[*,*,3]))*Radius/(chi*!pi) ;; gas scale height

    tdc[*,*,ncolStep+7]=tdc[*,*,4]*tdc[*,*,4]*tdc[*,*,4]*tdc[*,*,4]*Radius*cmperkpc* $
        mdotext0/(chi*chi*vphiR*1d5*tdc[*,*,3] * MSol)  $
        *  (sign(tdc[*,*,48-1] - fH2Min*1.0001)+1.000001)/2;;;  2d Jeans Mass in Solar Masses in the SF region

    evArray[1,*] = evArray[1,*]*2*!pi*Radius/(vphiR*.977813952) ;; kpc/(km/s) = .9778 Gyr -- convert times to Ga.

    ;;; Replace Toomre Mass with 2d Jeans Mass
    FOR i=0,n_elements(evArray[7,*])-1 DO BEGIN
    	evArray[7,i] = MAX(tdc[i,*,ncolStep+7])
    ENDFOR
    ;;; Replace Disk Mass with its Dimensional Version
    evArray[5,*] = evArray[5,*] * mdotext0*Radius*cmperkpc/(vphiR*1d5*MSol) ;; g/s*kpc*cm/kpc / (km/s*cm/km*g/msol) = msol
    evArray[10,*] = evArray[10,*] * mdotext0 /(MSol/speryear) 
    evArray[8,*] = evArray[8,*] * mdotext0 /(MSol/speryear) 
    evArray[3,*] = evArray[3,*] * (mdotext0/(MSol/speryear)) * (2*!pi*Radius*cmperkpc/(vphiR*1d5))

    tdc[*,*,ncolStep+8] = tdc[*,*,0]*Radius ;; r
    tdc[*,*,ncolStep+9] = tdc[*,*,3]*mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*Msol)  ;; Sigma_g -- in  (g/s) * (cm/pc)*(cm/pc) / (km/s * kpc * cm/kpc * (cm/km) (g/Msun)) = Msun/pc^2
    tdc[*,*,ncolStep+24-1]=alog10(tdc[*,*,ncolStep+9])
    tdc[*,*,ncolStep+10] = tdc[*,*,4]*vphiR ;; sigma_g
    tdc[*,*,ncolStep+11] = tdc[*,*,5]*mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) ;; Sigma_* --  solar masses / pc^2
    tdc[*,*,ncolStep+25-1]=alog10(tdc[*,*,ncolStep+11])
    tdc[*,*,ncolStep+12] = tdc[*,*,6]*vphiR ;; sigma_*
    ;; g/s cm/pc cm/pc pc/kpc pc/kpc cm/s s/yr / (km/s kpc cm/kpc cm/km g/Msun cm)
    ;; Msun /kpc^2/yr
    tdc[*,*,ncolStep+13] = tdc[*,*,28] * mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) *1d6 *(vphiR*1d5*speryear)/(2*!pi*Radius*cmperkpc) ;; col SFR-  solar masses/kpc^2/yr
    tdc[*,*,ncolStep+14] = tdc[*,*,51] * mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) ;; cumulative SF
    tdc[*,*,ncolStep+15] = -1.0*tdc[*,*,34] * vphiR ;; vr_*
    tdc[*,*,ncolStep+16] = -1.0*tdc[*,*,36] * vphiR ;; vr_gas
	
    ;;ang. mom. of.. 
    tdc[*,*,ncolstep+28-1]= tdc[*,*,3]*tdc[*,*,0] ;; gas
    tdc[*,*,ncolstep+29-1]= tdc[*,*,5]*tdc[*,*,0] ;; stars
    tdc[*,*,ncolstep+30-1]= tdc[*,*,51]*tdc[*,*,0]*MLF ;; cumulative outflow

    ;; compute the average age of stars as a function of radius and time
    ;; compute both the average age at z=0, and the age in situ.
    FOR i=0,NPassive DO BEGIN
	    stellarAgesIS[*,i] = stellarAges[*,i] - (evArray[1,n_elements(evArray[1,*])-1]-evArray[1,*]) * 1.0e9
    ENDFOR
    FOR xi=0,nx-1 DO BEGIN ;; loop over space
	FOR zi=0,n_elements(tdc[*,0,0])-1 DO BEGIN ;; loop over time
  	    tdc[zi,xi,ncolstep+31-1] = $
		TOTAL(starsHyperCube[zi,*,xi,0] * stellarAges[zi,*])/TOTAL(starsHyperCube[zi,*,xi,0])
	    tdc[zi,xi,ncolstep+32-1] = $
		TOTAL(starsHyperCube[zi,*,xi,0] * stellarAgesIS[zi,*])/TOTAL(starsHyperCube[zi,*,xi,0])
	ENDFOR
    ENDFOR


    tdc[*,*,ncolstep+33-1] = -tdc[*,*,1] ;; - tau ( for purposes of log plotting)
    tdc[*,*,ncolstep+34-1] = -tdc[*,*,2] ;; - tau'    '''

    tdc[0,*,ncolstep+14-1] = tdc[1,*,ncolstep+14-1] ;; replace the 0th time step (all zeroes) with the first.
    fac = 1.0;; fac = MLF+Rf
    tdc[*,*,ncolstep+35-1] = tdc[*,*,ncolstep+10-1]*pcperkpc*pcperkpc/(fac*tdc[*,*,ncolstep+14-1]);; depletion time = col/coldtSF --- Msol/pc^2 / (Msol/kpc^2/yr) * (1000 pc/kpc)^2
    tdc[*,*,ncolstep+36-1] = tdc[*,*,ncolstep+9-1]*kmperkpc/((tdc[*,*,ncolstep+17-1])*speryear);; viscous time = r / v_r --- kpc/(km/s)
    tdc[*,*,ncolstep+37-1] = tdc[*,*,ncolstep+35-1]*tdc[*,*,ncolstep+17-1]*speryear/(tdc[*,*,ncolstep+9-1]*kmperkpc);; depletion time / viscous time
    ; 1/t_visc
    tdc[*,*,ncolstep+38-1] = tdc[*,*,ncolstep+17-1]*speryear / (tdc[*,*,ncolstep+9-1]*kmperkpc) ; year^-1
    x25s = dblarr(n_elements(tdc[*,0,0]))
    x25s_2 = x25s[*]
    x25s_3 = x25s[*]
    x25s_4 = x25s[*]
    colTranses = x25s[*]
    FOR ti=0, n_elements(tdc[*,0,0])-1 DO BEGIN
        ; Locations where fH2 cross 0.5, i.e. (fH2(x_i)-.5)*(fH2(x_{i+1})-.5)<=0
        wh = WHERE((tdc[ti,0:(nx-2),48-1]-0.5)*(tdc[ti,1:(nx-1),48-1]-0.5) LE 0.0, count)
        ;; If no such locations exist, just take a guess
        IF(count LE 0) THEN BEGIN 
            ;colTrans = 6.75558 ; KMT value for Zsun
            ;x25_4 = 1.0
            maxfH2 = MAX(tdc[ti,*,48-1],loc)  
            colTrans = tdc[ti,loc,ncolstep+10-1] * 0.5 / maxfH2
            x25_4 = tdc[ti,loc,0]*1.0/0.449659
        ENDIF
        IF(count GE 1) THEN BEGIN
            colTrans = sqrt(tdc[ti,wh[count-1],ncolstep+10-1] * tdc[ti,wh[count-1]+1,ncolstep+10-1])
            x25_4 = sqrt(tdc[ti,wh[count-1],0]*tdc[ti,wh[count-1]+1,0]) * 1.0/0.449659
        ENDIF
        colTranses[ti] = colTrans
        dummy = MIN(abs(tdc[ti,*,ncolstep+12-1] - 6.619),index)
        x25 = tdc[ti,index,0] ; location of 6.619 solar masses /pc^2 in stellar mass
        dummy = MIN(abs(tdc[ti,*,ncolstep+14-1] - 6.74d-4),index)
        x25_2= tdc[ti,index,0]  ; -19.5
        dummy = MIN(abs(tdc[ti,*,ncolstep+14-1] - 1.19d-3),index)
        x25_3 = tdc[ti,index,0] ; -18.9
        ; dimensionless ratio of real value over predicted value by Bigiel & Blitz
        ;tdc[ti,*,ncolstep+39-1] = tdc[ti,*,ncolstep+10-1]/(2.1*colTrans*exp(-1.65*tdc[ti,*,0]/x25))
        ;tdc[ti,*,ncolstep+40-1] = tdc[ti,*,ncolstep+10-1]/(2.1*colTrans*exp(-1.65*tdc[ti,*,0]/x25_2))
        ;tdc[ti,*,ncolstep+41-1] = tdc[ti,*,ncolstep+10-1]/(2.1*colTrans*exp(-1.65*tdc[ti,*,0]/x25_3))
        tdc[ti,*,ncolstep+44-1] = tdc[ti,*,ncolstep+10-1]/(2.1*colTrans*exp(-1.65*tdc[ti,*,0]/x25_4))
        
        x25s[ti] = x25
        x25s_2[ti] = x25_2
        x25s_3[ti] = x25_3
        x25s_4[ti] = x25_4
    ENDFOR

    ; accretion rate & timescale
    ;; md0/(vphiR Radius) / (2piR/vphiR) = md0/(2pi R^2)
    tdc[*,*,ncolstep+42-1] = tdc[*,*,30-1] * md0 / (2.*!pi * Radius*Radius * 1.0d6) ; Msun/yr/pc^2
    tdc[*,*,ncolstep+43-1] = tdc[*,*,ncolstep+10-1]/tdc[*,*,ncolstep+42-1] ;; yr

    ; H2 depletion timescale = fH2 * depletion time
    tdc[*,*,ncolstep+39-1] = tdc[*,*,48-1]*tdc[*,*,ncolstep+35-1]


    ;; This one I think will be slightly more useful: just (mu+f_R)*colSF / colAccr
    tdc[*,*,ncolstep+41-1] = (mlf+Rf)*tdc[*,*,ncolstep+14-1] *1.0e-6 / tdc[*,*,ncolstep+42-1]


    tdc[*,*,ncolStep+18-1] = ZstNum[*,0,*,0]/ZstDen[*,0,*,0];; Z_*
    tdc[*,*,ncolStep+19-1] = -1.0*tdc[*,*,2]*mdotext0*speryear /(MSol*uu*(1+bbeta)) ; mdot (msol/yr) - centered, GI only. tdc[*,*,39-1] is mdot (MSun/yr) incl MRI

    dlnx=-alog(tdc[0,0,0])/float(n_elements(tdc[0,*,0])-1)

    ; dcoldt_transport / col_accr
    tdc[*,*,ncolStep+45-1] = (tdc[*,*,8-1]*(md0/(Radius*2.0*!pi*Radius*1.0d6)) - tdc[*,*,ncolstep+42-1]) / tdc[*,*,ncolstep+42-1] + tdc[*,*,ncolstep+41-1]
    ;; dcoldt_tr / (|dcoldt_tr| + |dcoldt_accr| + |dcoldt_SF| )
    tdc[*,*,ncolstep+46-1] = tdc[*,*,ncolstep+45-1]/((abs(tdc[*,*,ncolstep+45-1])+1+tdc[*,*,ncolstep+41-1]) )

    ; ratio of Sigma to Sigma_eq (low column density)
    ; Msun/pc^2   / (3 pi colAccr(Msun/pc^2/yr)^2 * Qlim * sigmath*10^5(cm/s) * radius*cmperkpc (cm) / (32 epsff^2 * fH2^2 * vphi*10^5(cm/s) * sqrt(2(b+1)) * G*gperMsun (cm^3/s^2 Msun) *(fR+mlf)^2))^(1/3)
    tdc[*,*,ncolstep+40-1] = tdc[*,*,ncolstep+10-1]/((3.0*!pi*tdc[*,*,ncolstep+42-1]^2*Qlim*sigmath * tdc[*,*,ncolstep+9-1]*cmperkpcpers2peryear2*cm2perpc2perMsol / (32.0*eps_ff^2 * tdc[*,*,48-1] * tdc[*,*,48-1] * tdc[*,*,16-1]*vphiR * sqrt(2.0*(bbeta+1)) * G*(Rf+mlf)^2))^(1./3.)) 

    ; viscous timescale = integral (gas mass interior to r) / mdot(r)
    tdc[*,*,ncolstep+47-1] = TOTAL(2.0*!pi*x*Radius*x*Radius*sinh(dlnx)*tdc[*,*,ncolstep+10-1]*1.0d6,2,/cumulative)/tdc[*,*,39-1]

    ; SF timescale = integral (gas mass interior to r) / integral(SF interior to r)
    tdc[*,*,ncolstep+48-1] = TOTAL(2.0*!pi*x*Radius*x*Radius*sinh(dlnx)*tdc[*,*,ncolstep+10-1]*1.0d6,2,/cumulative)/TOTAL(2.0*!pi*x*Radius*x*Radius*sinh(dlnx)*tdc[*,*,ncolstep+14-1],2,/cumulative)

    ;; dcoldt_accr / all
    tdc[*,*,ncolstep+49-1] = tdc[*,*,ncolstep+46-1]/tdc[*,*,ncolstep+45-1]

    ;; dcoldt_SF / all
    tdc[*,*,ncolstep+50-1] = tdc[*,*,ncolstep+49-1] *tdc[*,*,ncolstep+41-1]

    ; dZ/dt dilution term
    tdc[*,*,ncolstep+51-1] = tdc[*,*,30-1] * (ZIGM - tdc[*,*,22-1])/tdc[*,*,4-1]
    ; dZ/dt sf term
    yREC=.054
    tdc[*,*,ncolstep+52-1] = yREC*(1.0-Rf)*zetaREC*tdc[*,*,29-1]/tdc[*,*,4-1]

    ;; Are we in equilibrium? -1 to 1, equilibrium ~ 0
    tdc[*,*,ncolstep+53-1] = tdc[*,*,8-1]/(tdc[*,*,29-1]*(MLF+Rf) + tdc[*,*,30-1] + ABS(tdc[*,*,8-1] -tdc[*,*,29-1]*(MLF+Rf) - tdc[*,*,30-1]))

    ;; Same quantity as above but normalized to Accretion Rate
    tdc[*,*,ncolstep+54-1] = tdc[*,*,8-1]/(tdc[*,*,30-1])

    ;; Sigma/SigmaCrit
    tdc[*,*,ncolstep+55-1] = (2.0/3.0) * (tdc[*,*,ncolstep+11-1]/tdc[*,*,24-1]) * (fixedQ/sigmath)

    ;; velocity dispersion / analytic
    NN = fixedQ * G * evArray[20-1,*]*Msol/speryear / (6.0*eta*sigmath^3*1.0d15)
    ;xiIn = xmin *exp(-dlnx) / (accScaleLength / Radius)
;    xiOut = TOTAL( tdc[*,*,ncolstep+10-1] * 2 * !pi * tdc[*,*,0]*tdc[*,*,0]*Radius*Radius*sinh(dlnx/2.0)*2.0 * 1.0d6, 2) *MSol *G *fixedQ / (3.0*vphiR*1.0d5*sigmath*1.0d5 * cmperkpc * accScaleLength)
    xiOut=NN*0
    xiIn=NN*0
    FOR iii=0,n_elements(xiOut)-1 DO BEGIN
        wh = WHERE( tdc[iii,*,12-1] LE fixedQ ,ct)
        IF(ct GT 0) THEN xiOut[iii] = MAX( tdc[iii,wh,ncolstep+9-1] ) ELSE xiOut[iii]=-max(tdc[iii,*,ncolstep+9-1])
        IF(ct GT 0) THEN xiIn[iii] = MIN( tdc[iii,wh,ncolstep+9-1] ) ELSE xiIn[iii] = -xmin*exp(-dlnx)*Radius/accScaleLength
    ENDFOR
    xiOut=xiOut/accScaleLength
    xiIn = xiIn/accScaleLength
    xiOutArr= tdc[*,*,0]*0
    xiInArr = tdc[*,*,0]*0
    NNarr=tdc[*,*,0]*0
    xi = tdc[*,*,ncolstep+9-1] / accScaleLength
    FOR iii=0,n_elements(xiOut)-1 DO xiOutArr[iii,*] = xiOut[iii]
    FOR iii=0,n_elements(xiIn)-1 DO xiInArr[iii,*] = xiIn[iii]
    FOR iii=0,n_elements(NN)-1 DO NNarr[iii,*] = NN[iii]
;    arg = (NNarr*(exp(-xiOutArr)*(1.0-xiIn/xi)*(2+xiOutArr/xi)/(xiOutArr-xiIn) - exp(-xiIn)*(2+xiIn)/(xiOutArr-xiIn) - exp(-xi)*(1.0+2/xi)))
    dXi = xiOutArr - xiInArr
    wh=WHERE(dXi EQ 0, ct)
    IF(ct GT 0) THEN dXi[wh]=.01
    arg = NNarr*(-exp(-xi)*(1.0+2.0/xi) + 1.0/dXi * (-exp(-xiInArr)*(2.0+xiInArr) + exp(-xiOutArr)*(2+xiOutArr)) + 1.0/(xi*dXi)*(exp(-xiInArr)*(2.0+xiInArr)*xiOutArr - exp(-xiOutArr)*(2.0+xiOutArr)*xiInArr))
    tdc[*,*,ncolstep+56-1] = tdc[*,*,ncolstep+11-1]/ $
        ( .5*(sign(xiOutArr-xi)+1.0)*.5*(sign(xi-xiInArr)+1.0)*(sigmath*sigmath*((arg*arg)^(1.0/3.0) + 1))^(1.0/2.0) + .5*(sign(xi-xiOutArr)+1.0)*sigmath + .5*(sign(xiInArr-xi)+1.0)*sigmath)



;    thetotal = tdc[*,*,ncolstep+50-1] + tdc[*,*,ncolstep+49-1] + abs(tdc[*,*,ncolstep+46-1])
	;also = tdc[*,*,ncolstep+46-1]*((tdc[*,*,ncolstep+41-1]+1)/tdc[*,*,ncolstep+45-1] + 1)

	nt = n_elements(tdc[*,0,0]) ;; nx already defined
	initialStellarMass = dblarr(nt,nx)
	initialGasMass = dblarr(nt,nx)
	FOR n=0, nt-1 DO BEGIN
		initialStellarMass[n,*] = tdc[0,*,ncolstep+11]
		initialGasMass[n,*] = tdc[0,*,ncolstep+9]
	ENDFOR

	;;; 	tdc    index mapping (indexed from 1)
	; 1- x,     2- tau,    3- tau',     4- S,        5- s,   6- S_*,    7- s_*
	; 8- dS/dT, 9- ds/dT, 10- 0, 11- ds_*/dT, 12- Q,  13- dZDiskdtAdv,    14- MdotMRI(Msun/yr)
	;15- beta,   16- uu,    17- f_g,     18- q        19- Toomre Length, 20- Toomre Mass
	;21- dZ/dt,22- Z,     23- Q_*      24- Q_g      25- Q_Raf 26- Q_WS,27- Q_RW, 
	;28- Q(q) ,29- dS/dT_SF, 30- colAccr   31- dQ/dS 32- dQ/ds 33- dQ/dS_err 
	;34- dQ/ds_err 35- y, 36- torqueErr 37- vrg, 38- Cumulative Stars Out, 
	;39- dimensional Mdot, 40- dsigdtTrans
	;41- dsigdtDdx
	;42- dsigdtHeat , 43- ddx(tau'), 44- s', 
	;45- dsigdtCool, 46- dZDiskdtDiff, 47- alpha, 
	;48- f_H2, 49- cuTorqueErr, 50- cuTorqueErr2,
	;51- tau''(meas),52- cuSF
        ;53- dcoldtIncoming 54- dcoldtOutgoing

	;	ncolstep+ ... (still indexed from 1)
	;1- col timescale, 2- sig timescale, 3- col* timescale, 4- sig* timescale,
	;5- phi, 6- stellar scale height, 7- gas scale height, 8- 2d Jeans Mass
	;9- radius, 10- Sigma_g, 11- sigma_g, 12- Sigma_*, 13- sigma_*, 14- col SFR
	;15- cuSFR, 16- -vr_*, 17- -vr_g, 18- Z_*, 19- mdot,
	;20- gas in - gas out - (Rf+mu)*SF
	;21- st in - st out + Rf*SF
	;22- Gas fractional error in mass
	;23- Stellar fractional error in mass
	;24- log Sigma_g, 25- log Sigma_*
	;26- 0, 27- 0 
	;28- ang. mom of gas, 29- ang. mom of stars
	;30- ang. mom of cumulative outflows
    ;; 33,34 = -tau, -tau'
    ;; 35,36,37: depletion time, viscous time, dep/visc
    ;; 38,39 v_rg/r=1/tvisc, tdepH2
    ;; 40,41 Sigma/SigmaEQ, (mu+Rf)*colSFR/colAccr
    ;; 42- Accretion col dens 43- accretion timescale, 44- BB4
    ;; 45- dcoldt_transport/ dcoldt_accr
    ;; 46- (dcoldt_transport+dcoldt_SF)/dcoldt_accr
    ;; 47 - integrated t_visc,  48 - t_SF within r
    ;; 49 - dcoldt_accr/all,  50 - dcoldt_SF/all
    ;; 51- dZDiskdt_dil, 52- dZDiskdt_sf -- these two are dimensionless
    ;; 53- equilibrium number, 54- dcoldt/dcoldt_accr, 55- Sigma/SigmaCrit
    ;; 56 - sigma/sigma_an

	;	ncolstep + npostprocess +... (still indexed from 1)
	; 1- S_*0,   2- s_*0,  3- Z_*0, 4- ScHeight (kpc), 5- Q_i ( initial population of stars )
	; stv*i+1- S_*,  stv*i+2- s_*i,  stv*i+3- Z_*i   --- i=0 thru NPassive+1


	;; Normalize all metallicities to solar and take the log10.
;	convert = [22,ncolstep+18,ncolstep+npostprocess+indgen(NPassive+1)*STVars+3, $
;		   ncolstep+npostprocess+indgen(NPassive+1)*STVars+7,  $
;		   ncolstep+npostprocess+indgen(NPassive+1)*STVars+8]-1
        convert = [22,ncolstep+18]-1
	tdc[*,*,convert] = alog10(tdc[*,*,convert]/0.02 + 1.0d-8)
	FOR j=0,n_elements(convert)-1 DO BEGIN
		convSquare = tdc[*,*,convert[j]] 
		c=where(convSquare LT -2,ct)
		IF(ct NE 0) THEN convSquare[c] = -2
		tdc[*,*,convert[j]] = convsquare
	ENDFOR

	starsHyperCube[*,*,*,3]=alog10(starsHyperCube[*,*,*,3]/0.02+1.0d-8)
	starsHyperCube[*,*,*,7:8]=alog10(starsHyperCube[*,*,*,7:8]/0.02+1.0d-8)
        starsHyperCubeA[*,*,*,3]=alog10(starsHyperCubeA[*,*,*,3]/.02+1.0d-8)
        starsHyperCubeA[*,*,*,7:8]=alog10(starsHyperCubeA[*,*,*,7:8]/.02+1.0d-8)

	convCube=starsHyperCube[*,*,*,3]
	c=where(convCube LT -2,ct)
	IF(ct NE 0) THEN convCube[c] = -2
	starsHyperCube[*,*,*,3]=convCube

	convHCube=starsHyperCube[*,*,*,7:8]
	c=where(convHCube LT -2,ct)
	IF(ct NE 0) THEN convHCube[c] = -2
	starsHyperCube[*,*,*,7:8]=convHCube


    wh = where(tdc NE tdc,ct)
    IF(ct GT 0) THEN stop,"NaN found in data array."


	model = {name:name,nx:nx,tmax:tmax,maxstep:nstepmax, $
		radius:Radius,vphiR:vphiR,eta:ETA,epsff:EPS_ff,fg0:FG0, $
		tauHeat:tauHeat,cosmology:CosmologyOn,analyticQ:analyticQ, $
		tol:TOL,comment:comment,evArray:evArray,dataCube:tdc,$
		stellarData:starsHyperCube,stellarAges:stellarAges,$
		convergence:conv,ncolev:ncolev,ncolstep:ncolstep, $
		NPostProcess:NPostProcess,NPassive:NPassive,STVars:STVars, $
		MLF:MLF,mdotext0:md0,Rf:Rf,whichAccHistory:whichAccHistory,$;}
		fixedQ:fixedQ,diskScaleLength:scaleLength,Qlim:Qlim,dlnx:dlnx,$
		mh0:mh0, minStsig:minStsig,zquench:zquench,xmin:xmin,$
		zrelax:zrelax, zetaREC:zetaREC, deltaOmega:deltaOmega, $
        NoutputsNominal:NoutputsNominal,  ND08attempts:ND08attempts, $
        x25s_2:x25s_2, x25s_3:x25s_3, x25s:x25s, x25s_4:x25s_4, $
        colTranses:colTranses, NN:NN, sigmath:sigmath,  $
        alphaAccProf:alphaAccProf, accScaleLength:accScaleLength, $
        rotCurveTurnover:b}
	
	RETURN,model ;; end of the readOutput function
END


