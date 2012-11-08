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
    "-tau","-tau'","depletion time","viscous time","tDep/tVisc"] ;; 33..34

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
  innerPowerLaw=ExtractCmtFlt(lunCom)
  softening=ExtractCmtFlt(lunCom)
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
  md0 = ExtractCmtFlt(lunCom)
  ND08attempts = ExtractCmtL(lunCom)

  ;; Now we need to rescale some values for mh0!=1e12
 ; MLF = MLF * (Mh0/1.0e12) ^ (-1.0/3.0);
 ; vphiR = vphiR * (Mh0/1.0e12)^  (1.0/3.0)
	

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

  ;; Some useful constants
  MSol = 1.98892d33 ;; solar masses in grams
  G=6.673d-8 ;; newton's constant in cgs
  speryear = 31556926d
  cmperkpc = 3.08568025d21	
  cmperpc =  3.08568025d18

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
		READU,lunStars,age
		stellarAges[timeOutputs,i]=age
		READU,lunStars,contOfAgeBin
		starsHyperCube[timeOutputs,i,0:nx-1,0:(STVars-1)]=contOfAgeBin
;;;	ENDIF
      ENDFOR
      FOR ii=0, stsSizeA-1 DO BEGIN
;;;	IF(EOF(lunStA) NE 1) THEN BEGIN
		READU,lunStarsA,age
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
    u=x/sqrt(b*b+x*x)
    beta=b*b/(b*b+x*x)

    FOR i=0,NPassive DO BEGIN
	shc[*,i,*,STVars-4] = shc[*,i,*,1]*dataCube[*,*,0]*Radius/u;; scale height
	shc[*,i,*,STVars-3] = shc[*,i,*,1]*sqrt(2*(beta+1))*u/(!pi*chi*dataCube[*,*,0]*shc[*,i,*,0]) ;; Q_i
	shc[*,i,*,STVars-2] = shc[*,i,*,2] - shc[*,i,*,3] ;; Z - sigma_Z
	shc[*,i,*,STVars-1] = shc[*,i,*,2] + shc[*,i,*,3] ;; Z + sigma_Z
    ENDFOR
    FOR ii=0,NActive DO BEGIN
	shcA[*,ii,*,STVars-4] = shcA[*,ii,*,1]*dataCube[*,*,0]*Radius/u;; scale height
	shcA[*,ii,*,STVars-3] = shcA[*,ii,*,1]*sqrt(2*(beta+1))*u/(!pi*chi*dataCube[*,*,0]*shcA[*,ii,*,0]) ;; Q_i
	shcA[*,ii,*,STVars-2] = shcA[*,ii,*,2] - shcA[*,ii,*,3] ;; Z - sigma_Z
	shcA[*,ii,*,STVars-1] = shcA[*,ii,*,2] + shcA[*,ii,*,3] ;; Z + sigma_Z
    ENDFOR
    starsHyperCube=(temporary(shc))[*,*,*,*]
    starsHyperCubeA=(temporary(shcA))[*,*,*,*]
    stellarAges=(temporary(stellarAges))[0:timeOutputs-1,0:NABp1-1] ;; trim down the stellarAges array


;    OPENR,lunConv,(name+"_convergence.dat"),/GET_LUN
    conv=dblarr(ncolconv)
;    IF(MAX(evArray[1,*]) GT 25.2) THEN READF,lunConv,conv
;    FREE_LUN,lunConv

;;   Change the reported errors in the numerical derivatives to absolute values
;    errs = dataCube[*,*,60:63]
;    errs=abs(errs)
;    errs[where(errs LT 1.e-21)] = 1.e-21
;    dataCube[*,*,60:63] = errs

    NPostProcess= 37
    tdc=dblarr(n_elements(dataCube[*,0,0]),n_elements(dataCube[0,*,0]),n_elements(dataCube[0,0,*])+NPostProcess+(NPassive+1)*STVars)
    tdc[*,*,0:(ncolStep-1)] = (temporary(dataCube))[*,*,*]
    ZstNum = dblarr(n_elements(starsHyperCubeA[*,0,0,0]),1,n_elements(starsHyperCubeA[0,0,*,0]),1)
    ZstDen = dblarr(n_elements(starsHyperCubeA[*,0,0,0]),1,n_elements(starsHyperCubeA[0,0,*,0]),1)
    FOR i=0,NPassive DO BEGIN
;	tdc[*,*,ncolstep+npostprocess+i*STVars:ncolstep+npostprocess+i*STVars-1] = starsHyperCube[*,i,*,*]
	tdc[*,*,ncolstep+npostprocess+i*STVars] = starsHyperCube[*,i,*,0] * mdotext0/(vphiR*1d5*Radius*cmperkpc) * cmperpc*cmperpc/MSol ;; g/cm^2 * (cm/pc)^2 *msol/g
	tdc[*,*,ncolstep+npostprocess+i*STVars+1] = starsHyperCube[*,i,*,1] * vphiR
;	ZstNum[*,0,*,0] += starsHyperCube[*,i,*,2]*starsHyperCube[*,i,*,0]
;	ZstDen[*,0,*,0] += starsHyperCube[*,i,*,0]
	tdc[*,*,ncolstep+npostprocess+i*STVars+2:ncolstep+npostprocess+(i+1)*STVars-1] = starsHyperCube[*,i,*,2:STVars-1]
    ENDFOR
    FOR ii=0,NActive DO BEGIN
	ZstNum[*,0,*,0] += starsHyperCubeA[*,ii,*,2]*starsHyperCube[*,ii,*,0]
	ZstDen[*,0,*,0] += starsHyperCubeA[*,ii,*,0]
    ENDFOR



    tdc[0,*,39-1] = tdc[1,*,39-1] ;; currently the initial value is zero.
    evArray[9-1,0] = evArray[9-1,1]

    tdc[*,*,ncolStep]=  abs(tdc[*,*,3])/(abs(tdc[*,*,7]) + .00001)
    tdc[*,*,ncolStep+1]=abs(tdc[*,*,4])/(abs(tdc[*,*,8]) + .00001)
    tdc[*,*,ncolStep+2]=abs(tdc[*,*,5])/(abs(tdc[*,*,9]) + .00001)
    tdc[*,*,ncolStep+3]=abs(tdc[*,*,6])/(abs(tdc[*,*,10]) + .00001)
    tdc[*,*,ncolStep+4]=abs(tdc[*,*,4])/(abs(tdc[*,*,6])) ;; gas / stellar temp
    ;; scale heights: 
    tdc[*,*,ncolStep+5]=tdc[*,*,6]*(tdc[*,*,6]/(tdc[*,*,5]+tdc[*,*,3]))*Radius/(chi*!pi) ;; stellar scale height
    tdc[*,*,ncolStep+6]=tdc[*,*,4]*(tdc[*,*,4]/(tdc[*,*,5]*tdc[*,*,4]/tdc[*,*,6]+tdc[*,*,3]))*Radius/(chi*!pi) ;; gas scale height

    tdc[*,*,ncolStep+7]=tdc[*,*,4]*tdc[*,*,4]*tdc[*,*,4]*tdc[*,*,4]*Radius*cmperkpc*mdotext0/(chi*chi*vphiR*1d5*tdc[*,*,3] * MSol)  ;;;  2d Jeans Mass in Solar Masses

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
    tdc[*,*,ncolStep+9] = tdc[*,*,3]*mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*Msol)  ;; Sigma_g
    tdc[*,*,ncolStep+24-1]=alog10(tdc[*,*,ncolStep+9])
    tdc[*,*,ncolStep+10] = tdc[*,*,4]*vphiR ;; sigma_g
    tdc[*,*,ncolStep+11] = tdc[*,*,5]*mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) ;; Sigma_* --  solar masses / pc^2
    tdc[*,*,ncolStep+25-1]=alog10(tdc[*,*,ncolStep+11])
    tdc[*,*,ncolStep+12] = tdc[*,*,6]*vphiR ;; sigma_*
    tdc[*,*,ncolStep+13] = tdc[*,*,28] * mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) *1d6 *(vphiR*1d5*speryear)/(2*!pi*Radius*cmperkpc) ;; col SFR-  solar masses/kpc^2/yr
    tdc[*,*,ncolStep+14] = tdc[*,*,51] * mdotext0*cmperpc*cmperpc/(vphiR*Radius*cmperkpc*1d5*MSol) ;; cumulative SF
    tdc[*,*,ncolStep+15] = -1.0*tdc[*,*,34] * vphiR ;; vr_*
    tdc[*,*,ncolStep+16] = -1.0*tdc[*,*,36] * vphiR ;; vr_gas
    tdc[*,*,ncolStep+26-1] = alog10(tdc[*,*,ncolStep+16])
    tdc[*,*,ncolStep+27-1] = alog10(tdc[*,*,ncolStep+15])	
	
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


    speryear=31556926.0
    kmperkpc=3.08568025d16
    pcperkpc=1d3

    tdc[*,*,ncolstep+33-1] = -tdc[*,*,1] ;; - tau ( for purposes of log plotting)
    tdc[*,*,ncolstep+34-1] = -tdc[*,*,2] ;; - tau'    '''

    tdc[0,*,ncolstep+14-1] = tdc[1,*,ncolstep+14-1] ;; replace the 0th time step (all zeroes) with the first.
    fac = 1.0;; fac = MLF+Rf
    tdc[*,*,ncolstep+35-1] = tdc[*,*,ncolstep+10-1]*pcperkpc*pcperkpc/(fac*tdc[*,*,ncolstep+14-1]);; depletion time = col/coldtSF --- Msol/pc^2 / (Msol/kpc^2/yr) * (1000 pc/kpc)^2
    tdc[*,*,ncolstep+36-1] = tdc[*,*,ncolstep+9-1]*kmperkpc/((tdc[*,*,ncolstep+17-1]+.00001)*speryear);; viscous time = r / v_r --- kpc/(km/s)
    tdc[*,*,ncolstep+37-1] = tdc[*,*,ncolstep+35-1]/tdc[*,*,ncolstep+36-1];; depletion time / viscous time
    

    tdc[*,*,ncolStep+17] = ZstNum[*,0,*,0]/ZstDen[*,0,*,0];; Z_*
    tdc[*,*,ncolStep+18] = -1.0*tdc[*,*,2]*mdotext0*speryear /(MSol*u*(1+beta)) ; mdot (msol/yr)

    dlnx=-alog(tdc[0,0,0])/float(n_elements(tdc[0,*,0])-1)

;;;;;	tdc[*,*,40-1]=tdc[*,*,52-1] * 2*!pi*tdc[*,*,0]*tdc[*,*,0]*dlnx*md0*radius/vphiR * kmperkpc/speryear ;; cumulative star formation in a given cell - solar masses
	nt = n_elements(tdc[*,0,0]) ;; nx already defined
	initialStellarMass = dblarr(nt,nx)
	initialGasMass = dblarr(nt,nx)
	FOR n=0, nt-1 DO BEGIN
		initialStellarMass[n,*] = tdc[0,*,ncolstep+11]
		initialGasMass[n,*] = tdc[0,*,ncolstep+9]
	ENDFOR
;;;;;	tdc[*,*,41-1]=(tdc[*,*,ncolstep+11]-initialStellarMass)  * 2.0*!pi*tdc[*,*,ncolstep+8]*dlnx*tdc[*,*,ncolstep+8] * pcperkpc*pcperkpc
;;;;;	tdc[*,*,42-1]=(tdc[*,*,ncolstep+9]-InitialGasMass) * 2.0*!pi*tdc[*,*,ncolstep+8]*dlnx*tdc[*,*,ncolstep+8] * pcperkpc*pcperkpc

;;;;;	tdc[*,0:nx-2,45-1]=tdc[*,1:nx-1,38-1]
;;;;;	tdc[*,0:nx-2,46-1]=tdc[*,1:nx-1,39-1]
;;;;;	tdc[*,nx-1,45-1] = evArray[17-1,*]*0 ;; stars (no stars migrating through the boundary)
;;;;;	tdc[*,nx-1,46-1] = evArray[17-1,*] ;; gas
;;;;;	tdc[*,*,ncolstep+20-1] = tdc[*,*,46-1] - tdc[*,*,39-1] - (Rf+MLF)*tdc[*,*,40-1] ;gas
;;;;;;	tdc[*,*,ncolstep+21-1] = tdc[*,*,45-1] - tdc[*,*,38-1] + Rf*tdc[*,*,40-1] ;star
	
;;;;;	fracErrGas=dblarr(nt,nx)
;;;;;	fracErrSt =dblarr(nt,nx)
;;;;;;	tdc[*,*,ncolstep+22-1] 
;;;;;	fracErrGas = abs((tdc[*,*,ncolstep+20-1] - tdc[*,*,42-1])/(ABS(tdc[*,*,42-1])+.0000001))
;;;;;;	tdc[*,*,ncolstep+23-1] 
;;;;;	fracErrSt= abs((tdc[*,*,ncolstep+21-1] - tdc[*,*,41-1])/(ABS(tdc[*,*,41-1])+.0000001))

;;;;;	largeErr = where(ABS(fracErrGas - .5) GE .5 ,ctr)
;;;;;	IF(ctr NE 0 ) THEN fracErrGas[largeErr]=1
;;;;;	tdc[*,*,ncolstep+22-1]=fracErrGas
;;;;;	largeErr = where(ABS(fracErrSt - .5) GE .5,ctr)
;;;;;	IF(ctr NE 0) THEN fracErrSt[largeErr]=1
;;;;;	tdc[*,*,ncolstep+23-1]=fracErrSt

	;;; 	tdc    index mapping (indexed from 1)
	; 1- x,     2- tau,    3- tau',     4- S,        5- s,   6- S_*,    7- s_*
	; 8- dS/dT, 9- ds/dT, 10- dS_*/dT, 11- ds_*/dT, 12- Q,  13- h0,    14- h1
	;15- h2,   16- H,     17- f_g,     18- q        19- Toomre Length, 20- Toomre Mass
	;21- dZ/dt,22- Z,     23- Q_*      24- Q_g      25- Q_Raf 26- Q_WS,27- Q_RW, 
	;28- Q(q) ,29- dS/dT_SF, 30- tau''(torqueEq)    31- dQ/dS 32- dQ/ds 33- dQ/dS_err 
	;34- dQ/ds_err 35- y, 36- torqueErr 37- vrg, 38- Cumulative Stars Out, 
	;39- Cumulative Gas Out, 40- dimensional Cumulative SF, 
	;41- dimensional change in st mass in a cell, 
	;42- dimensional change in gas mass in a cell, 43- ddx(tau'), 44- s', 
	;45- Cumulative Stars In, 46- Cumulative Gas in, 47- alpha, 
	;48- f_H2, 49- cuTorqueErr, 50- cuTorqueErr2,
	;51- tau''(meas),52- cuSF

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
	;26- log vr_g, 27- log vr_*
	;28- ang. mom of gas, 29- ang. mom of stars
	;30- ang. mom of cumulative outflows
        ;; 33,34 = -tau, -tau'
        ;; 35,36,37: depletion time, viscous time, dep/visc

	;	ncolstep + npostprocess +... (still indexed from 1)
	; 1- S_*0,   2- s_*0,  3- Z_*0, 4- ScHeight (kpc), 5- Q_i ( initial population of stars )
	; stv*i+1- S_*,  stv*i+2- s_*i,  stv*i+3- Z_*i   --- i=0 thru NPassive+1

	;; Normalize all metallicities to solar and take the log10.
	convert = [22,ncolstep+18,ncolstep+npostprocess+indgen(NPassive+1)*STVars+3, $
		   ncolstep+npostprocess+indgen(NPassive+1)*STVars+7,  $
		   ncolstep+npostprocess+indgen(NPassive+1)*STVars+8]-1
	tdc[*,*,convert] = alog10(tdc[*,*,convert]/0.02)
	FOR j=0,n_elements(convert)-1 DO BEGIN
		convSquare = tdc[*,*,convert[j]] 
		c=where(convSquare LT -2,ct)
		IF(ct NE 0) THEN convSquare[c] = -2
		tdc[*,*,convert[j]] = convsquare
	ENDFOR

	starsHyperCube[*,*,*,2]=alog10(starsHyperCube[*,*,*,2]/0.02)
	starsHyperCube[*,*,*,6:7]=alog10(starsHyperCube[*,*,*,6:7]/0.02)

	convCube=starsHyperCube[*,*,*,2]
	c=where(convCube LT -2,ct)
	IF(ct NE 0) THEN convCube[c] = -2
	starsHyperCube[*,*,*,2]=convCube

	convHCube=starsHyperCube[*,*,*,6:7]
	c=where(convHCube LT -2,ct)
	IF(ct NE 0) THEN convHCube[c] = -2
	starsHyperCube[*,*,*,6:7]=convHCube

;	success=(file_lines(name+"_stde.txt") EQ 0)

	model = {name:name,nx:nx,tmax:tmax,maxstep:nstepmax, $
		radius:Radius,vphiR:vphiR,eta:ETA,epsff:EPS_ff,fg0:FG0, $
		tauHeat:tauHeat,cosmology:CosmologyOn,analyticQ:analyticQ, $
		tol:TOL,comment:comment,evArray:evArray,dataCube:tdc,$
		stellarData:starsHyperCube,stellarAges:stellarAges,$
		convergence:conv,ncolev:ncolev,ncolstep:ncolstep, $
		NPostProcess:NPostProcess,NPassive:NPassive,STVars:STVars, $
		MLF:MLF,mdotext0:md0,Rf:Rf,whichAccHistory:whichAccHistory,$;}
		fixedQ:fixedQ,diskScaleLength:scaleLength,Qlim:Qlim,dlnx:dlnx,$
		mh0:mh0, minStsig:minStsig,zquench:zquench,xmin:xmin }
	
	RETURN,model ;; end of the readOutput function
END

FUNCTION readmodels,nameList
        model=readOutput(nameList[0])
        modelList=ptrarr(n_elements(nameList), /allocate_heap)
        *(modelList[0])=model
        FOR i=1, n_elements(nameList)-1 DO *(modelList[i])=readOutput(nameList[i])
        RETURN,modelList
END
