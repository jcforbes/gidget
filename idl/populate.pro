FUNCTION findindex,xbounds,ybounds,nx,ny,x,y
  IF(xbounds[0] GE xbounds[1]) THEN message,"findindex in populate.pro: Invalid xbounds"
  IF(ybounds[0] GE ybounds[1]) THEN message,"findindex in populate.pro: Invalid ybounds"

  xv = x
  yv = y

  if(x LT xbounds[0]) THEN xv=xbounds[0]*1.01
  if(x GT xbounds[1]) THEN xv=xbounds[1]*.99
  if(y LT ybounds[0]) THEN yv=ybounds[0]*1.01
  if(y GT ybounds[1]) THEN yv=ybounds[1]*.99
  i = FIX( (nx-1.0) * alog10(xv/xbounds[0]) / alog10(xbounds[1]/xbounds[0])  + 0.5)
  j = FIX( (ny-1.0) * alog10(yv/ybounds[0]) / alog10(ybounds[1]/ybounds[0])  + 0.5)

  RETURN,[i,j]
END

FUNCTION weight,weights,mh0
  index = WHERE(ABS(weights[0,*]-alog10(mh0)) LT .00001, count)
  if(count EQ 1) THEN RETURN,weights[1,index]
  MESSAGE,("Could not find the weight for Mh0 = "+string(mh0))
  RETURN,0.0 ;; never get here
END

FUNCTION GetAxis, ybounds,Ny
  y = dblarr(Ny)
  FOR j=0,Ny-1 DO y[j]=ybounds[0]*(ybounds[1]/ybounds[0])^(float(j)/(Ny-1.0))
  RETURN,y
END

FUNCTION ReadWeights,filename
  OPENR, WL, 'nM_z0.dat', /GET_LUN
  logmh12 = 0d
  theWeight = 0d
  weights = dblarr(2,100000) ; a very large array
  ctr=0
  WHILE (NOT EOF(WL)) DO BEGIN
    READF, WL, logmh12, theWeight
    weights[0,ctr] = logmh12+12.0
    weights[1,ctr] = theWeight
    ctr=ctr+1
  ENDWHILE
  CLOSE,WL
  RETURN,weights[*,0:ctr-1]
END

FUNCTION ReadFromDatFile,filename,N
  OPENR, lun, filename, /GET_LUN
  data = dblarr(13,N)
  READF, lun, data
  FREE_LUN, lun
  RETURN, data
END

FUNCTION ReadFromRawFiles,filename,namelist
  data = dblarr(13,n_elements(namelist))
  FOR i=0,n_elements(namelist)-1 DO BEGIN ;; loop over every valid model
    IF(i MOD 50 EQ 0) THEN PRINT,"Reading in model #",i," of ",n_elements(namelist)
    modelPtr=ReadModels([nameList[i]])
    model = (*(modelPtr[0]))
    modelinfo = diskstats(model,z=0.0) ;; get a bit of info about the model at z=0
    data[*,i] = [model.mh0,modelinfo[14],modelinfo[13],modelinfo[16],modelinfo[15],modelinfo[10],modelinfo[17],modelinfo[0],modelinfo[1],modelinfo[2],modelinfo[3],modelinfo[8],modelinfo[9]]
    ptr_free,modelPtr[0]
  ENDFOR ; end loop over every valid model

  OPENW,lun, filename, /get_lun
  PRINTF,lun,data
  FREE_LUN,lun

  RETURN, data
END

;; Generate a grid (e.g. M* vs SFR) and populate it with a ton of models
PRO populate,expName,N,sv
  compile_opt idl2
  set_plot,'x'
  DEVICE,DECOMPOSED=0
  DEVICE,RETAIN=2
  SETCT,1,0,0

  nameList = ListValidModels(expName,N)

  expName2 = expName+'_'+strcompress(string(MIN([n_elements(nameList),N])),/remove) ;; e.g. rk5_113
  
  xbounds = [1.0d5,1.0d11] ; M* : MSol
  ySFRbounds = [1.0d-5,1.0d2] ; sfr: MSol/yr
  yZbounds = [.1,3] ;; Z/Zsolar (linear, not log)
  yFgbounds = [.01,1.5] ;; gas fraction
  yfH2bounds = [.01,1.5] ;; fH2 fraction
  yEffbounds = [1.0e-3,1.0e-1] ;; efficiency of turning baryons into stars

  Nx=500
  Ny=500

  exists=FILE_TEST(expname2+"_0d.dat")
  IF(exists) THEN theData= ReadFromDatFile(expname2+"_0d.dat",n_elements(namelist))$
   ELSE theData = ReadFromRawFiles(expname2+"_0d.dat",namelist)
  
  

  histo_mst_sfr = dblarr(Nx,Ny)+1.0
  histo_mst_Z = dblarr(Nx,Ny)+1.0
  histo_mst_fg = dblarr(Nx,Ny)+1.0
  histo_mst_fH2 = dblarr(Nx,Ny)+1.0
  histo_mst_eff = dblarr(Nx,Ny)+1.0


  weights= ReadWeights('nM_z0.dat')
;;  weights *= 1.0 / MIN(weights) ;; ensure that all weights are >= 1 (probably not necessary)


  FOR i=0,n_elements(nameList)-1 DO BEGIN ;; loop over every valid model
    ;; M*, SFR, Z, Fg, fH2, efficiency
    indSFR = findindex(xbounds,ySFRbounds,Nx,Ny,theData[1,i],theData[2,i])
    indZ = findindex(xbounds,yZbounds,Nx,Ny,theData[1,i],theData[3,i])
    indFg = findindex(xbounds,yFgbounds,Nx,Ny,theData[1,i],theData[4,i])
    indfH2 = findindex(xbounds,yfH2bounds,Nx,Ny,theData[1,i],theData[5,i])
    indEff = findindex(xbounds,yEffbounds,Nx,Ny,theData[1,i],theData[6,i])

    ;; weight = dN/dMh12 or 1 or 1/(# of successful runs in this Mh0 bin)?
    theWeight = weight(weights,theData[0,i])
    histo_mst_sfr[indSFR[0],indSFR[1]] += theWeight ; weight(weights,model.mh0)
    histo_mst_Z[indZ[0],indZ[1]] += theWeight
    histo_mst_fg[indFg[0],indFg[1]] += theWeight
    histo_mst_fH2[indfH2[0],indfH2[1]] += theWeight
    histo_mst_eff[indEff[0],indEff[1]] += theWeight

  ENDFOR ; end loop over every valid model

  x=GetAxis(xbounds,Nx);dblarr(Nx)
  NL=10
;  window,0
  FIGUREInit,(expName2+'_Mst_SFR'),sv,1,1
  CONTOUR,alog10(histo_mst_sfr),x,GetAxis(ySFRbounds,Ny),/xlog,/ylog,COLOR=0,BACKGROUND=255,XTITLE="M*",YTITLE="SFR",nlevels=NL,XSTYLE=1,YSTYLE=1
  FigureClean,(expName2+'_Mst_SFR'),sv

  FIGUREINIT,(expName2+'_Mst_Z'),sv,1,1
  CONTOUR,alog10(histo_mst_Z),x,GetAxis(yZbounds,Ny),/xlog,/ylog,COLOR=0,BACKGROUND=255,XTITLE="M*",YTITLE="Z/Zsol",nlevels=NL,XSTYLE=1,YSTYLE=1
  FigureClean,(expName2+'_Mst_Z'),sv

  FigureInit,(expName2+'_Mst_fg'),sv,1,1
  CONTOUR,alog10(histo_mst_fg),x,GetAxis(yFgbounds,Ny),/xlog,/ylog,COLOR=0,BACKGROUND=255,XTITLE="M*",YTITLE="f_g",nlevels=NL,XSTYLE=1,YSTYLE=1
  FigureClean,(expName2+'_Mst_fg'),sv

  FigureInit,(expName2+'_Mst_fH2'),sv,1,1
  CONTOUR,alog10(histo_mst_fH2),x,GetAxis(yfH2bounds,Ny),/xlog,/ylog,COLOR=0,BACKGROUND=255,XTITLE="M*",YTITLE="f_H2",nlevels=NL,XSTYLE=1,YSTYLE=1
  FigureClean,(expName2+'_Mst_fH2'),sv

  FigureInit,(expName2+'_Mst_eff'),sv,1,1
  CONTOUR,alog10(histo_mst_eff),x,GetAxis(yEffbounds,Ny),/xlog,/ylog,COLOR=0,BACKGROUND=255,XTITLE="M*",YTITLE="SF eff",nlevels=NL,XSTYLE=1,YSTYLE=1
  FigureClean,(expName2+'_Mst_eff'),sv


END


