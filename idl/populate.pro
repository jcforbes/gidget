FUNCTION findindex(xbounds,ybounds,nx,ny,x,y)
  IF(xbounds[0] GE xbounds[1]) THEN message,"findindex in populate.pro: Invalid xbounds"
  IF(ybounds[0] GE ybounds[1]) THEN message,"findindex in populate.pro: Invalid ybounds"

  if(x LT xbounds[0]) THEN xv=xbounds[0]
  if(x GT xbounds[1]) THEN xv=xbounds[1]
  if(y LT ybounds[0]) THEN yv=ybounds[0]
  if(y GT ybounds[1]) THEN yv=ybounds[1]
  i = FIX( (nx-1.0) * log(xv/xbound[0]) / log(xbound[1]/xbound[0]) )
  j = FIX( (ny-1.0) * log(yv/ybound[0]) / log(ybound[1]/ybound[0]) )

  RETURN,[i,j]
END


;; Generate a grid (e.g. M* vs SFR) and populate it with a ton of models
PRO populate,expName,N,sv
  compile_opt idl2
  set_plot,'x'
  DEVICE,DECOMPOSED=0
  DEVICE,RETAIN=2

  nameList = ListValidModels(expName,N)

  modelPtr = ReadModels([nameList[0]])
  model = (*(modelPtr[0]))

  ctr = 0
  time = model.evArray[1,*]

  xbounds = [1.0d7,1.0d11] ; M* : MSol
  ybounds = [1.0d-5,1.0d2] ; sfr: MSol/yr
  histogram = dblarr(Nx,Ny)

  FOR i=0,n_elements(nameList)-1 DO BEGIN ;; loop over every valid model
    IF(i MOD 50 EQ 0) THEN PRINT,"Processing model #",i," of ",n_elements(nameList)
    IF(i NE 0) THEN modelPtr=ReadModels([nameList[i]])
    IF(i NE 0) THEN model = (*(modelPtr[0]))

    modelinfo = diskstats(model,0.0) ;; get a bit of info about the model at z=0
    ind = findindex(xbounds,ybounds,Nx,Ny,modelinfo[14],modelinfo[13])
    
    histogram[ind[0],ind[1]] += weight 

  ENDFOR ; end loop over every valid model

  x=dblarr(Nx)
  y=dblarr(Ny)
  FOR i=0,Nx-1 DO x[i]=xbounds[0]*(xbounds[1]/xbounds[0])^(float(i)/(Nx-1.0))
  FOR j=0,Ny-1 DO y[j]=ybounds[0]*(ybounds[1]/ybounds[1])^(float(j)/(Ny-1.0))

  CONTOUR,histogram,x,y,/xlog,/ylog

END
