FUNCTION OP,j
  IF(j EQ 0) THEN ret=0 ELSE ret=1
  RETURN,ret
END

PRO MakeComparisonMovie,modelList,Xind,ind,npx,npy,fname,xranges,ranges,sv,cmt,leg
  !p.multi=[0,npx,npy]
  IF(sv EQ 3) THEN cg=1 ELSE cg=0
  IF(cg EQ 1) THEN cs = 1 ELSE cs=2
  nxwin=1024
  nywin=nxwin

  fn=""
  IF(n_elements(modelList) LT 6) THEN $
    FOR i=0,n_elements(modelList)-1 DO $
      fn+= (ExtractSlash((modelList[i]),1) + "_") $
   ELSE fn+= ExtractSlash((modelList[0]),0)+"_"
  fn+=fname
  dn="movie_"+fn
  set_plot,'x'
  if(sv EQ 3 ) THEN set_plot,'z' 
  if(cg NE 1) THEN WINDOW,1,XSIZE=nxwin,YSIZE=nywin,RETAIN=2
  if(sv EQ 3) THEN cgDisplay,nxwin,nywin
  IF(sv EQ 1) THEN mpeg_id=MPEG_OPEN([nxwin,nywin],FILENAME=(fn+".mpg"))
  IF (sv EQ 2 || sv EQ 3) THEN BEGIN
    FILE_MKDIR,dn
  ENDIF
  count=0

  maxTime=0.
  whichmodel=0
  nts=1
  FOR i=0,n_elements(modelList)-1 DO BEGIN
    maxTime=MAX([maxTime,MAX((*(modelList[i])).evArray[1,*])],replaced)
    IF(replaced EQ 1) THEN BEGIN
      whichmodel=i
      nts=n_elements( (*(modelList[i])).evArray[1,*])
    ENDIF
  ENDFOR
  FOR n=1,nts DO BEGIN
    theTime=(*(modelList[whichmodel])).evArray[1,n-1]
    theRedshift=(*(modelList[whichmodel])).evArray[9,n-1]
    ls=intarr(n_elements(modelList))                
    ti=intarr(n_elements(modelList))
    FOR i=0,n_elements(modelList)-1 DO BEGIN
      nts2 = n_elements((*(modelList[i])).evArray[1,*])
;      ti[i] = MIN([n,nts])-1
      ti[i] = nearest("time",1,(*(modelList[i])),theTime)
      IF(nts2-1 EQ ti[i])THEN ls[i] = 1 ;; make the line dotted if that run has stopped
    ENDFOR

    count=count+1
    ;; make the plots for this frame
    ns=0.
    ns2=0.
    FOR i=1,n_elements(ind) DO BEGIN
      IF(ind[i-1] GT 0) THEN BEGIN
                               
        FOR j=0,n_elements(modelList)-1 DO PlotVsXind,*(modelList[j]),ti[j],Xind[i-1],ind[i-1],OP(j),1,XRANGE=xranges[*,i-1], $
          YRANGE=ranges[*,i-1],CHARSIZE=cs,COLOR=j,LINESTYLE=ls[j];;,THICK=3,YSTYLE=1
        FOR j=0,n_elements(modelList)-1 DO PlotVsXind,*(modelList[j]),0,Xind[i-1],ind[i-1],1,1,COLOR=j,LINESTYLE=2;;,THICK=3
        IF i EQ 1 THEN BEGIN
          XYOUTS, .1,2./float(npy)-.008,/NORMAL,("z = "+strcompress(string(theRedshift),/remove)+" "+ CMT),COLOR=0,CHARSIZE=.8*cs,WIDTH=ns
          XYOUTS, .1,1./float(npy)-.008,/NORMAL,("T = "+strcompress(string(theTime),/remove)+" Ga "+ CMT),COLOR=0,CHARSIZE=.8*cs,WIDTH=ns
          ns2=.1+ns
          FOR k=0,n_elements(leg)-1 DO BEGIN
            XYOUTS, ns2,1./float(npy)-.008,/NORMAL," "+leg[k],COLOR=k,CHARSIZE=.7*cs,WIDTH=ns
            ns2=ns2+ns
          ENDFOR
        ENDIF
      ENDIF ELSE BEGIN
        IF(ind[i-1] EQ -1) THEN dQdqContour,*(modelList[0]),ti[0]
        IF(ind[i-1] EQ -2) THEN dQdqContour,*(modelList[1]),ti[1]
        IF(ind[i-1] EQ -3) THEN Raf3,*(modelList[0]),ti[0]
        IF(ind[i-1] EQ -4) THEN Raf3,*(modelList[1]),ti[1]
      END
    ENDFOR
    IF(sv EQ 1) THEN MPEG_PUT,mpeg_id,WINDOW=1,FRAME=count,/ORDER
    IF(sv EQ 2) THEN WRITE_PNG, dn+'/frame_' + STRING(count, FORMAT="(I4.4)") + '.png', TVRD(TRUE=1)
    IF(sv EQ 3) THEN dummy=cgSnapshot(filename=(dn+'/frame_'+STRING(count,FORMAT="(I4.4)")),/png,/true,/nodialog)
    IF(sv EQ 0) THEN wait,.05
  ENDFOR
  IF(sv EQ 1) THEN MPEG_SAVE, mpeg_id
  IF(sv EQ 1) THEN MPEG_CLOSE, mpeg_id
  IF(sv EQ 2 || sv EQ 3) THEN spawn,("rm "+ fn+".mpg")
  IF(sv EQ 2 || sv EQ 3) THEN spawn,("ffmpeg -f image2 -qscale 4 -i "+dn+"/frame_%04d.png " + fn + ".mpg")
  IF(sv EQ 2 || sv EQ 3) THEN spawn,("rm -rf " + dn)
        
  !p.multi=[0,1,1] ; clean up
  set_plot,'x'
END


PRO diskAnalyzeCompare,modelList,sv,cmt,leg
	IF(n_elements(modelList) EQ 0 || n_elements(sv) EQ 0 || n_elements(cmt) EQ 0 || n_elements(leg) EQ 0) THEN BEGIN
		message,("Usage: diskAnalyzeCompare,['<name1>','<name2>',...],sv,comment,['<legend1>','<legend2>',...] "+ $
                  "-- sv=0->don't save sv=1->save as mpg sv=2->save png images to be constructed into a movie by ffmpeg,"+ $
                  " sv=3->same as 2 but don't plot in an X window.")
	ENDIF

	SETCT,0,0,0
	SETCT,1,0,0
	FOR j=0, n_elements(modelList)-1 DO PRINT,"Model ",j," ends at redshift ",$
		(*(modelList[j])).evArray[9,n_elements((*(modelList[j])).evArray[9,*])-1]," with ", $
		n_elements((*(modelList[j])).evArray[0,*])," outputs."

	ind=[2,3,43,4,5,24,6,7,23]  ;;;tau,tau',tau'',S,s,Qg,S*,s*,Qs
	Xind=ind*0+1
	xr = getranges5(modelList,Xind)
	xr[0,*]*=.5
        xr[1,*]*=2
	MakeComparisonMovie,modelList,(ind*0+1),ind,3,3,"stt",xr,getranges5(modelList,ind),sv,CMT,leg

	ind=[62,63,64,65,66,67,68,69]
	Xind=ind*0+61
	MakeComparisonMovie,modelList,Xind,ind,2,4,"dim",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg
	
	ncs=(*(modelList[0])).ncolstep
	npp=(*(modelList[0])).npostprocess
	stv=(*(modelList[0])).STVars
	np=(*(modelList[0])).NPassive
	ind=[8,ncs+1,9,ncs+2,10,ncs+3,11,ncs+4,21,29] ;; Sdot,sdot,S*dot,s*dot,dZ/dt,SFR
	Xind=ind*0+1
;	MakeComparisonMovie,modelList,(ind*0+1),ind,2,5,"ddt",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg

	ind=[45,60,18,47,17,57,48,22,37,35,58,59] ;; Toomre Length, Toomre Mass,q,alpha, f_g, ratio of sig/sig_st, fH2, Z, vrg, yy
	Xind=ind*0+61
	MakeComparisonMovie,modelList,Xind,ind,3,4,"phv",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg

	ind=[50,25,26,27,  36,12,31,32,  13,14,15,16] ;; CumulativeTorqueErr,Q_R,Q_WS,Q_RW,  torqueErr,Q,dQdS,dQds, h0,h1,h2,H
	xind=ind*0+1
;	MakeComparisonMovie,modelList,(ind*0+1),ind,4,3,"dg",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg

;	ind=[41,ncs+21,45,ncs+23,   42,ncs+20,46, ncs+22 ]
	ind=[41,42,   ncs+21,ncs+20,   45,46,   ncs+23,ncs+22]
	xind=ind*0+61
;	MakeComparisonMovie,modelList,Xind,ind,2,4,"massBudget",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg


        ind = [ncs+33, ncs+34,12,44,   8,9,43,11,  13,14,15,16] ;; -tau, -tau', Q, sig',   Sdot, sdot, tau'', s*dot,    h0,h1,h2,H
        xind=ind*0+1
        xr=getranges5(modelList,Xind)
	xr[0,*]*=.5
	xr[1,*]*=2
        MakeComparisonMovie,modelList,Xind,ind,4,3,"dg",xr,getranges5(modelList,ind),sv,CMT,leg

	ind=[ncs+npp+stv*0 + 1, ncs+npp+stv*1 + 1, ncs+npp+stv*3 + 1, ncs+npp+stv*np + 1,$
		ncs+npp+stv*0 +2,ncs+npp+stv*1+ 2, ncs+npp+stv*3 + 2, ncs+npp+stv*np + 2,$
		ncs+npp+stv*0 +3,ncs+npp+stv*1+ 3, ncs+npp+stv*3 + 3, ncs+npp+stv*np + 3]
	xind=ind*0+61
	ranges=dblarr(2,12)
	ranges[0,0:3]=.1
	ranges[1,0:3]=1000
	ranges[0,4:7]=10
	ranges[1,4:7]=100
	ranges[0,8:11]=-1.5
	ranges[1,8:11]=.5
	MakeComparisonMovie,modelList,Xind,ind,4,3,"stars",getranges5(modelList,Xind),getranges5(modelList,ind),sv,CMT,leg

;;	ind=[4,6,17]
;;	xind=ind*0+1
;;	MakeComparisonMovie,modelList,(ind*0+1),ind,1,3,"ColDen",[[.05,1],[.05,1],[.05,1]],[[1,10000],[1,10000],[0,.6]],sv,CMT,leg

;;	ind=[39,41,53]
;;	xind=ind*0+1
;;	MakeComparisonMovie,modelList,(ind*0+1),ind,1,3,"VelDis",[[.05,1],[.05,1],[.05,1]],[[5,220],[5,220],[0,1]],sv,CMT,leg

	

;	makeplots,cmt,sv,modelList=modelList,legend=leg

	FOR j=0, n_elements(modelList)-1 DO ptr_free,modelList[j]
END
