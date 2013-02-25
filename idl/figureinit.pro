PRO newWindow,nx,ny
        window,(!D.WINDOW+1),XSIZE=512*nx,YSIZE=512*ny,RETAIN=2
END

;; Two complementary functions which respectively close 
;; or open a window/png file/eps file for figures.


PRO figureClean,nm,sv
	IF(sv EQ 1) THEN WRITE_PNG,nm+'.png',TVRD(true=1)
	IF(sv EQ 1) THEN wdelete
	IF(sv EQ 2) THEN BEGIN
		DEVICE,/CLOSE
		SET_PLOT,'x'
	ENDIF
	IF(sv EQ 4) THEN BEGIN
		dummy = cgSnapshot(filename=nm,/png,/true,/nodialog)
		SET_PLOT,'x'
	ENDIF
END


PRO figureInit,nm,sv,nx,ny,sm=sm
        ; PRINT, "Initializng figure ",nm,"with sv,nx,ny: ",sv,nx,ny
    IF(n_elements(sm) EQ 0) THEN sm=0
	IF(sv EQ 1 || sv EQ 0) THEN newWindow,nx,ny
	
  	!p.thick=2.5
	!x.thick=2.5
	!y.thick=2.5
	!z.thick=2.5



	IF(sv EQ 4) THEN BEGIN
		set_plot,'z'
		cgDisplay,4096*nx,4096*ny
	ENDIF
	

	IF(sv EQ 2) THEN BEGIN
        pageInfo = PSWINDOW(aspectratio=float(ny)/nx,margin=.15,yfudge=.25)
		Set_Plot,'ps'
        IF(sm EQ 0) THEN $
        DEVICE, _Extra=pageInfo,encapsulated=eps,/helvetica,/isolatin1,bits_per_pixel=8,/color,filename=nm+'.eps',language_level=2
		!p.font = 0
        IF(sm EQ 1) THEN $
		DEVICE,encapsulated=eps,/helvetica,/isolatin1,bits_per_pixel=8,/color,filename=nm+'.eps',xsize=8.89*nx,ysize=8.89*ny
		;DEVICE,encapsulated=eps,/helvetica,/isolatin1,bits_per_pixel=8,/color,filename=nm+'.eps',xoffset=1.0,yoffset=1.0

		
	ENDIF
    !p.background=255
END
