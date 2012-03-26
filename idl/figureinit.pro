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
END


PRO figureInit,nm,sv,nx,ny
	IF(sv EQ 1 || sv EQ 0) THEN newWindow,nx,ny
	!p.thick=1.5
	!x.thick=1.5
	!y.thick=1.5
	!z.thick=1.5
	IF(sv EQ 2) THEN BEGIN
		Set_Plot,'ps'
		!p.font = 0
		DEVICE,encapsulated=eps,/helvetica,/isolatin1,bits_per_pixel=8,/color,filename=nm+'.eps',xsize=8.89*nx,ysize=8.89*ny
		
	ENDIF
END