PRO newWindow,nx,ny
	window,(!D.WINDOW+1),XSIZE=512*nx,YSIZE=512*ny,RETAIN=2
END

;; Choose which color table to use.
;; which=1 <=> each color is quite distinct from the others
;; which=2 <=> a gradual progression from black to red to blue
;; which=3 <=> a gradual progression within a certain set of colors, determined by "color"
;;                  n refers to the number of gradations in this color
PRO setct,which,n,color
	IF(which EQ 0) THEN BEGIN
		compile_opt idl2
		set_plot,'x'
		DEVICE,DECOMPOSED=0
		DEVICE,RETAIN=2
	ENDIF
	
	;;; Construct a color table
	r=BYTE(INDGEN(256))
	g=r
	b=r ;; at this point, everything is gray

	nn=fix(n)
	IF(nn GT 255) THEN nn=254

	IF(which EQ 2) THEN BEGIN


		r[0:13]=byte([  0,255,215,185,150,125,100, 50, 50,  0,  0,  0,  0,  0])
		g[0:13]=byte([  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0])
		b[0:13]=byte([  0,  0,  0,  0,  0,  0,  0,  0, 50,  0,100,150,200,255])
	
		TVLCT,r,g,b	
	ENDIF
	IF(which EQ 3) THEN BEGIN
;;		IF(n GT 256) THEN n=256
		;;; 0 -> blue, 1->red
		IF(color EQ 0 || color EQ 3 || color EQ 4 || color EQ 6) THEN wb=1 ELSE wb=0
		IF(color EQ 1 || color EQ 3 || color EQ 5 || color EQ 6) THEN wr=1 ELSE wr=0
		IF(color EQ 2 || color GE 4) THEN wg=1 ELSE wg=0
		r[0:nn-1]=byte((indgen(nn)+1)*fix(256.0/float(nn+1))*wr)
		g[0:nn-1]=byte((indgen(nn)+1)*fix(256.0/float(nn+1))*wg)
		b[0:nn-1]=byte((indgen(nn)+1)*fix(256.0/float(nn+1))*wb)
	ENDIF
	IF(which EQ 1) THEN BEGIN
		r[0:10]=byte([  0,255,  0,255,175,  6,116,  0,  7,255,255]) ; black,red,blue,orange,purple,cyan,maroon,dark green,yellow,magenta
		g[0:10]=byte([  0,  0,  0,125,  2,192,  5,255,112,255,  9])
		b[0:10]=byte([  0,  0,255,  0,175,175,  5,  0,  7,  0,178])
	
	ENDIF

	IF(which EQ 4) THEN BEGIN
		r=BYTE(RANDOMU(color,256)*256)
		g=BYTE(RANDOMU(color+1,256)*256)
		b=BYTE(RANDOMU(color+2,256)*256)
	ENDIF

	IF(which EQ 5) THEN BEGIN
		;; n_elements(styles),j
		;; rr,gg,bb=50
		;; rr->256
		;; gg->256
		;; rr->50
		;; bb->256
		;; gg->50
		;; rr,gg->256
		;; rr,gg,bb->50
		
		mmm = 70
		sp = 255-mmm
;		j=color MOD (7*sp)
		j = ( fix(float(color * 6.0 * sp) / float(n)  ) ) MOD (6*sp)
		
		ct = (j) MOD sp
		rr=mmm
		gg=mmm
		bb=mmm
		IF(j LE sp) THEN rr=mmm+ct
		IF(j GT sp AND j LE 2*sp) THEN BEGIN
			rr=255
			gg=mmm+ct
		ENDIF
		IF(j GT 2*sp AND j LE 3*sp) THEN BEGIN
			gg=255
			rr=255 - ct
		ENDIF
		IF(j GT 3*sp AND j LE 4*sp) THEN BEGIN
			bb=mmm+ct
			gg=255
		ENDIF
		IF(j GT 4*sp AND j LE 5*sp) THEN BEGIN
			bb=255
			gg=255-ct
		ENDIF
		IF(j GT 5*sp AND j LE 6*sp) THEN BEGIN
			bb=255
;			gg=mmm+ct
			rr=mmm+ct
		ENDIF
		IF(j GT 6*sp AND j LE 7*sp) THEN BEGIN
			bb=255 -ct
			gg=mmm +ct
			rr=255 
		ENDIF
		ind=j MOD 256
		r[*]= rr
		g[*]= gg
		b[*]= bb
		;stop
	ENDIF

	IF(which NE 5) THEN BEGIN
		r[0]=byte(0)
		g[0]=byte(0)
		b[0]=byte(0)
		r[255]=byte(255)
		g[255]=byte(255)
		b[255]=byte(255)
	ENDIF

;	stop,which
	TVLCT,r,g,b
END
