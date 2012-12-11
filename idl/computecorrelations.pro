
;; Given an array of mdot's vs. other variables, compute the cross-correlations between mdot and the others
;; and also the autocorrelation of each variable.
PRO ComputeCorrelations,vsmdot,colors,time,labels,names,name,sv=sv,nt0=nt0,thicknesses=thicknesses,logarithms=logarithms
    IF(n_elements(nt0) EQ 0) THEN nt0 = 1
    IF(n_elements(sv) EQ 0 ) THEN sv=4
    IF(n_elements(thicknesses) EQ 0) THEN thicknesses=intarr(n_elements(vsMdot[0,0,0,*]))+1 ; set each model's thickness to 1
    IF(n_elements(logarithms) EQ 0) THEN logarithms=intarr(n_elements(vsMdot[0,0,*,0])) ; set each variable to be not logarithmic
    cs = 2
    ct = 1
    IF(sv EQ 4) THEN cs = 3
    IF(sv EQ 4) THEN ct = 3
    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)*3
    ls=0
    IF(n_elements(vsMdot[0,0,0,*]) GT 50) THEN ls =1 

    ;; Construct some cross-correlations
    nt = n_elements(vsMdot[*,0,0,0])
    ;;; (# of delay times to test) x (# of cross-correlations) x (# of models)
    crossCorrelations = dblarr(nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    autoCorrelations = dblarr(nt-1-nt0,n_elements(vsMdot[0,0,*,0]),n_elements(vsMdot[0,0,0,*]))
    ;;; loop over models
    FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
        ;; loop over variables
        FOR k=0, n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
            laggingVar = vsMdot[nt0:(nt-1),0,k,j]
            IF(logarithms[k] EQ 1) THEN laggingVar = alog10(laggingVar)
            mdot = vsMdot[nt0:(nt-1),0,0,j]
            IF(logarithms[0] EQ 1) THEN mdot = alog10(mdot)
            crossCorrelations[*,k,j] = c_correlate(mdot,laggingVar,indgen(nt-1-nt0))
            autoCorrelations[*,k,j] = c_correlate(laggingVar,laggingVar,indgen(nt-1-nt0))
        ENDFOR

    ENDFOR

    ;; loop over variables
    FOR k=0,n_elements(vsMdot[0,0,*,0])-1 DO BEGIN
        setct,1,0,2
        FIGUREINIT,(name+"_cc_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="Lag Time (Ga)",YTITLE="x-corr: "+labels[0]+" with "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
            OPLOT, time[1:(nt-1-nt0)],crossCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
        ENDFOR
        FIGURECLEAN,(name+"_cc_"+names[k]),sv

        FIGUREINIT,(name+"_ac_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[.8*MIN(time[2:nt-2]),MAX(time)],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="Lag Time (Ga)",YTITLE="autocorr: "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
             OPLOT, time[1:(nt-1-nt0)],autoCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
        ENDFOR
        FIGURECLEAN,(name+"_ac_"+names[k]),sv

        FIGUREINIT,(name+"_xcc_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[-1.0,1.0],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="mdot autocorrelation",YTITLE="x-corr: "+labels[0]+" with "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
            OPLOT, crossCorrelations[*,0,j] , crossCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
        ENDFOR
        FIGURECLEAN,(name+"_xcc_"+names[k]),sv

        FIGUREINIT,(name+"_xac_"+names[k]),sv,2,2
        PLOT,[0],[0], COLOR=0,BACKGROUND=255, XRANGE=[-1.0,1.0],YRANGE=[-1.0,1.0],XSTYLE=1,YSTYLE=1,THICK=1,XTITLE="mdot autocorrelation",YTITLE="autocorr: "+labels[k],XLOG=1,CHARSIZE=cs,CHARTHICK=ct
        FOR j=0, n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
            OPLOT, crossCorrelations[*,0,j] , autoCorrelations[*,k,j],COLOR=colors[0,j],THICK=thicknesses[j],LINESTYLE=ls
        ENDFOR
        FIGURECLEAN,(name+"_xac_"+names[k]),sv


    ENDFOR
    IF(sv EQ 4) THEN thicknesses = temporary(thicknesses)/3.0
END


