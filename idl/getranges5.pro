FUNCTION GetRanges,model,ind
        ni=n_elements(ind)
        ranges=dblarr(2,ni)
        FOR i=0,ni-1 DO ranges[*,i]=[MIN(model.dataCube[*,*,ind[i]-1]),MAX(model.dataCube[*,*,ind[i]-1])]
        RETURN, ranges
END

FUNCTION GetRanges2,model1,model2,ind
        ni=n_elements(ind)
        ranges=dblarr(2,ni)
        FOR i=0,ni-1 DO BEGIN
                IF(ind[i] GE 0) THEN ranges[*,i]=[MIN([MIN(model1.dataCube[*,*,ind[i]-1]),MIN(model2.dataCube[*,*,ind[i]-1])]),  MAX([MAX(model1.dataCube[*,*,ind[i]-1]),MAX(model2.dataCube[*,*,ind[i]-1])])] ELSE ranges[*,i] = [0.,1.]
        ENDFOR
        RETURN, ranges
END

FUNCTION GetRanges4,model1,model2,ind
        ni=n_elements(ind)
        ranges=dblarr(2,ni)
        ranges2=GetRanges2(model1,model2,ind)
        FOR i=0,ni-1 DO BEGIN
                IF(ind[i] GE 0) THEN BEGIN
                        md1=0d
                        RESISTANT_Mean,(model1.dataCube[*,*,ind[i]-1]),5.,md1
                        md2=0d
                        RESISTANT_Mean,(model2.dataCube[*,*,ind[i]-1]),5.,md2

;                       md1=median(model1.dataCube[*,*,ind[i]-1])
;                       md2=median(model2.dataCube[*,*,ind[i]-1])
                        sg1= robust_sigma(model1.dataCube[*,*,ind[i]-1])
                        sg2= robust_sigma(model2.dataCube[*,*,ind[i]-1])
                        mult=8
                        ranges[*,i]=[MIN([md1-mult*sg1,md2-mult*sg2]),MAX([md1+mult*sg1,md2+mult*sg2])]
                        IF(ranges[1,i]-ranges[0,i] NE 0.) THEN rat = (ranges2[1,i]-ranges2[0,i])/(ranges[1,i]-ranges[0,i]) ELSE rat=0. ;; if there is no variance, then use the larger range.
                        IF(rat LT 15) THEN ranges[*,i]=ranges2[*,i]
                        IF(logvar(ind[i],model1) EQ 1) THEN ranges[*,i] = ranges2[*,i]
                ENDIF ELSE ranges[*,i]=[0.,1.]

        ENDFOR
        RETURN,ranges
END
FUNCTION GetRanges3,model1,model2,ind
        ni=n_elements(ind)
        m1t=n_elements(model1.dataCube[*,0,0])
        m2t=n_elements(model2.dataCube[*,0,0])
        ranges=dblarr(2,ni)
        ranges2=GetRanges2(model1,model2,ind)
        FOR i=0,ni-1 DO BEGIN
                a= MIN([Min(model1.dataCube[m1t-1,*,ind[i]-1]),Min(model2.dataCube[m2t-1,*,ind[i]-1])]) 
                b= MIN([Min(model1.dataCube[*,*,ind[i]-1]),Min(model2.dataCube[*,*,ind[i]-1])])
;               b= MIN([Min(model1.dataCube[0,*,ind[i]-1]),Min(model2.dataCube[0,*,ind[i]-1])])

                IF(b GT a) THEN BEGIN
                        ;; swap a and b so that a is always the least larger value (usually this will mean the less extreme value, since the minima will be negative)
                        tempA = a
                        a=b
                        b=tempA
                ENDIF
                IF(a*b GE 0.) THEN lower=(abs(a)^.8)*(abs(b)^.2)*(FIX(a GT 0) - FIX(a LT 0)) ELSE lower = FIX(a LT 0)*a + FIX(b LT 0)*b
                c= MAX([Max(model1.dataCube[m1t-1,*,ind[i]-1]),Max(model2.dataCube[m2t-1,*,ind[i]-1])]) 
                d= MAX([Max(model1.dataCube[*,*,ind[i]-1]),Max(model2.dataCube[*,*,ind[i]-1])])
;               d= MAX([Max(model1.dataCube[0,*,ind[i]-1]),Max(model2.dataCube[0,*,ind[i]-1])])

                IF (c GT d) THEN BEGIN
                                tempC= c
                                c=d
                                d=tempC
                ENDIF
                IF(c*d GE 0.) THEN upper=(abs(c)^.8)*(abs(d)^.2)*(FIX(c GT 0) - FIX(c LT 0)) ELSE upper = FIX(c GT 0)*c + FIX(d GT 0)*d

;               IF lower LT 0 THEN lower*=2.
;               IF lower GT 0 THEN lower=lower/2.
;               IF upper LT 0 THEN upper=upper/2.
;               IF upper GT 0 THEN upper=upper*2.
                ranges[*,i]=[lower,upper]
                
                IF ((c-a)/(d-b) LT 10. && (c-a)/(d-b) GT (1./10.)) THEN BEGIN
                        ranges[*,i]=ranges2[*,i]
                ENDIF

        ENDFOR
                

        RETURN,ranges
END


FUNCTION getranges5,modelList,ind
        ni=n_elements(ind)
        nm=n_elements(modelList)
        ranges=dblarr(2,ni)
        ranges2=dblarr(2,ni)


        FOR i=0,ni-1 DO BEGIN
	        log = logvar(ind[i],(*(modelList[0])))
                IF(ind[i] GE 0) THEN BEGIN
                        max=-1.e30
                        min=1.e30
			;; take the max so far and compare it to the MAX from this model
                        FOR j=0,nm-1 DO max=MAX([max,MAX((*(modelList[j])).dataCube[*,*,ind[i]-1])])
			;; take the min so far and compare it to the MIN from this model
                        IF(log EQ 0) THEN FOR j=0,nm-1 DO min=MIN([min,MIN((*(modelList[j])).dataCube[*,*,ind[i]-1])])
			IF(log EQ 1) THEN BEGIN
				FOR j=0,nm-1 DO BEGIN 
                                        arr = (*(modelList[j])).dataCube[*,*,ind[i]-1]
					valInd = WHERE(arr GT 0.0, ct)
					IF(ct GT 0) THEN min=MIN([min,MIN(arr[valInd])])
;					minTry=MIN([min,MIN((*(modelList[j])).dataCube[*,*,ind[i]-1])])
;					IF(minTry GT 0.0) THEN min=minTry
				ENDFOR
			ENDIF
			IF(log EQ 1) THEN IF(max LE 0) THEN max=1.0
			IF(log EQ 1 AND min GE max) THEN min=.001*max
			IF(log EQ 1 and max/min GT 1000000.0) THEN min=max/1000000.0
			ranges2[*,i] = [min,max]
                END ELSE ranges2[*,i]=[0.,1.]
        ENDFOR
        
;	stop,ind
	;;; just return the min and max...
	RETURN,ranges2

        FOR i=0,ni-1 DO BEGIN
                IF(ind[i] GE 0) THEN BEGIN
                        md=dblarr(nm)
                        sg=dblarr(nm)
                        tempMd=0.
                        FOR j=0,nm-1 DO BEGIN 
                                RESISTANT_Mean,((*(modelList[j])).dataCube[*,*,ind[i]-1]),5.,tempMd
                                md[j]=tempMD
                                sg[j]=robust_sigma((*(modelList[j])).dataCube[*,*,ind[i]-1])
                        ENDFOR
                        mult=4.5
                        maximum=MAX(md+mult*sg)
                        minimum=MIN(md-mult*sg)
                        ranges[*,i]=[minimum,maximum]
                        IF(ranges[1,i]-ranges[0,i] NE 0.) THEN rat = (ranges2[1,i]-ranges2[0,i])/(ranges[1,i]-ranges[0,i]) ELSE rat=0. ;; if there is no variance then use the larger range.
                        IF(rat LT 4) THEN ranges[*,i]=ranges2[*,i] 
                        IF(rat GE 4 && ranges2[1,i] LT ranges[1,i] ) THEN ranges[1,i]=ranges2[1,i]
                        IF(rat GE 4 && ranges2[0,i] GT ranges[0,i] ) THEN ranges[0,i]=ranges2[0,i]
                        IF(logvar(ind[i],*(modelList[0])) EQ 1) THEN ranges[*,i]=ranges2[*,i]
                ENDIF ELSE ranges[*,i]=[0.,1.]
        ENDFOR

        ;; don't do anything fancy,actually
        ranges=ranges2

        FOR i=0,ni-1 DO BEGIN
                IF(logvar(ind[i],*(modelList[0])) EQ 1 && ranges[0,i] LE 0) THEN BEGIN
                        min = 1.d50
                        FOR j=0,nm-1 DO BEGIN
                                a= (*(modelList[j])).dataCube[*,*,ind[i]-1]
                                IF(n_elements(where(a GT 0)) NE 1) THEN min=MIN([min,MIN(a[where(a GT 0)])])
                        ENDFOR
                        ranges[0,i]=min
                ENDIF
        ENDFOR

        RETURN,ranges
END


