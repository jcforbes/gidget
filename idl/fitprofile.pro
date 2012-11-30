;; Here we aim to fit a 1D profile (presumably some column density) w/ a exponential + sersic profile


;; functions to fit to stellar column density profiles
PRO gfunct, X, p, F, pder
        A=p[0]
        IF(A LT 0) THEN A=0
        B=p[1]
        IF(B LT 0) THEN B=0
        R0=p[2]
        R1=p[3]
        n=p[4]
        ;; sersic + exponential
        F = A*exp(-x/R0) + B*exp(-(x/R1)^(1.0/n))
        IF(n_params() GE 4) THEN BEGIN
                pder=dblarr(n_elements(x),n_elements(p))
                pder[*,0]= exp(-x/R0)
                pder[*,1]=exp(-(x/R1)^(1.0/n))
                pder[*,2]=A*exp(-x/R0)*x/(R0^2)
                pder[*,3]=B*exp(-(x/R1)^(1.0/n))*x*(x/R1)^(1.0/n-1.0)/ (n*R1*R1)
                pder[*,4]=B*exp(-(x/R1)^(1.0/n))*(x/R1)^(1.0/n)*alog(x/R1)/(n*n)
                
        ENDIF
END

PRO gfunct2, X, p, F, pder        
	A=p[0]
	R0=p[1]
        F=A*exp(-x/R0)
        IF(n_params() GE 4) THEN BEGIN
                pder=dblarr(n_elements(x),n_elements(p))
                pder[*,0] = exp(-x/R0)
                pder[*,1] = A*exp(-x/R0)*x/(R0*R0)
        ENDIF
END






PRO fitprofile,r,prof,p,chi2,single=single,colfit=colfit
    IF(n_elements(single) EQ 0) THEN single=0
    
    weights = r[*]

    ; guess:
    IF(single NE 0) THEN p = [prof[0],3.5] ELSE p = [.5*prof[0],.5*prof[0],3.5,1.0,3.0]

    ; fit
    funcName='gfunct'
    IF(single NE 0) THEN funcName='gfunct2'
    colfit = CURVEFIT(r,prof,weights,p,FUNCTION_NAME=funcName,ITMAX=500,chisq=chi2)


END
