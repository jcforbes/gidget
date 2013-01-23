;; Here we aim to fit a 1D profile (presumably some column density) w/ a exponential + sersic profile


;; functions to fit to stellar column density profiles
FUNCTION gfunct, X, p; , F, pder
        A=p[0]
        IF(A LT 0) THEN A=0
        B=p[1]
        IF(B LT 0) THEN B=0
        R0=p[2]
        R1=p[3]
        n=p[4]
        ;; sersic + exponential
        return, A*exp(-x/R0) + B*exp(-(x/R1)^(1.0/n))
;        IF(n_params() GE 4) THEN BEGIN
;                pder=dblarr(n_elements(x),n_elements(p))
;                pder[*,0]= exp(-x/R0)
;                pder[*,1]=exp(-(x/R1)^(1.0/n))
;                pder[*,2]=A*exp(-x/R0)*x/(R0^2)
;                pder[*,3]=B*exp(-(x/R1)^(1.0/n))*x*(x/R1)^(1.0/n-1.0)/ (n*R1*R1)
;                pder[*,4]=B*exp(-(x/R1)^(1.0/n))*(x/R1)^(1.0/n)*alog(x/R1)/(n*n)
;                
;        ENDIF
END

FUNCTION gfunct2, X, p; , F, pder        
	A=p[0]
	R0=p[1]
    return,A*exp(-x/R0)
;        IF(n_params() GE 4) THEN BEGIN
;                pder=dblarr(n_elements(x),n_elements(p))
;                pder[*,0] = exp(-x/R0)
;                pder[*,1] = A*exp(-x/R0)*x/(R0*R0)
;        ENDIF
END


; fit the radial profile prof to an exp or exp+sersic profile.
; chi2 is defined in the usual way, though beware that the values of the weights are fairly ad-hoc
; single tells the function whether to fit just an exp (for single=1) or the exp+sersic
; colfit stores the result of the fit
; fake forces the function to not be fit by mpfitfun, and just returns the initial guess.
FUNCTION mpfitprofile,r,prof,chi2,single=single,colfit=colfit,fake=fake
    IF(n_elements(single) EQ 0) THEN single=0
    IF(n_elements(fake) EQ 0) THEN fake=0
    ;    weights = r[*] ;    avg = exp(average(alog(prof)))
    average = avg(prof)
    weights = intarr(n_elements(prof))+1.0/(.4*average)^2
;    errs = 0.1*sqrt(average*prof)
    errs = .01*prof

    IF (single EQ 1) THEN np = 2 ELSE np=5
    ; require that all parameters be positive:
    pi = replicate({fixed:0, limited:[1,0], limits:[0.D,0.D]},np)

    maxlength=300.D 
    minsersic=0.01D
    minlength=0.01D ; ten parsecs
    IF(single EQ 1) THEN BEGIN
        pi[1].limited[1]=1
        pi[1].limits[1]=maxlength
    ENDIF ELSE BEGIN
        pi[2].limited[1]=1
        pi[3].limited[1]=1
        pi[4].limited[1]=1
        pi[2].limits[1]=maxlength
        pi[3].limits[1]=maxlength
        pi[4].limits[1]=10.D
        pi[2].limits[0]=minlength
        pi[3].limits[0]=minlength
        pi[4].limits[0]=minsersic
    END

    ;; a guess for the length scale:
    length = -(r[n_elements(r)-1] - r[0])/alog(prof[n_elements(prof)-1]/prof[0])
    IF(length GT maxlength) THEN length=maxlength
    IF(length LT minlength) THEN length=maxlength ;; not a typo. If the guess for length is negative, then the best we can do is guess a flat profile, i.e. we want a large length, not a small one.
    ; guess:
    IF(single NE 0) THEN p = [prof[0],length] ELSE p = [.7*prof[0],.3*prof[0],length,length,1.0]

    ; fit
    funcName='gfunct'
    IF(single NE 0) THEN funcName='gfunct2'
;    colfit = CURVEFIT(r,prof,weights,p,FUNCTION_NAME=funcName,ITMAX=200,chisq=chi2,TOL=.1)
    IF(fake EQ 1) THEN BEGIN
        result=p[*]
        bestnorm=1
        dof=1
        status=1
    ENDIF ELSE result = MPFITFUN( funcName, r, prof, errs, p, PARINFO=pi, QUIET=1, STATUS=status, GTOL=1D-5, FTOL=1D-5,BESTNORM=bestnorm,DOF=dof)
    IF(single EQ 1) THEN colfit=gfunct2(r,result) ELSE colfit=gfunct(r,result)
    ;PRINT,"mpfitprofile status: ", status
    chi2 = bestnorm/dof
    IF(status LE 0) THEN stop,"status: ",status

    return,result


END
