; The is a procedure for convenient placement of all labels, mapping, etc.
; We want to keep on file:
;  - axis labels for each element of dataCube, vsMdot,...
;  - abbreviated names which will appear in the names of plotfiles for each variable
;  - whether each variable is better expressed as a logarithm
;  - tex-ified versions of the labels
;  - the element of the appropriate array where this variable resides

FUNCTION GetLabels, keys=keys, names=names, tex=tex, labels=labels, log=log, offset=offset, base=base, ranges=ranges


  wrtXyt=['R (kpc)','Gas Column Density (Msun/pc^2)','Velocity Dispersion (km/s)','Stellar Column Density (Msun/pc^2)','Stellar Velocity Dispersion (km/s)','Gas Fraction','[Z/Zsun]','SFR Column Density (Msun/yr/kpc^2)','Q','Gas Q','Stellar Q','Molecular Fraction','Age At z=0 (Ga)','Age (Ga)','Depletion Time (yr)','Viscous Time (Ga)','Depletion/Viscous Time','Mass Flux Through Disk (Msun/yr)',"Gas v_r (km/s)","Star v_r (km/s)","Inv Viscous Time (yr^-1)","tdepH2","Sigma/SigmaEQ","(mu+Rf)*colSFR/colAccr","col timescale (yr)","Accr (Msun/pc^2/yr)","Accretion timescale","Ratio with Universal Prof 4","beta","v/vcirc","transPerAccr","transPerAll","IntegralTVisc","tSF","Positive Mdot (Msun/yr)","accrPerAll","sfPerAll","dsigdt_tr","dsigdt_ddx","dsigdtHeat","dsigdtCool","dZdtAdv","dZdtDiff","dZdtDil","dZdtSF","dcoldtIn","dcoldtOut","Balance","dcoldtPerAccr","colPerCrit","sigmaPerAn"]
  wrtXytex=['r (kpc)','$\Sigma\ (M_\odot/pc^2)$','$\sigma$ (km/s)','$\Sigma_*\ (M_\odot/pc^2)$','$\sigma_{*,r}$ (km/s)','$f_g$','$[Z/Z_\odot]$','$\dot{\Sigma}_*^{SF}\ (M_\odot/yr/kpc^2)$','$Q$','$Q_g$','$Q_*$','$f_{H_2}$','Age At $z=0$ (Ga)','Age (Ga)','$t_{dep}$ (yr)','$t_{visc}$ (yr)','$t_{dep}/t_{visc}$','$\dot{M}_{GI}+\dot{M}_{MRI}\ (M_\odot/yr)$',"$v_r$ (km/s)","$v_{r,*}$ (km/s)","$t_{visc}^{-1}\ (yr^{-1})$","$t_{dep,H_2}$ (yr)","$\Sigma/\Sigma_{eq}$","$(\mu+f_R)\dot{\Sigma}^{SF}/\dot{\Sigma}_{cos}$","$\Sigma/\dot{\Sigma}$ (yr)","$\dot{\Sigma}_{cos}\ (M_\odot/pc^2/yr)$","$\Sigma/\dot{\Sigma}_{cos}$","$\Sigma/\Sigma_{UP}$","$\beta$","$v/v_{circ}$","$\dot{\Sigma}_{tr}/\dot{\Sigma}_{acc}$","$\dot{\Sigma}_{tr}/(|\dot{\Sigma}_{tr}|+\dot{\Sigma}_{cos}+\dot{\Sigma}^{SF})$","Integrated $t_{visc}$","$\int^r M_g/\int^r SFR$","$\dot{M}>0 (M_\odot/yr)$","$\dot{\Sigma}_{acc}/(|\dot{\Sigma}_{tr}|+\dot{\Sigma}_{cos}+\dot{\Sigma}^{SF})$","$\dot{\Sigma}^{SF}/(|\dot{\Sigma}_{tr}|+\dot{\Sigma}_{cos}+\dot{\Sigma}^{SF})$","$\partial\sigma/\partial T_{tr}$","$\partial\sigma/\partial T_{ddx}$","$\partial\sigma/\partial T_{visc}$","$\partial\sigma/\partial T_{Cool}$","$\partial Z/\partial T_{Adv}$","$\partial Z/\partial T_{Diff}$","$\partial Z/\partial T_{Dil}$","$\partial Z/\partial T_{SF}$","$\dot{\Sigma}_{in}$","$\dot{\Sigma}_{out}$","$\dot{\Sigma}/A$","$\dot{\Sigma}/\dot{\Sigma}_{accr}$","$\Sigma/\Sigma_{crit}$","$\sigma/\sigma_{an}$"]
  wrtXyy=[9,10,11,12,13, $
      17,22,14,12,24, $
      23,48,31,32,35, $
      36,37,39,17,16, $
      38,39,40,41, 1, $
      42,43,44,15,16, $
      45,46,47,48,39, $
      49,50,40,41,42, $
      45,13,46,51,52, $
      53,54,53,54,55, $
      56] ;; index
  wrtXyp=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1] ;; offset? 
  wrtXyn=[ $
  'r', 'col','sig','colst','sigst', $   ; 5
  'fg','Z','colsfr','Q','Qg', $ ;10
  'Qst','fH2','ageAtz0','age','tdep', $ ; 15
  'tvisc','tdepOverTvisc','mdotDisk','vrg','vst', $ ;20
  'tviscinv','tdepH2','ratcoleq','sfrPerAccr','tcol', $ ;25
  'accr','taccr','BB4','beta','uu', $ ; 30
  'transPerAccr','transPerAll','intTVisc','tSF','mdotGtr0', $ ; 35
  "accrPerAll","sfPerAll","dsigTr","dsigDdx","dsigVisc", $ ;40
  "dsigCool","dZdtAdv","dZdtDiff","dZdtDil","dZdtSF", $ ; 45
  "dcoldtIn","dcoldtOut","equilibrium","dcoldtPerAccr","colPerCrit", $ ; 50
  "sigPerAn"]

  wrtXyr=dblarr(2,n_elements(wrtXyn))
  ;; these are sometimes forsaken in favor of the output of simpleRanges, which tends to do a decent job
  wrtXyr[*,*] = [ $
      [-5,45], [.5,150], [5,50], [.5,2000], [1,200], $
      [-.1,1.1],[-1.1,1.1],[1.0d-5,10.0],[1,10], [1,10], $
      [1,10],[0,1],[0,13.7d9],[0,13.7d9],[1.0d8,1.0d11], $
      [1.0d8,1.0d11],[1.0d-3,1.0d3],[-5,5],[-10,50],[-1,1], $
      [0,1.0d-8], [1.0d8,3.0d9], [.01,100],[.001,1000],[1.0d8,1.0d11], $
      [1.0d-9,1.0d-5], [1.0d8,1.0d11],[.03,30],[-.1,1.1],[-.1,1.1], $
      [-5,50],[-1,1],[1.0d8,1.0d11],[1.0d8,1.0d11],[1.0d-5,5], $
      [-.1,1.1],[-.1,1.1], [-1,1], [-1,1], [-1,1], $
      [-1,1] , [-1,1], [-1,1], [-1,1], [-1,1], $
      [-1,1], [-1,1], [-1,1], [-1,1], [.001,10], $
      [.01,100] ]


;  wrtXyr=[[0,20],[.1,1000],[5,300],[.1,1000],[5,300],[0,1],[-1,1.0],$
;    [1d-6,10],[1.5,5.0],[.5,50],[.5,50],[0,1],[1d9,14d9],[1d5,14d9],[1.0e8,1.0e10],[1.0e8,1.0e10],[.01,100]]
  wrtXyl=[1,1,0,1,1, $
      0,0,1,1,1, $
      1,0,0,0,1, $
      1,1,0,0,0, $
      0,0,1,1,1, $
      1,1,1,0,0, $
      0,0,1,1,1, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,1, $
      0] ;; log plot

  IF(n_elements(keys) GT n_elements(wrtXyt) OR n_elements(keys) EQ 0 ) THEN keys=indgen(n_elements(wrtXyt))


  labels=wrtXyt[keys]
  base=wrtXyy[keys]
  offset=wrtXyp[keys]
  names=wrtXyn[keys]
  log=wrtXyl[keys] 
  ranges=wrtXyr[*,keys]
  tex=wrtXytex[keys]

END




