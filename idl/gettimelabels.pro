FUNCTION GetTimeLabels, keys=keys, names=names, tex=tex, labels=labels, log=log, offset=offset, base=base, ranges=ranges, strt=strt


       vsMdotIndices = [ $
            20,11,20, 2, 3, $
             4, 8, 9,10,11, $ ; 10
            12,13,50,19,18, $ 
            51, 6,20,21,52, $ ; 20
            16, 5,23,24,25, $
            26,27,28,29,30, $ ; 30
            31,32,33,34,35, $
            36,37,38,39,40, $ ; 40
            41,42,43,44,45, $
            46,47,48,53,54, $ ; 50
            55,56,57,58,49, $
            59,60,61]
        vsMdotOffsets = [ $
             0, 0, 1, 1, 1, $
             1, 1, 1, 1, 1, $ ; 10
             1, 1, 1, 1, 1, $
             1, 1, 0, 1, 1, $ ; 20
             1, 0, 1, 1, 1, $
             1, 1, 1, 1, 1, $ ; 30
             1, 1, 1, 1, 1, $
             1, 1, 1, 1, 1, $ ; 40
             1, 1, 1, 1, 1, $
             1, 1, 1, 1, 1, $ ; 50
             1, 1, 1, 1, 1, $
             1, 1, 1]

  vsMdotLabels = [ $ 
      "Mdot (Msun/yr)", $
      "DiskSFR (Msun/yr)", $
      "Stellar Mass (Msun)", $
      "Peak Radius of Column Density (kpc)", $
      "Radius of HI transition (kpc)", $ ;; 5
      "Column Density at r=8 kpc (Msun/pc^2)", $
      "Specific SFR (yr^-1)", $
      "Bulge Gas Mass (Msun)", $
      "Bulge Stellar Mass (Msun)", $
      "H_2 Fraction", $             ; 10 
      "Mdot Gas into Bulge (Msun/yr)", $
      "Mdot Stars into Bulge (Msun/yr)", $
      "Bulge to Total Ratio", $
      "SFR Including Gas Flux Into Bulge (Msun/yr)", $
      "efficiency (Mstar / (f_b M_h))", $ ;; 15
      "Z","Gas Fraction in SF Region", $
      "Mdot (Msun/yr)", $
      "Gas Fraction in Optical Region", $
      "Z in Optical Region", $
      "Gas Fraction", $             ;; 20
      "Z Bulge", $
      "Average inward gas radial velocity (km/s)", $
      "Inner Edge of GI Region (kpc)", $
      "Outer Edge of GI Region (kpc)", $
      "Average velocity dispersion (km/s)", $  ; 25
      "Average radial velocity for non-outward velocities only (km/s)", $
      "Molecular Depletion Time (Ga)", $
      "Gas Sc Length (kpc)", $
      "St Sc Length (kpc)", $
      "sersic", $                      ; 30
      "BT Meas", $
      "gas chi^2", $
      "st chi^2", $
      "Half SFR radius (kpc)", $
      "Half gas radius (kpc)", $         ;35 
      "Half stellar radius (kpc)", $
      "Central Density (Msun/pc^2)", $
      "specific J_*", $
      "specific J_g", $
      "specific J_out", $               ; 40
      "gas Z decr (dex)", $
      "st Z decr (dex)", $
      "mol gas fraction", $
      "Depletion Time for all gas (yr)", $
      "Z*", $                           ; 45
      "Stellar Age (Gyr)", $
      "col decr (dex)", $
      "r25 1", $
      "r25 2", $
      "r25 3", $                         ; 50
      "transition col density", $
      "Halo Mass (Msun)", $
      "r25", $
      "col decr avg", $
      "sigMax", $                       ; 55
      "maxSigma", $
      "sigRatio" $
      ]
  vsMdotTexLabels = ["$\dot{M} (M_\odot/yr)$","Disk SFR ($M_\odot$/yr)","Stellar Mass ($M_\odot$)","Peak Radius (kpc)","Radius of HI trans. (kpc)","$\Sigma$ at $r=8$ kpc ($M_\odot/pc^2$)","sSFR ($yr^{-1}$)","Bulge Gas Mass ($M_\odot$)","Bulge Stellar Mass ($M_\odot$)","$H_2$ Fraction","$\dot{M}_g$ into Bulge ($M_\odot$/yr)","$\dot{M}_*$ into Bulge ($M_\odot$/yr)","Bulge:Total Ratio","SFR + $\dot{M}_{g,\rightarrow\mathrm{bulge}}$ ($M_\odot$/yr)","$\epsilon$ ($M_*$ / ($f_b M_h$))","$Z=M_Z/M_g$","$f_g$ in SF Region","$\dot{M}$ ($M_\odot$/yr)","$f_g$ in Optical Region","$Z=M_Z/M_g$ in Optical Region","$f_g$","$Z_\mathrm{bulge}=M_Z/M_*$","Average $-v_r$ (km/s)","Innermost GI Region (kpc)","Outermost GI Region (kpc)","$\langle\sigma\rangle$ (km/s)","Average $-v_r\ge 0$ (km/s)","$t_{dep} = M_g f_{H_2} / \dot{M}^{SF}$ (Ga)","Gas Sc Length (kpc)","St Sc Length (kpc)","sersic","BT Meas","gas $\chi^2$","st $\chi^2$","$r_{\mathrm{SFR},\frac12}$ (kpc)","$r_{g,\frac12}$ (kpc)","$r_{*,\frac12}$","$\Sigma_0\ (M_\odot/pc^2)$","$J_*/M_*$","$J_g/M_g$","$J_{\mathrm{out}}/M_\mathrm{out}$","$\max \Delta [Z_g]$","$\max \Delta [Z_*]$","$f_{g,\mathrm{mol}}$","$t_{dep} = M_g/\dot{M}_*$ (yr)","$Z_*=M_{Z,*}/M_*$","Stellar Age (Gyr)","$\log \Sigma(r_\mathrm{peak})/\Sigma(r_0)$ (dex)","$r_{25}$","$r_{25,2}$","$r_{25,3}$","$\Sigma_\mathrm{trans}(M_\odot/pc^2)$","$M_h (M_\odot)$","$r_{25} (kpc)$","col decr avg",'$\sigma_{max}$','$\max(\sigma)$','$\max{\sigma}/\sigma_{max}$']


  vsMdotNames  = [ $
      "mdot","DiskSFR","Mst","rPeak","rHI", $
      "colsol","sSFR","BulgeGasMass","BulgeStMass","fH2", $ ; 10
      "mdotBulgeGas","mdotBulgeStars","BT","SFRplusMdotIn","efficiency", $
      "Z","fgInSF","mdot","fgL","ZL", $ ; 20
      "fg","ZBulge","vrAvg","rQin","rQout", $
      "sigAvg","vrgGtr0","tdepAvg","GasScLength","StScLength", $ ; 30
      "sersic","BTMeas","gaschi2","stchi2","rHalfSFR", $
      "rHalfGas","rHalfSt","CentralDensity","spJSt","spJg", $ ; 40
      "spJout","DZg","DZst","fgmol","tdepAll", $
      "Zst","stAge","dCol","r25","r25b", $  ; 50
      "r25c","colTrans","Mh","r25d",'dcolAvg', $
      'sigMax','maxSig','sigRat']
  vsMdotRanges = [ $
      [.1,1000],[.1,1000],[1d9,1d11],[.001,40],[.001,40], $
      [.1,30],[1.0d-11,1.0d-9],[1.0d8,1.0d10],[1.0d8,1.0d10],[0,1], $ ; 10
      [1.0d-3,10],[1.0d-3,1],[0,1],[.1,100],[0,1], $
      [0,0],[0,0],[0,0],[0,0],[0,0], $ ; 20
      [0,0],[0,0],[0,0],[0,0],[0,0], $
      [0,0],[0,0],[0,0],[0,0],[0,0], $ ;; 30
      [0,0],[0,0],[0,0],[0,0],[0,0], $
      [0,0],[0,0],[0,0],[0,0],[0,0], $ ; 40
      [0,0],[0,0],[0,0],[0,0],[0,0], $
      [0,0],[0,0],[0,0],[0,0],[0,0], $ ; 50
      [0,0],[0,0],[0,0],[0,0],[0,0], $
      [0,0],[0,0],[0,0]  ]
  vsMdotStrt =  [ $
      1,1,0,0,0, $
      0,0,0,0,0, $
      1,1,0,1,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0,0,0, $
      0,0,0]
  vsMdotToLog = [ $
      1,1,1,1,0, $
      1,1,1,1,0, $ ; 10
      1,1,0,1,1, $
      1,1,1,0,1, $ ; 20
      1,1,0,1,1, $
      1,1,1,1,1, $ ; 30
      0,0,1,1,1, $
      1,1,1,1,1, $ ; 40
      1,0,0,1,1, $
      1,0,0,1,1, $ ; 50
      1,0,1,0,0, $
      0,0,0] 



  IF(n_elements(keys) EQ 0) THEN keys=indgen(n_elements(vsMdotLabels))
  names=vsMdotNames[keys]
  tex=vsMdotTexLabels[keys]
  labels=vsMdotLabels[keys]
  log=vsMdotToLog[keys]
  offset=vsMdotOffsets[keys]
  base=vsMdotIndices[keys]
  ranges=vsMdotRanges[keys,*]
  strt=vsMdotStrt[keys]




END

