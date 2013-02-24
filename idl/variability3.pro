

;; given an array whose first column indexes time, second column is irrelevant, third column
;;  indexes different data sets fourth column indexes models, plus corresponding time
;; axis labels, colors, and styles, make a series of simple plots
PRO simplePlot,data,time,labels,names,colors,name,sv,styles=styles
  IF(n_elements(styles) EQ 0) THEN styles = intarr(n_elements(data[0,0,0,*]))
  FOR k=0,n_elements(labels)-1 DO BEGIN;; loop over y-axis variables to be plotted against time
    set_plot,'x'
    figureInit,name+"_"+names[k],sv,2,2
    
    PLOT,[MAX(time)],[MAX(data[*,0,k,*])],COLOR=0, BACKGROUND=255,THICK=2,XTITLE="time",YTITLE=labels[k],XSTYLE=1,YSTYLE=1
    FOR j=0,n_elements(data[0,0,0,*])-1 DO BEGIN
      ;PLOT,time,data[*,0,k,j],COLOR=colors[0,j],BACKGROUND=255,linestyle=styles[j],XTITLE="time",YTITLE=labels[k],XSTYLE=1,YSTYLE=1
      OPLOT,time,data[*,0,k,j],COLOR=colors[0,j],linestyle=styles[j]
    ENDFOR

    figureClean,name+"_"+names[k]
  ENDFOR
END

FUNCTION getElements,nvardim,i
  telements=(unpack(i,nVarDim))
  srtd=nVarDim[reverse(sort(nvardim))]
  elements=telements[reverse(sort(nvardim))]
  IF(n_elements(elements) LT 3) THEN BEGIN
    qelements=intarr(3)
    qelements[0:n_elements(elements)-1]=elements[0:n_elements(elements)-1]
    elements=qelements
  ENDIF

  RETURN,elements
END


FUNCTION vsMstarTransform,original,NVS
    vsMstarLabels = original[*]
    vsMstarLabels[0] = vsMstarLabels[2]
    vsMstarLabels[2:(NVS-2)] = vsMstarLabels[3:(NVS-1)]
    vsMstarLabels=vsMstarLabels[0:(NVS-2)]
    RETURN,vsMstarLabels
END




;;; Meant to do an analysis of many runs at once, particularly when those 
;;; correspond to different accr. histories
;;; expName is an experiment name, like one would specify in exper.py. This 
;;;   code assumes that you are 
;;;     1) in the gidget/analysis directory, and 2) you have a bunch of 
;;;     runs in a folder with names of the form
;;;      expName/expName<identifiers>_<filenames>, 
;;;     where identifiers can be anything as long as it is unique to each run
;;;     and <filenames> encompasses all the files output by the code, e.g. 
;;;     evolution.dat, radial.dat, etc.
;;; keys is an array which lets you specify which variable to plot from the 
;;;   list in wrtXyt. This is useful for very large experiments, where even 
;;;   storing the modest-size arrays here which take a limited amount of
;;;   data for each run can become cumbersome. By default the code will
;;;   analyze all the variables
;;; N controls how many runs in the experiment to analyze. For large number 
;;;   runs, sometimes it's useful to just look at a few runs. If N is larger 
;;;   than the number of runs in the experiment, all the runs are analyzed.
;;; sv- behaves like other instances of sv; 0 to just show the plots on the 
;;;   screen, 2 to save them as png's and subsequently movies. 4 saves 
;;;   movies in a format conducive to showing on a low-resolution projector.
PRO variability3,expNames,N=N,sv=sv,keys=keys,notations=notations
  IF(n_elements(N) EQ 0) THEN N=1000000
  IF(n_elements(sv) EQ 0) THEN sv=3
  IF(n_elements(notations) EQ 0) THEN notations=expNames[*]

  compile_opt idl2
  set_plot,'x'
  DEVICE,DECOMPOSED=0
  DEVICE,RETAIN=2

  proftimes=dblarr(5)
  proftimes[0]=systime(1)

  nameList2=ListValidModels(expNames[0],N)
  modelcounter = [0,  n_elements(nameList2) ]
  FOR expInd=1,n_elements(expNames)-1 DO BEGIN
      nameList3 = ListValidModels(expNames[expInd],N)
      nameList2 = [nameList2,nameList3]
      modelcounter = [modelcounter, n_elements(namelist2)]
  ENDFOR

  ctr=0

  IF(n_elements(keys) EQ 0) THEN keys=indgen(100)+1
  tk=[0,keys[*]]

  ;; variables to plot- title, column index, log plot, name, y-range
  dummy = GetLabels(keys=tk,tex=wrtXytex, labels=wrtXyt, names=wrtXyn, log=wrtXyl, base=wrtXyy, offset=wrtXyp)

  dummy = where( keys EQ 28, n1)
  IF(n1 EQ 2) THEN n1=1




  ;; now loop over every model in the list - 
  ;; all we're doing here is reading in each model then freeing 
  ;; the associated memory
  ;; so that we can iterate over a very large number. For each model, 
  ;; record the data specified in the arrays above (e.g. yy)
  modelPtr=ReadModels([nameList2[0]])
  model = (*(modelPtr[0]))
  theNTS = model.NOutputsNominal+2

  whichFrames = indgen(theNTS) ;; plot every output
  whichFrames = indgen((theNTS+2)/3)*3 ; plot every 3rd output.
  whichFrames3 = [1,51,201] ;; roughly z=2,1,0

  time = model.evArray[1,*]  
  z = model.evArray[10-1,*]

  ;; (# timesteps) x (nx) x (# of columns) x (# models)
  theData= dblarr(theNTS,n_elements(model.dataCube[0,*,0]),n_elements(wrtXyy),n_elements(nameList2))

  NVS = 55
  ;; (# timesteps) x (# pts/model frame = 1) x (one thing vs another = 20) x (# models)
  vsMdot = dblarr(theNTS,1,NVS,n_elements(nameList2))
  nameList4 = strarr(n_elements(nameList2))

  ; so far we have set up two data structures, namely one for quantities which vary in radius & time, 
  ; and one for quantities which vary only in time. Now in analogy to the first structure, we will 
  ; construct an array for stellar properties at redshift zero which vary with radius & age.
  byAge = dblarr(model.npassive+1,n_elements(model.dataCube[0,*,0]),5+1,n_elements(nameList2))
  vsAge = dblarr(3,model.npassive+1,5+1,n_elements(nameList2))

  proftimes[1]=systime(1)

;  FOR i=0,n_elements(nameList2)-1 DO BEGIN ;; loop over every valid model
  ctr=0
  modelcounter2 = modelcounter*0

  FOR expInd=1, n_elements(expNames) DO BEGIN
    low = modelcounter[expInd-1]
    high= modelcounter[expInd]-1
    ctr2=0
    FOR i=low,high DO BEGIN
 
       IF(i MOD 50 EQ 0) THEN PRINT,"Processing model #",i," of ",n_elements(nameList2)

       ;; the number of colors necessary is srtd[0], 
       ;; which color set to use is elements[2]
       IF(i NE 0) THEN modelPtr=ReadModels([nameList2[i]])
       IF(i NE 0) THEN model = (*(modelPtr[0]))
       zi=nearest("redshift",1,model,0.0)
       offsets=[0,model.ncolstep,model.ncolstep+model.npostprocess]
       IF(n_elements(model.evArray[1,*]) EQ theNTS) THEN BEGIN
         time = model.evArray[1,*]  
         z = model.evArray[10-1,*]
       ENDIF
       x = model.dataCube[0,*,0]
       r = x[*]*model.Radius

       ;; if this model has the appropriate number of time steps..
       ;;PRINT,"Number of time steps",n_elements(model.dataCube[*,0,0]),n_elements(theData[*,0,0,0])
       IF(n_elements(model.dataCube[*,0,0]) EQ theNTS) THEN BEGIN
         FOR k=0,n_elements(wrtXyy)-1 DO BEGIN	
           theData[*,*,k,ctr] = model.dataCube[*,*,wrtXyy[k]+offsets[wrtXyp[k]]-1]
         ENDFOR

         positionIndices=[nearest("position",1,model,1.0/model.Radius), $
            nearest("position",1,model,8.0/model.Radius), $
            nearest("position",1,model,15.0/model.Radius)]
;         IF(model.npassive GT n_elements(byAge[*,0,0,0])) THEN BEGIN
;            byAgeT=temporary(byAge)[*,*,*,*]
;            byAge=dblarr(model.npassive+1,n_elements(model.dataCube[0,*,0]),5+1,n_elements(namelist2))
;            byAge[0:n_elements(byAgeT[*,0,0,0])-1,*,*,*] = temporary(byAgeT)[*,*,*,*]
;         ENDIF
         FOR popI=0, model.npassive DO BEGIN
            byAge[popI, *, 0, ctr] = model.dataCube[0,*,0]*model.Radius
            byAge[popI, *, 1:5, ctr] = model.stellarData[theNTS-1, popI, *, 0:4]
            FOR m=0,n_elements(positionIndices)-1 DO BEGIN
                vsAge[m,popI,1:5, ctr] = model.stellarData[theNTS-1,popI, positionIndices[m], 0:4]
            ENDFOR
         ENDFOR
         ages = model.stellarAges[theNTS-1,*]
         FOR m=0,n_elements(positionIndices)-1 DO vsAge[m,*,0,ctr] = ages[*]

         vsMdot[*,0,0,ctr] = model.evArray[20-1,*] ;; Cosmological accretion rate.
         vsMdot[*,0,1,ctr] = model.evArray[10,*] ;; SFR [in the disk!]

         FOR zi=0, n_elements(model.evArray[9,*])-1 DO BEGIN ; loop over redshift
           modelInfo = diskStats(model,z=model.evArray[9,zi])
           vsMdot[zi,0,3,ctr] = modelInfo[1] ;; peak radius of column density (kpc)
           vsMdot[zi,0,4,ctr] = modelInfo[2] ;; edge of SF region (kpc)
           vsMdot[zi,0,5,ctr] = modelInfo[3] ;; col @ r=8kpc (solar masses/pc^2)
           vsMdot[zi,0,6,ctr] = modelInfo[7] ;; sSFR (yr^-1)
           vsMdot[zi,0,7,ctr] = modelInfo[8] ;; gas mass in bulge (MSol)
           vsMdot[zi,0,8,ctr] = modelInfo[9] ;; stellar mass in bulge
           vsMdot[zi,0,9,ctr] = modelInfo[10] ;; total f_H2
           vsMdot[zi,0,10,ctr] = modelInfo[11] ;; mdot bulge gas
           vsMdot[zi,0,11,ctr] = modelInfo[12] ;; mdot bulge stars
           vsMdot[zi,0,2,ctr] = modelInfo[19] ;; stellar mass in disk + bulge
           vsMdot[zi,0,12,ctr] = modelInfo[8]/modelInfo[19] ;; B:T
           vsMdot[zi,0,13,ctr] = modelInfo[18] ;; SFR + mdotbulgeGas
           vsMdot[zi,0,14,ctr] = modelInfo[17] ;; efficiency
           vsMdot[zi,0,15,ctr] = modelInfo[16]*.02 ; M_Z/M_gas ;alog10(modelInfo[16]) ;; gas metallicity
           vsMdot[zi,0,16,ctr] = modelInfo[5] ;; fG in SF region.
           vsMdot[zi,0,17,ctr] = model.evArray[20-1,zi] ;; Cosmological accretion rate again   (modelInfo[18] - modelInfo[11])/(modelInfo[18]+1.0d-10) ;; fraction of SF in the bulge
           vsMdot[zi,0,18,ctr] = modelInfo[20] ;; fgL
           vsMdot[zi,0,19,ctr] = modelinfo[21]*.02 ;alog10(modelInfo[21]) ;; ZL
           vsMdot[zi,0,20,ctr] = modelInfo[15] ;; fg
	       vsMdot[zi,0,21,ctr] = model.evArray[5-1,zi] ;M_Z,bulge/M_bulge alog10(model.evArray[5-1,zi]/.02) ;; ZBulge
	       vsMdot[zi,0,22,ctr] = modelInfo[23-1];; vrAvg (mass-averaged, km/s)
	       vsMdot[zi,0,23,ctr] = modelInfo[24-1];; GIin (kpc)
	       vsMdot[zi,0,24,ctr] = modelInfo[25-1];; Giout (kpc)
       	   vsMdot[zi,0,25,ctr] = modelInfo[26-1];; velocity dispersion (mass-averaged, km/s)
       	   vsMdot[zi,0,26,ctr] = modelInfo[27-1];; vrAvgGt0
	       vsMdot[zi,0,27,ctr] = modelInfo[28-1];; tdep = gas mass / SFR
           vsMdot[zi,0,28,ctr] = modelInfo[29-1];; gas scale length (kpc)
           vsMdot[zi,0,29,ctr] = modelInfo[30-1];; stellar scale length (kpc)
           vsMdot[zi,0,30,ctr] = modelInfo[31-1];; stellar sersic index
           vsMdot[zi,0,31,ctr] = modelInfo[32-1];; measured BT
           vsMdot[zi,0,32,ctr] = modelInfo[33-1];; chi^2 gas
           vsMdot[zi,0,33,ctr] = modelInfo[34-1];; chi^2 stellar
           vsMdot[zi,0,34,ctr] = modelInfo[35-1];; sfrHalfRadius (kpc)
           vsMdot[zi,0,35,ctr] = modelInfo[36-1];; gasHalfRadius (kpc)
           vsMdot[zi,0,36,ctr] = modelInfo[37-1];; stHalfRadius (kpc)
           vsMdot[zi,0,37,ctr] = modelInfo[38-1];; central density (Msun/pc^2)
           vsMdot[zi,0,38,ctr] = modelInfo[39-1];; specific stellar J (kpc km/s)
           vsMdot[zi,0,39,ctr] = modelInfo[40-1];; specific gas J (kpc km/s)
           vsMdot[zi,0,40,ctr] = modelInfo[41-1];; specific outflow J (kpc km/s)
           vsMdot[zi,0,41:47,ctr] = modelInfo[42-1:48-1];gasDZ(dex),stDZ(dex),fgmol,tdep[all](yr),stZ[],stAge(Gyr),dCol(dex)
           vsMdot[zi,0,48:51,ctr] = [model.x25s[zi]*model.Radius,model.x25s_2[zi]*model.Radius, $
                                    model.x25s_3[zi]*model.Radius,model.colTranses[zi]]
           vsMdot[zi,0,53,ctr] = model.x25s_4[zi]*model.Radius

           vsMdot[zi,0,52,ctr] = modelInfo[19] /(.18 * modelInfo[17]) ;; = M_h
           vsMdot[zi,0,54,ctr] = modelInfo[49-1] ; col decr (dex) averaged over inner 10 cells
         ENDFOR ;; end loop over redshift

         nameList4[ctr] = nameList2[i]

         ctr2 = ctr2+1
         ctr = ctr+1
       ENDIF ELSE BEGIN
          print,"Warning: discarding model i: ",i," named ",model.name," because it has ",n_elements(model.dataCube[*,0,0])," time steps instead of ",n_elements(theData[*,0,0,0])

       END ;; end check of whether this model has the correct # of time steps.

       IF(i EQ 0) THEN print,model.name,diskStats(model)
       ptr_free,modelPtr[0]
    ENDFOR ;; end loop over specific models (i).
    modelcounter2[expInd] = ctr
  ENDFOR ;; end loop over expNames


;['mdot','DiskSFR','Mst','rPeak','rHI','colsol','sSFR','BulgeGasMass','BulgeStMass','fH2','mdotBulgeGas','mdotBulgeStars',"BT","SFRplusMdot_b_gas","efficiency","Z","f_gInSF","mdot","fgL","ZL","fg","ZBulge","vrgAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg"]

  vsMdotLabels = ["Mdot (Msun/yr)","DiskSFR (Msun/yr)","Stellar Mass (Msun)","Peak Radius of Column Density (kpc)","Radius of HI transition (kpc)","Column Density at r=8 kpc (Msun/pc^2)","Specific SFR (yr^-1)","Bulge Gas Mass (Msun)","Bulge Stellar Mass (Msun)","H_2 Fraction","Mdot Gas into Bulge (Msun/yr)","Mdot Stars into Bulge (Msun/yr)","Bulge to Total Ratio","SFR Including Gas Flux Into Bulge (Msun/yr)","efficiency (Mstar / (f_b M_h))","Z","Gas Fraction in SF Region","Mdot (Msun/yr)","Gas Fraction in Optical Region","Z in Optical Region","Gas Fraction","Z Bulge","Average inward gas radial velocity (km/s)","Inner Edge of GI Region (kpc)","Outer Edge of GI Region (kpc)","Average velocity dispersion (km/s)","Average radial velocity for non-outward velocities only (km/s)","Molecular Depletion Time (Ga)","Gas Sc Length (kpc)","St Sc Length (kpc)","sersic","BT Meas","gas chi^2","st chi^2","Half SFR radius (kpc)","Half gas radius (kpc)","Half stellar radius (kpc)","Central Density (Msun/pc^2)","specific J_*","specific J_g","specific J_out","gas Z decr (dex)","st Z decr (dex)","mol gas fraction","Depletion Time for all gas (yr)","Z*","Stellar Age (Gyr)","col decr (dex)","r25","r25 2","r25 3","transition col density","Halo Mass (Msun)","r25 4","col decr avg"]
  vsMdotTexLabels = ["$\dot{M} (M_\odot/yr)$","Disk SFR ($M_\odot$/yr)","Stellar Mass ($M_\odot$)","Peak Radius (kpc)","Radius of HI trans. (kpc)","$\Sigma$ at $r=8$ kpc ($M_\odot/pc^2$)","sSFR ($yr^{-1}$)","Bulge Gas Mass ($M_\odot$)","Bulge Stellar Mass ($M_\odot$)","$H_2$ Fraction","$\dot{M}_g$ into Bulge ($M_\odot$/yr)","$\dot{M}_*$ into Bulge ($M_\odot$/yr)","Bulge:Total Ratio","SFR + $\dot{M}_{g,\rightarrow\mathrm{bulge}}$ ($M_\odot$/yr)","$\epsilon$ ($M_*$ / ($f_b M_h$))","$Z=M_Z/M_g$","$f_g$ in SF Region","$\dot{M}$ ($M_\odot$/yr)","$f_g$ in Optical Region","$Z=M_Z/M_g$ in Optical Region","$f_g$","$Z_\mathrm{bulge}=M_Z/M_*$","Average $-v_r$ (km/s)","Innermost GI Region (kpc)","Outermost GI Region (kpc)","$\langle\sigma\rangle$ (km/s)","Average $-v_r\ge 0$ (km/s)","$t_{dep} = M_g f_{H_2} / \dot{M}^{SF}$ (Ga)","Gas Sc Length (kpc)","St Sc Length (kpc)","sersic","BT Meas","gas $\chi^2$","st $\chi^2$","$r_{\mathrm{SFR},\frac12}$ (kpc)","$r_{g,\frac12}$ (kpc)","$r_{*,\frac12}$","$\Sigma_0\ (M_\odot/pc^2)$","$J_*/M_*$","$J_g/M_g$","$J_{\mathrm{out}}/M_\mathrm{out}$","$\max \Delta [Z_g]$","$\max \Delta [Z_*]$","$f_{g,\mathrm{mol}}$","$t_{dep} = M_g/\dot{M}_*$ (yr)","$Z_*=M_{Z,*}/M_*$","Stellar Age (Gyr)","$\log \Sigma(r_\mathrm{peak})/\Sigma(r_0)$ (dex)","$r_{25}$","$r_{25,2}$","$r_{25,3}$","$\Sigma_\mathrm{trans}$","$M_h (M_\odot)$","$r_{25,4} (kpc)$","col decr avg"]

  vsMdotNames  = ["mdot","DiskSFR","Mst","rPeak","rHI","colsol","sSFR","BulgeGasMass","BulgeStMass","fH2","mdotBulgeGas","mdotBulgeStars","BT","SFRplusMdotIn","efficiency","Z","fgInSF","mdot","fgL","ZL","fg","ZBulge","vrAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg","GasScLength","StScLength","sersic","BTMeas","gaschi2","stchi2","rHalfSFR","rHalfGas","rHalfSt","CentralDensity","spJSt","spJg","spJout","DZg","DZst","fgmol","tdepAll","Zst","stAge","dCol","r25","r25b","r25c","colTrans","Mh","r25d",'dcolAvg']
  vsMdotStrt =  [1,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  vsMdotToLog = [1,1,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1,0] 

  byAgeLabels=["Radius (kpc)","Column Density","r Vel Dis","z Vel Dis","Z","var Z"]
  byAgeTexLabels=["Radius (kpc)","$\Sigma_{*,i}$","$\sigma_{*,r}$","$\sigma_{*,z}$","Z","var Z"]
  byAgeNames=["r","stCol","sigr","sigz","stZ","stVZ"]
  byAgeToLog=[1,1,1,1,0,0]

  vsAgeLabels=byAgeLabels[*]
  vsAgeTexLabels=byAgeTexLabels[*]
  vsAgeNames=byAgeNames[*]
  vsAgeToLog=byAgeToLog[*]
  vsAgeLabels[0]="Ages (Ga)"
  vsAgeTexLabels[0]="Ages (Ga)"
  vsAgeNames[0]="ages"

  modelcounter = modelcounter2[*]
  nameList2 = nameList4[0:ctr-1]

  proftimes[2]=systime(1)

  theData = (temporary(theData))[*,*,*,0:ctr-1]
  colors = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))
  vsMdot = (temporary(vsMdot))[*,*,*,0:ctr-1]
  byAge = (temporary(byAge))[*,*,*,0:ctr-1]
  vsAge = (temporary(vsAge))[*,*,*,0:ctr-1]

  thePercentiles = [0.025,.16,.5,.84,.975]
  nper = n_elements(thePercentiles)
  intervals = dblarr(n_elements(theData[*,0,0,0]),  n_elements(theData[0,*,0,0]), n_elements(theData[0,0,*,0]), nper*n_elements(expNames))  
  ;; having accumulated theData, let's do some analysis!
;  sortdData = theData * 0.0
;  sortdData[*,*,0,*] = theData[*,*,0,*] ;; don't sort radii..
  ;sortdVsMdot = vsMdot*0.0
  FOR expInd=1, n_elements(expNames) DO BEGIN
    low = modelcounter[expInd-1]
    high= modelcounter[expInd]-1
    ;; for all time & space, set the x-values of intervals = appropriate x-values from theData
    FOR perInd = 0, nper -1 DO $
        intervals[*,*,0,(expInd-1)*nper+perInd] = theData[*,*,0,low]

    FOR j=0,n_elements(theData[*,0,0,0])-1 DO BEGIN ;; loop over time
      FOR k=1,n_elements(wrtXyy)-1 DO BEGIN ;; loop over variables
        FOR xx=0,n_elements(theData[0,*,0,0])-1 DO BEGIN ;; loop over x-coord
    	  ;;; Sort the models, within each experiment
          intervals[j,xx,k,(expInd-1)*nper:(expInd)*nper-1] = percentiles(theData[j,xx,k,low:high],thePercentiles,wrtXyl[k])
	    ENDFOR ;; end loop over x-coord 
      ENDFOR ;; end loop over variables
    ENDFOR ;; end loop over time
  ENDFOR ;; end loop over experiments



  proftimes[3]=systime(1)

  FOR expInd=1, n_elements(expNames) DO BEGIN
    low=modelcounter[expInd-1]
    high=modelcounter[expInd]-1
    colors[*,low:high] = expInd-1
  ENDFOR


  nmodels = n_elements(vsMdot[0,0,0,*])
  ncolors = 5
  IF(n_elements(expNames) EQ 1) THEN BEGIN
    zi0 = n_elements(vsMdot[*,0,0,0])-1
    FOR zi = 0, zi0 DO BEGIN
      FOR nc=0, ncolors-1 DO BEGIN
        ; sort according to z=0 sSFR - feel free to change this!
        ;; examples include structural parameters, e.g. CentralDensity, rPeak, rHalfSFR,..
        ;; also metallicity, initial mass,...
        colors[zi,GetIthBin(nc, vsMdot[zi0,0,6,*], ncolors)] = nc
        
      ENDFOR
    ENDFOR
  ENDIF

  thicknesses=intarr(nmodels)+1
  Mh0 = vsMdot[theNTS-1,0,2,*]/vsMdot[theNTS-1,0,14,*]
  FOR nc=0, ncolors-1 DO BEGIN
      thicknesses[GetIthBin(nc,Mh0,ncolors)] = 12.0 + .5*float(nc) ;; thicker lines have higher Mh0.
  ENDFOR


  expName2 = expNames[0]
  IF(n_elements(expNames) LT 5) THEN BEGIN
      FOR i=1, N_elements(expNames)-1 DO BEGIN
          expName2 = expName2 + '_' + expNames[i];; +'_'+strcompress(string(MIN([n_elements(nameList2),N])),/remove) ;; e.g. rk5_113
      ENDFOR
  ENDIF ELSE expName2 +='_' + strcompress(string(n_elements(expNames)),/remove)
  expName2 = expName2 + '_' + strcompress(string(MIN([n_elements(nameList2),N*n_elements(expNames)])),/remove) ;; e.g. rk5_rk6_113

  unsRanges = simpleRanges(intervals[2:n_elements(time)-1,*,*,*],wrtXyl)
  unsRanges[0,[1,3]] = 0.1 ;; manually set minimum col and colst to 1 Msun/pc^2
  unsRanges[0,7] =  1.0d-6 ;; manually set min SFR col density to 10^-5 Msun/yr/kpc^2
  unsRanges[1,[8,9,10]] = 50.0 ;; manually set max Q,Qgas,Qst
  unsRanges[0,17] = -2.0
  unsRanges[1,17] = 4.0
  
  unsRangesLin = unsRanges[*,*]
  unsRangesLin[0,0] = -2.0
  unsRangesLin[1,0] = 42.0

      linestyles = intarr(n_elements(expNames)*nper)+2;nper-1
      midpoints = where(indgen(n_elements(expNames)*nper) MOD nper EQ nper/2)
      linestyles[midpoints] = 0
      intervalsColors0 = indgen(n_elements(expNames)*nper)/nper
      intervalsColors = intarr(n_elements(time),n_elements(intervalsColors0))
      FOR tt=0,n_elements(time)-1 DO intervalsColors[tt,*] = intervalsColors0 


    unsThicknesses=colors*0+4
    unsThicknesses[0]=12
	


  ;; average mdot over a number of time steps nt
  nt = 40
  avgdVsMdot = vsMdot[*,*,*,*]

  FOR ti=0,n_elements(time)-1 DO BEGIN
    FOR mm=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
      avgdVsMdot[ti,0,0,mm] = TOTAL(vsMdot[MAX([0,ti-nt]):ti,0,0,mm]) / (ti - MAX([0,ti-nt]) + 1)
    ENDFOR
  ENDFOR



     vsAvgMdotLabels = vsMdotLabels[*]
     vsAvgMdotLabels[0] = "Average Mdot Over Past 40 Outputs (Msun/yr)"
     vsAvgMdotTexLabels = vsMdotTexLabels[*]
     vsAvgMdotTexLabels[0] = "$\dot{M}_{40 \mathrm{outputs}} (M_\odot/yr)$"
     vsAvgMdotNames = vsMdotNames[*]
     vsAvgMdotNames[0] = "avgMdot"



    vsSpecificSFR = vsMdot[*,*,*,*]
    vsSpecificSFR[*,*,0,*] = vsMdot[*,*,6,*]
    vssSFRLabels = vsMdotLabels[*]
    vssSFRLabels[0] = vsMdotLabels[6]
    vssSFRTexLabels = vsMdotTexLabels[*]
    vssSFRTexLabels[0] =vsMdotTexLabels[6]
    vssSFRNames = vsMdotNames[*]
    vssSFRNames[0] = vsMdotNames[6]
    vssSFRToLog = vsMdotToLog[*]
    vssSFRToLog[0] = vsMdotToLog[6]


  ;; Plot these quantities vs. the star formation rate.
  vsSFR = vsMdot[*,*,1:n_elements(vsMdot[0,0,*,0])-1,*] ;; cut out the mdot slice of the array

    vsSFRLabels = vsMdotLabels[1:(n_elements(vsMdotLabels)-1)]
    vsSFRTexLabels = vsMdotTexLabels[1:n_elements(vsMdotTexLabels)-1]
    vsSFRNames = vsMdotNames[1:(n_elements(vsMdotNames)-1)]
    vsSFRToLog = vsMdotToLog[1:(n_elements(vsMdotToLog)-1)]
    vsSFRstrt= vsMdotstrt[1:(n_elements(vsMdotStrt)-1)]



  ;; Now do a series of plots vs. M*
  vsMstar = vsMdot[*,*,*,*]
  vsMstar[*,*,0,*] = vsMstar[*,*,2,*] ;replace Mdot with M* as the zeroeth column of the data.
  vsMstar[*,*,2:(NVS-2),*] =vsMstar[*,*,3:(NVS-1),*] ;move over the other columns so that we don't produce an M* vs M* plot
  vsMstar=vsMstar[*,*,0:(NVS-2),*] ; eliminate the extraneous last column.

  vsMstarLabels = vsMstarTransform(vsMdotLabels,NVS)
  vsMstarTexLabels = vsMstarTransform(vsMdotTexLabels,NVS)
  vsMstarNames = vsMstarTransform(vsMdotNames,NVS)
  vsMstarStrt = vsMstarTransform(vsMdotStrt,NVS)
  vsMstarToLog = vsMstarTransform(vsMdotToLog,NVS)



  vsTimeLabels = vsMdotLabels[*]
  vsTimeLabels[0] = "Time Since z=2 (Ga)"
  vsTimeTexLabels = vsMdotTexLabels[*]
  vsTimeTexLabels[0] = "Time since $z=2$ (Ga)"
  vsTimeNames = vsMdotNames[*]
  vsTimeToLog = vsMdotToLog[*]
  vsTimeToLog[0] = 0

  nmodels = n_elements(vsMdot[0,0,0,*])
  vsTime = dblarr(1,n_elements(vsMdot[*,0,0,0]),n_elements(vsMdot[0,0,*,0]),nmodels)
  FOR tt=0, n_elements(vsMdot[*,0,0,0])-1 DO BEGIN
    vsTime[0,tt,*,*] = vsMdot[tt,0,*,*]
    vsTime[0,tt,0,*] = time[tt]
  ENDFOR

  vsTimeStyles=colors*0
  vsTimeColors=colors[*,*]
  IF(n_elements(vsTime[0,0,0,*]) GT 50) THEN vsTimeStyles=temporary(vsTimeStyles)+1

  svSinglePlot=2 ;; 2=eps, 4=png

  setct,3,n_elements(vsMdot[0,0,0,*]),0
;  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2,sv=4,nt0=2,thicknesses=thicknesses
  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2+"_log",sv=svSinglePlot,nt0=2,thicknesses=thicknesses,logarithms=vsMdotToLog,normalize=(n_elements(expNames) NE n_elements(vsMdot[0,0,0,*]))



  ;;;; MOVIES ;;;;;;
;;  simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervals",sv,axislabels=wrtXyt,whichFrames=whichFrames,ranges=unsRanges
;;  IF(n_elements(vsSpecificSFR[0,0,0,*]) LE 50) THEN simpleMovie, vsSpecificSFR,z,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFR",sv,PSYM=1,prev=1,axisLabels=vssSFRLabels,taillength=2,whichFrames=whichFrames
;;  IF(n_elements(vsSpecificSFR[0,0,0,*]) GT 50) THEN simpleMovie, vsSpecificSFR,z,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFRDist",sv,PSYM=1,prev=0,axisLabels=vssSFRLabels,NIndVarBins=5,whichFrames=whichFrames
;;  IF(n_elements(vsMstar[0,0,0,*]) LE 50) THEN simpleMovie, vsMstar, z,vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstar",sv,PSYM=1,prev=0,axisLabels=vsMstarLabels,whichFrames=whichFrames
;;  IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, z,vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstarDist",sv,PSYM=1,prev=0,axisLabels=vsMstarLabels,NIndVarBins=5,percentileList=[.17,.5,.83],whichFrames=whichFrames
;;  simpleMovie, avgdVsMdot,z,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot",sv,PSYM=1,prev=1,axisLabels=vsAvgMdotLabels,whichFrames=whichFrames,thicknesses=unsThicknesses
;;  simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsorted",sv,axislabels=wrtXyt,whichFrames=whichFrames,ranges=unsRanges
;;;    simpleMovie, vsSFR, z,vsSFRNames, colors, colors*0, vsSFRtoLog,expName2+"_vsSFR",sv,vsSFRstrt,PSYM=1,prev=1,axisLabels=vsSFRLabels
;;;    simpleMovie, vsMdot,z,vsMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsmdot",sv,vsMdotStrt,PSYM=1,prev=1,axisLabels=vsMdotLabels


   ;;;;; STANDALONE PNG's ;;;;;
   ;vs radius. First two are log radius, second two are identical with linear radius
  TheNumberOfModels = n_elements(vsMstar[0,0,0,*])

  whichFramesz0 = [n_elements(z)-1]
  whichFramesz2 = [1]
  whichFrames5 = [1,20,51,105,201]; a guess for z=2,1.5,1,.5,0

  setct,1,n_elements(linestyles),2

  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervalsLogR", $
          5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRanges,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsortedLogR", $
      5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRanges,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_intervals", $
          5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRanges2,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsorted", $
      5,axislabels=wrtXyt,whichFrames=whichFrames3,ranges=unsRanges2,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex


  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervalsLogRz0", $
          5,axislabels=wrtXyt,whichFrames=whichFramesz0,ranges=unsRanges,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsortedLogRz0", $
      5,axislabels=wrtXyt,whichFrames=whichFramesz0,ranges=unsRanges,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_intervalsz0", $
          5,axislabels=wrtXyt,whichFrames=whichFramesz0,ranges=unsRanges2,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedz0", $
      5,axislabels=wrtXyt,whichFrames=whichFramesz0,ranges=unsRanges2,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex

    IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervalsLogRz2", $
          5,axislabels=wrtXyt,whichFrames=whichFramesz2,ranges=unsRanges,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsortedLogRz2", $
      5,axislabels=wrtXyt,whichFrames=whichFramesz2,ranges=unsRanges,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_intervalsz2", $
          5,axislabels=wrtXyt,whichFrames=whichFramesz2,ranges=unsRanges2,horizontal=1,thicknesses=linestyles[*]+12, $
          svSinglePlot=svSinglePlot,texLabels=wrtXytex
  simpleMovie,theData,z,wrtXyn,colors,colors*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedz2", $
      5,axislabels=wrtXyt,whichFrames=whichFramesz2,ranges=unsRanges2,horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=wrtXytex


  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsTime,[z[n_elements(time)-1]],vsTimeNames,vsTimecolors,linestyles,vsTimeToLog,expName2+"_vstimeDist", $
          5,axisLabels=vsTimeLabels,whichFrames=[0],NIndVarBins=20,thicknesses=unsThicknesses,$
          svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels,percentileList=[.025,.16,.5,.84,.975]
  simpleMovie, vsTime,[z[n_elements(time)-1]],vsTimeNames,vsTimecolors,vsTimeStyles,vsTimeToLog,expName2+"_vstime", $
      5,axisLabels=vsTimeLabels,whichFrames=[0],thicknesses=unsThicknesses,svSinglePlot=svSinglePlot,texLabels=vsTimeTexLabels
;  simpleMovie, avgdVsMdot,z,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot", $
;      5,PSYM=1,axisLabels=vsAvgMdotLabels,horizontal=1,whichFrames=whichFrames3,thicknesses=unsThicknesses
  simpleMovie, byAge,ages/1.0d9,byAgeNames,colors,colors*0,byAgeToLog,expName2+"_byAgeLogR", $
      5,axisLabels=byAgeLabels,whichFrames=[0,1,5,10],horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=byAgeTexLabels,timeText="age = "
  simpleMovie, byAge,ages/1.0d9,byAgeNames,colors,colors*0,[0,byAgeToLog[1:n_elements(byAgeToLog)-1]],expName2+"_byAge", $
      5,axisLabels=byAgeLabels,whichFrames=[0,1,5,10],horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=byAgeTexLabels,timeText="age = "
  simpleMovie, vsAge,r[positionIndices],vsAgeNames,colors,colors*0,vsAgeToLog,expName2+"_vsAge", $
      5,axisLabels=vsAgeLabels,whichFrames=[0,1,2],horizontal=1,thicknesses=unsThicknesses, $
      svSinglePlot=svSinglePlot,texLabels=vsAgeTexLabels,timeText="r = "
  IF(TheNumberOfModels GT 10) THEN $
      simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarDistH", $
          5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201], $
          NIndVarBins=5,horizontal=1,thicknesses=unsThicknesses,svSinglePlot=svSinglePlot
    


  ;; We've made some plots where each panel is a specific time and each line is a specific model. Maybe we can do the converse.
  ;; Right now unsorted ~ intervals ~ (time) (pos) (variable) (model).
  ;; The thing simpleMovie really does is ~ (pane) (pos) (plot) (line)
  IF(TheNumberOfModels LE 3) THEN BEGIN
    singlePane = dblarr(TheNumberOfModels,n_elements(theData[0,*,0,0]),n_elements(theData[0,0,*,0]),5)
    FOR j=0,TheNumberOfModels-1 DO BEGIN
      FOR ii=0,5-1 DO BEGIN
	    singlePane[j,*,*,ii] = theData[whichFrames5[ii],*,*,j]
      ENDFOR
    ENDFOR
    spColors= indgen(5)
    FOR j=1,TheNumberOfModels-1 DO spColors=[[spColors],[indgen(5)]]
    simpleMovie,singlePane,intarr(thenumberofmodels),wrtXyn,transpose(spcolors), $
        transpose(spcolors)*0,wrtXyl,expName2+"_unsortedSPlogR", $
    	5,axisLabels=wrtXyt,whichFrames=indgen(TheNumberOfModels),ranges=unsRanges, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations
    simpleMovie,singlePane,intarr(thenumberofmodels),wrtXyn,transpose(spcolors), $
        transpose(spcolors)*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_unsortedSP", $
     	5,axisLabels=wrtXyt,whichFrames=indgen(TheNumberOfModels),ranges=unsRanges2, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations

    ; each pane is a model, separate plots for z.
    fractionsPanePerModel = dblarr(TheNumberOfModels, n_elements(theData[0,*,0,0]),  n_elements(whichFrames3), 4)
    fractionsPanePerZ = dblarr(n_elements(whichFrames3), n_elements(theData[0,*,0,0]), TheNumberOfModels, 4) ; more like the standard thing
    FOR j=0,ThenumberOfModels-1 DO BEGIN
        FOR ii=0,n_elements(whichFrames3)-1 DO BEGIN
	        fractionsPanePerModel[j, *, ii, 1] = abs(theData[whichFrames3[ii], *, 31 ,j]) ; trPerAll
        	fractionsPanePerModel[j, *, ii, 2] = fractionsPanePerModel[j,*,ii,1] + theData[whichFrames3[ii], *, 36, j ]; tr+sf perAll
        	fractionsPanePerModel[j, *, ii, 3] = fractionsPanePerModel[j,*,ii,2] + theData[whichFrames3[ii], *, 35, j]

    		fractionsPanePerZ[ii,*,j,1]=fractionsPanePerModel[j,*,ii,1]
	    	fractionsPanePerZ[ii,*,j,2]=fractionsPanePerModel[j,*,ii,2]
		    fractionsPanePerZ[ii,*,j,3]=fractionsPanePerModel[j,*,ii,3]
    	ENDFOR
    ENDFOR
  ENDIF



 IF(TheNumberOfModels GT 10) THEN BEGIN
    singlePane = dblarr(n_elements(expNames),n_elements(theData[0,*,0,0]),n_elements(theData[0,0,*,0]),5)
    FOR j=0,n_elements(expNames)-1 DO BEGIN
      FOR ii=0,5-1 DO BEGIN
	    singlePane[j,*,*,ii] = intervals[whichFrames5[ii],*,*,j*nper+nper/2]
      ENDFOR
    ENDFOR
    spColors= indgen(5)
    FOR j=1,n_elements(expNames)-1 DO spColors=[[spColors],[indgen(5)]]
    simpleMovie,singlePane,intarr(n_elements(expNames)),wrtXyn,transpose(spColors), $
        transpose(spColors)*0,wrtXyl,expName2+"_intervalsSPlogR", $
    	5,axisLabels=wrtXyt,whichFrames=indgen(n_elements(expNames)),ranges=unsRanges, $
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations
    simpleMovie,singlePane,intarr(n_elements(expNames)),wrtXyn,transpose(spColors), $
        transpose(spColors)*0,[0,wrtXyl[1:n_elements(wrtXyl)-1]],expName2+"_intervalsSP", $
	    5,axisLabels=wrtXyt,whichFrames=indgen(n_elements(expNames)),ranges=unsRanges2,$
        horizontal=1,svSinglePlot=svSinglePlot,texLabels=wrtXytex,replacementText=notations

  ENDIF


  ;; Now make some standalone plots.
  ; roughly redshifts 2,1.5,1,.5,0 are [1,20,51,104,201]
;  simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstar",5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201]
;  IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarDist",5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201],NIndVarBins=5
;  simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarH",5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201],horizontal=1
  

	
  proftimes[4]=systime(1)

  PRINT,proftimes[1]-proftimes[0]," seconds to check model validities and initialize"
  PRINT,proftimes[2]-proftimes[1]," seconds to read in data"
  PRINT,proftimes[3]-proftimes[2]," seconds to sort the data"
  PRINT,proftimes[4]-proftimes[3]," seconds to generate movies"

  STOP

END
