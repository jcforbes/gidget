;; This is the same as variability2 (forked August 18, 2012), except that I'll
;; attempt to add support for an arbitrary number of data arrays.



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
;;;   data for each run can become cumbersome. If you input an empty array,
;;;   by default the code will analyze all the variables
;;; N controls how many runs in the experiment to analyze. For large number 
;;;   runs, sometimes it's useful to just look at a few runs. If N is larger 
;;;   than the number of runs in the experiment, all the runs are analyzed.
;;; sv- behaves like other instances of sv; 0 to just show the plots on the 
;;;   screen, 2 to save them as png's and subsequently movies. 4 saves 
;;;   movies in a format conducive to showing on a low-resolution projector.
PRO variability3,expNames,keys,N,sv
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

  ;; variables to plot- title, column index, log plot, name, y-range
  wrtXyt=['R (kpc)','Gas Column Density (Msun/pc^2)','Velocity Dispersion (km/s)','Stellar Column Density (Msun/pc^2)','Stellar Velocity Dispersion (km/s)','Gas Fraction','[Z/Zsun]','SFR Column Density (Msun/yr/kpc^2)','Q','Gas Q','Stellar Q','Molecular Fraction','Age At z=0 (Ga)','Age (Ga)','Depletion Time (Ga)','Viscous Time (Ga)','Depletion/Viscous Time','Mass Flux Through Disk (Msun/yr)']
  wrtXyy=[9,10,11,12,13,17,22,14,12,24,23,48,31,32,35,36,37,39] ;; index
  wrtXyp=[1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1,0] ;; offset? 
  wrtXyn=['r', 'col','sig','colst','sigst','fg','Z','colsfr','Q','Qg','Qst','fH2','ageAtz0','age','tdep','tvisc','tdepOverTvisc','mdotDisk']
;  wrtXyr=[[0,20],[.1,1000],[5,300],[.1,1000],[5,300],[0,1],[-1,1.0],$
;    [1d-6,10],[1.5,5.0],[.5,50],[.5,50],[0,1],[1d9,14d9],[1d5,14d9],[1.0e8,1.0e10],[1.0e8,1.0e10],[.01,100]]
  wrtXyl=[0,1,1,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1] ;; log plot 18

  IF(n_elements(keys) EQ 0) THEN keys=indgen(n_elements(wrtXyt)-1)+1

  dummy = where( keys EQ 19, n1)
  IF(n1 EQ 2) THEN n1=1

  ;; only plot the variables present in 'keys'	
  tk = [[0],keys]
  wrtXyt=wrtXyt[tk]
  wrtXyy=wrtXyy[tk]
  wrtXyp=wrtXyp[tk]
  wrtXyn=wrtXyn[tk]
  ;wrtXyr=wrtXyr[*,tk]
  wrtXyl=wrtXyl[tk]  

  ;; now loop over every model in the list - 
  ;; all we're doing here is reading in each model then freeing 
  ;; the associated memory
  ;; so that we can iterate over a very large number. For each model, 
  ;; record the data specified in the arrays above (e.g. yy)
  modelPtr=ReadModels([nameList2[0]])
  model = (*(modelPtr[0]))
  theNTS = model.NOutputsNominal+2

  time = model.evArray[1,*]  
  z = model.evArray[10-1,*]
  

  ;; (# timesteps) x (nx) x (# of columns) x (# models)
  theData= dblarr(theNTS,n_elements(model.dataCube[0,*,0]),n_elements(wrtXyy),n_elements(nameList2))

  NVS = 41
  ;; (# timesteps) x (# pts/model frame = 1) x (one thing vs another = 20) x (# models)
  vsMdot = dblarr(theNTS,1,NVS,n_elements(nameList2))
  nameList4 = strarr(n_elements(nameList2))


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

       ;; if this model has the appropriate number of time steps..
       ;;PRINT,"Number of time steps",n_elements(model.dataCube[*,0,0]),n_elements(theData[*,0,0,0])
       ;IF(n_elements(model.dataCube[*,0,0]) EQ n_elements(theData[*,0,0,0])) THEN BEGIN
       IF(n_elements(model.dataCube[*,0,0]) EQ theNTS) THEN BEGIN
         FOR k=0,n_elements(wrtXyy)-1 DO BEGIN	
           theData[*,*,k,ctr] = model.dataCube[*,*,wrtXyy[k]+offsets[wrtXyp[k]]-1]
         ENDFOR

;         vsMdot[*,0,0,ctr] =  model.dataCube[*,model.nx-1,39-1]  ;; external accretion rate
         
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
           vsMdot[zi,0,15,ctr] = alog10(modelInfo[16]) ;; metallicity
           vsMdot[zi,0,16,ctr] = modelInfo[5] ;; fG in SF region.
           vsMdot[zi,0,17,ctr] = model.evArray[20-1,zi] ;; (modelInfo[18] - modelInfo[11])/(modelInfo[18]+1.0d-10) ;; fraction of SF in the bulge
           vsMdot[zi,0,18,ctr] = modelInfo[20] ;; fgL
           vsMdot[zi,0,19,ctr] = alog10(modelInfo[21]) ;; ZL
           vsMdot[zi,0,20,ctr] = modelInfo[15] ;; fg
	   vsMdot[zi,0,21,ctr] = alog10(model.evArray[5-1,zi]/.02) ;; ZBulge
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

  vsMdotLabels = ["Mdot (Msun/yr)","DiskSFR (Msun/yr)","Stellar Mass (Msun)","Peak Radius of Column Density (kpc)","Radius of HI transition (kpc)","Column Density at r=8 kpc (Msun/pc^2)","Specific SFR (yr^-1)","Bulge Gas Mass (Msun)","Bulge Stellar Mass (Msun)","H_2 Fraction","Mdot Gas into Bulge (Msun/yr)","Mdot Stars into Bulge (Msun/yr)","Bulge to Total Ratio","SFR Including Gas Flux Into Bulge (Msun/yr)","efficiency (Mstar / (f_b M_h))","[Z/Zsun]","Gas Fraction in SF Region","Mdot (Msun/yr)","Gas Fraction in Optical Region","[Z/Zsun] in Optical Region","Gas Fraction","[Z Bulge/Zsun]","Average inward gas radial velocity (km/s)","Inner Edge of GI Region (kpc)","Outer Edge of GI Region (kpc)","Average velocity dispersion (km/s)","Average radial velocity for non-outward velocities only (km/s)","Depletion Time (Ga)","Gas Sc Length (kpc)","St Sc Length (kpc)","sersic","BT Meas","gas chi^2","st chi^2","Half SFR radius (kpc)","Half gas radius (kpc)","Half stellar radius (kpc)","Central Density (Msun/pc^2)","specific J_*","specific J_g","specific J_out"]
  vsMdotTexLabels = ["$\dot{M} (M_\odot/yr)$","Disk SFR ($M_\odot$/yr)","Stellar Mass ($M_\odot$)","Peak Radius (kpc)","Radius of HI trans. (kpc)","$\Sigma$ at $r=8$ kpc ($M_\odot/pc^2$)","sSFR ($yr^{-1}$)","Bulge Gas Mass ($M_\odot$)","Bulge Stellar Mass ($M_\odot$)","$H_2$ Fraction","$\dot{M}_g$ into Bulge ($M_\odot$/yr)","$\dot{M}_*$ into Bulge ($M_\odot$/yr)","Bulge:Total Ratio","SFR + $\dot{M}_{g,\rightarrow\mathrm{bulge}}$ ($M_\odot$/yr)","$\epsilon$ ($M_*$ / ($f_b M_h$))","$[Z/Z_\odot]$","$f_g$ in SF Region","$\dot{M}$ ($M_\odot$/yr)","$f_g$ in Optical Region","$[Z/Z_\odot]$ in Optical Region","$f_g$","$[Z_\mathrm{bulge}/Z_\odot]$","Average $-v_r$ (km/s)","Innermost GI Region (kpc)","Outermost GI Region (kpc)","$\langle\sigma\rangle$ (km/s)","Average $-v_r\ge 0$ (km/s)","$t_{dep} = M_g f_{H_2} / \dot{M}^{SF}$ (Ga)","Gas Sc Length (kpc)","St Sc Length (kpc)","sersic","BT Meas","gas $\chi^2$","st $\chi^2$","$r_{\mathrm{SFR},\frac12}$ (kpc)","$r_{g,\frac12}$ (kpc)","$r_{*,\frac12}$","$\Sigma_0\ (M_\odot/pc^2)$","$J_*/M_*$","$J_g/M_g$","$J_{\mathrm{out}}/M_\mathrm{out}$"]

  vsMdotNames  = ["mdot","DiskSFR","Mst","rPeak","rHI","colsol","sSFR","BulgeGasMass","BulgeStMass","fH2","mdotBulgeGas","mdotBulgeStars","BT","SFRplusMdotIn","efficiency","Z","fgInSF","mdot","fgL","ZL","fg","ZBulge","vrAvg","rQin","rQout","sigAvg","vrgGtr0","tdepAvg","GasScLength","StScLength","sersic","BTMeas","gaschi2","stchi2","rHalfSFR","rHalfGas","rHalfSt","CentralDensity","spJSt","spJg","spJout"]
  vsMdotStrt =  [1,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  vsMdotToLog = [1,1,1,1,0,1,1,0,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1] 

  modelcounter = modelcounter2[*]
  nameList2 = nameList4[*]

  proftimes[2]=systime(1)

  theData = (temporary(theData))[*,*,*,0:ctr-1]
  colors = intarr(n_elements(theData[*,0,0,0]),n_elements(nameList2))

  intervals = dblarr(n_elements(theData[*,0,0,0]),  n_elements(theData[0,*,0,0]), n_elements(theData[0,0,*,0]), 3*n_elements(expNames))  
  ;; having accumulated theData, let's do some analysis!
;  sortdData = theData * 0.0
;  sortdData[*,*,0,*] = theData[*,*,0,*] ;; don't sort radii..
  ;sortdVsMdot = vsMdot*0.0
  FOR expInd=1, n_elements(expNames) DO BEGIN
    low = modelcounter[expInd-1]
    high= modelcounter[expInd]-1
    FOR j=0,n_elements(theData[*,0,0,0])-1 DO BEGIN ;; loop over time
      FOR k=1,n_elements(wrtXyy)-1 DO BEGIN ;; loop over variables
        FOR xx=0,n_elements(theData[0,*,0,0])-1 DO BEGIN ;; loop over x-coord
	  ;;; Sort the models, within each experiment
          intervals[j,xx,k,(expInd-1)*3:(expInd)*3-1] = percentiles(theData[j,xx,k,low:high],[.05,.5,.95],wrtXlog[k])

;	  srtd = sort(theData[j,xx,k,low:high])
;          sortdData[j,xx,k,low:high] = theData[j,xx,k,low+srtd]
          IF(xx EQ 4 AND j EQ n_elements(theData[*,0,0,0])-1) THEN BEGIN
          ENDIF
	ENDFOR ;; end loop over x-coord 
      ENDFOR ;; end loop over variables
    ENDFOR ;; end loop over time
  ENDFOR ;; end loop over experiments



  proftimes[3]=systime(1)

  FOR expInd=1, n_elements(expNames) DO BEGIN
    low=modelcounter[expInd-1]
    high=modelcounter[expInd]-1
    delta = high-low
    intervals[*,*,*, (expInd-1)*3+0]=sortdData[*,*,*, fix(low+.025*delta)] ;; fix(low+.025*(high-low))
    intervals[*,*,*, (expInd-1)*3+1]=sortdData[*,*,*, fix(low+.5*delta)]
    intervals[*,*,*, (expInd-1)*3+2]=sortdData[*,*,*, fix(low+.975*delta)]
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
      thicknesses[GetIthBin(nc,Mh0,ncolors)] = 1.0 + .5*float(nc) ;; thicker lines have higher Mh0.
  ENDFOR




  expName2 = expNames[0]
  IF(n_elements(expNames) LT 5) THEN BEGIN
      FOR i=1, N_elements(expNames)-1 DO BEGIN
          expName2 = expName2 + '_' + expNames[i];; +'_'+strcompress(string(MIN([n_elements(nameList2),N])),/remove) ;; e.g. rk5_113
      ENDFOR
  ENDIF ELSE expName2 +='_' + strcompress(string(n_elements(expNames)),/remove)
  expName2 = expName2 + '_' + strcompress(string(MIN([n_elements(nameList2),N*n_elements(expNames)])),/remove) ;; e.g. rk5_rk6_113



  IF(n_elements(nameList2) GT 50) THEN BEGIN
      setct,1,0,0
      linestyles = intarr(n_elements(expNames)*3)+2
      midpoints = where(indgen(n_elements(expNames)*3) MOD 3 EQ 1)
      linestyles[midpoints] = 0
      intervalsColors0 = indgen(n_elements(expNames)*3)/3
      intervalsColors = intarr(n_elements(time),n_elements(intervalsColors0))
      FOR tt=0,n_elements(time)-1 DO intervalsColors[tt,*] = intervalsColors0 
      simpleMovie,intervals,z,wrtXyn,intervalsColors,linestyles,wrtXyl,expName2+"_intervals",sv,axislabels=wrtXyt
  ENDIF



  setct,3,n_elements(theData[0,0,0,*]),0
;  IF(n_elements(nameList2) LT 105) THEN BEGIN 	
    simpleMovie,theData,z,wrtXyn,colors,colors*0,wrtXyl,expName2+"_unsorted",sv,axislabels=wrtXyt
;  ENDIF
	

  setct,3,n_elements(vsMdot[0,0,0,*]),0
  IF(n1 EQ 1) THEN BEGIN
;    simpleMovie, vsMdot,z,vsMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsmdot",sv,vsMdotStrt,PSYM=1,prev=1,axisLabels=vsMdotLabels
  ENDIF

  ;; average mdot over a number of time steps nt
  nt = 40
  setct,3,n_elements(vsMdot[0,0,0,*]),0
  avgdVsMdot = vsMdot[*,*,*,*]

  FOR ti=0,n_elements(time)-1 DO BEGIN
    FOR mm=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
      avgdVsMdot[ti,0,0,mm] = TOTAL(vsMdot[MAX([0,ti-nt]):ti,0,0,mm]) / (ti - MAX([0,ti-nt]) + 1)
    ENDFOR
  ENDFOR
  IF(n1 EQ 1) THEN BEGIN
     vsAvgMdotLabels = vsMdotLabels[*]
     vsAvgMdotLabels[0] = "Average Mdot Over Past 40 Outputs (Msun/yr)"
     vsAvgMdotTexLabels = vsMdotTexLabels[*]
     vsAvgMdotTexLabels[0] = "$\dot{M}_{40 \mathrm{outputs}} (M_\odot/yr)$"
     vsAvgMdotNames = vsMdotNames[*]
     vsAvgMdotNames[0] = "avgMdot"
     simpleMovie, avgdVsMdot,z,vsAvgMdotNames,colors,colors*0,vsMdotToLog,expName2+"_vsavgmdot",sv,PSYM=1,prev=1,axisLabels=vsAvgMdotLabels
  ENDIF


  IF(n1 EQ 1) THEN BEGIN
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
    simpleMovie, vsSpecificSFR,z,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFR",sv,PSYM=1,prev=1,axisLabels=vssSFRLabels,taillength=2
    IF(n_elements(vsSpecificSFR[0,0,0,*]) GT 50) THEN simpleMovie, vsSpecificSFR,z,vssSFRNames,colors,colors*0,vssSFRToLog,expName2+"_vsSpSFRDist",sv,PSYM=3,prev=0,axisLabels=vssSFRLabels,NIndVarBins=5

  ENDIF    


  ;; Plot these quantities vs. the star formation rate.
  setct,3,n_elements(vsMdot[0,0,0,*]),0
  vsSFR = vsMdot[*,*,1:n_elements(vsMdot[0,0,*,0])-1,*] ;; cut out the mdot slice of the array
  IF(n1 EQ 1) THEN BEGIN 

    vsSFRLabels = vsMdotLabels[1:(n_elements(vsMdotLabels)-1)]
    vsSFRTexLabels = vsMdotTexLabels[1:n_elements(vsMdotTexLabels)-1]
    vsSFRNames = vsMdotNames[1:(n_elements(vsMdotNames)-1)]
    vsSFRToLog = vsMdotToLog[1:(n_elements(vsMdotToLog)-1)]
    vsSFRstrt= vsMdotstrt[1:(n_elements(vsMdotStrt)-1)]

;    simpleMovie, vsSFR, z,vsSFRNames, colors, colors*0, vsSFRtoLog,expName2+"_vsSFR",sv,vsSFRstrt,PSYM=1,prev=1,axisLabels=vsSFRLabels
  ENDIF


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


  setct,3,n_elements(vsMstar[0,0,0,*]),0
  IF(n1 EQ 1) THEN BEGIN
    simpleMovie, vsMstar, z,vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstar",sv,PSYM=1,prev=0,axisLabels=vsMstarLabels
    IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, z,vsMstarNames, colors, colors*0, vsMstarToLog,expName2+"_vsMstarDist",sv,PSYM=3,prev=0,axisLabels=vsMstarLabels,NIndVarBins=5,percentiles=[.17,.5,.83]

  ENDIF



  setct,3,n_elements(vsMdot[0,0,0,*]),0
  vsTime = vsMdot[*,*,*,*]
  vsTimeLabels = vsMdotLabels[*]
  vsTimeLabels[0] = "Time Since z=2 (Ga)"
  vsTimeTexLabels = vsMdotTexLabels[*]
  vsTimeTexLabels[0] = "Time since $z=2$ (Ga)"
  vsTimeNames = vsMdotNames[*]
  vsTimeToLog = vsMdotToLog[*]
  vsTimeToLog[0] = 0
  vsTimeStrt =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]

  ;; for each model, set the independent variable equal to the time: 
  FOR j=0, n_elements(vsTime[0,0,0,*]) -1 DO BEGIN
    vsTime[*,0,0,j] = time[*]
  ENDFOR
  IF(n1 EQ 1) THEN BEGIN
    simpleMovie, vsTime,z,vsTimeNames,colors,colors*0,vsTimeToLog,expName2+"_vstime",5,PSYM=3,prev=1,taillength=1000,axisLabels=vsTimeLabels,whichFrames=[n_elements(z)-1]
  ENDIF


  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2,sv=4,nt0=2,thicknesses=thicknesses
  ComputeCorrelations,vsmdot,colors,time,vsMdotLabels,vsMdotNames,expName2+"_log",sv=4,nt0=2,thicknesses=thicknesses,logarithms=vsMdotToLog

  ;; Now make some standalone plots.
  ; roughly redshifts 2,1.5,1,.5,0 are [1,20,51,104,201]
  simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstar",5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201]
  IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarDist",5,PSYM=3,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201],NIndVarBins=5
  simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarH",5,PSYM=1,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201],horizontal=1
  IF(n_elements(vsMstar[0,0,0,*]) GT 50) THEN simpleMovie, vsMstar, z, vsMstarNames, colors, colors*0, vsMstarToLog, expName2+"_vsMstarDistH",5,PSYM=3,prev=0,axisLabels=vsMstarLabels,texLabels=vsMstarTexLabels,whichFrames=[1,51,201],NIndVarBins=5,horizontal=1
  

  Mh = vsMdot[*,0,2,*]/(.18*vsMdot[*,0,14,*])
  zn = n_elements(Mh[*,0,0,0])
	
  ;; Plot the accretion rates.
  ;SETCT,3,n_elements(vsMdot[0,0,0,*]),0
  FIGUREInit,(expName2+'_accRates'),1,3,3
  yr=[MIN(vsMdot[zn-1,0,0,*])*.5,MAX(vsMdot[0,0,0,*])*2];[.01,1000];[MIN(sortdVsMdot[*,0,0,*]),MAX(sortdVsmdot[*,0,0,*])]
  Plot,[time[0]],[vsMdot[0,0,0,0]],COLOR=0,BACKGROUND=255,XRANGE=[MIN(time),MAX(time)],YRANGE=yr,YLOG=1,XSTYLE=1,XTITLE="Time",YTITLE="Accretion Rate",CHARSIZE=2,THICK=2,CHARTHICK=1.5
  FOR m=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
    ;setct,5,n_elements(vsMdot[0,0,0,*]),m
    OPLOT,time,vsMdot[*,0,0,m],COLOR=colors[0,m],THICK=thicknesses[m]
	
  ENDFOR
  FIGURECLEAN,(expName2+'_accRates'),1



  ;; Plot the halo mass
  ;SETCT,3,n_elements(vsMdot[0,0,0,*]),0
  FIGUREInit,(expName2+'_Mh'),1,3,3
  ;Plot,time,Mh[*,0,0,0],COLOR=0,BACKGROUND=255,YRANGE=[Min(Mh[zn-1,0,0,*])*.1,Max(Mh[zn-1,0,0,*])*2],YLOG=1,XSTYLE=1,XTITLE="Time",YTITLE="Halo Mass",CHARSIZE=3,THICK=thicknesses[0],CHARTHICK=2
  PLOT,[time[0]],[Mh[0,0,0,0]],COLOR=0,BACKGROUND=255,XRANGE=[MIN(time),MAX(time)],YRANGE=[MIN(Mh[zn-1,0,0,*])*.1,MAX(Mh[zn-1,0,0,*])*2],YLOG=1,XSTYLE=1,YSTYLE=1,XTITLE="Time (Ga)",YTITLE="Halo Mass",CHARSIZE=2,THICK=2,CHARTHICK=1.5
  FOR m=0,n_elements(vsMdot[0,0,0,*])-1 DO BEGIN
    ;setct,5,n_elements(vsMdot[0,0,0,*]),m
    OPLOT,time,Mh[*,0,0,m],COLOR=colors[0,m],THICK=thicknesses[m]
	
  ENDFOR
  FIGURECLEAN,(expName2+'_Mh'),1



;  FIGUREInit,(expName2+'_z0AccRates'),1,1,1
;  cgHistoplot, alog10(vsMdot[n_elements(vsMdot[*,0,0,0])-1,0,0,*]), BinSize=.1
;  FIGURECLEAN,(expName2+'_z0AccRates'),1
	
  proftimes[4]=systime(1)

  PRINT,proftimes[1]-proftimes[0]," seconds to check model validities and initialize"
  PRINT,proftimes[2]-proftimes[1]," seconds to read in data"
  PRINT,proftimes[3]-proftimes[2]," seconds to sort the data"
  PRINT,proftimes[4]-proftimes[3]," seconds to generate movies"

  STOP

END
