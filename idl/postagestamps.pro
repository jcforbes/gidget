;; We'd like to produce some postage stamps.
;; Use the results of variability3, and stitch together
;; N x M set of plots at a particular time.
;; I think this will just involve constructing the 
;; appropriate matrix to send to simpleMovie

FUNCTION GetStampData,fiducialName,comparisonName,whichVariables,whichRedshifts,offsetVars,comparisonNumbers,third
    expNames = [fiducialName, comparisonName+'a', comparisonName+'b']
    IF(third) THEN expNames=[fiducialName,comparisonName+'a', comparisonName+'b', comparisonName+'c']

    ;; ntime x nx x nvar x nModels
    nModelsMax = 100
    theData = dblarr(n_elements(whichRedshifts),200,n_elements(whichVariables),nModelsMax)
    comparisonNumbers=intarr(n_elements(expNames))

    FOR i=0, n_elements(expNames)-1 DO BEGIN
        PRINT,"Loading experiment ",expNames[i]
        nameList = ListValidModels(expNames[i],nModelsMax)
        comparisonNumbers[i] = n_elements(nameList)
        FOR j=0, n_elements(nameList)-1 DO BEGIN
            modelPtr = ReadModels([nameList[j]])
            model = (*(modelPtr[0]))
            offsetModel = 0
            IF(i NE 0) THEN offsetModel = TOTAL(comparisonNumbers[0:i-1])
            theData[*,*,*,j+offsetModel] = model.dataCube[whichRedshifts,*,whichVariables+offsetVars*model.ncolstep-1]
            ptr_free,modelPtr[0]
        ENDFOR

    ENDFOR


    theData = theData[*,*,*,0:TOTAL(comparisonNumbers)-1]

    RETURN, theData

END



PRO postagestamps

    compile_opt idl2
    set_plot,'x'
    device,decomposed=0
    device,retain=2

    columns = 6
    rows = 4
    setct,1,0,0
    fid = 'rg01'
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2
    whichRedshifts=[201]
    zs=['0']
    ; Radius, Col, Sig, Z, sfrPerAccr, ageAtz0, equilibrium, colPerCrit, BB4, Q
    dummy = GetLabels(keys=[0,1,2,6,23,12,47,49,27,8], names=names, labels=labels, tex=texLabels, base=vars, log=logs, offset=offsets,ranges=ranges)

        ranges[0,0] = -.7
        ranges[1,0] = 39.99
        ranges[0,1] = 2.0 ; min column dens of 2 solar masses/yr
        ranges[1,3] = 0.6 ; max metallicity [Z/Zsun]

    stampLabels=[ $
        'racc','kz','kzh','QGI','rIC', $
        'QGIt','mri','mu','eta','rb', $
        'fg0','vcir','epf','xi','Qlim', $
        'fR','fRt','nRC','fcl','bet0', $
        'eff','efft','phi0','fH2m','Tgas', $
        'tSC','ZIGM','Mh0','Mh0+','racct', $
        'aAcc','yrec']
    texStampLabels=[ $
        '$r_{acc}$','$\kappa_Z$','$\kappa_Z(H)$','$Q_{GI}$','$r_{IC}$', $
        '$Q_{GI}(\tau)$','$\alpha_{MRI}$','$\mu$','$\eta$','$r_b$',$ 
        '$f_{g,0}$','$v_{circ}$','$\epsilon_\mathrm{ff}$','$\xi$', '$Q_{lim}$', $
        '$f_R$','$f_R(t)$','$n$','$f_{cool}$','$\beta_0$',$
        '$\epsilon_{in}$','$\epsilon_{in}(t)$','$\phi_0$','$f_{H_2,min}$','$T_{gas}$', $
        '$t_{SC}$','$Z_{IGM}$','$M_{h,0}$','$M_{h,0}(\gamma)$','$r_{acc}(t)$', $
        '$\alpha_{r}$','$y$']
    comparisonExperiments= $
        ['rg02','rg03','rg04','rg05','rg34', $
         'rg06','rg07','rg08','rg09','rg10', $
         'rg11','rg12','rg14','rg15','rg16', $
         'rg18','rg19','rg24','rg25','rg26', $
         'rg27','rg28','rg29','rg30','rg31', $
         'rg33','rg38','rg39','rg40','rg53', $
         'rg54','rg56']
     lowRange= $
         ['1-6.6 kpc', '0.1-1x fid.', '0.1-1 x fid', '1.3-1.9', '1.0-3.2', $
          '1.3-1.9','0-.009','0.1-0.4','0.5-1.5','0-2.5', $
          '0.1-0.45','160-215','.0031-.0095','-0.3:-0.1','1.8-2.4', $
          '0.4-0.5','0.4-0.5','1.0-1.9','.20-.55','0 - .45', $
          '.003-.2','.003-.2','.5-.95','.003-.029','100-6900', $
          '1-1.9', '1/100 - 1/10', '1e10-1e12','1e10-1e12','1-6.6 kpc', $
          '0-.3','.05-.052']
     hiRange= $
         ['9.9-20.8 kpc', '1-10 x fid.', '1-10 x fid.', '2.1-3.0', '3.8-20.8', $
          '2.1-3.0','.011-0.5','0.6-1.5','1.5-4.5','3.5-10', $
          '.55-.7','225-250','.0105-.031','0.1-1.0','2.6-3.0', $
          '0.6-0.7','0.6-0.7','2.1-3.0','.65-1.0','.55-.95', $
          '.31-.46','.31-.46','1.05-2','.031-.3','7100-30000', $
          '2.1-4', '1/10-1','1e12-1e13','1e12-1e13','9.9-20.8 kpc', $
          '.35-1','.056-.07']
      exRange = REPLICATE('',n_elements(hiRange))
      exRange[19] = '-0.5 - 0'
;    keys= [ $
;        0,1,3,4,6,7, $
;        8,9,10,11,12,13, $
;        14,15,17,18,19,20, $
;        22,23,24,26,27,28 $
;        ]
    keys = [ $
        30, 19, 12, 26,  3, 20, $
        10, 17, 23, 13,  6, 27, $
        18,  9, 25,  1,  8, 15, $
        22, 11,  7, 31, 24,  0  $
        ]
    stampLabels = stampLabels[keys]
    texStampLabels=texStampLabels[keys]
    comparisonExperiments=comparisonExperiments[keys]
    lowRange=lowRange[keys]
    hiRange=hiRange[keys]
    exRange=exRange[keys]

    ;; collect data from the disk.
    stampList = ptrarr(n_elements(stampLabels), /allocate_heap)
    FOR k=0, n_elements(stampLabels)-1 DO BEGIN
        theData = GetStampData(fid,comparisonExperiments[k],vars,whichRedshifts,offsets,numberOfModels,(exRange[k] NE ''))
    	stamp = {fid:fid, comp:comparisonExperiments[k], $
                 theData:theData,numberOfModels:numberOfModels}
        *(stampList[k]) = stamp
     ENDFOR

    FOR a=0, 100 DO BEGIN

	     makeThePostageStampPlots, stampList, whichRedshifts,vars,zs,th, $
             labels,stampLabels,texLabels,texStampLabels, $
             lowRange, hiRange, exRange, $
             names,columns,rows,svSinglePlot,cs,chth,logs,ranges
         stop
    ENDFOR

END

