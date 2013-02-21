;; We'd like to produce some postage stamps.
;; Use the results of variability3, and stitch together
;; N x M set of plots at a particular time.
;; I think this will just involve constructing the 
;; appropriate matrix to send to simpleMovie

FUNCTION GetStampData,fiducialName,comparisonName,whichVariables,whichRedshifts,offsetVars,comparisonNumbers
    expNames = [fiducialName, comparisonName+"a", comparisonName+"b"]

    ;; ntime x nx x nvar x nModels
    nModelsMax = 100
    theData = dblarr(n_elements(whichRedshifts),200,n_elements(whichVariables),nModelsMax)
    comparisonNumbers=intarr(n_elements(expNames))

    FOR i=0, n_elements(expNames)-1 DO BEGIN
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

    columns = 4
    rows = 6
    setct,1,0,0
    fid = "rw01"
    th=2.0
    chth=1
    cs=.3
    svSinglePlot=2
    whichRedshifts=[1,51,201]
    zs=['2','1','0']
    names=["r","col","Z"]
    labels=["r","col","Z"]
    texLabels=["$r$ (kpc)","$\Sigma (M_\odot/pc^2)$","$\log_{10}(Z/Z_\odot)$"]
    vars = [9,10,22]
    logs = [0,1,0]
    offsets=[1,1,0]
;    stampLabels = ["j","kappa_Z","kappa_z(S)","Q_GI", $
;	"Q_GIdyn","alphaMRI","mu","eta", $
;	"b","fg0","vcirc","rgauss"]
;    texStampLabels=["$r_{acc}=r_{IC}$","$\kappa_Z$","$\kappa_Z(S)$","$Q_{GI}$", $
;	"$Q_{GI,dyn}$","$\alpha_{MRI}$","$\mu$","$\eta$", $
;	"b","$f_{g,0}$","$v_{circ}$","$r_{gauss}$"]
;    comparisonExperiments=["rv02","rv03","rv04","rv05", $
;	"rv06","rv07","rv08","rv09", $
;	"rv10","rv11","rv12","rv13"]

    stampLabels=[ $
        "j","kz","kzh","QGI","rIC", $
        "QGIt","mri","mu","eta","rb", $
        "fg0","vcir","epf","zet","Qlim", $
        "fR","fRt","nRC","fcl","bet0", $
        "eff","efft","phi0","fH2m","Tgas", $
        "tSC"]
    texStampLabels=[ $
        "$\lambda$","$\kappa_Z$","$\kappa_Z(H)$","$Q_{GI}$","$r_{IC}$", $
        "$Q_{GI}(\tau)$","$\alpha_{MRI}$","$\mu$","$\eta$","$r_b$",$ 
        "$f_{g,0}$","$v_{circ}$","$\epsilon_\mathrm{ff}$","$\zeta$", "$Q_{lim}$", $
        "$f_R$","$f_R(t)$","$n$","$f_{cool}$","$\beta_0$",$
        "$\epsilon_{in}$","$\epsilon_{in}(t)$","$\phi_0$","$f_{H_2,min}$","$T_{gas}$", $
        "$t_{SC}$"]
    comparisonExperiments= $
        ["rw02","rw03","rw04","rw05","rw34", $
         "rw06","rw07","rw08","rw09","rw10", $
         "rw11","rw12","rw14","rw15","rw16", $
         "rw18","rw19","rw24","rw25","rw26", $
         "rw27","rw28","rw29","rw30","rw31", $
         "rw33"]
;    keys = indgen(12)
    keys= [0,1,3,4, 6,7,8,9, 10,11,12,13, 14,15,17,18, 19,20,22,23, 24,26,2,5]
    stampLabels = stampLabels[keys]
    texStampLabels=texStampLabels[keys]
    comparisonExperiments=comparisonExperiments[keys]

    ;; collect data from the disk.
    stampList = ptrarr(n_elements(stampLabels), /allocate_heap)
    FOR k=0, n_elements(stampLabels)-1 DO BEGIN
        theData = GetStampData(fid,comparisonExperiments[k],vars,whichRedshifts,offsets,numberOfModels)
    	stamp = {fid:fid, comp:comparisonExperiments[k], $
                 theData:theData,numberOfModels:numberOfModels}
        *(stampList[k]) = stamp
     ENDFOR


     ;; now that we have the data, let's make the plots.
     ; loop over variables
     ranges=dblarr(2,n_elements(vars))
     ranges[*,0]=[-5,45]
     ranges[*,1]=[.5,140]
     ranges[*,2]=[-1.1,.7]
     FOR j=0,n_elements(whichRedshifts)-1 DO BEGIN
	    FOR i=1,n_elements(vars)-1 DO BEGIN
            filename = "stamps_z"+zs[j]+"_"+names[i]
            FigureInit,filename,svSinglePlot,columns,rows
            modmultiplot,[columns,rows],mxtitle=labels[0],myTitle=labels[i],mxTitSize=cs*2,myTitSize=cs*2;,xgap=.001,ygap=.001
            !p.noerase=0
            FOR k=0, n_elements(stampLabels)-1 DO BEGIN
                IF(k NE 0) THEN !p.noerase = 1
                theData = (*(stampList[k])).theData
                numberOfModels = (*(stampList[k])).numberOfModels
                
                PLOT,[0],[0],/nodata,COLOR=0,BACKGROUND=255,XRANGE=ranges[*,0],YRANGE=ranges[*,i], $
                    XSTYLE=1,YSTYLE=1,ylog=logs[i],CHARSIZE=cs,CHARTHICK=chth,THICK=th,XTHICK=th, $
                    YTHICK=th,LINESTYLE=0;,ticklen=.05
                FOR exper=0,n_elements(numberOfModels)-1 DO BEGIN
                    offsetModel = 0
		            tthh = th
                    IF(exper NE 0) THEN tthh=1.0
                    IF(exper NE 0) THEN offsetModel = TOTAL(numberOfModels[0:exper-1])
                    FOR m=0,numberOfModels[exper]-1 DO OPLOT, theData[j,*,0,m+offsetModel],theData[j,*,i,m+offsetModel],COLOR=exper,THICK=tthh
		            gap =0.05; 0.2
                    xoff = gap + (1.0-2.29*gap) * (.9/double(columns) + double((k MOD columns))/double(columns))
                    yoff = 1.0 - 1.22*gap - (1.0-2.29*gap)* double(k/columns)/double(rows)
        		    XYOUTS,xoff,yoff,stampLabels[k],/normal,color=0,charsize=cs,charthick=chth
                ENDFOR

                modmultiplot
            ENDFOR
            FigureClean,filename,svSinglePlot
    	    latexify,(filename+".eps"),[labels,stampLabels],[texLabels,texStampLabels],[replicate(2.0,n_elements(labels)),replicate(1.0,n_elements(stampLabels))],tempname=("TMP_"+filename);height=8.89*rows,width=8.89*columns,tempname=("TMP_"+filename)
            modmultiplot,/default
	    spawn,"ps2pdf -dEPSCrop "+filename+".eps"
        ENDFOR
     ENDFOR
    
    set_plot,'x'
END
