
FUNCTION GetIthBin, index, THELIST, NBINS, subset=SUBSET
  IF(n_elements(subset) EQ 0) THEN subset = lindgen(n_elements(thelist))

  nmodels = N_ELEMENTS(subset)

  st = SORT( thelist([subset]))

  wh = where( st GE  index*   (float(nmodels)/float(NBins)) and $
              st LT (index+1)*(float(nmodels)/float(NBins)))
  ll = fix(float(index)*float(nmodels)/float(NBins))
  hh = fix(float(index+1)*(float(nmodels)/(float(NBins))))-1
  IF(hh LT 0) THEN hh=0
  IF(hh LT ll) THEN hh=ll


  RETURN,st[ll:hh]


END
