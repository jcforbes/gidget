;; Given a 4-d array of the form discussed in the procedures below,
;; generate a 2-d array of the form [[min,max], [min,max], ...] for each variable
;; where min and max should encompass the entire data range (for all times and models)
FUNCTION simpleranges,data,wrtxlog
  ranges=dblarr(2,n_elements(data[0,0,*,0]))
  FOR k=0,n_elements(data[0,0,*,0])-1 DO BEGIN ;;; loop over variables
    ;; do the simplest thing: look at the min and the max in the data.
    ranges[*,k] = [MIN(data[*,*,k,*]),MAX(data[*,*,k,*])]

    ;; If the variable is supposed to be plotted logarithmically but the data 
    ;; isn't plottable that way (i.e. <=0 in some places), only plot the valid data.
    IF(wrtXlog[k] EQ 1) THEN BEGIN
      flat = data[*,*,k,*]
      ind = WHERE(flat GT 0.0, ct)
      IF(ct NE 0) THEN ranges[0,k] = MIN(flat[ind]) ELSE ranges[0,k]=.1
      IF(ct NE 0) THEN med = median(flat[ind]) ELSE med =1
      IF(ct EQ 0) THEN ranges[*,k] = [.1,1]
      fac = 1.0e7
      ;; Set a maximum dynamic range so that a few bad apples don't spoil the readability of the plot
      IF(ranges[1,k] / ranges[0,k] GT fac) THEN BEGIN
          ; If we're going to limit ourselves to some dynamic range, we have to decide 
          ; to anchor to the minimum or maximum. We try to anchor to the one containing the median value..
    	  IF(med GT ranges[1,k]/fac) THEN ranges[0,k] = ranges[1,k]/fac ELSE IF(med LT ranges[0,k]*fac) THEN ranges[1,k]=ranges[0,k]*fac ELSE ranges[0:1,k] = [med/sqrt(fac), med*sqrt(fac)]
      ENDIF

      minDynRange = 40.0
      IF(ranges[1,k] / ranges[0,k] LT minDynRange) THEN BEGIN
            IF(med GT sqrt(ranges[1,k]*ranges[0,k])) THEN ranges[1,k] = ranges[0,k]*minDynRange ELSE ranges[0,k]=ranges[1,k]/minDynRange
      ENDIF

    ENDIF

    ;; Give us a little breathing room: add 10% in each direction for each axis.
    IF(wrtXlog[k] EQ 1) THEN BEGIN
      dexRange = alog10(ranges[1,k]/ranges[0,k])
      ranges[1,k] = ranges[1,k]*(10.0 ^ (dexRange/10.0))
      ranges[0,k] = ranges[0,k]/(10.0 ^ (dexRange/10.0))

    ENDIF
    IF(wrtXlog[k] EQ 0) THEN BEGIN
      linRange = ranges[1,k]-ranges[0,k]
      ranges[1,k] = ranges[1,k] + linRange/10.0
      ranges[0,k] = ranges[0,k] - linRange/10.0
    ENDIF

 
  ENDFOR
  return,ranges
END


