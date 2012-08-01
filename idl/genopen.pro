;; This is a generalization of the OPENR routine which can handle symlinks..
PRO genopen,unit,filename,_EXTRA=EXTRA_KEYWORDS

   thePath = FILE_READLINK(filename,/allow_nonsymlink)
   IF(thePath EQ "") THEN OPENR, unit, filename, /GET_LUN ELSE OPENR, unit,thePath,/GET_LUN

END


