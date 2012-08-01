;; taken from the internet; counts number of lines in a file
;; from http://physics.nyu.edu/grierlab/idl_html_help/F12.html#wp894303
FUNCTION filelinesLoc, filename 
;   thePath = FILE_READLINK(filename,/allow_nonsymlink)
;   IF(thePath EQ "") THEN OPENR, unit, filename, /GET_LUN ELSE OPENR, unit,thePath,/GET_LUN
   genOpen,unit,filename,/get_lun
   str = '' 
   count = 0ll 
   WHILE ~ EOF(unit) DO BEGIN 
      READF, unit, str 
      count = count + 1 
   ENDWHILE 
   FREE_LUN, unit 
   RETURN, count 
END

;; Is this run successful?
FUNCTION success, name
   dummy2=0
   dummy=filelinesLoc(name+"_stde.txt") EQ 0
   RETURN,dummy
END
