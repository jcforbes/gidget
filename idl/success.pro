;; taken from the internet; counts number of lines in a file
;; from http://physics.nyu.edu/grierlab/idl_html_help/F12.html#wp894303
FUNCTION file_lines, filename 
   OPENR, unit, filename, /GET_LUN 
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
   RETURN,file_lines(name+"_stde.txt") EQ 0
END
