;+
;PURPOSE
;	to pad a string with zeros
;SYNTAX
;	result=padzero(arr, decimal)
;INPUTS
;	arr: an array of strings
;	decimal: to what power of ten do you want to pad [default 1]
;OUPUTS
;	result: the padded string
;EXAMPLE
;	arr='4'
;       result=padzero(arr, 1)
;	print, result
;	04
;
;Written by R. da Silva, UCSC, 4-7-2010
;-


function padzero, arr, decimal
if not keyword_set(decimal) then decimal=1

for i=0, decimal-1 do begin
wh=where(float(arr) LT 10^(decimal-i), ct)
if ct NE 0 then arr[wh]='0'+arr[wh]
endfor

return, arr
end
