;+
;PURPOSE
;	to edit the bounding box of a ps file
;SYNTAX
;	ps_mod, filename[, x0=x0, y0=y0, x1=x1, y1=y1, $
;	       dy0=dy0, dy1=dy1, dx0=dx0, dx1=dx1, /high]
;INPUTS
;	filename: name of file
;	/high: set for highest quality PDFs
;  Bounding Box Keywords
;	x0: x position of lower left hand corner
;	y0: y position of lower left hand corner
;	x1: x position of upper right hand corner
;	y1: y position of upper_right hand corner
;       dx0: offset for x position of lower left hand corner
;       dy0: offset for y position of lower left hand corner
;       dx1: offest for x position of upper right hand corner
;       dy1: offest for y position of upper_right hand corner
;	/skipsize: use to delete line defining paper size
;  
;Written by R. da Silva, UCSC, 9-23-10
;-

pro ps_mod, filename, x0=x0, y0=y0, x1=x1, y1=y1, $
	dy0=dy0, dy1=dy1, dx0=dx0, dx1=dx1, skipsize=skipsize, $
	coords=coords, high=high
fluff='x1123123213x'

len=strlen(filename)
if strmatch(strmid(filename, len-2, 2), 'ps') EQ 0 then begin
  splog, 'FILE MUST BE PS or EPS'
endif else begin
if strmatch(strmid(filename, len-3,3), 'eps') then iseps=1 else iseps=0
pslen=(iseps?3:2)
psname=(iseps?'.eps':'.ps')

filename2=strmid(filename, 0, len-1-pslen)+fluff+psname
lin=''
openr, lun, filename, /get_lun
openw, lun2, filename2, /get_lun
;read first line
readf, lun, lin
printf, lun2, lin
;read bounding box line
readf, lun, lin
while strmatch(lin, '%%BoundingBox*') EQ 0 do readf, lun, lin
;edit it accordingly
split2=strsplit(lin, ':', /extract)
split3=strsplit(split2[1], /extract)
split3=long(split3)
if keyword_set(x0) then split3[0]=x0
if keyword_set(y0) then split3[1]=y0
if keyword_set(x1) then split3[2]=x1
if keyword_set(y1) then split3[3]=y1

if keyword_set(dx0) then split3[0]+=dx0
if keyword_set(dy0) then split3[1]+=dy0
if keyword_set(dx1) then split3[2]+=dx1
if keyword_set(dy1) then split3[3]+=dy1

newline=split2[0]+': '+rstring(split3[0])+' '+rstring(split3[1])+' '+$
	rstring(split3[2])+' '+rstring(split3[3])
coords=split3
;print, newline
;print the new line
printf, lun2, newline

;read rest of lines
while ~EOF(lun) do begin
  readf, lun, lin
  if strmatch(lin, '%%DocumentPaperSizes*') then continue
  if strmatch(lin, '%%EndProlog*') AND keyword_set(high) then begin
	printf, lun2, 'systemdict /setdistillerparams known {'
	printf, lun2, '<< /AutoFilterColorImages false /ColorImageFilter /FlateEncode >> setdistillerparams'
	printf, lun2, '} if'

  endif

  printf, lun2, lin
endwhile
close, lun2
close, lun
free_lun, lun, lun2

spawn, 'mv '+filename2+' '+filename

endelse
end
