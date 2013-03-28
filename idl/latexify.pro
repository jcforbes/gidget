;;;; Modified by JF Feb 5, 2013 - added tempname option so 
;;;; that multiple instances of latexify can be run at the
;;;; same time.
;+
;PURPOSE
;	to use PSfrag to do latex syntax correctly
;SYNTAX
;	latexify, filename, tags, tex, scale
;INPUTS
;	filename: name of eps file
;	tags: name of dummy placeholder
;	tex: the tex you want to replace it
;	scale: the scale you want the latex to be [default to 1]
;NOTES:
;	requires the following latex packages properly installed:
;		geometry, graphicx, psfrag
;	requires you to use font=0
;
;	Follows from discussion by Saurav
;		http://www.astrobetter.com/idl-psfrag/
;	Writes and deletes files called x144temp.*
;EXAMPLE
;
;   SET_PLOT, 'PS'
;   !P.FONT=0
;   DEVICE, FILENAME='figure.eps'
;   plot,findgen(10),xtitle='xtitle'
;   DEVICE, /CLOSE
;
;   latexify, 'figure.eps', 'xtitle', '$M_{\odot}$'
;	
;Written by R. da Silva, UCSC, 1-4-11	
;-
pro latexify, filename, tag, tex, scale, outname=outname, example=example,$
    	height=height, width=width, full=full, tempname=tempname

if keyword_set(example) then begin
       SET_PLOT, 'PS'
       !P.FONT=0
       DEVICE, FILENAME='figure.eps'
          plot,findgen(10),xtitle='xtitle'
       DEVICE, /CLOSE
       
       filename='figure.eps'
       tag='xtitle'
       tex='Mass ($M_{\odot}$)'
       outname='out.eps'
       scale=1
endif
;eps_bounding_box,filename, coord=coord
ps_mod,filename, coord=coord
xsi=coord[2]-coord[0]
ysi=coord[3]-coord[1]



if ~keyword_set(width) then width=11.9
if ~keyword_set(height) then height=width*float(ysi)/xsi
if ~keyword_set(tempname) then tempname='x144temp'

if ~keyword_set(scale) then scale=replicate(1, n_elements(tag))
scale=rstring(scale)

if keyword_set(outname) EQ 0 then begin
	 outname=tempname+'2.eps'
	noname=1
endif
openw, lun, tempname+'.tex', /get_lun
printf, lun, '\documentclass{article}'
printf, lun, '\usepackage{geometry, graphicx, psfrag}'
printf, lun,'\pagestyle{empty}'
printf, lun,'\geometry{paperwidth='+rstring(width+0.2)+'cm,'+$
                    'paperheight='+rstring(height+0.2)+'cm,margin=0pt}'

printf, lun,'\begin{document}'
for i=0, n_elements(tag)-1 do $
    	printf, lun,'\psfrag{'+tag[i]+'}[c][]['+scale[i]+']{'+tex[i]+'}'
printf, lun,'\includegraphics[width='+$
    	rstring(width)+'cm,height='+rstring(height)+'cm]{'+filename+'}'
printf, lun,'\end{document}'
close, lun
free_lun, lun

spawn, 'latex '+tempname+'.tex', ls0
spawn, 'dvips -o '+outname+' '+tempname+'.dvi', ls1
;spawn, 'dvips -o x144temp.ps x144temp.dvi', ls1
;spawn, 'ps2epsi x144temp.ps x144temp.epsi, ls2
;spawn, "perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/' <"+$
;	 'x144temp.epsi > '+outname, ls3

;garbage cleanup
if keyword_set(noname) then begin
	spawn, 'mv '+tempname+'2.eps '+filename
	outname=filename
endif
;; This line seems like it could go very wrong, so I'll 
; manually specify the files to be removed- JF
;spawn, 'rm -f '+tempname+'*'
spawn, 'rm -f '+tempname+'.eps'
spawn, 'rm -f '+tempname+'.aux'
spawn, 'rm -f '+tempname+'.dvi'
spawn, 'rm -f '+tempname+'.blg'
spawn, 'rm -f '+tempname+'.log'
spawn, 'rm -f '+tempname+'.out'
spawn, 'rm -f '+tempname+'.bbl'
;spawn, 'rm -f '+tempname+'.tex'


ps_mod, outname, /skip

;if keyword_set(full) then ps_mod, outname, y0=600, x1=360
if keyword_set(full) then ps_mod, outname, y0=835-360*float(ysi)/xsi, x1=360

end
