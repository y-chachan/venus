;; make an image from volimage

wav = 1000.0                      ; wav to observe at, this is currently irrelevant

;; get model file
fs = 'volimage.dat'

;; radius and simple P(r) array
nr = 100
rs = linarrw(min=0.01,max=3,n=nr)
bigpr = FLTARR(nr)+1.
save,filename='temppr.vwz.xdr',bigpr,wav,rs
   
;; set up and make the image, most of this is ignored and re-set below
ptr1   = createptr1(modelfile=fs,skiplines=0,fluxfile='temppr.vwz.xdr',band=0, $
                          n1=1,n2=1,dxy=1,thz1=0.0,thx=0.0,thz2=0.0)

;; now fiddle values to make it the right kind of observation
(*ptr1).viewinfo.obstype = 4    ; aitoff

;; satellite position				  A 	  B
(*ptr1).viewinfo.asat = 1.0    ; a   0.96   1.042
(*ptr1).viewinfo.esat = 0.     ; ecc  0.005   0.04
(*ptr1).viewinfo.isat = 0.     ; inc  0.12   0.29
(*ptr1).viewinfo.omsat = 0.     ; asc node  214   336
(*ptr1).viewinfo.wbarsat = 0.   ; wbar (asc+peri)  308   123 or 483
(*ptr1).viewinfo.longsat = 0.  ; longitude

;; viewing info
(*ptr1).viewinfo.elongmin = 0.   ; minimum elongation
(*ptr1).viewinfo.elongmax = 360. ; maximum elongation

(*ptr1).viewinfo.nstep = 200.    ; integration steps
(*ptr1).viewinfo.n1 = 500.       ; x res
(*ptr1).viewinfo.n2 = 500.       ; y res
(*ptr1).viewinfo.band = 0.       ; observation band (set to 0)

;; restrict area for lat vs. long if needed
 (*ptr1).viewinfo.tsmin = (*ptr1).viewinfo.longsat + 180 - 70    ; longitude min
 (*ptr1).viewinfo.tsmax = (*ptr1).viewinfo.longsat + 180 + 70   ; longitude max
 (*ptr1).viewinfo.psmin = -40.   ; latitude min
 (*ptr1).viewinfo.psmax = 40.    ; latitude max

tempwav1  = readmodel(ptr1) 
losint,ptr1,/noprint,/noplot,CALLC=1
im = (*(*ptr1).ptr_intop).inu

imdisp,im
save,im,filename='res_disk_sv.sav'

end
