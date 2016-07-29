;; make a scattered light image from a volimage file, using losint

pro ring_image,filein,fileout,a,e,i,asc,wbar,meanlong,nlon=nlon,nlat=nlat,lon1=lon1,lon2=lon2,lat1=lat1,lat2=lat2,nstep=nstep,sat=sat

  IF ~ KEYWORD_SET(nstep) then nstep = 200
          
  ;; get model file
  filein = filein

  ;; radius and simple P(r) array, STEREO observe scattered light
  nr = 100
  rs = linarrw(min=0.01,max=3,n=nr)
  bigpr = FLTARR(nr)+1.
  wav = 1000.0                  ; wav to observe at, this is currently irrelevant
  save,filename='temppr.vwz.xdr',bigpr,wav,rs
   
  ;; set up and make the image, most of this is ignored and re-set below
  ptr1   = createptr1(modelfile=filein,skiplines=0,fluxfile='temppr.vwz.xdr',band=0, $
                      n1=1,n2=1,dxy=1,thz1=0.0,thx=0.0,thz2=0.0)

  ;; now fiddle values to make it the right kind of observation
  (*ptr1).viewinfo.obstype = 4  ; aitoff

  ;; satellite position
  (*ptr1).viewinfo.asat = a           ; a
  (*ptr1).viewinfo.esat = e           ; ecc
  (*ptr1).viewinfo.isat = i           ; inc
  (*ptr1).viewinfo.omsat = asc        ; asc node
  (*ptr1).viewinfo.wbarsat = wbar ; wbar (asc+peri)
  (*ptr1).viewinfo.longsat = meanlong ; mean longitude
                 
  ;; viewing info
  (*ptr1).viewinfo.elongmin = 0.   ; minimum elongation
  (*ptr1).viewinfo.elongmax = 360. ; maximum elongation

  (*ptr1).viewinfo.nstep = 200  ; integration steps
  (*ptr1).viewinfo.n1 = nlon    ; x res
  (*ptr1).viewinfo.n2 = nlat    ; y res
  (*ptr1).viewinfo.band = 0.    ; observation band (set to 0)

  ;; restrict area for lat vs. long if needed
  (*ptr1).viewinfo.tsmin = lon1 ; longitude min
  (*ptr1).viewinfo.tsmax = lon2 ; longitude max
  (*ptr1).viewinfo.psmin = lat1 ; latitude min
  (*ptr1).viewinfo.psmax = lat2 ; latitude max

  tempwav1  = readmodel(ptr1) 
  losint,ptr1,/noprint,/noplot,/CALLC,scatter=0 ; do scattered light with g=0
  im = (*(*ptr1).ptr_intop).inu
  lons = (*(*ptr1).ptr_intop).arr1
  lats = (*(*ptr1).ptr_intop).arr2
  sunlat = (*(*ptr1).ptr_intop).sunlat
  sunlong = (*(*ptr1).ptr_intop).sunlong

  ;; disp,im,xvec=lons,yvec=lats
  ;; imdisp,im

  save,file=fileout,im,lons,lats,sunlong,sunlat,a,e,i,asc,wbar,meanlong
  print,'success'
end
