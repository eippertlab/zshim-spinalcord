PRO gaussfilt,array,filter_width
;gaussian filters the 2-d + time or 3-d + time array

if n_params() eq 0 then begin 
	print,'gaussfilt,array,filter_width'
	print,'2D or 3D filtering of image time series'
	print,'CAUTION: arrays have to be of type float !'
	print,'filter width given in Voxels (2D or 3D array)'
	return
endif

if filter_width(0) eq 0.0 then return

sz = SIZE(array)
numdims = N_ELEMENTS(filter_width)
IF numdims LT 2 OR numdims GT 3 THEN WIDG_ERROR,$
   "Wrong number of dimensions."
IF sz(0) LT numdims THEN WIDG_ERROR,$
   "Dimensions of array less then # dimensions."


; Make the convolution kernel first
; Go out two sigma's on either side & cut it
num_filt = FIX(4. * filter_width + 1.) > 3
evennum = (num_filt/2.) EQ (num_filt/2)
num_filt = num_filt + evennum ; make odd number
start = num_filt/2
IF numdims EQ 3 THEN BEGIN ; need 3-d filtering
     filt = FLTARR(num_filt[0],num_filt[1],num_filt[2])
     FOR i=-start[0],start[0] DO $
     FOR j=-start[1],start[1] DO $
     FOR k=-start[2],start[2] DO BEGIN
     	xnum = FLOAT(i)/filter_width[0]
     	ynum = FLOAT(j)/filter_width[1]
     	znum = FLOAT(k)/filter_width[2]
     	filt(i+start[0],j+start[1],k+start[2]) = xnum^2+ynum^2+znum^2
     ENDFOR
ENDIF ELSE BEGIN ; 2-d filtering
     filt = FLTARR(num_filt[0],num_filt[1])
     FOR i=-start[0],start[0] DO $
     FOR j=-start[1],start[1] DO BEGIN
     	xnum = FLOAT(i)/filter_width[0]
     	ynum = FLOAT(j)/filter_width[1]
     	filt(i+start[0],j+start[1]) = xnum^2+ynum^2
     ENDFOR
ENDELSE
filt = EXP(-filt/2.)

; Make sure size is OK.
sizeoffilt = SIZE(filt)
sizeofarr = SIZE(array)
CASE numdims OF
   3 : index = [1,2,3]
   2 : index = [1,2]
ENDCASE
sizeoffilt = sizeoffilt[index]
sizeofarr = sizeofarr[index]
junk = WHERE(sizeoffilt GE sizeofarr,toobig)
IF toobig NE 0 THEN MESSAGE,"Filter size too big."

; Now do the filtering!

IF sz(0) LT 4 THEN ntimes = 1 ELSE ntimes = sz(4)
IF sz(0) LT 3 THEN nslices = 1 ELSE nslices = sz[3]

FOR m=0,ntimes-1 DO BEGIN
   IF numdims EQ 3 THEN BEGIN
   	array(*,*,*,m) = CONVOL(array(*,*,*,m),filt,TOTAL(filt),/EDGE_TRUNCATE)
   ENDIF ELSE FOR k=0,nslices-1 DO BEGIN
	array[*,*,k,m] = CONVOL(array[*,*,k,m],filt,TOTAL(filt),/EDGE_TRUNCATE)
   ENDFOR
ENDFOR
END




