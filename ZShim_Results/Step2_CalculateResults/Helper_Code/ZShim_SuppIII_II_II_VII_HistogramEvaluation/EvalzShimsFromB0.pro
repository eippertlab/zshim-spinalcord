; AUTHOR:
; Toralf Mildner, MPI for Human Cognitive and Brain Sciences, Stephanstr. 1a, 04103 Leipzig, Germany
; mildner@cbs.mpg.de
; NOTICE: the procedures 'readdata.pro' and 'writeniftimap.pro' use the IDL dynamic link module (DLM) 'NIFTI' which was kindly provided by Ron Kneusel
; this DLM can be provided on request

; read the path for the output data, the file name of the B0 map to be evaluated, ...
; the file name of the spinal cord mask, and a flag determining whether the GRE (vendor) or
; the CBSFL3D (in-house) fieldmap is provided
outputpath=''
B0mapName=''
T2mapName=''
useFM_GRE=''
read,PROMPT='outputpath: ',outputpath
read,PROMPT='B0mapName: ',B0mapName
read,PROMPT='T2mapName: ',T2mapName
read,PROMPT='useFM_GRE: ',useFM_GRE

; gyromagnetic ratio of protons in rad/s/mT
GAMMA_RAD = 2.675222099e5
Geff_max = 0.21 ; mT/m
NR_Zref = 21
nSurrounds = 2
; echo spacing in gre_fieldmapping
dTE = 2.46e-3

; read and store the geometry information of the input B0 map in order to
; copy it later to the output file
creategeometryinfonifti,B0mapName,geometryInfoB0

; read the nifti file of the B0 map which was obtained by SPM Dicom-to-Nifti conversion
readdata,B0mapName,b0map_temp
NX = (size(b0map_temp))(1)
NY = (size(b0map_temp))(2)
NZ = (size(b0map_temp))(3)

; scale the B0 map to the unit of Hz
if useFM_GRE eq 'false' then begin
	b0map = float(b0map_temp)*0.097 - 200.0
endif else begin
	b0map = float(b0map_temp)/4096.0/2.0/dTE
endelse

; create the output array of the B0 gradient and set the width of the 2D Gaussian filter to 1x1 voxels
B0_gradz = fltarr(NX,NY,NZ)
widthGauss = [1.0,1.0]

; calculate the B0 gradient using the function gradient.pro and filter it using gaussfilt.pro in a loop over sagittal slices
for x = 0, NX - 1 do begin
	 inputMapB0 = REFORM(b0map(x,*,*),NY,NZ)
	 gradMap = gradient(inputMapB0,/VECTOR)
	 map2Filt = gradMap(*,*,1)
	 gaussfilt,map2Filt,widthGauss
	 B0_gradz(x,*,*) = map2Filt(*,*)
endfor

; scale the B0 gradient map to the unit of mT/m and save it to the output path
B0_gradz = B0_gradz*2.0*!Pi/GAMMA_RAD	; mT/mm
B0_gradz = -1e3*B0_gradz	  	; mT/m
writeniftimap,outputpath + '/B0_gradz.nii', B0_gradz, geometryInfoB0

; read the nifti file of the spinal cord mask
readdata,T2mapName,maskSC
help,maskSC

; create and open the output file for the z-shim indices
openw,1,outputpath + '/zShimsMean_vs_Histogram.dat'

; loop for calculating of z-shim indices (shrunk to the actual slab of the EPI slices in z-direction (ranging from fieldmap sub-slice 30 to 149 (zero-based)
; the step size of this loop is 5 (i.e., it goes through the sub-slices of the fieldmap in steps corresponding to each EPI slice)
for i = 30, 150 - 1, 5 do begin
	; get the mask and B0 gradient values inside the 5 sub-slices corresponding to one EPI slices
	msk = maskSC(*,*,i:i + 4)
	slice = B0_gradz(*,*,i:i + 4)

	; create the roi containing all voxels inside the spinal cord of the EPI sub-slice
	roi = where(msk ne 0.0)

	; calculate the mean B0 gradient inside the spinal cord ROI
	valsB0gradzMean = mean(slice(roi))

	; calculate the corresponding z-shim index
	zShimIndexMean = round((valsB0gradzMean + Geff_max)/Geff_max*(NR_Zref - 1)*0.5) + 1

	; create a histogram of B0 gradients inside the spinal cord ROI with a binsize of 0.01 mT/m
	his = histogram(slice(roi),BINSIZE=0.01)

	; smooth the histogram using the IDL function smooth
	his = smooth(his,n_elements(his)/20)

	; get the array of x (B0 gradient) values in the histogram
	aa=min(slice(roi))
	bb=max(slice(roi))
	xx = aa+(bb-aa)/n_elements(his)*indgen(n_elements(his))

	; sort the y values of the histogram (yields sorted indices beginning from the smallest)
	sort_his = SORT(his)

	; create and fill an array containing the intensity averages over the center +/- 2 neighboring points of each of the three most frequent bins
	vAverages = fltarr(3)
	for base = 0, 2 do begin
		range = sort_his(n_elements(his) - base - 1) - nSurrounds + indgen(2*nSurrounds + 1)
		validRange = where(range ge 0 and range lt n_elements(his))
		integral = mean(his(range(validRange)))
		vAverages(base) = integral/n_elements(validRange)
	endfor

	; sort this array
	idxMax = SORT(vAverages)

	; get the index of the main peak in the histogram from the index of the maximum (idxMax(2)) of the array of intensity averages
	idxSelectedPeak = sort_his(n_elements(his) - idxMax(2) - 1)

	; perform a weighted summation of B0 gradient bins (x values or array 'xx') over a range of +/- 10 points around the main peak
	; use a scaled version of the histogram (scaled_his) by setting the intensity of the main peak to the value of 1 (or 100%)
	xxNorm = 0.0
	xxWeighted = 0.0
	scaled_his = 1.0*his/his(idxSelectedPeak)
	for ofs = -10, 10 do begin
		if (idxSelectedPeak + ofs) ge 0 and (idxSelectedPeak + ofs) lt n_elements(his) then begin

			; the summation range is additionally shrunk by the condition that the intensity values should be larger than 25% of the main peak
			if scaled_his(idxSelectedPeak + ofs) gt 0.25 then begin
				xxWeighted = xxWeighted + scaled_his(idxSelectedPeak + ofs)*xx(idxSelectedPeak + ofs)
				xxNorm = xxNorm + scaled_his(idxSelectedPeak + ofs)
			endif
		endif
	endfor

	; calculate the B0 gradient resulting from this shrunk histogram summation by normalization
	valsB0gradzHist = xxWeighted/xxNorm

	; calculate the z-shim index of the histogram-based evaluation
	zShimIndexHist = round((valsB0gradzHist + Geff_max)/Geff_max*(NR_Zref - 1)*0.5) + 1

	; the next part is for the case that you want to save the histogram and show the z-shim values for a certain EPI slice number
	if i eq -1 then begin
		print,'valsB0gradzMean:',valsB0gradzMean
		print,'valsB0gradzHist:',valsB0gradzHist
		print,'range B0 grad in Z: ',bb - aa
		val_histo=fltarr(2,n_elements(xx))
		val_histo(*)=0.0
		val_histo(0,*)=transpose(xx)
		val_histo(1,*)=transpose(his)
		openw,2,outputpath + '/histogramSlice.dat'
		printf,2,val_histo
		close,2
	endif

	; append the mean and histogram-based z-shim values of the current EPI slice to the output file
	printf,1,FORMAT='(I2.1,2(%"\t",F7.2))',1 + floor(i/5) - 6,zShimIndexMean,zShimIndexHist
endfor

close,1


end
