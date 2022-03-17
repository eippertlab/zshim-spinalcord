PRO ReadData,filename,data,scan,reconr

; ***  Bruker oder Dicom und einlesen *** 
IF N_PARAMS() GT 2 THEN BEGIN
	IF FILE_TEST(filename+'/subject') EQ 1 THEN BEGIN
		IF N_PARAMS() EQ 3 THEN ReadBrukerData,filename,data,scan
		IF N_PARAMS() EQ 4 THEN ReadBrukerData,filename,data,scan,reconr
	ENDIF ELSE BEGIN
		IF N_PARAMS() EQ 3 THEN ReadDicomData,filename,data,scan
	ENDELSE
ENDIF

;*** Datei-Endung unterscheiden

A=STRSPLIT(filename,'.', /EXTRACT)
file_ext=A(n_elements(A) - 1)

if file_ext eq 'nii' or file_ext eq 'gz' then begin
	ReadNiftiData,filename,data
endif

if file_ext eq 'v' then begin
	ReadVistaData,filename,1,temp
	B=size(temp)
	Dim=B(0)

	if Dim eq 4 then begin
		NX = B(1)
		NY = B(2)
		NR = B(3) 
		NS = B(4)

		ReadVistaHeader,filename,header
		ParseVistaHeader,header,1,'repn',val_repn
	
		if val_repn eq 'ubyte' then begin
    			data=bytarr(NX,NY,NS,NR)
		endif
	
		if val_repn eq 'short' then begin 
        		data=intarr(NX,NY,NS,NR)
		endif
	
		if val_repn eq 'float' then begin
			data=fltarr(NX,NY,NS,NR)
		endif
	
		for rep=0, NR-1 do begin
			for slice=0, NS-1 do begin
				data(*,*,slice,rep)=temp(*,*,rep,slice)
			endfor
		endfor
	endif else data = temp
endif

end
