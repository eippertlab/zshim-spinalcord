PRO WriteNiftiMap,filename,data,geometryInfo
; save float data in a nifti file

nim = nifti_new()
ec = nifti_set(nim,data)

if n_params() eq 3 then setgeometryinfonifti,nim,geometryInfo

ec = nifti_write(nim,filename)
ec = nifti_free(nim)

end
