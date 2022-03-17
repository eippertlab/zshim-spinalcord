PRO ReadNiftiData,filename,data

nim = nifti_read(filename)

data = nifti_get(nim, /COPY)

end
