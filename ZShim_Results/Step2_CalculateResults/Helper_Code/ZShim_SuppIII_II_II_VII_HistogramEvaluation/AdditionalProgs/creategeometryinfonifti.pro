PRO CreateGeometryInfoNifti,filename,Info

Info = {dim:INTARR(8),pixdim:FLTARR(8),xyz_units:2,qf_code:1,sf_code:2,qto_xyz:FLTARR(4,4),qto_ijk:FLTARR(4,4),sto_xyz:FLTARR(4,4),sto_ijk:FLTARR(4,4),$
quatern_b:0.0,quatern_c:0.0,quatern_d:0.0,qoffset_x:0.0,qoffset_y:0.0,qoffset_z:0.0,qfac:1.0}

dlm_load,'nifti'
nim = nifti_new()
nim = nifti_read(filename)

Info.dim = nifti_getProperty(nim,'dim')
Info.pixdim = nifti_getProperty(nim,'pixdim')
Info.xyz_units = nifti_getProperty(nim,'xyz_units')
Info.qf_code = nifti_getProperty(nim,'qform_code')
Info.sf_code = nifti_getProperty(nim,'sform_code')
Info.qto_xyz = nifti_getProperty(nim,'qto_xyz')
Info.qto_ijk = nifti_getProperty(nim,'qto_ijk')
Info.sto_xyz = nifti_getProperty(nim,'sto_xyz')
Info.sto_ijk = nifti_getProperty(nim,'sto_ijk')
Info.quatern_b = nifti_getProperty(nim,'quatern_b')
Info.quatern_c = nifti_getProperty(nim,'quatern_c')
Info.quatern_d = nifti_getProperty(nim,'quatern_d')
Info.qoffset_x = nifti_getProperty(nim,'qoffset_x')
Info.qoffset_y = nifti_getProperty(nim,'qoffset_y')
Info.qoffset_z = nifti_getProperty(nim,'qoffset_z')
Info.qfac = nifti_getProperty(nim,'qfac')

end
