PRO SetGeometryInfoNifti,nim,geometryInfo

ec = nifti_setProperty(nim,'NX',geometryInfo.dim(1))
ec = nifti_setProperty(nim,'NY',geometryInfo.dim(2))
ec = nifti_setProperty(nim,'NZ',geometryInfo.dim(3))
; ec = nifti_setProperty(nim,'NT',geometryInfo.dim(4)) ; nifti library fills this field depending on matrix dimension 4
ec = nifti_setProperty(nim,'DX',geometryInfo.pixdim(1))
ec = nifti_setProperty(nim,'DY',geometryInfo.pixdim(2))
ec = nifti_setProperty(nim,'DZ',geometryInfo.pixdim(3))
ec = nifti_setProperty(nim,'DT',geometryInfo.pixdim(4))
ec = nifti_setProperty(nim,'xyz_units',long(geometryInfo.xyz_units))
ec = nifti_setProperty(nim,'qform_code',long(geometryInfo.qf_code))
ec = nifti_setProperty(nim,'sform_code',long(geometryInfo.sf_code))
ec = nifti_setProperty(nim,'qto_xyz',geometryInfo.qto_xyz)
ec = nifti_setProperty(nim,'qto_ijk',geometryInfo.qto_ijk)
ec = nifti_setProperty(nim,'sto_xyz',geometryInfo.sto_xyz)
ec = nifti_setProperty(nim,'sto_ijk',geometryInfo.sto_ijk)
ec = nifti_setProperty(nim,'quatern_b',geometryInfo.quatern_b)
ec = nifti_setProperty(nim,'quatern_c',geometryInfo.quatern_c)
ec = nifti_setProperty(nim,'quatern_d',geometryInfo.quatern_d)
ec = nifti_setProperty(nim,'qoffset_x',geometryInfo.qoffset_x)
ec = nifti_setProperty(nim,'qoffset_y',geometryInfo.qoffset_y)
ec = nifti_setProperty(nim,'qoffset_z',geometryInfo.qoffset_z)
ec = nifti_setProperty(nim,'qfac',geometryInfo.qfac)

end
