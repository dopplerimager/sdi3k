;--This procedure finds all SDI netCDF spectra files (with the specified extensions) in "log_dir" and makes new copies
;  in "stripped_dir" - where the new copies are the same as the originals, but with the background sky images and any
;  analysis results removed. The original files are then moved to an archive directory.
;
;  The purpose of this routine is to prepare smaller copies of the data that can be beamed in from remote field sites
;  with limited bandwidth links.


;-------------------------------------------------------------------
;   This routine creates an empty destination netCDF file and defines
;   the required dimensions and variables for SDI data. Then, it
;   copies the "header" info from the corresponding source netCDF
;   file:
pro istrip_create_header, src, dest

       print, 'Creating file: ' + dest
       sdi3k_read_netcdf_data, src, metadata=mm, /close
       rem = ncdf_open (src)
       ncdf_varget,  rem, ncdf_varid(rem, 'Zone_Radii'),   radii,   offset = 0, count=mm.rings+1
       ncdf_varget,  rem, ncdf_varid(rem, 'Zone_Sectors'), sectors, offset = 0, count=mm.rings
       ncdf_close, rem

				ncdid = ncdf_create(dest, /clobber)

			;\\ Create some dimensions
				chan_dim_id = ncdf_dimdef(ncdid, 'Channel', mm.scan_channels)
				zone_dim_id = ncdf_dimdef(ncdid, 'Zone',    mm.nzones)
				time_dim_id = ncdf_dimdef(ncdid, 'Time',    /unlimited)
				rid 		= ncdf_dimdef(ncdid, 'Rings', 	mm.rings+1)
				rid2 		= ncdf_dimdef(ncdid, 'Rings2', 	mm.rings)
;---------------The following two dimensions are not used in the realtime local copy
;               files, but are needed for the netCDF reader to work correctly:
				xdim_id = ncdf_dimdef(ncdid, 'XDim',  mm.columns)
				ydim_id = ncdf_dimdef(ncdid, 'YDim',  mm.rows)

				bdate = bin_date(systime(/ut))
				date = string(bdate(2)) + '/' + string(bdate(1)) + '/' + string(bdate(0))


			;\\ Create the global attributes
				ncdf_attput, ncdid, /global, 'Start_Date_UT',dt_tm_mak(js2jd(0d)+1, mm.start_time, format='0d$/0n$/Y$'),  	       /char
				ncdf_attput, ncdid, /global, 'Site',      	 mm.site,      /char
				ncdf_attput, ncdid, /global, 'Site_code', 	 mm.site_code, /char
				ncdf_attput, ncdid, /global, 'Latitude',  	 mm.latitude,  /float
				ncdf_attput, ncdid, /global, 'Longitude', 	 mm.longitude, /float
				ncdf_attput, ncdid, /global, 'Operator',  	 mm.operator,  /char
				ncdf_attput, ncdid, /global, 'Comment',   	 mm.comment,   /char
				ncdf_attput, ncdid, /global, 'Software',  	 'sdi3k_ncdf_filestrip',  /char

			;\\ Create the variables
				id = ncdf_vardef  (ncdid, 'Start_Time',      time_dim_id, /long)
		       	id = ncdf_vardef  (ncdid, 'End_Time',        time_dim_id, /long)
    		   	id = ncdf_vardef  (ncdid, 'Number_Scans',    time_dim_id, /short)
       			id = ncdf_vardef  (ncdid, 'X_Center',        time_dim_id, /float)
      			id = ncdf_vardef  (ncdid, 'Y_Center',        time_dim_id, /float)
      			id = ncdf_vardef  (ncdid, 'Cam_Temp',        time_dim_id, /float)
      			id = ncdf_vardef  (ncdid, 'Cam_Gain',        time_dim_id, /short)
      			id = ncdf_vardef  (ncdid, 'Cam_Exptime',     time_dim_id, /float)
      			id = ncdf_vardef  (ncdid, 'X_Bin',             	  /short)
      			id = ncdf_vardef  (ncdid, 'Y_Bin',                /short)
     			id = ncdf_vardef  (ncdid, 'Gap',                  /float)
       			id = ncdf_vardef  (ncdid, 'Nm_Per_Step',     time_dim_id, /float)
       			id = ncdf_vardef  (ncdid, 'Scan_Channels',        /short)
       			id = ncdf_vardef  (ncdid, 'Gap_Refractive_Index', /float)
       			id = ncdf_vardef  (ncdid, 'Zone_Radii',      rid, /float)
       			id = ncdf_vardef  (ncdid, 'Zone_Sectors',    rid2, /byte)
				id = ncdf_vardef  (ncdid, 'Wavelength',      	  /float)
				id = ncdf_vardef  (ncdid, 'Leg1_Start_Volt',  time_dim_id, /short)
				id = ncdf_vardef  (ncdid, 'Leg2_Start_Volt',  time_dim_id, /short)
				id = ncdf_vardef  (ncdid, 'Leg3_Start_Volt',  time_dim_id, /short)
				id = ncdf_vardef  (ncdid, 'Leg1_Offset', 	  	  /float)
				id = ncdf_vardef  (ncdid, 'Leg2_Offset', 	  	  /float)
				id = ncdf_vardef  (ncdid, 'Leg3_Offset', 	  	  /float)
				id = ncdf_vardef  (ncdid, 'Spectra', [zone_dim_id, chan_dim_id, time_dim_id], /long)

			;\\ Write the attributes
				ncdf_attput, ncdid, ncdf_varid(ncdid, 'Start_Time'),           'Units', 'Julian seconds', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'End_Time'),             'Units', 'Julian seconds', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Number_Scans'),         'Units', 'Etalon scans', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Zone_Radii'),           'Units', 'Zone ring radii percent fov', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Zone_Sectors'),         'Units', 'Sectors per ring', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'X_Center'),             'Units', 'Image pixel number', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Y_Center'),             'Units', 'Image pixel number', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Cam_Temp'),             'Units', 'Degrees', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Cam_Gain'),             'Units', 'Dimensionless', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Cam_Exptime'),          'Units', 'Seconds', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'X_Bin'), 	           'Units', 'Image x binning in pixels', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Y_Bin'),     	       'Units', 'Image y binning in pixels', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Gap'),                  'Units', 'mm', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Nm_Per_Step'),          'Units', 'nm', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Scan_Channels'),        'Units', 'Etalon steps per interference order', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Gap_Refractive_Index'), 'Units', 'Dimensionless', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Wavelength'),   		   'Units', 'nm', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg1_Start_Volt'),      'Units', 'Digital voltage', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg2_Start_Volt'),      'Units', 'Digital voltage', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg3_Start_Volt'),      'Units', 'Digital voltage', /char
       			ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg1_Offset'),   	   'Units', 'Dimensionless', /char
				ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg2_Offset'),   	   'Units', 'Dimensionless', /char
				ncdf_attput, ncdid, ncdf_varid(ncdid, 'Leg3_Offset'),   	   'Units', 'Dimensionless', /char
                ncdf_attput, ncdid, ncdf_varid(ncdid, 'Spectra'),              'Units', 'Camera digital units', /char
                ncdf_attput, ncdid, ncdf_varid(ncdid, 'Accumulated_Image'),    'Units', 'Camera digital units', /char

				ncdf_control, ncdid, /endef

;---------------Write the static variables:
				ncdf_varput, ncdid, ncdf_varid(ncdid, 'X_Bin'),	                mm.xbin
      			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Y_Bin'),                 mm.ybin
      			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Gap'),                   mm.gap_mm
      			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Scan_Channels'),         mm.scan_channels
       			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Gap_Refractive_Index'),  mm.gap_refractive_index
       			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Wavelength'),  		    mm.wavelength_nm
       			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Zone_Radii'),            radii
       			ncdf_varput, ncdid, ncdf_varid(ncdid, 'Zone_Sectors'),          sectors

				ncdf_control,ncdid, /sync
				ncdf_close,  ncdid

end

;-------------------------------------------------------------------
;   This routine appends the most recent data to the dstal copies
;   of the current SDI data files:
pro istrip_copy_records, srcfile, dstfile, count
    count = 0

    src   = ncdf_open (srcfile)
    dst   = ncdf_open (dstfile, /write)

    ncdf_diminq, src, ncdf_dimid(src, 'Time'), dummy, maxsrc
    ncdf_diminq, dst, ncdf_dimid(dst, 'Time'), dummy, maxdst

    if maxsrc gt maxdst then begin
       ncdf_diminq, src, ncdf_dimid(src, 'Zone'),    dummy,  nz
       ncdf_diminq, src, ncdf_dimid(src, 'Channel'), dummy,  nchan
       spectra = ulonarr(nz, nchan)
       for record=0, maxsrc-1 do begin
           print, 'Appending record ' + strcompress(string(record), /remove_all) + ' to ' + dstfile
;----------Read the exposure data from the SDI source file:
           ncdf_varget,  src, ncdf_varid(src, 'Spectra'),        spectra, offset=[0,0,record], $
                         count=[n_elements(spectra(*,0)), n_elements(spectra(0,*)), 1]
           ncdf_varget,  src, ncdf_varid(src, 'Start_Time'),     stime, offset=[record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'End_Time'),       etime, offset=[record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Number_Scans'),   scanz, offset=[record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'X_Center'),       x_center, 			offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Y_Center'),       y_center, 			offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Nm_Per_Step'),    nm_per_step, 		offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Leg1_Start_Volt'),leg1_start_volt, 	offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Leg2_Start_Volt'),leg2_start_volt, 	offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Leg3_Start_Volt'),leg3_start_volt, 	offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Cam_Temp'),       cam_temp, 			offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Cam_Gain'),       cam_gain, 			offset = [record], count=[1]
           ncdf_varget,  src, ncdf_varid(src, 'Cam_Exptime'), 	 cam_exptime, 		offset = [record], count=[1]

;----------Write the exposure to the dstal destination file:
           ncdf_varput,  dst, ncdf_varid(dst, 'Spectra'),        spectra, offset=[0,0,record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Start_Time'),     stime,   offset=[record]
           ncdf_varput,  dst, ncdf_varid(dst, 'End_Time'),       etime,   offset=[record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Number_Scans'),   scanz,   offset=[record]

           ncdf_varput,  dst, ncdf_varid(dst, 'X_Center'),       x_center, 			offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Y_Center'),       y_center, 			offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Nm_Per_Step'),    nm_per_step, 		offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Leg1_Start_Volt'),leg1_start_volt, 	offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Leg2_Start_Volt'),leg2_start_volt, 	offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Leg3_Start_Volt'),leg3_start_volt, 	offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Cam_Temp'),       cam_temp, 			offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Cam_Gain'),       cam_gain, 			offset = [record]
           ncdf_varput,  dst, ncdf_varid(dst, 'Cam_Exptime'), 	 cam_exptime, 		offset = [record]
           count = count + 1
           wait, 0.01
       endfor
    endif

	ncdf_control,src, /sync
    ncdf_close, src
    ncdf_close, dst
    wait, 0.5

end


pro ncdf_image_strip, fname, stripped_dir
    finfo    = mc_fileparse(fname)
    destname = stripped_dir + finfo.name_only + '_stripped' + finfo.extension
    istrip_create_header, fname, destname
    istrip_copy_records,  fname, destname, count
end

;------------------------------------------------------------
;---This is the main entry point for sdi3k_ncdf_filestrip:

pro sdi3k_ncdf_filestrip, log_dir, archive_dir, stripped_dir, filter=filter

if not(keyword_set(filter)) then filter = ['*.pf', '*.nc', '*.sky', '*.las']

flis = file_search(log_dir + filter)
for j=0,n_elements(flis)-1 do begin
    ncdf_image_strip, flis(j), stripped_dir
    file_move, flis(j), archive_dir
endfor
end