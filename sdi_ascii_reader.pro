
pro sdi_ascii_parse_structure, struc, text

	tags = tag_names(struc)
	for t = 0, n_elements(tags) - 1 do begin

		split = strsplit(text[t], ':', /extract)

	 	type = size(struc.(t), /type)

		if (type ge 1 and type le 5) then begin
			struc.(t) = fix(split[1], type=type)
			continue
		endif

		if (type gt 12 and type le 15) then begin
			struc.(t) = fix(split[1], type=type)
			continue
		endif

		if (type eq 7) then begin
			struc.(t) = split[1]
			continue
		endif

		print, 'BAD TYPE'
		stop
	endfor

end


;\\ Parse a key : value pair
function sdi_ascii_parse_field, text
	split = strsplit(text, ':', /extract)
	return, {key:split[0], val:split[1]}
end

;\\ Parse a numeric (float) array that spans multiple lines
function sdi_ascii_parse_array, text, loc, n_values

	vals = 0
	arr = fltarr(n_values)

	while vals lt n_values do begin
		line = text[loc]
		split = float(strsplit(line, ' ', /extract))
		arr[vals:vals+n_elements(split)-1] = split
		vals += n_elements(split)
		loc++
	endwhile

	return, arr
end

pro sdi_ascii_parse_header, section, text, output

	hdr = {site:'',$
           latitude:0.0, $
           longitude:0.0, $
           date_ut:'', $
           records:0L, $
       	   wavelength_in_nm:0.0}

	sdi_ascii_parse_structure, hdr, text[section.text_start:section.text_end]

	output = create_struct(output, section.name, hdr)

end

pro sdi_ascii_parse_fov, section, text, output

	;\\ Custom parsing here, dealing with variable length arrays

	loc = section.text_start
	half_angle = (sdi_ascii_parse_field(text[loc])).val
	loc ++
	oval_angle = (sdi_ascii_parse_field(text[loc])).val
	loc ++
	rotation_angle = (sdi_ascii_parse_field(text[loc])).val
	loc ++
	zones = (sdi_ascii_parse_field(text[loc])).val
	loc ++
	zone_radii = (sdi_ascii_parse_field(text[loc])).val
	loc ++
	zone_sectors = (sdi_ascii_parse_field(text[loc])).val
	loc ++

	zones = fix(zones)
	zone_radii = fix(strsplit(zone_radii, ' ', /extract))
	zone_sectors = fix(strsplit(zone_sectors, ' ', /extract))

	loc ++ ;\\ now at line after 'Zenith Angles'
	zenith_angles = sdi_ascii_parse_array(text, loc, zones)

	loc ++ ;\\ now at line after 'Azimuth Angles'
	azimuth_angles = sdi_ascii_parse_array(text, loc, zones)

	loc ++ ;\\ now at line after 'Zone Longitudes'
	zone_longitudes = sdi_ascii_parse_array(text, loc, zones)

	loc ++ ;\\ now at line after 'Zone Latitudes'
	zone_latitudes = sdi_ascii_parse_array(text, loc, zones)

	fov = {fov_half_angle:half_angle, $
		   auroral_oval_angle:oval_angle, $
		   rotation_from_oval:rotation_angle, $
		   viewing_directions:zones, $
		   zone_radii:zone_radii, $
		   zone_sectors:zone_sectors, $
		   zenith_angles:zenith_angles, $
		   azimuth_angles:azimuth_angles, $
		   zone_longitudes:zone_longitudes, $
		   zone_latitudes:zone_latitudes}

	output = create_struct(output, section.name, fov)

end

pro sdi_ascii_parse_allsky_tmp_int, section, text, output, timeseries

	recs = output.header.records

	data = replicate({temperature:0.0, $
					  sigma_t:0.0, $
					  intensity:0.0, $
					  sigma_inten:0.0}, recs)

	parse_times = 0
	if (size(timeseries, /type) eq 0) then begin
		timeseries = fltarr(recs, 2)
		parse_times = 1
	endif

	for i = 0, recs - 1 do begin
		split = strsplit(text[i+1], ' ', /extract)

		if parse_times eq 1 then begin
			timeseries[i,0] = total(float(strsplit(split[0], ':', /extract)) * [1.0, 1./60., 1./3600.])
			timeseries[i,1] = total(float(strsplit(split[1], ':', /extract)) * [1.0, 1./60., 1./3600.])
		endif

		data[i] = {temperature:float(split[2]), $
				   sigma_t:float(split[3]), $
				   intensity:float(split[4]), $
				   sigma_inten:float(split[5])}
	endfor

	output = create_struct(output, section.name, data)
end

pro sdi_ascii_parse_local_winds, section, text, output, timeseries

	recs = output.header.records

	data = replicate({zonal_wind:0.0, $
					  sigma_zon:0.0, $
					  merid_wind:0.0, $
					  sigma_mer:0.0, $
					  vertical_wind:0.0, $
					  sigma_vz:0.0}, recs)

	parse_times = 0
	if (size(timeseries, /type) eq 0) then begin
		timeseries = fltarr(recs, 2)
		parse_times = 1
	endif

	for i = 0, recs - 1 do begin
		split = strsplit(text[i+1], ' ', /extract)

		if parse_times eq 1 then begin
			timeseries[i,0] = total(float(strsplit(split[0], ':', /extract)) * [1.0, 1./60., 1./3600.])
			timeseries[i,1] = total(float(strsplit(split[1], ':', /extract)) * [1.0, 1./60., 1./3600.])
		endif

		data[i] = {zonal_wind:float(split[2]), $
				   sigma_zon:float(split[3]), $
				   merid_wind:float(split[4]), $
				   sigma_mer:float(split[5]), $
				   vertical_wind:float(split[6]), $
				   sigma_vz:float(split[7]) }
	endfor

	output = create_struct(output, section.name, data)
end


pro sdi_ascii_parse_local_gradients, section, text, output, timeseries

	recs = output.header.records

	data = replicate({dudx:0.0, $
					  dudy:0.0, $
					  dvdx:0.0, $
					  dvdy:0.0, $
					  vorticity:0.0, $
					  divergence:0.0}, recs)

	parse_times = 0
	if (size(timeseries, /type) eq 0) then begin
		timeseries = fltarr(recs, 2)
		parse_times = 1
	endif

	for i = 0, recs - 1 do begin
		split = strsplit(text[i+1], ' ', /extract)

		if parse_times eq 1 then begin
			timeseries[i,0] = total(float(strsplit(split[0], ':', /extract)) * [1.0, 1./60., 1./3600.])
			timeseries[i,1] = total(float(strsplit(split[1], ':', /extract)) * [1.0, 1./60., 1./3600.])
		endif

		data[i] = {dudx:float(split[2]), $
				   dudy:float(split[3]), $
				   dvdx:float(split[4]), $
				   dvdy:float(split[5]), $
				   vorticity:float(split[6]), $
				   divergence:float(split[7]) }
	endfor

	output = create_struct(output, section.name, data)
end


pro sdi_ascii_parse_skymap, section, text, output

	recs = output.header.records
	zons = output.fov.viewing_directions

	data = fltarr(recs, zons)

	loc = section.text_start + 1

	for i = 0, recs - 1 do begin
		data[i,*] = sdi_ascii_parse_array(text, loc, zons)
		if strmatch(text[loc], '>>>>>> End Section*') eq 1 and i lt recs - 1 then begin
			print, 'RECORDS # SKYMAP MISMATCH - DUPLICATING FIRST RECORD'
			data[1:*,*] = data[0:recs-2,*]
			break
		endif else begin
			loc ++
		endelse
	endfor

	output = create_struct(output, section.name, data)
end


pro sdi_ascii_parse_section, section, text, output, timeseries

	print, 'PARSING SECTION ' + section.name

	;\\ Is this a skymap section?
	split = strsplit(section.name, '_', /extract)
	n = n_elements(split)
	if (n gt 1) then begin

		if (strupcase(split[n-1]) eq 'SKYMAP') then begin
			sdi_ascii_parse_skymap, section, text, output
			return
		endif

	endif

	sectext = text[section.text_start:section.text_end]
	case strupcase(section.name) of
		'ALLSKY_TMP_INT':  sdi_ascii_parse_allsky_tmp_int, section, sectext, output, timeseries
		'LOCAL_GEO_WINDS': sdi_ascii_parse_local_winds, section, sectext, output, timeseries
		'LOCAL_MAG_WINDS': sdi_ascii_parse_local_winds, section, sectext, output, timeseries
		'WIND_GRADIENTS':  sdi_ascii_parse_local_gradients, section, sectext, output, timeseries
		else: print, 'UNKNOWN SECTION: ' + section.name
	endcase

end


function sdi_ascii_reader, filename

	output = {filename:filename}

	if (file_test(filename) eq 0) then return, output

	openr, hnd, filename, /get
	strtext = strarr(file_lines(filename))
	readf, hnd, strtext
	free_lun, hnd

	start_marker = '>>>>>> Begin Section'
	end_marker = '>>>>>> End Section'

	for lineIdx = 0L, n_elements(strtext) - 1 do begin

		line = strtext[lineIdx]

		if strmatch(line, start_marker + '*') eq 1 then begin

			split = strsplit(line, ' ', /extract)
			sec_number = fix((strsplit(split[3], ':', /extract))[0])
			sec_name = strmid(split[4], 1, strlen(split[4]) - 2)
			sec_start = lineIdx + 1

		endif else begin

			if strmatch(line, end_marker + '*') eq 1 then begin

				newSection = {text_start: sec_start, $
							  text_end: lineIdx, $
							  number: sec_number, $
							  name: sec_name }

				sec_start = 0
				sec_number = -1
				sec_name = ''

				if (size(sections, /type) eq 0) then begin
					sections = [newSection]
				endif else begin
					sections = [sections, newSection]
				endelse
			endif

		endelse

	endfor ;\\ line loop

	;\\ Parse the header and FOV sections first, these are required
	hdr = where(strupcase(sections.name) eq 'HEADER', n_hdr)
	fov = where(strupcase(sections.name) eq 'FOV', n_fov)
	if (n_hdr ne 1) then begin
		print, 'NO HEADER INFORMATION'
		return, output
	endif
	if (n_fov ne 1) then begin
		print, 'NO FIELD-OF-VIEW INFORMATION'
		return, output
	endif

	sdi_ascii_parse_header, sections[hdr[0]], strtext, output
	sdi_ascii_parse_fov, sections[fov[0]], strtext, output

	;\\ Parse remining sections
	parse = where(strupcase(sections.name) ne 'HEADER' and $
			      strupcase(sections.name) ne 'FOV', n_parse)

	for i = 0, n_parse - 1 do begin
		sdi_ascii_parse_section, sections[parse[i]], strtext, output, timeseries
	endfor

	if (size(timeseries, /type) ne 0) then begin

		max_0 = 0.
		max_1 = 0.
		for i = 0, n_elements(timeseries[*,0]) - 1 do begin
			if timeseries[i,0] lt max_0 then timeseries[i,0] += 24
			if timeseries[i,1] lt max_1 then timeseries[i,1] += 24

			if timeseries[i,0] gt max_0 then max_0 = timeseries[i,0]
			if timeseries[i,1] gt max_1 then max_1 = timeseries[i,1]
		endfor

		month_num = ['january', 'february', 'march', 'april', $
					 'may', 'june', 'july', 'august', 'september', $
					 'october', 'november', 'december' ]

		recs = output.header.records
		date = strsplit(output.header.date_ut, '-', /extract)

		month = replicate(where(month_num eq strlowcase(date[1]), nm), recs)
		year = replicate(fix(date[2]), recs)
		day = replicate(fix(date[0]), recs)

		output = create_struct(output, 'time', {year:year, $
												month:month, $
												day:day, $
												ut_begin:timeseries[*,0], $
												ut_end:timeseries[*,1]})
	endif

	return, output

end
