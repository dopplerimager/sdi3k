;---This procedure prints out one record skymap data over multiple lines, with 15 items per line.
;   Optionally, you can prefix it with aheader line, which would typically contain a time stamp.
pro sdi3k_texport_printarr, fout, headline, array, format, noheader=noheader
    if not(keyword_set(noheader)) then printf, fout, headline, format='(a25)'
    j = 0
    while j lt n_elements(array) do begin
        subarr = array(j: (j + 14) < (n_elements(array) - 1))
        printf, fout, string(subarr, format='(' + string(n_elements(subarr)) + format + ')' ), format='(a)'
        j = j + 15
    endwhile
end

;---This procedure adds the "begin" and "end" lines to a data section. It also increments the section number.
pro sdi3k_texport_section, fout, scount, keystring, description, endsec=endsec
    if keyword_set(endsec) then printf, fout, '>>>>>> End Section ' + strcompress(string(scount, format='(i6)'), /remove_all)
    if keyword_set(endsec) then return

    scount = scount + 1
    printf, fout, '>>>>>> Begin Section ' +  strcompress(string(scount, format='(i2.2)'), /remove_all) + $
                  ': [' + strupcase(keystring) + '] -- ' + description
end

;-----------------------------------------------------------------------------------
;
;  Main entry point starts here:

pro sdi3k_ascii_export, setup = setup, files = flis, outpath = outpath, skip_existing=skip_existing

       stp = {export_allsky: 1, $
              export_skymaps: 1, $
              export_spectra: 0, $
              apply_smoothing: 1, $
              time_smoothing: 1.1, $
              space_smoothing: 0.09}
    if not(keyword_set(setup)) then begin
        setup = stp
        obj_edt, setup
    endif else if size(setup, /tname) ne 'STRUCT' then setup = stp

    if not(keyword_set(flis)) then begin
       drive = get_drive()
;      flis = drive + '\Users\sdi3000\Data\poker\PKR 2010_011_Poker_630nm_Red_Sky_Date_01_11.nc'
;      flis = 'D:\users\SDI3000\Data\Poker\PKR 2010_034_Poker_630nm_Red_Sky_Date_02_03.nc'
;      goto, skip_filesel
       fpath = '\users\SDI3000\Data'
       xx = alldisk_findfiles(fpath)
       flis = dialog_pickfile(filter="*.nc", path=xx.(0), title='Select SDI netCDF data files:', /multiple)
    endif
    if not(keyword_set(skip_existing)) then skip_existing = 0

;---Read in the SDI data from a netCDF file. In this case data are returned in the following structures:
;   metadata: mm
;   level-1 spectral fit results: spekfits
;   fitted 2D vector wind skymaps: winds
;   all-sky average wind data: windpars
    for thisfile=0,n_elements(flis)-1 do begin
    ncfile = flis(thisfile)

;---Check if we need to skip because the output file has already been made:
    fspec   = mc_fileparse(ncfile, /lowercase)
    if not(keyword_set(outpath)) then outpath = fspec.path
    outfile = outpath + fspec.name_only + '.txt'
    if file_test(outfile) and skip_existing then goto, SKIP_FILE
    sdi3k_read_netcdf_data,  ncfile, $
                             metadata=mm, spex=spex, winds=winds, spekfits=spekfits, windpars=windpars, $
                             zone_centers=zone_centers, /preprocess

;---Optional smoothing for skymap data:
    if not(setup.apply_smoothing) then goto, NO_SMOOTHING
    posarr = spekfits.velocity
    print, 'Time smoothing winds...'
    sdi3k_timesmooth_fits,  posarr, setup.time_smoothing, mm
    print, 'Space smoothing winds...'
    sdi3k_spacesmooth_fits, posarr, setup.space_smoothing, mm, zone_centers
    spekfits.velocity = posarr
    if mm.maxrec gt 3 then spekfits.velocity = spekfits.velocity - total(spekfits(1:mm.maxrec-2).velocity(0))/n_elements(spekfits(1:mm.maxrec-2).velocity(0))
    tprarr = spekfits.temperature
    print, 'Time smoothing temperatures...'
    sdi3k_timesmooth_fits,  tprarr, setup.time_smoothing, mm
    print, 'Space smoothing temperatures...'
    sdi3k_spacesmooth_fits, tprarr, setup.space_smoothing, mm, zone_centers
    spekfits.temperature = tprarr

NO_SMOOTHING:

;---Now open the output text file:
    openw, fout, outfile, /get_lun
    scount = 0

;---Remaining code outputs each data section:
    sdi3k_texport_section, fout, scount, 'Header', 'Header information'
    printf, fout, '    Site: ', mm.site, format='(a25, a)'
    printf, fout, '    Latitude: ', mm.latitude, format='(a25, f8.4)'
    printf, fout, '    Longitude: ', mm.longitude, format='(a25, f9.4)'
    printf, fout, '    Date UT: ', dt_tm_mk(js2jd(0d)+1, (winds(0).start_time + winds(0).end_time)/2, format='0d$-N$-Y$'), format='(a25, a)'
    printf, fout, '    Records: ' , mm.maxrec, format='(a25, i4.4)'
    printf, fout, '    Wavelength in nm: ', mm.wavelength_nm, format='(a25, f5.1)'
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'FOV', 'Field of view'
    printf, fout, '    FOV Half-Angle: ', mm.sky_fov_deg, format='(a25, f4.1)'
    printf, fout, '    Auroral Oval Angle: ', mm.oval_angle, format='(a25, f5.1)'
    printf, fout, '    Rotation from Oval: ', mm.rotation_from_oval, format='(a25, f5.1)'
    printf, fout, '    Viewing Directions: ', mm.nzones, format='(a25, i3.3)'
    zsel = where(mm.zone_radii gt 0., nz)
    printf, fout, '    Zone Radii: ', strtrim(string(mm.zone_radii(zsel), format='(' + string(nz) + 'i3.2)'), 2), format='(a25, a)'
    zsel = where(mm.zone_radii gt 0., nz)
    printf, fout, '    Sectors: ', strtrim(string(mm.zone_sectors(zsel), format='(' + string(nz) + 'i3.2)'), 2), format='(a25, a)'
    sdi3k_texport_printarr, fout, '    Zenith Angles: ', winds(0).zeniths, 'f12.1'
    sdi3k_texport_printarr, fout, '    Azimuth Angles: ', winds(0).azimuths, 'f12.1'
    sdi3k_texport_printarr, fout, '    Zone Longitudes: ', winds(0).zone_longitudes, 'f12.4'
    sdi3k_texport_printarr, fout, '    Zone Latitudes: ', winds(0).zone_latitudes, 'f12.4'
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    if not(setup.export_allsky) then goto, NO_ALLSKY
    sdi3k_texport_section, fout, scount, 'ALLSKY_TMP_INT', 'All-sky average temperature in K and intensity in arbirary units'
    printf, fout, 'Begin Time', 'End Time', $
                  'Temperature', 'Sigma T', $
                  'Intensity', 'Sigma Inten', format='(6a15)'
    for j=0, n_elements(spekfits) - 1 do begin
        printf, fout, dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      median(spekfits(j).temperature), stddev(spekfits(j).temperature), $
                      median(spekfits(j).intensity), stddev(spekfits(j).intensity), $
                      format='(2a15, 2f15.1, 2e15.4)'


    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'LOCAL_GEO_WINDS', 'Winds in m/s at the station location, aligned GEOGRAPHICALLY'
    printf, fout, 'Begin Time', 'End Time', $
                  'Zonal Wind', 'Sigma Zon', $
                  'Merid Wind', 'Sigma Mer', $
                  'Vertical Wind', 'Sigma Vz', format='(8a15)'
    for j=0, n_elements(spekfits) - 1 do begin
        printf, fout, dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      windpars(j).geo_zonal_wind, windpars(j).siggeozon, $
                      windpars(j).geo_meridional_wind, windpars(j).siggeomer, $
                      windpars(j).vertical_wind, spekfits(j).sigma_velocity(0), $
                      format='(2a15, 6f15.1)'
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'LOCAL_MAG_WINDS', 'Winds in m/s at the station location, aligned GEOMAGNETICALLY'
    printf, fout, 'Begin Time', 'End Time', $
                  'Zonal Wind', 'Sigma Zon', $
                  'Merid Wind', 'Sigma Mer', $
                  'Vertical Wind', 'Sigma Vz', format='(8a15)'
    for j=0, n_elements(spekfits) - 1 do begin
        printf, fout, dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      windpars(j).mag_zonal_wind, windpars(j).sigmagzon, $
                      windpars(j).mag_meridional_wind, windpars(j).sigmagmer, $
                      windpars(j).vertical_wind, spekfits(j).sigma_velocity(0), $
                      format='(2a15, 6f15.1)'
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

NO_ALLSKY:
    if not(setup.export_skymaps) then goto, NO_SKYMAPS

    sdi3k_texport_section, fout, scount, 'TEMP_SKYMAP', 'Temperatures in K at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).temperature, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'SIGMA_TEMP_SKYMAP', 'Temperature uncertainty in K at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).sigma_temperature, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'LOS_WIND_SKYMAP', 'Line-of_sight winds in m/s at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).velocity, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

sdi3k_texport_section, fout, scount, 'SIGMA_LOS_WIND_SKYMAP', 'Line-of_sight wind uncertainty in m/s at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).sigma_velocity, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''


    brg = mm.oval_angle
    sdi3k_texport_section, fout, scount, 'GEO_ZONAL_WIND_SKYMAP', 'Zonal winds in m/s at each viewing location, aligned GEOGRAPHICALLY'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', winds(j).zonal_wind*cos(!dtor*brg) + winds(j).meridional_wind*sin(!dtor*brg), 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'GEO_MERID_WIND_SKYMAP', 'Meridional winds in m/s at each viewing location, aligned GEOGRAPHICALLY'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', winds(j).meridional_wind*cos(!dtor*brg) - winds(j).zonal_wind*sin(!dtor*brg), 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'MAG_ZONAL_WIND_SKYMAP', 'Zonal winds in m/s at each viewing location, aligned GEOMAGNETICALLY'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', winds(j).zonal_wind, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'MAG_MERID_WIND_SKYMAP', 'Meridional winds in m/s at each viewing location, aligned GEOMAGNETICALLY'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, winds(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, winds(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', winds(j).meridional_wind, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'INTENSITY_SKYMAP', 'Intensity in arbitrary units at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, spekfits(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, spekfits(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).intensity/100. < 1e6, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'SNR_SKYMAP', 'Spectral signal-to-noise ratio at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, spekfits(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, spekfits(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).signal2noise < 1e6, 'f8.1', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

    sdi3k_texport_section, fout, scount, 'CHI_SQUARED_SKYMAP', 'Spectral chi-squared ratio at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, spekfits(j).start_time, format='h$:m$:s$'), $
                      ' to ',    dt_tm_mk(js2jd(0d)+1, spekfits(j).end_time,   format='h$:m$:s$'), $
                      format='(a, a8, a, a8)'
        sdi3k_texport_printarr, fout, ' ', spekfits(j).chi_squared, 'f8.2', /noheader
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

NO_SKYMAPS:
    if not(setup.export_spectra) then goto, NO_SPECTRA
;---Export the individual spectra:
    sdi3k_texport_section, fout, scount, 'SPECTRA', 'Spectra in arbitrary units at each viewing location'
    for j=1, n_elements(spekfits) - 1 do begin
        for k=0, n_elements(spex(j).spectra(*,0))-1 do begin
            printf, fout, 'Times ', dt_tm_mk(js2jd(0d)+1, spex(j).start_time, format='h$:m$:s$'), $
                          ' to ',    dt_tm_mk(js2jd(0d)+1, spex(j).end_time,   format='h$:m$:s$'), $
                          ', summed over ' + string(spex(j).scans, format='(i3.3)') + ' scans', $
                          '. Record Number: ' + string(j, format='(i4.4)'), $
                          '. Zone Number: ' + string(k, format='(i3.3)'), $
                          format='(a, a8, a, a8, a, a, a)'
            sdi3k_texport_printarr, fout, ' ', spex(j).spectra(k,*), 'i12.10', /noheader
        endfor
    endfor
    sdi3k_texport_section, fout, scount, /endsec
    printf, fout, ''

NO_SPECTRA:
;---Close the text file and end:
    close, fout
    free_lun, fout
SKIP_FILE:
    endfor
end
