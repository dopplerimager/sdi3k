
pro sdi3k_batch_windfitz, skyfile, drift_mode=drift_mode

if not(keyword_set(drift_mode)) then drift_mode = 'data'

sdi3k_read_netcdf_data, skyfile, metadata=mmsky, zonemap=zonemap, zone_centers=zone_centers, zone_edges=zone_edges, spekfits=spekfits
if mmsky.start_time eq mmsky.end_time then return
;---Determine the wavelength:
    doing_sodium = 0
    doing_red    = 0
    doing_green  = 0
    doing_OH     = 0
    if abs(mmsky.wavelength_nm - 589.0) lt 5. then begin
       lamda = '5890'
       doing_sodium = 1
    endif
    if abs(mmsky.wavelength_nm - 557.7) lt 5. then begin
       lamda = '5577'
       doing_green = 1
    endif
    if abs(mmsky.wavelength_nm - 630.03) lt 5. then begin
       lamda = '6300'
       doing_red = 1
    endif
    if abs(mmsky.wavelength_nm - 843.0) lt 5. then begin
       lamda = '8430'
       doing_oh = 1
    endif
    height = 120.
    if doing_red    then height = 240.
    if doing_green  then height = 120.
    if doing_sodium then height = 90.
    if doing_oh     then height = 97.

    wind_settings = {time_smoothing: 1.4, $
                    space_smoothing: 0.08, $
                    dvdx_assumption: 'dv_dx=0', $
                          algorithm: 'Fourier_Fit', $
                     assumed_height: height, $
                           geometry: 'none'}
;   wind_settings.dvdx_assumption = 'Minimum Gradients'    ; Added MC Jan 26, 2011. Doesn't work very well, but.
    dvdx_zero = wind_settings.dvdx_assumption ne 'dv_dx=1_over_epsilon_times_dv_dt'
;    dvdx_zero = 0

;---Reduced smoothing for 2008 data!
    if mmsky.year gt 2007 then begin
       wind_settings.time_smoothing  = 0.9
       wind_settings.space_smoothing = 0.06
    endif
;---Extra smoothing for pre-1998 data!
    if mmsky.year lt 1998 then begin
       wind_settings.time_smoothing  = 1.6
       wind_settings.space_smoothing = 0.1
    endif
;---Less smoothing for HAARP winds:
    if mmsky.site_code eq 'HRP' then begin
       if doing_red then begin
          wind_settings.time_smoothing  = 0.6
          wind_settings.space_smoothing = 0.05
       endif
       if doing_green then begin
          wind_settings.time_smoothing  = 0.5
          wind_settings.space_smoothing = 0.04
       endif
    endif
       if doing_oh then begin
          wind_settings.time_smoothing  = 1.3
          wind_settings.space_smoothing = 0.08
       endif

    wind_offset = spekfits(0).velocity*0.
    if doing_green then sdi3k_get_wind_offset, getenv('SDI_GREEN_ZERO_VELOCITY_FILE'), wind_offset, mmsky
    if doing_red   then sdi3k_get_wind_offset, getenv('SDI_RED_ZERO_VELOCITY_FILE'),   wind_offset, mmsky
    if doing_oh    then sdi3k_get_wind_offset, getenv('SDI_OH_ZERO_VELOCITY_FILE'),    wind_offset, mmsky
    snrarr = fltarr(n_elements(spekfits))
    chiarr = fltarr(n_elements(spekfits))
    for j=0,n_elements(spekfits) - 1 do begin
        spekfits(j).velocity = spekfits(j).velocity - wind_offset
        chiarr(j) = mean(spekfits(j).chi_squared)
        snrarr(j) = mean(spekfits(j).signal2noise)
    endfor

;---Replace any spectral fits with really high chi-squareds with nearest good record:
    chilim = 1.9
    if abs(mmsky.wavelength_nm - 557.7) lt 1. then chilim = 5.
    posarr = spekfits.velocity
    goods = where(spekfits.chi_squared lt chilim and spekfits.signal2noise gt 70., ngg)
    if ngg le 0 then return
    bads  = where(spekfits.chi_squared ge chilim or spekfits.signal2noise le 70., nn)
    for j=0, nn-1 do begin
        distz = abs(bads(j) - goods)
        best  = where(distz eq min(distz))
        best  = best(0)
        posarr(bads(j)) = posarr(goods(best))
    endfor
    spekfits.velocity = posarr

;     spekfits.velocity = (spekfits.velocity + 2.25*mmsky.scan_channels) mod mmsky.scan_channels



    data_based_drift = strupcase(drift_mode) eq 'DATA'

    sdi3k_drift_correct, spekfits, mmsky, /force, data_based=data_based_drift, insfile=drift_mode ;########
    sdi3k_remove_radial_residual, mmsky, spekfits, parname='VELOCITY'
;    stop
    spekfits.velocity = mmsky.channels_to_velocity*spekfits.velocity

    posarr = spekfits.velocity
    sdi3k_timesmooth_fits,  posarr, wind_settings.time_smoothing, mmsky
    pos2arr = posarr


    sdi3k_spacesmooth_fits, posarr,      wind_settings.space_smoothing, mmsky, zone_centers
    sdi3k_spacesmooth_fits, pos2arr, 1.5*wind_settings.space_smoothing, mmsky, zone_centers
    spekfits.velocity = posarr
    spekfits.velocity(0) = reform(pos2arr(0,*))

    nobs = n_elements(spekfits)
    if nobs gt 2 then spekfits.velocity = spekfits.velocity - total(spekfits(1:nobs-2).velocity(0))/n_elements(spekfits(1:nobs-2).velocity(0))
    sdi3k_fit_wind, spekfits, mmsky, dvdx_zero=dvdx_zero, windfit, wind_settings, zone_centers
    ncid = sdi3k_nc_get_ncid(skyfile, write_allowed=1)
    sdi3k_write_windfitpars, ncid, mmsky, windfit, wind_settings, wind_offset
    ncid = sdi3k_nc_get_ncid(skyfile, write_allowed=0)

;---The following line forces the winds to be calculated relative to GEOGRAPHIC north:
;    mmsky.rotation_from_oval = mmsky.rotation_from_oval + mmsky.oval_angle
;    sdi3k_fit_wind, spekfits, mmsky, dvdx_zero=dvdx_zero, windfit, wind_settings, zone_centers
;    ncid = sdi3k_nc_get_ncid(skyfile, write_allowed=1)
;    sdi3k_write_windfitpars, ncid, mmsky, windfit, wind_settings, /geofit
;    ncid = sdi3k_nc_get_ncid(skyfile, write_allowed=0)

 end





