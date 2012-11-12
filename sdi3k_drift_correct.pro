pro sdi3k_drift_correct, spekfits, mmsky, force=fdc, data_based=dbase, insfile=insfile

;---Find out if we actually need to do drift corrections:
    if not keyword_set(fdc) then begin
       mcchoice, 'Correct for wavelength drift?', ['Yes, Laser Based', 'Yes, Data Based','No'], choice
       fdc   = strpos(choice.name, 'Yes')  ge 0
       dbase = strpos(choice.name, 'Data') ge 0
    endif else fdc = 1
    if not fdc then return

    nz   = mmsky.nzones
    nobs = n_elements(spekfits)
    if keyword_set(dbase) then goto, data_corr

;---Find out the drift calibration file name:
       if keyword_set(insfile) then drfile = insfile else $
       drfile = dialog_pickfile(file=fitfile, $
                                filter='*_CAL_*.nc', $
                                group=widx, title='Select a file of sky spectra: ', $
                                path=host.operation.logging.log_directory)

;---Build an array of drift results:
    sdi3k_read_netcdf_data, drfile, metadata=mmdrft, spekfits=driftarr

;---Get drift times and sky times:
    dt = (driftarr.start_time + driftarr.end_time)/2.
    st = (spekfits.start_time  + spekfits.end_time)/2.

    ndrft  = n_elements(driftarr)
    drvals = fltarr(ndrft)

;---Find where we have good drift spectra:
    nn     = 0
    goodrex= where(driftarr.signal2noise(0) gt 5000., nn)
    if nn gt 2 then begin
;------Find any bad drift spectra:
       badrex  = where(driftarr.signal2noise(0) le 5000., nn)
;------Replace bad drift spectra with the most recent good one:
       if nn gt 0 then begin
          badrex = badrex(sort(driftarr(badrex).record))
          for jj=0,nn-1 do begin
              kk=badrex(jj)
              if kk gt 0 then driftarr(kk) = driftarr(kk-1)
          endfor
       endif

       caldmode = 'byzone'
;       caldmode = 'allsky'
       if caldmode eq 'allsky' then begin
          for j=0,ndrft-1 do begin
              driftarr(j).velocity = median(driftarr(j).velocity(0:42))
          endfor
       endif

       for zidx=0,nz-1 do begin
;---------get a time series of drift values:
          drvals = driftarr.velocity(zidx)
          drvals = drvals - drvals(0)
;---------Interpolate the calibration drift data onto the sky times, and subtract it from the sky data:
          drift = interpol(drvals, dt, st, /lsquad)*mmdrft.wavelength_nm/mmsky.wavelength_nm
          drift = mc_im_sm(drift, 3, /nomedian)
          for j=0,nobs-1 do begin
              if strupcase(getenv('sdi3k_LASER_DRIFTCORR')) ne 'NO' then spekfits(j).velocity(zidx) = spekfits(j).velocity(zidx) - drift(j)
          endfor
       endfor
    endif
;    goto, zero_mean ;##############


;---Do a data-based drift correction:
data_corr:

;stop
if n_elements(spekfits) lt 10 then goto, zero_mean

    deltol = 25.0

    chilim = 1.9
    if abs(mmsky.wavelength_nm - 557.7) lt 1. then chilim = 5.
    bads  = where(total(spekfits.chi_squared, 1)/mmsky.nzones ge chilim, nb)
    velfix = spekfits.velocity
    if nb gt 0 then begin
       goods = where(total(spekfits.chi_squared, 1)/mmsky.nzones lt chilim)
       for j=0,mmsky.nzones - 1 do begin
           repair = INTERPOL( spekfits(goods).velocity(j), (spekfits(goods).start_time  + spekfits(goods).end_time)/2., $
                             (spekfits.start_time  + spekfits.end_time)/2.)
           velfix(j, bads) = repair(bads)
       endfor
       spekfits.velocity = velfix
    endif
;---This is an updated version of the old Inuvik "etalon jump" fixer:
    for zz=0, mmsky.nzones - 1 do begin
        for j = 1,nobs-1 do begin
;       lodx = (j - 3) > 0
;       hidx = (j + 1) < (nobs-1)
;       before = median(aparr(lodx:(j-1)>0 ))
;       after  = mean(aparr((j)<(nobs-1):(hidx)<(nobs-1)))

;       if abs(aparr(j) - before) gt deltol then begin
;          corr = mean(aparr((j+1)<(nobs-1):(hidx + 1)<(nobs-1))) - before
;          aparr(j:*) = aparr(j:*) - corr
;          spekfits(j:*).velocity = spekfits(j:*).velocity - corr
;       endif
;       if abs(aparr(j) - aparr(j-1)) gt deltol then begin
       if abs(spekfits(j).velocity(zz) - spekfits(j-1).velocity(zz)) gt deltol then begin
;       stop
          sgn = spekfits(j).velocity(zz) - spekfits(j-1).velocity(zz)
          corr = -128.
          if sgn gt 0 then corr = 128.

;          aparr(j:*) = aparr(j:*) - corr

              spekfits(j:*).velocity(zz) = spekfits(j:*).velocity(zz) - corr

           endif
       endfor
    endfor
jump_jumpfix:



    sdi3k_zone_angles, mmsky, sky_fov, rad, theta, ridx
;    nobs  = n_elements(spekfits.velocity(0))
;    aparr = fltarr(nobs)
;    for j=0L,nobs-1 do begin
;        zdat     = spekfits(j).velocity
;        zdat     = zdat(sort(zdat))
;        aparr(j) = (20*total(zdat(0.33*nz:0.67*nz))/n_elements(zdat(0.33*nz:0.67*nz)) + 2*spekfits(j).velocity(0))/22
;        aparr(j) =  total(zdat(0.10*nz:0.90*nz))/n_elements(zdat(0.10*nz:0.90*nz))
;    endfor



;---Actually, we'll only use peak positions derived from the inner rings:
    velarr = spekfits.velocity
    velarr = median(velarr, 5)
    inner_rings = where(ridx gt 0 and ridx le 3)
;    aparr = (7*total(spekfits.velocity(inner_rings), 1)/n_elements(inner_rings) + 1.*spekfits.velocity(0))/8.
    aparr = total(velarr(inner_rings, *), 1)/n_elements(inner_rings)
;goto, jump_jumpfix



;stop
;---Check for data at start and end of the night that may be lasers rather than skies. Replace aparr entries that look suspect:
    for j=min([nobs-2,10]),0,-1 do begin
        if median(spekfits(j).temperature) lt 110. then begin
           aparr(j) =aparr(j+1)
        endif
    endfor
    for j=max([nobs-12,1]),nobs-1 do begin
        if median(spekfits(j).temperature) lt 110. then begin
           aparr(j) =aparr(j-1)
        endif
    endfor

;---Taper the first and last few points to average of nearby data:
    if nobs gt 8 then begin
       for j=3,0,-1 do begin
              aparr(j) = total(aparr(j+1:j+3))/3
       endfor
       for j=nobs-4,nobs-1 do begin
              aparr(j) = total(aparr(j-3:j-1))/3
       endfor
    endif

;---Check for bad chi-squared values:
    chiarr = fltarr(nobs)
    for j=0,nobs-1 do begin
        chiarr(j) = mean(spekfits(j).chi_squared)
    endfor

    chilim = 1.9
    if abs(mmsky.wavelength_nm - 557.7) lt 1. then chilim = 5.

    goods = where(spekfits.chi_squared(0) lt chilim and spekfits.signal2noise(0) gt 120. and chiarr lt chilim)
    if keyword_set(dbase) or nobs gt 1 then begin
       if n_elements(goods) gt 10 then begin
          st = (spekfits.start_time  + spekfits.end_time)/2.
          aparr = aparr(goods)
          mcpoly_filter, st(goods), aparr, /lowpass
          drft = interpol(aparr, st(goods), st)
          for j=0l,nobs-1 do begin
              spekfits(j).velocity = spekfits(j).velocity - drft(j)
          endfor
       endif
    endif


;---Force zero-median vertical velocity over the whole night:
zero_mean:
;    if n_elements(spekfits) gt 3 then spekfits.velocity = spekfits.velocity - median(smooth(spekfits.velocity(0), 5, /edge)) ;#################
    if n_elements(spekfits) gt 3 then spekfits.velocity = spekfits.velocity - median(spekfits.velocity(0))
end

