pro sdi3k_spacetime_interpol, obs_coords, out_coords, obsarr, outarr, settings, weighting=weighting, progress=progress

    if keyword_set(progress) then begin
       progressBar = Obj_New("SHOWPROGRESS", message='Percent Completion', title='Interpolating...')
       progressBar->Start
    endif


    deg2km = 111.12
    outarr = replicate(obsarr(0), n_elements(out_coords))
    for j=0L,n_elements(outarr)-1 do begin
        dx    = (obs_coords.lons - out_coords(j).lons)*deg2km*cos(!dtor*out_coords(j).lats)
        dy    = (obs_coords.lats - out_coords(j).lats)*deg2km
        dt    = (obs_coords.js   - out_coords(j).js)*settings.charspeed_kmps
        distz = sqrt(dx^2 + dy^2 + dt^2)
        if keyword_set(weighting) then begin
           if strupcase(weighting) eq '1/D' then weights =1/(distz + settings.dist_par_km)
        endif else weights = exp(-(distz/settings.dist_par_km)^2)
        maxwgt = max(weights)
        use = where(weights gt maxwgt/1000.)

        twgt = total(weights(use))
        for k=0,n_tags(outarr(0)) - 1 do outarr(j).(k) = total(obsarr(use).(k)*weights(use))/twgt
        if keyword_set(progress) then progressbar->update, 100.*j/float(n_elements(outarr)-1)
        wait, 0.00005
    endfor
    if keyword_set(progress) then begin
       progressBar->Destroy
       Obj_Destroy, progressBar
    endif

end