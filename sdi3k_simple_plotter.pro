pro sdi3k_simple_plotter, filename, plotarr, canvas_size=canvas_size, what_plot=what_plot, $
                          time_smoothing=tsm, space_smoothing=ssm, yrange=yrange, $
                          quality_filter=quality_filter, lambda=lambda, nzones=nzones, psplot=psplot, $
                          notitle=notitle, cadence=cadence, time_format=time_format
    if not(keyword_set(cadence)) then cadence=1
    if not(keyword_set(time_format)) then time_format = 'y$doy$ h$:m$'
    if not(keyword_set(what_plot)) then begin
       plotsel = ['LOSWIND',  'TEMPERATURE',     'INTENSITY',  'RGBMAPS', $
                  'WINDMAPS', 'XY_TEMPERATURE',  'MEAN_WINDS', 'WIND_GRADIENTS', 'DIVOR']
       mcchoice, 'Plot Type: ', plotsel, choice, $
                  heading = {text: 'Choose a Plot Type', font: 'Helvetica*Bold*Proof*30'}
       what_plot = choice.name

    endif
    what_plot = strupcase(what_plot)

    whoami, dir, file
    pathinfo = mc_fileparse(dir + file)
    drive = str_sep(pathinfo(0).path, '\')
    drive = drive(0)

   pathinfo = mc_fileparse(filename)
   local_path = pathinfo(0).path
   flats    = findfile(local_path + "Wind_flat_field*.sav")
   mcchoice, 'Wind flat field?', [flats, 'None'], choice, help='This flat field correction will be applied regardless of wavelelength'
   if choice.name ne 'None' then flatfield = flats(choice.index) else flatfield = 'None'
   if flatfield ne 'None' then restore, flatfield

    firstspek = 1
    firstwind = 1
    firstwpar = 1
    for j=0,n_elements(filename)-1 do begin
        sdi3k_read_netcdf_data, filename(j), metadata=mm, winds=wndz, spekfits=spkftz, windpars=wndprz, zonemap=zonemap, zone_centers=zone_centers, /close
        skip = 0
        if size(mm, /tname) ne 'STRUCT' then goto, not_this_file
        if keyword_set(lambda) then skip = abs(mm.wavelength_nm - lambda) gt 2.
        if keyword_set(nzones) then skip = skip or mm.nzones ne nzones
        if skip then goto, not_this_file
        if flatfield ne 'None' then for j=0, n_elements(spkftz) -1 do spkftz(j).velocity = spkftz(j).velocity - wind_flat_field.wind_offset
        if firstspek and n_elements(spkftz) gt 1 then begin
           spekfits  = spkftz
           firstspek = 0
        endif else if n_elements(spkftz) gt 1 then  spekfits = [spekfits, spkftz]
        if firstwind and n_elements(wndz) gt 1 then begin
           winds     = wndz
           firstwind = 0
        endif else if n_elements(wndz) gt 1 then  winds = [winds, wndz]
        if firstwpar and n_elements(wndprz) gt 1 then begin
           windpars  = wndprz
           firstwind = 0
        endif else if n_elements(wndprz) gt 1 then  windpars = [windpars, wndprz]
 not_this_file:
    endfor

;---Determine the wavelength:
    doing_sodium = 0
    doing_red    = 0
    doing_green  = 0
    doing_oh     = 0
    if abs(mm.wavelength_nm - 589.0) lt 5. then begin
       lamda = '5890'
       doing_sodium = 1
       lamstring = '_sodium'
    endif
    if abs(mm.wavelength_nm - 557.7) lt 5. then begin
       lamda = '5577'
       doing_green = 1
       lamstring = '_green'
    endif
    if abs(mm.wavelength_nm - 630.0) lt 5. then begin
       lamda = '6300'
       doing_red = 1
       lamstring = '_red'
    endif
    if abs(mm.wavelength_nm - 843.0) lt 5. then begin
       lamda = '8430'
       doing_oh = 1
       lamstring = '_OH'
    endif

;---Apply any zero-velocity offset correction maps that have been selected:
    if doing_green and getenv('SDI_GREEN_ZERO_VELOCITY_FILE') ne '' then begin
       restore, getenv('SDI_GREEN_ZERO_VELOCITY_FILE')
       print, 'Using vzero map: ', getenv('SDI_GREEN_ZERO_VELOCITY_FILE')
       for j=0,n_elements(spekfits) - 1 do begin
           spekfits(j).velocity = spekfits(j).velocity - wind_offset
       endfor
    endif
    if doing_red and getenv('SDI_RED_ZERO_VELOCITY_FILE') ne '' then begin
       restore, getenv('SDI_RED_ZERO_VELOCITY_FILE')
       print, 'Using vzero map: ', getenv('SDI_RED_ZERO_VELOCITY_FILE')
       for j=0,n_elements(spekfits) - 1 do begin
           spekfits(j).velocity = spekfits(j).velocity - wind_offset
       endfor
    endif

    sdi3k_drift_correct, spekfits, mm, /force, /data
    spekfits.velocity = spekfits.velocity*mm.channels_to_velocity

    while !d.window ge 0 do wdelete, !d.window

;---Initialize various data:
    load_pal, culz

    year      = strcompress(string(mm.year),             /remove_all)
    doy       = strcompress(string(mm.start_day_ut, format='(i3.3)'),     /remove_all)

;---Build the time information arrays:
    tcen   = (spekfits.start_time + spekfits.end_time)/2
    tlist  = dt_tm_mk(js2jd(0d)+1, tcen, format=time_format)

    mcchoice, 'Start Time: ', tlist, choice, $
               heading = {text: 'Start at What Time?', font: 'Helvetica*Bold*Proof*30'}
    jlo = choice.index
    mcchoice, 'End Time: ', tlist, choice, $
               heading = {text: 'End at What Time?', font: 'Helvetica*Bold*Proof*30'}
    jhi = choice.index


    marg = 16 < n_elements(resarr)/2 - 1
    if not(doing_sodium) then rex = indgen(mm.maxrec) else begin
       tsbrite = fltarr(mm.maxrec)
       tsbgnd  = fltarr(mm.maxrec)
       for j=0,mm.maxrec-1 do begin
           tsbrite(j) = median(spekfits(j).intensity)
           tsbgnd(j)  = median(spekfits(j).background)
       endfor
;      rex = where(tsbrite gt 2. and tsbrite lt 100. and tsbgnd lt 1000.)
       rex = where(tsbrite gt 100. and tsbrite lt 5e6 and tsbgnd lt 1e7)
;       rex = [indgen(marg), indgen(marg) + n_elements(resarr) - marg - 1]
    endelse
    if n_elements(rex) lt 2 then goto, PLOTZ_DONE
    sdi3k_remove_radial_residual, mm, spekfits, parname='VELOCITY',    recsel=rex
    sdi3k_remove_radial_residual, mm, spekfits, parname='TEMPERATURE', recsel=rex, /zero_mean
    sdi3k_remove_radial_residual, mm, spekfits, parname='INTENSITY',   recsel=rex, /multiplicative

;---Setup wind, briteness and temperature scaling:
    britearr   = reform(spekfits.intensity)
    intord     = sort(britearr)
    britescale = [0., 1.4*britearr(intord(0.97*n_elements(intord)))]
    medtemp = median(spekfits.temperature)
    medtemp = 100*fix(medtemp/100)
    if not(keyword_set(yrange)) then begin
;       if doing_red    then tprscale = [medtemp - 200. > 0, medtemp + 200.]
       if doing_red    then tprscale = [650., 800.]
       if doing_red    then wndscale = [-330., 330.]
;       if doing_green  then tprscale = [medtemp - 250. > 0, medtemp + 250.]
       if doing_green  then tprscale = [200, 550]
       if doing_green  then wndscale = [-220., 220.]
       if doing_sodium then tprscale = [medtemp - 350. > 0, medtemp + 350.]
       if doing_OH     then tprscale = [150, 300.]
    endif else begin
       tprscale = yrange
       wndscale = yrange
    endelse
    f107 = 70.
    if medtemp gt 850  then f107 = 120.
    if medtemp gt 950  then f107 = 150.
    if medtemp gt 1050 then f107 = 200.
    alt = 240.
    if abs(mm.wavelength_nm - 630.03) gt 10. then alt = 120.

STAGE_SMOOTHING:
;---Now do some smoothing prior to making the maps:
    if not(keyword_set(tsm)) then begin
       if doing_red    then tsm = 1.
       if doing_green  then tsm = 0.6
       if doing_sodium then tsm = 1.5
       if doing_red and mm.year lt 2007   then tsm = 1.5
    endif
    if not(keyword_set(ssm)) then begin
       if doing_red    then ssm = 0.07
       if doing_green  then ssm = 0.03
       if doing_sodium then ssm = 0.03
       if doing_green and mm.site_code eq 'HRP' then begin
          tsm = 0.4
          ssm = 0.02
       endif
    endif

    posarr = spekfits.velocity
    print, 'Time smoothing winds...'
    sdi3k_timesmooth_fits,  posarr, tsm, mm
    print, 'Space smoothing winds...'
    sdi3k_spacesmooth_fits, posarr, ssm, mm, zone_centers
    spekfits.velocity = posarr
    if mm.maxrec gt 3 then spekfits.velocity = spekfits.velocity - total(spekfits(1:mm.maxrec-2).velocity(0))/n_elements(spekfits(1:mm.maxrec-2).velocity(0))

    if not(keyword_set(tsm)) then begin
       if doing_red    then tsm = 1.2
       if doing_green  then tsm = 0.8
       if doing_sodium then tsm = 1.2
       if doing_red and mm.year lt 2007   then tsm = 1.8
    endif
    if not(keyword_set(ssm)) then begin
       if doing_red    then ssm = 0.07
       if doing_green  then ssm = 0.05
       if doing_sodium then ssm = 0.07
       if doing_green and mm.site_code eq 'HRP' then begin
          tsm = 0.6
          ssm = 0.05
       endif
       if doing_red and mm.year lt 2007   then ssm = 0.1
    endif

    tprarr = spekfits.temperature
    print, 'Time smoothing temperatures...'
    sdi3k_timesmooth_fits,  tprarr, tsm, mm
    print, 'Space smoothing temperatures...'
    sdi3k_spacesmooth_fits, tprarr, ssm, mm, zone_centers
    spekfits.temperature = tprarr

;---Set the MSIS and HWM parameters:
    f107 = 70.
    if medtemp gt 850  then f107 = 120.
    if medtemp gt 950  then f107 = 150.
    if medtemp gt 1050 then f107 = 200.
    alt = 240.
    if abs(mm.wavelength_nm - 630.03) gt 10. then alt = 120.

    msis  = {tsplot_msis,   plot_msis: 0, $
                               f10pt7: f107, $
                                   ap: 5., $
                          msis_height: alt}
    hwm   =                 {plot_hwm: 1, $
                               f10pt7: f107, $
                                   ap: 5., $
                           hwm_height: alt}


    xsize    = canvas_size(0)
    ysize    = canvas_size(1)
    if keyword_set(psplot) then begin
       set_plot, 'PS'
       device, /landscape, xsize=25, ysize=15
       device, bits_per_pixel=8, /color, /encapsulated
       device, filename=dialog_pickfile(path='C:\Users\Conde\Main\Poker_SDI\', filter='*.eps')
       !p.charsize = 0.5
       note_size = 0.4
       !p.position = [0.18, 0.13,0.95,0.9]
    endif else begin
       if not(keyword_set(canvas_size)) then canvas_size = [1300, 1000]
       while !d.window ge 0 do wdelete, !d.window
       window, xsize=xsize, ysize=ysize
    endelse

;---Set the bitmap size for skymaps:
    thumsize = 0.7*sqrt(xsize*ysize/(jhi - jlo + 1))
    geo   = {xsize:  xsize, ysize: ysize}
    scale = {yrange: [-150., 150.], auto_scale: 0}
    skymap_settings = {scale: scale, parameter: 'Velocity', black_bgnd: 1, geometry: geo, records: [jlo, jhi]}

    lamlab  = ' !4k!3=' + strcompress(string(mm.wavelength_nm, format='(f12.1)'), /remove_all) + ' nm'
    title = mm.site + ': ' + dt_tm_mk(js2jd(0d)+1, mm.start_time(0), format='d$-n$-Y$,')  + lamlab

    if what_plot eq 'LOSWIND'            then goto, LOSWIND
    if what_plot eq 'TEMPERATURE'        then goto, TEMPERATURE
    if what_plot eq 'INTENSITY'          then goto, INTENSITY
    if what_plot eq 'RGBMAPS'            then goto, RGBMAPS
    if what_plot eq 'WINDMAPS'           then goto, WINDMAPS
    if what_plot eq 'XY_TEMPERATURE'     then goto, XY_TEMPERATURE
    if what_plot eq 'MEAN_WINDS'         then goto, MEAN_WINDS
    if what_plot eq 'WIND_GRADIENTS'     then goto, WIND_GRADIENTS
    if what_plot eq 'DIVOR'              then goto, DIVOR


;---Temperature Time Series:
XY_TEMPERATURE:
    scale = {time_range: [0., 1.08], yrange: tprscale, auto_scale: 0}
    geo   = {xsize:  xsize, ysize: ysize}
    style =   {tsplot_style, charsize: 3.0, $
                            charthick: 3, $
                       line_thickness: 3, $
                       axis_thickness: 4, $
                    menu_configurable: 1, $
                        user_editable: [0,1,2,3]}
    if keyword_set(psplot) then style.charsize=1.5
    npts = n_elements(tcen)
    nz = mm.nzones
    zz = strarr(nz)
    for j=0,nz-1 do zz(j) = 'Zone ' + strcompress(string(j), /remove_all)
    first = 1
    repeat begin
           mcchoice, 'Select zone usage:', ['All', 'Zenith', 'Median', zz, 'Done'], choice
           if choice.name ne 'Done' then begin
              if first eq 1 then begin
                 zsel   = choice.name
                 pastel = 0
              endif else begin
                 zsel = [zsel, choice.name]
                 if choice.name eq 'All' then pp = 0.5 else pp = 0
                 pastel = [pastel, pp]
              endelse
              first = 0
           endif
    endrep until choice.name eq 'Done'

    idxdir = drive + "\users\conde\main\indices\"
    get_solterr_indices, [min(tcen) - 100L*86400, max(tcen) + 100L*86400], idxdir, indices
;    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, alt, indices, msis_pts, /progress
;    epoints = {time: msis_pts.time, y:msis_pts.tz, psym: 6, style: 0, color: culz.rose, thick: 1, symsize: 0.4}

    tsplot_settings = {scale: scale, parameter: 'Temperature', zones: zsel, black_bgnd: 1, geometry: geo, records: [jlo, jhi], msis: msis, hwm: hwm, style: style}
;    sdi3k_tseries_plotter, tlist, tcen, mm, spekfits, tsplot_settings, culz, quality_filter=quality_filter, pastel=pastel, epoints=epoints

    sdi3k_tseries_plotter, tlist, tcen, mm, spekfits, tsplot_settings, culz, quality_filter=quality_filter, pastel=pastel, /notitle
    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, 240, indices, msis_pts
    if n_elements(msis_pts) lt 2 then goto, PLOTZ_DONE
    oplot, msis_pts.time, msis_pts.tz, linestyle=0, color=culz.red, thick=5
    xyouts, msis_pts(npts-1).time, msis_pts(npts-1).tz , ' 240 km', align=0., charthick=2, charsize=0.7*style.charsize, color=culz.red
    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, 200, indices, msis_pts
    oplot, msis_pts.time, msis_pts.tz, linestyle=2, color=culz.red, thick=3
    xyouts, msis_pts(npts-1).time, msis_pts(npts-1).tz, ' 200 km', align=0., charthick=2, charsize=0.7*style.charsize, color=culz.red
    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, 170, indices, msis_pts
    oplot, msis_pts.time, msis_pts.tz, linestyle=2, color=culz.red, thick=3
    xyouts, msis_pts(npts-1).time, msis_pts(npts-1).tz, ' 170 km', align=0., charthick=2, charsize=0.7*style.charsize, color=culz.red
    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, 150, indices, msis_pts
    oplot, msis_pts.time, msis_pts.tz, linestyle=2, color=culz.red, thick=3
    xyouts, msis_pts(npts-1).time, msis_pts(npts-1).tz, ' 150 km', align=0., charthick=2, charsize=0.7*style.charsize, color=culz.red
    get_msis, tcen, {lon: mm.longitude, lat: mm.latitude}, 140, indices, msis_pts
    oplot, msis_pts.time, msis_pts.tz, linestyle=2, color=culz.red, thick=3
    xyouts, msis_pts(npts-1).time, msis_pts(npts-1).tz, ' 140 km', align=0., charthick=2, charsize=0.7*style.charsize, color=culz.red
    goto, PLOTZ_DONE

MEAN_WINDS:
    mc_npanel_plot,  layout, yinfo, /setup
    layout.position = [0.15, 0.14, 0.975, 0.93]
    layout.panels = 2
    layout.time_axis =1
    layout.xrange = [tcen(jlo), tcen(jhi)]
    layout.title  = title
    layout.ypad(0:layout.panels-2) = 0.01 + fltarr(layout.panels-1)
    layout.panel_rgb_background = {factor: 0.05, color: culz.cyan, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}
    layout.plot_panel_background = 1

    yinfo.title = 'Magnetic Zonal!C !CWind [m/s]'
    yinfo.range = wndscale
    layout.charscale = 1.25

    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).mag_zonal_wind}
    mc_npanel_plot,  layout, yinfo, wnddata

    layout.panel_rgb_background = {factor: 0.05, color: culz.yellow, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}
    yinfo.title = 'Magnetic Meridional!C !CWind [m/s]'
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).mag_meridional_wind}
    mc_npanel_plot,  layout, yinfo, wnddata

    goto, PLOTZ_DONE

WIND_GRADIENTS:
    mc_npanel_plot,  layout, yinfo, /setup
    erase, color=culz.white
    layout.position = [0.125, 0.14, 0.975, 0.93]
    layout.panels = 4
    layout.time_axis =1
    layout.xrange = [tcen(jlo), tcen(jhi)]
    layout.title  = title
    layout.erase = 0
    layout.ypad(0:layout.panels-2) = 0.01 + fltarr(layout.panels-1)
    layout.panel_rgb_background = {factor: 0.05, color: culz.cyan, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}
    layout.plot_panel_background = 1
    layout.charscale = 1.05


    yinfo.range = [-1.1, 1.1]
    if keyword_set(yrange) then yinfo.range = yrange
    yinfo.charsize = 1.5

    yinfo.title = 'du/dx!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 1
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*winds(jlo:jhi).dudx(j)}
        mc_npanel_plot,  layout, yinfo, wnddata, panel=0
    endfor
    load_pal, culz
    yinfo.thickness = 3
    yinfo.line_color = culz.black
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).du_dx}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=0

    yinfo.title = 'du/dy!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 1
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*winds(jlo:jhi).dudy(j)}
        mc_npanel_plot,  layout, yinfo, wnddata, panel=1
    endfor
    load_pal, culz
    yinfo.thickness = 3
    yinfo.line_color = culz.black
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).du_dy}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=1

    layout.panel_rgb_background = {factor: 0.05, color: culz.yellow, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}

    yinfo.title = 'dv/dx!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 1
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*winds(jlo:jhi).dvdx(j)}
        mc_npanel_plot,  layout, yinfo, wnddata, panel=2
    endfor
    load_pal, culz
    yinfo.thickness = 3
    yinfo.line_color = culz.black
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).dv_dx}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=2

    yinfo.title = 'dv/dy!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 1
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*winds(jlo:jhi).dvdy(j)}
        mc_npanel_plot,  layout, yinfo, wnddata, panel=3
    endfor
    load_pal, culz
    yinfo.thickness = 3
    yinfo.line_color = culz.black
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).dv_dy}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=3

    goto, PLOTZ_DONE

DIVOR:
    mc_npanel_plot,  layout, yinfo, /setup
    layout.position = [0.18, 0.15, 0.88, 0.92]
    layout.charscale = 1.5
    layout.charthick = 4
    erase, color=culz.white
    layout.panels = 2
    layout.time_axis =1
    layout.xrange = [tcen(jlo), tcen(jhi)]
    layout.title  = title
    layout.erase = 0
    layout.ypad(0:layout.panels-2) = 0.01 + fltarr(layout.panels-1)
    layout.panel_rgb_background = {factor: 0.05, color: culz.cyan, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}
    layout.plot_panel_background = 1

    yinfo.range = [-1.4, 1.4]
    if keyword_set(yrange) then yinfo.range = yrange
    yinfo.charsize = 1.8
;    yinfo.charthick = 4

    yinfo.title = 'Divergence!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 2
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*(winds(jlo:jhi).dudx(j) + winds(jlo:jhi).dvdy(j))}
        yinfo.legend.show = 1
        yinfo.legend.charsize = 1.2
        yinfo.legend.charthick = 3
        yinfo.legend.n_items = mm.rings
        yinfo.legend.item = j -1
        yinfo.legend.text = 'Ring ' + strcompress(string(j, format='(i2)'), /remove_all)
        yinfo.legend.color = yinfo.line_color
        mc_npanel_plot,  layout, yinfo, wnddata, panel=0
        layout.plot_panel_background = 0
    endfor
    load_pal, culz
    yinfo.thickness = 4
    yinfo.line_color = culz.black
    yinfo.legend.item = mm.rings - 1
    yinfo.legend.text = 'All Sky'
    yinfo.legend.color = yinfo.line_color
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).divergence}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=0
    yinfo.legend.show = 0

    layout.panel_rgb_background = {factor: 0.05, color: culz.yellow, $
                                  hilite_factor: 0.18, hilite_color: culz.yellow, xhilite: [-999d, -999d]}
    layout.plot_panel_background = 1
    yinfo.title = 'Vorticity!C !C[1000 x s!U-1!N]'
    yinfo.symsize = 0.2
    pastel_palette, factor=0.8
    for j=1,mm.rings-1 do begin
        yinfo.thickness = 2
        yinfo.line_color = culz.imgmin + float(j)*(culz.imgmax - culz.imgmin -2)/float(mm.rings - 1)
        wnddata = {x: tcen(jlo:jhi), y: 1000*(winds(jlo:jhi).dvdx(j) - winds(jlo:jhi).dudy(j))}
        yinfo.legend.show = 1
        yinfo.legend.charsize = 1.2
        yinfo.legend.charthick = 3
        yinfo.legend.n_items = mm.rings
        yinfo.legend.item = j -1
        yinfo.legend.text = 'Ring ' + strcompress(string(j, format='(i2)'), /remove_all)
        yinfo.legend.color = yinfo.line_color
        mc_npanel_plot,  layout, yinfo, wnddata, panel=1
        layout.plot_panel_background = 0
    endfor
    load_pal, culz
    yinfo.thickness = 4
    yinfo.line_color = culz.black
    yinfo.legend.item = mm.rings - 1
    yinfo.legend.text = 'All Sky'
    yinfo.legend.color = yinfo.line_color
    wnddata = {x: tcen(jlo:jhi), y: windpars(jlo:jhi).vorticity}
    mc_npanel_plot,  layout, yinfo, wnddata, panel=1
    yinfo.legend.show = 0

    goto, PLOTZ_DONE


;---Map LOS winds:
LOSWIND:
    if doing_red    then scale = {yrange: [-300., 300.], auto_scale: 0}
    if doing_sodium then scale = {yrange: [-150., 150.], auto_scale: 0}
    if doing_green  then scale = {yrange: [-150., 150.], auto_scale: 0}
    if keyword_set(yrange) then scale.yrange = yrange
    skymap_settings.scale = scale

    sdi3k_sky_mapper, tlist, tcen, mm, spekfits, skymap_settings, culz, zonemap, cadence=cadence
    load_pal, culz
    goto, PLOTZ_DONE

;---Map temperatures:
TEMPERATURE:
    if doing_green and mm.site_code eq 'HRP' then tprscale = [150, 550.]
    skymap_settings.parameter = 'Temperature'
    skymap_settings.scale.yrange = tprscale
    sdi3k_sky_mapper, tlist, tcen, mm, spekfits, skymap_settings, culz, zonemap, cadence=cadence
    goto, PLOTZ_DONE

;---Map intensities:
INTENSITY:
    skymap_settings.parameter = 'Intensity'
    skymap_settings.scale.yrange = britescale
    sdi3k_sky_mapper, tlist, tcen, mm, spekfits, skymap_settings, culz, zonemap, cadence=cadence
    goto, PLOTZ_DONE

STAGE_RGBMAPS:
RGBMAPS:
;---Set the bitmap size for RGB skymaps:
    geo   = {xsize:  xsize, ysize: ysize}

;---Make temperature/Intensity RGB skymaps:
    yrange = [[tprscale], [britescale], [0., 4e15]]
    rgbmap_scale = {rgbmap_scale, auto_scale: 0, $
                                      yrange: yrange, $
                           menu_configurable: 1, $
                               user_editable: [0,1]}
    rgbmap_geom  = {rgbmap_geom, viewing_from_above: 0, $
                              radius_maps_to_distance: 0, $
                                       north_rotation: 0, $
                                    menu_configurable: 1, $
                                        user_editable: [0,1,2]}
    rgbmap_settings = {scale: rgbmap_scale, parameter: ['Temperature', 'Intensity', 'Background'], black_bgnd: 1, map_view: rgbmap_geom, geometry: geo, records: [jlo, jhi]}

    sdi3k_sky_rgbmap, tlist, tcen, mm, spekfits, rgbmap_settings, zonemap, culz, /no_purple
    goto, PLOTZ_DONE

STAGE_WINDMAPS:
WINDMAPS:
    if n_elements(winds) lt 3 then return

;---Setup the wind mapping:

    if doing_green  then scale = {yrange: 150., auto_scale: 0, rbscale: tprscale, gscale: britescale, pscale: [0, 9e9]}
    if doing_sodium then scale = {yrange: 200., auto_scale: 0, rbscale: tprscale, gscale: britescale, pscale: [0, 9e9]}
    if doing_red    then scale = {yrange: 150., auto_scale: 0, rbscale: tprscale, gscale: britescale, pscale: [0, 9e9]}
    if doing_green  then load_pal, culz, idl=[8,  0]
    if doing_sodium then load_pal, culz, idl=[14, 0]
    if doing_red    then load_pal, culz, idl=[3,  0]

    mcchoice, 'Plot Cadence: ', string(indgen(5)+1), choice, $
                  heading = {text: 'Choose a Plot Cadence', font: 'Helvetica*Bold*Proof*30'}
    cadence  = choice.index + 1
    geo   = {xsize:  xsize, ysize: ysize}
    pp  = 'Map'
    oo  = 'Magnetic Noon at Top'
    oo  = 'Magnetic North at Top'
    windmap_settings = {scale: scale, perspective: pp, orientation: oo, black_bgnd: 1, geometry: geo, records: [jlo, jhi]}
;    sdi3k_read_netcdf_data, filename, metadata=mm_dummy, images=images, cadence=cadence
;    images = images(jlo:jhi)
    sdi3k_wind_mapper, tlist, tcen, mm, winds, windmap_settings, culz, spekfits, zonemap, images=images, cadence=cadence

PLOTZ_DONE:
    if keyword_set(psplot) then begin
       device, /close
       set_plot, 'WIN'
       ps_run = 0
    endif
    plotarr = tvrd(/true)

end





