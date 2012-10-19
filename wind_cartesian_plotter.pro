data_path = 'd:\users\sdi3000\data'
if getenv('computername') eq 'VERTEX'   then data_path = 'G:\'
if getenv('computername') eq 'FLYWHEEL' then data_path = 'D:\users\SDI3000\Data\Spectra\'
ncfile = dialog_pickfile(path=data_path, get_path=data_path, /read, /must_exist)
;ncfile = 'G:\SDI3000\Data\Spectra\PKR 2008_041_Poker_630nm_Red_Sky_Date_02_10.nc'

sdi3k_read_netcdf_data, ncfile, metadata=mm, winds=winds
timelis = dt_tm_mk(js2jd(0d)+1, (winds.start_time + winds.end_time)/2, format='h$:m$:s$')
mcchoice, 'Start Time: ', timelis, choice, $
           heading = {text: 'Start Plot at What Time?', font: 'Helvetica*Bold*Proof*30'}
jlo = choice.index
mcchoice, 'End Time: ', timelis, choice, $
           heading = {text: 'End Plot at What Time?', font: 'Helvetica*Bold*Proof*30'}
jhi = choice.index
sdi3k_read_netcdf_data, ncfile, metadata=mm, winds=winds, spekfits=spekfits, windpars=windpars, range=[jlo,jhi], /preprocess
timlimz = [(min(winds.start_time) + min(winds.end_time))/2, (max(winds.start_time) + max(winds.end_time))/2]
deltime = timlimz(1) - timlimz(0)

load_pal, culz, prop=0.7;, idl=[3,0]

wscales = [10, 15, 20, 30, 50,  80, 100, 120, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1200, 1500, 2000]
mcchoice, 'Wind scale: ', string(wscales, format='(i4)'), choice, $
           heading = {text: 'Scale factor for wind plot?', font: 'Helvetica*Bold*Proof*30'}
wrange = 1.11*wscales(choice.index)
;wrange = 435.
vscale = 0.5*[-wrange, wrange]
vscale = [-87, 87]
divor_range = [-1.1, 1.1]
divor_range = [-0.4, 0.4]

;----Determine whether we're plotting magnetic or geographic wind components:
     mcchoice, 'Plot Coordinates: ', ['Magnetic', 'Geographic'], coords, $
                heading = {text: 'Magnetic or Geographic Wind Components?', font: 'Helvetica*Bold*Proof*30'}
     zonal  = windpars.mag_zonal_wind
     merid  = windpars.mag_meridional_wind
     sigzon = windpars.sigmagzon
     sigmer = windpars.sigmagmer
     if coords.index eq 1 then begin
        zonal  = windpars.geo_zonal_wind
        merid  = windpars.geo_meridional_wind
        sigzon = windpars.siggeozon
        sigmer = windpars.siggeomer
           brg = mm.oval_angle
   geozon_wind = winds.zonal_wind*cos(!dtor*brg)      + winds.meridional_wind*sin(!dtor*brg)
   geomer_wind = winds.meridional_wind*cos(!dtor*brg) - winds.zonal_wind*sin(!dtor*brg)
   winds.zonal_wind = geozon_wind
   winds.meridional_wind = geomer_wind

     endif

;---Setup the HWM stuff:
    hwm_dll  = 'd:\users\conde\main\idl\hwm\nrlhwm93.dll'
    if not(file_test(hwm_dll)) then hwm_dll = 'c:\users\conde\main\idl\hwm\nrlhwm93.dll'
    nhwm     = 32
    hwm_vals = fltarr(2, 3, nhwm)
    tmx       = findgen(nhwm)*deltime/(nhwm-1) + timlimz(0)
    sec   = 0.
    lat   = float(mm.latitude)
    lon   = float(mm.longitude)
    if lon lt 0 then lon = lon+360.
    hwm   = {altitude: 240., F107: 100., ap: 15, plot_hwm: 0}
    obj_edt, hwm, title='HWM Parameters'
    f107  = hwm.f107
    alt   = hwm.altitude
    f107a = f107
    ap    = fltarr(7)
    ap(0) = hwm.ap
    magrot= !dtor*(mm.oval_angle + mm.rotation_from_oval)
    if coords.index eq 1 then magrot = 0.
    mass  = 48L
    t     = fltarr(2)
    d     = fltarr(8)
    delz  = 20.
    if alt lt 180. then delz = 10.
    if alt lt 115. then delz = 5.
    flags = fltarr(25) + 1; Control flags (see HWM-93 FORTRAN source code for usage)
    w     = fltarr(2)
    aplab    = strcompress(string(ap(0), format='(i8)'), /remove_all)
    f107lab  = strcompress(string(f107,  format='(i8)'), /remove_all)
    hwm_lab  = 'NRL-HWM93 model at ' + 'Ap=' + aplab + ', F10.7=' + f107lab

    for idz=-1,1  do begin
        for j=0,nhwm-1 do begin
            yyddd  = long(dt_tm_mk(js2jd(0d)+1, tmx(j), format='y$doy$'))
            js2ymds, tmx(j), yy, mmm, dd, ss
            ss     = float(ss)
            lst    = ss/3600. + lon/15.
            if lst lt 0  then lst = lst + 24.
            if lst gt 24 then lst = lst - 24.
            zv     = alt + idz*delz
            w      = fltarr(2)
            result = call_external(hwm_dll,'nrlhwm93', yyddd, ss, zv, lat, lon, lst, f107a, f107, ap, flags,w)
;            print, yyddd, ss, zv, lat, lon, lst, f107a, f107, ap, flags,w
            hwm_vals(1, idz+1, j) = w(1)*cos(magrot) - w(0)*sin(magrot)
            hwm_vals(0, idz+1, j) = w(0)*cos(magrot) + w(1)*sin(magrot)
;            hwm_vals(0, idz+1, j) = w(0)
;            hwm_vals(1, idz+1, j) = w(1)
            wait, 0.002
         endfor
     endfor
mcchoice, 'Plot as Postscript: ', ['No -- No, use Windows screen plotting', 'Yes -- Do use Postscript'], choice, $
           heading = {text: 'Select Type of Plot Output?', font: 'Helvetica*Bold*Proof*30'}
psplot = choice.index

;--Now start plotting:

    if psplot then begin
       set_plot, 'PS'
       device, /landscape, xsize=26, ysize=20
       device, bits_per_pixel=8, /color, /encapsulated
       device, filename=dialog_pickfile(path='d:\Users\Conde\Main\poker_sdi\', filter='*.eps')
       charscale = 0.6
    endif else begin
       canvas_size = [1400, 900]
       xsize    = canvas_size(0)
       ysize    = canvas_size(1)
       while !d.window ge 0 do wdelete, !d.window
       window, xsize=xsize, ysize=ysize
       charscale=1.
    endelse

erase, color=culz.white

lamlab  = '!4k!3=' + strcompress(string(mm.wavelength_nm, format='(f12.1)'), /remove_all) + ' nm'
posarr = spekfits.velocity
;sdi3k_timesmooth_fits,  posarr, 1.5, mm
;drftfit   = poly_fit((spekfits.start_time + spekfits.end_time)/2 - spekfits(0).start_time, posarr(0,*), 4, measure_errors=1.+spekfits.sigma_velocity(0), yfit=drift, /double)
;vz      = mm.channels_to_velocity*(spekfits.velocity(0) - drift)
vz      = spekfits.velocity(0)


deg_tik,  timlimz, ttvals, nttix, minor, minimum=5
xtn     = dt_tm_mk(js2jd(0d)+1, ttvals, format='h$:m$')

timlimz = timlimz + 0.05*[-deltime, deltime]

    ;---Draw background color fills:
    !p.position = [0.14, 0.09, 0.88, 0.93]
    plot, timlimz, [0., 1.], /noerase, color=culz.white, /xstyle
    wavetlo = ymds2js(2010, 01, 10, 11.0*3600D)
    wavethi = ymds2js(2010, 01, 10, 13.5*3600D)


    x_verts = [wavetlo, wavetlo, wavethi, wavethi, wavetlo]
    y_verts = [0, 1, 1, 0, 0]
    tvlct, r, g, b, /get
    pastel_palette, factor=0.10
    polyfill, x_verts, y_verts, color=culz.slate, /data

    tvlct, r, g, b

    !p.position = [0.14, 0.65, 0.88, 0.93]
    plot, timlimz, [-wrange, wrange], /noerase, color=culz.black, /nodata, /xstyle, /ystyle, $
          xticks=1, xtickname=[' ', ' '], xminor=1, xthick=3, ythick=3, charsize=2.*charscale, charthick=3, $
          title=mm.site + ': ' + dt_tm_mk(js2jd(0d)+1, timlimz(0), format='0d$-n$-Y$') + ', ' + lamlab, $
          ytitle=coords.name + ' Zonal!C!CWind [m s!U-1!N]'
          axis, xaxis=0, xtickname=replicate(' ', 30), color=culz.black, charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          axis, xaxis=1, xtickname=replicate(' ', 30), color=culz.black, charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          oplot, timlimz, [0., 0.], color=culz.ash, linestyle=2

          tvlct, rr, gg, bb, /get
          pastel_palette, factor=0.4
          for j=1,n_elements(winds(0).zonal_wind)-1 do begin
              oplot, (winds.start_time + winds.end_time)/2, winds.zonal_wind(j), $
                      color=culz.imgmin + float(j)*(culz.imgmax - culz.imgmin - 1)/mm.nzones, $
                      psym=6, symsize=0.1, thick=4
          endfor
          tvlct, rr, gg, bb

          mc_oploterr, (windpars.start_time + windpars.end_time)/2, zonal, sigzon, $
                        bar_color=culz.black, thick=4, symbol_color=culz.blue, line_color=culz.blue, psym=6, symsize=.3, errthick=3


        tm       = timlimz(0) + 0.97*deltime
        ym       = 0.8*wrange
        if hwm.plot_hwm then begin
        for idz=-1,1 do begin
            thk = 1 + (idz eq 0)
            oplot, tmx, hwm_vals(1, idz+1, *), thick=thk, linestyle=1, color=culz.orange
            zlab = strcompress(string(alt + idz*delz, format='(i8)'), /remove_all)
            t1   = timlimz(1) - 0.03*deltime
            sgn  = 0
            if idz ne 0 then sgn  = (hwm_vals(0, idz+1, nhwm-1) - hwm_vals(0, 1, nhwm-1))/abs(hwm_vals(0, idz+1, nhwm-1) - hwm_vals(0, 1, nhwm-1))
            y1   = hwm_vals(1, 1, nhwm-4) + abs(idz)*sgn*wrange*0.09; - 0.06*wrange
            xyouts, t1, y1, zlab +' km', charsize=0.9*charscale, color=culz.orange, align=1
        endfor
        xyouts, tm, ym, hwm_lab, charsize=1.2*charscale, color=culz.orange, align=1
        endif
    !p.position = [0.14, 0.37, 0.88, 0.65]
    plot, timlimz, [-wrange, wrange], /noerase, color=culz.black, /nodata, /xstyle, /ystyle, $
          xticks=1, xtickname=[' ', ' '], xminor=1, xthick=3, ythick=3, charsize=2*charscale, charthick=3, $
          ytitle=coords.name + ' Meridional!C!CWind [m s!U-1!N]'
          axis, xaxis=0, xtickname=replicate(' ', 30), color=culz.black, charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          axis, xaxis=1, xtickname=replicate(' ', 30), color=culz.black, charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          oplot, timlimz, [0., 0.], color=culz.ash, linestyle=2

          tvlct, rr, gg, bb, /get
          pastel_palette, factor=0.4
          for j=1,n_elements(winds(0).meridional_wind)-1 do begin
              oplot, (winds.start_time + winds.end_time)/2, winds.meridional_wind(j), $
                      color=culz.imgmin + float(j)*(culz.imgmax - culz.imgmin - 1)/mm.nzones, $
                      psym=6, symsize=0.1, thick=4
          endfor
          tvlct, rr, gg, bb

          mc_oploterr, (windpars.start_time + windpars.end_time)/2, merid, sigmer, $
                        bar_color=culz.black, thick=4, symbol_color=culz.red, line_color=culz.red, psym=6, symsize=.3, errthick=3

        tm       = timlimz(0) + 0.97*deltime
        ym       = 0.8*wrange
        if hwm.plot_hwm then begin
                for idz=-1,1 do begin
            thk = 1 + (idz eq 0)
            oplot, tmx, hwm_vals(0, idz+1, *), thick=thk, linestyle=1, color=culz.orange
            zlab = strcompress(string(alt + idz*delz, format='(i8)'), /remove_all)
            t1   = timlimz(1) - 0.03*deltime
            sgn  = 0
            if idz ne 0 then sgn  = (hwm_vals(0, idz+1, nhwm-1) - hwm_vals(0, 1, nhwm-1))/abs(hwm_vals(0, idz+1, nhwm-1) - hwm_vals(0, 1, nhwm-1))
            y1   = hwm_vals(0, 1, nhwm-4) + abs(idz)*sgn*wrange*0.09; - 0.06*wrange
            xyouts, t1, y1, zlab +' km', charsize=0.9*charscale, color=culz.orange, align=1
        endfor
        xyouts, tm, ym, hwm_lab, charsize=1.2*charscale, color=culz.orange, align=1
        endif

    !p.position = [0.14, 0.09, 0.88, 0.37]
     plot,timlimz, vscale, /noerase, color=culz.black, /nodata, /xstyle, ystyle=9, $
          xticks=1, xtickname=[' ', ' '], xminor=1, xthick=3, ythick=3, charsize=2*charscale, charthick=3, $
          ytitle='Vertical!C!CWind [m s!U-1!N]'
          axis, xaxis=0, xtickname=xtn, color=culz.black, xtitle='Time [Hours UT]', charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          axis, xaxis=1, xtickname=replicate(' ', 30), color=culz.black, charsize=2.*charscale, charthick=3, $
                xminor=minor, xticks=nttix, xticklen=0.02, xtickv=ttvals, xthick=4
          oplot, timlimz, [0., 0.], color=culz.ash, linestyle=2
;          mc_oploterr, (spekfits.start_time + spekfits.end_time)/2, vz, mm.channels_to_velocity*spekfits.sigma_velocity(0), $
;                        bar_color=culz.ash, thick=1, symbol_color=culz.black, line_color=culz.olive, /two, psym=6, symsize=.3
          mc_oploterr, (spekfits.start_time + spekfits.end_time)/2, vz, spekfits.sigma_velocity(0), $
                        bar_color=culz.ash, thick=4, symbol_color=culz.black, line_color=culz.olive, /two, psym=6, symsize=.3, errthick=3
    plot, timlimz, divor_range, /noerase, color=culz.white, /nodata, /xstyle, ystyle=9, $
          xticks=1, xtickname=[' ', ' '], xminor=1, xthick=3, ythick=3, charsize=2*charscale, charthick=3, $
          ytitle='Vertical!C!CWind [m s!U-1!N]'
    axis, yaxis=1, color=culz.rose, $
          ytitle='1000 x Divergence [s!U-1!N]!C!C1000 x Vorticity [s!U-1!N]', $
          charsize=2*charscale, charthick=3, ythick=3
    axis, yaxis=1, color=culz.slate, $
          ytitle='1000 x Divergence [s!U-1!N]', $
          charsize=2*charscale, charthick=3, ythick=3
    axis, yaxis=1, color=culz.black, $
          charsize=2*charscale, charthick=3, ythick=3
    oplot,(windpars.start_time + windpars.end_time)/2, windpars.divergence,    color=culz.slate, thick=4
    oplot,(windpars.start_time + windpars.end_time)/2, windpars.vorticity,     color=culz.rose,  thick=4
;    xyouts, timlimz(0) + 0.08*(timlimz(1) - timlimz(0)), -0.9, 'Vertical Wind', color=culz.olive, charsize=1.5, charthick=2
;    xyouts, timlimz(0) + 0.22*(timlimz(1) - timlimz(0)), -0.9, 'Divergence',    color=culz.slate, charsize=1.5, charthick=2
;    xyouts, timlimz(0) + 0.34*(timlimz(1) - timlimz(0)), -0.9, 'Vorticity',     color=culz.rose,  charsize=1.5, charthick=2

    plot, timlimz, vscale, /noerase, color=culz.black, /nodata, /xstyle, ystyle=9, $
          xticks=1, xtickname=[' ', ' '], xminor=1, xthick=3, ythick=3, charsize=2*charscale, charthick=3, $
          ytitle='Vertical!C!CWind [m s!U-1!N]', ytick_get=tkv
    axis, yaxis=0, color=culz.olive, $
          ytitle='Vertical!C!CWind [m s!U-1!N]', $
          charsize=2*charscale, charthick=3, ythick=3, ytickv=tkv, yminor=1, yticks=n_elements(tkv)+1
    axis, yaxis=0, color=culz.black, $
          charsize=2*charscale, charthick=3, ythick=3, ytickv=tkv, yminor=1, yticks=n_elements(tkv)+1

    xyouts, 0.89, 0.91, 'Note: Colored !Cbackground points !Cdenote wind !Cmeasurements in !Ceach zone. ' + $
                        'Error !Cbars denote !C1-!7r!3 standard !Cdeviations (across !Call zones) of !Czonal and !Cmeridional wind !Ccomponents at !Ceach observation !Ctime.', $
            color=culz.black, charsize=1.3*charscale, charthick=2, /normal
    empty
!p.position = 0
sdi3k_read_netcdf_data, ncfile, /close
    if psplot then begin
       device, /close
       set_plot, 'WIN'
    endif
end
