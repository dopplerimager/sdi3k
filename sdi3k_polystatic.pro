function sdi3k_2d_interp, valbefore, valafter, timebefore, timeafter, when
         return, valbefore + (valafter - valbefore)*(when - timebefore)/(timeafter - timebefore)
end

function overead_asc, fname, js=js, ascrot=ascrot, return_time=return_time
    if strpos(strupcase(fname), '.FITS') gt 0 then begin
       asc = readfits(fname, hdr)
       asc = asc - min(asc)
       asc = reverse(asc, 1)
       tt  = str_sep(hdr(40), "'")
       tt  = tt(1)
       yr  = fix(strmid(tt, 0, 4))
       mm  = fix(strmid(tt, 5, 2))
       dd  = fix(strmid(tt, 8, 2))
       ss  = 3600L*strmid(tt, 11, 2) + 60L*strmid(tt, 14, 2) + strmid(tt, 17, 2)
    endif else begin
      read_jpeg, fname, asc
      asc = reverse(asc, 2)
      hr  = strmid(fname, strlen(fname(0)) - 10, 2)
      mnt = strmid(fname, strlen(fname(0)) - 8, 2)
      sec = strmid(fname, strlen(fname(0)) - 6, 2)
      yr  = strmid(fname, strpos(fname, 'Year_') + 5, 4)
      mm  = strmid(fname, strpos(fname, 'Year_') + 10, 2)
      dd  = strmid(fname, strpos(fname, 'Year_') + 12, 2)
      ss  = 3600.* hr + 60.*mnt + sec
    endelse
    asc = long(asc)
    if keyword_set(ascrot) then asc = rot(asc, ascrot)
    js = ymds2js(yr, mm, dd, ss)
    if keyword_set(return_time) then return, js
    return, long(asc)
end

pro sid3k_get_asc_b4_after, flis, js, b4_aft
    b4_aft = -1
    lidx = [0, n_elements(flis)-1]
    lijs = [overead_asc(flis(lidx(0)), /return_time), overead_asc(flis(lidx(1)), /return_time)]
    if overead_asc(flis(lidx(0)), /return_time) gt js then return
    if overead_asc(flis(lidx(1)), /return_time) lt js then return

    while max(lidx) - min(lidx) gt 1 do begin
          midjs = overead_asc(flis(total(lidx)/2), /return_time)
          if midjs lt js then lidx(0) = total(lidx)/2 else lidx(1) = total(lidx)/2
    endwhile
    b4_aft = lidx
end

pro sdi3k_load_asc_fits, flis, js, image, ascrot
    sid3k_get_asc_b4_after, flis, js, b4_aft
    if n_elements(b4_aft) lt 2 then begin
       image = 0.*overead_asc(flis(0))
    endif else begin
       ib4   = overead_asc(flis(b4_aft(0)), js=jsb4,  ascrot=ascrot)
       iaft  = overead_asc(flis(b4_aft(1)), js=jsaft, ascrot=ascrot)
       image = sdi3k_2d_interp(ib4, iaft, jsb4, jsaft, js)
    endelse
end

pro sdi3k_asc_refpoints, flis, refpoints, culz, ascrot=ascrot
    if strpos(strupcase(flis(0)), '.FITS') gt 0 then $
       timelis = strmid(flis, strlen(flis(0)) - 15, 2) + '-' + strmid(flis, strlen(flis(0)) - 13, 2) + '-' + strmid(flis, strlen(flis(0)) - 11, 2) else $
       timelis = strmid(flis, strlen(flis(0)) - 10, 2) + '-' + strmid(flis, strlen(flis(0)) - 8, 2)  + '-' + strmid(flis, strlen(flis(0)) - 6, 2)
    subsidx = 10*(indgen(n_elements(timelis)/10))
    timelis = timelis(subsidx)
    mcchoice, 'Example Image: ', timelis, choice, $
               heading = {text: 'Choose an ASC image to use for scaling points.', font: 'Helvetica*Bold*Proof*30'}

       this_asc = overead_asc(flis(subsidx(choice.index)), ascrot=ascrot)
       window, 5, xsize=512, ysize=512
       xcen=256
       ycen=256

       tv, culz.greymin+bytscl(hist_equal(this_asc), top=(culz.greymax - culz.greymin))
       xyouts, xcen, ycen, "Click cursor on the LEFT horizon", $
     color=culz.green, /device, align=0.5, charthick=2, charsize=1.2
       cursor, x, y, 4, /device
       lox = x
       tv, culz.greymin+bytscl(hist_equal(this_asc), top=(culz.greymax - culz.greymin))
       xyouts, xcen, ycen, "Click cursor on the RIGHT horizon", $
     color=culz.orange, /device, align=0.5, charthick=2, charsize=1.2
       cursor, x, y, 4, /device
       hix = x
       tv, culz.greymin+bytscl(hist_equal(this_asc), top=(culz.greymax - culz.greymin))
       xyouts, xcen, ycen, "Click cursor on the BOTTOM horizon", $
     color=culz.green, /device, align=0.5, charthick=2, charsize=1.2
       cursor, x, y, 4, /device
       loy = y
       tv, culz.greymin+bytscl(hist_equal(this_asc), top=(culz.greymax - culz.greymin))
       xyouts, xcen, ycen, "Click cursor on the TOP horizon", $
     color=culz.orange, /device, align=0.5, charthick=2, charsize=1.2
       cursor, x, y, 4, /device
       hiy = y
       wdelete, 5

       zang = 360.*mc_dist(n_elements(this_asc(*,0)), n_elements(this_asc(0, *)), (lox+hix)/2, (loy+hiy)/2, x=xx, y=yy)/(lox + hix + loy + hiy)
       azi = atan(yy, xx)
       azi = 180 + rotate(reverse(azi, 1), 1)/!dtor
       rdist  = 100.*tan(!dtor*zang)
       xdist  = rdist*sin(!dtor*azi)
       ydist  = rdist*cos(!dtor*azi)
       useful = where(zang lt 85.)
       useord = sort(zang(useful))
       useful = useful(reverse(useord))
       dims   = [n_elements(this_asc(*,0)), n_elements(this_asc(0, *))]
       refpoints = {dims: dims, horizon: [lox, hix, loy, hiy], zang: zang, azimuth: azi, xdist: xdist, ydist: ydist, useful: useful}
end

pro sdi3k_get_asc_indices, xpix, ypix, mm, refpoints, mapvec, ascvec

    dummy   = mc_dist(xpix, ypix, 0., 0., x=xx, y=yy)
    lonlats = convert_coord(xx, yy, /device, /to_data)
    dx      = 111.12*(lonlats(0, *) - mm.longitude)*cos(!dtor*mm.latitude)
    dy      = 111.12*(lonlats(1, *) - mm.latitude)
    rads    = sqrt(dx^2 + dy^2)
    use     = where(rads lt 450.)
    zang    = atan(rads(use)/105.)/!dtor
    pixrad  = 0.25*total(refpoints.horizon)*zang/(90.)
    maxrad  = max(rads(use))
    azi     = atan(dy(use), dx(use))
    radhere = sqrt((dx(use))^2 + (dy(use))^2)
    xx      = pixrad*dx(use)/radhere + (refpoints.horizon(0) + refpoints.horizon(1))/2
    yy      = pixrad*dy(use)/radhere + (refpoints.horizon(2) + refpoints.horizon(3))/2
    mapvec  = use
    ascvec  = long(yy)*refpoints.dims(0) + xx
end

pro sdi3k_get_tmp_indices_old, xpix, ypix, mm, tmp_coords, tmapvec, twght

    if abs(mm.wavelength_nm - 630.) lt 1. then height = 240. else height=120.

    dummy   = mc_dist(xpix, ypix, 0., 0., x=xx, y=yy)
    lonlats = convert_coord(xx, yy, /device, /to_data)
    dx      = 111.12*(lonlats(0, *) - mm.longitude)*cos(!dtor*mm.latitude)
    dy      = 111.12*(lonlats(1, *) - mm.latitude)
    rads    = sqrt(dx^2 + dy^2)
    zang    = atan(rads/height)/!dtor
    use     = where(zang lt 75.)
    print, 'Building Temperature Map...'
    twght   = {idx: lonarr(8), weight: fltarr(8), totweight: 0.}
    twght   = replicate(twght, n_elements(use))
    progressBar = Obj_New("SHOWPROGRESS", message='Percent Completion', title='Indexing...')
    progressBar->Start

    for j=0L, n_elements(use) - 1 do begin
        dd = sqrt(((lonlats(0, use(j)) - tmp_coords.lons)*cos(!dtor*mm.latitude))^2 + (lonlats(1, use(j)) - tmp_coords.lats)^2)
        close = where(dd lt 3.0)
        best  = sort(dd(close))
        best  = close(best(0:7))
        twght(j).idx    = best
        twght(j).weight = exp(-(dd(best)/0.4)^2)
        twght(j).totweight = total(twght(j).weight)
        wait, 0.000001
        if j mod 1000 eq 0 then progressbar->update, 100.*j/float(n_elements(use))
    endfor
    tmapvec = use
    progressBar->Destroy
    Obj_Destroy, progressBar
end


pro sdi3k_get_tmp_indices, xpix, ypix, mm, tmp_coords, tmapvec, twght
    nsum = 10
    if abs(mm.wavelength_nm - 630.) lt 1. then height = 240. else height=120.

    dummy   = mc_dist(xpix, ypix, 0., 0., x=xx, y=yy)
    lonlats = convert_coord(xx, yy, /device, /to_data)
    dx      = 111.12*(lonlats(0, *) - mm.longitude)*cos(!dtor*mm.latitude)
    dy      = 111.12*(lonlats(1, *) - mm.latitude)
    rads    = sqrt(dx^2 + dy^2)
    zang    = atan(rads/height)/!dtor
    use     = where(zang lt 75.)
    print, 'Building Temperature Map...'
    twght   = {idx: lonarr(nsum), weight: fltarr(nsum), totweight: 0.}
    twght   = replicate(twght, n_elements(use))
    progressBar = Obj_New("SHOWPROGRESS", message='Percent Completion', title='Indexing...')
    progressBar->Start

    for j=0L, n_elements(use) - 1 do begin
        dd = sqrt(((lonlats(0, use(j)) - tmp_coords.lons)*cos(!dtor*mm.latitude))^2 + (lonlats(1, use(j)) - tmp_coords.lats)^2)
        best  = sort(dd)
        best  = best(0:nsum-1)
        nearest = 1.2*dd(best(0)) > 0.60
        twght(j).idx    = best
        twght(j).weight = exp(-(dd(best)/(nearest))^2)
        twght(j).totweight = total(twght(j).weight)
        wait, 0.000001
        if j mod 1000 eq 0 then progressbar->update, 100.*j/float(n_elements(use))
    endfor
    tmapvec = use
    progressBar->Destroy
    Obj_Destroy, progressBar
end

pro  sdi3k_map_windvex, xlons, ylats, zonwind, merwind, base_color, arrow_color, thick=thick
     if not(keyword_set(thick)) then thick=1
     for j=0L, n_elements(xlons)-1 do begin
         lat = ylats(j)
         lon = xlons(j)
         hh  = convert_coord(lon, lat, /data, /to_device)
         dlat = convert_coord(lon, lat + 1., /data, /to_device)
         dlon = convert_coord(lon + 1., lat, /data, /to_device)
         lathat = [dlat(0) - hh(0), dlat(1) - hh(1)]
         lonhat = [dlon(0) - hh(0), dlon(1) - hh(1)]/cos(!dtor*lat)
         plots, lon, lat, psym=1, symsize=0.25, thick=3, color=base_color
         xend = hh(0) + 0.008*(zonwind(j)*lonhat(0) + merwind(j)*lathat(0))
         yend = hh(1) + 0.008*(zonwind(j)*lonhat(1) + merwind(j)*lathat(1))
         arrow, hh(0), hh(1), xend, yend, color=arrow_color, thick=thick, hthick=thick, hsize=10
     endfor
end


;-------------------------------------------------------------------------------------------
;
;  Main program starts here:

;   goto, test_spot
cpu,TPOOL_MIN_ELTS =2000
load_pal, culz
interp_settings =  {charspeed_kmps: 0.005, dist_par_km: 4.}
addasc = 1
order = 1
xpix   = 900
ypix   = 900
xshrink = 1.
tstep = 600.
drive = 'c'
vertz = 0
nterms  = (vertz + 2)*total(indgen(order+1)+1)
used    = intarr(nterms) + 1
used(3) = 0
tmp_interp_settings =  {charspeed_kmps: 0.10, dist_par_km: 80.}
site_labels    = replicate({name: 'PKR', lon:    -147.430, lat: 65.1192, color: culz.wheat, thick: 5, psym: 6, symsize: 2, charsize: 2, charthick: 3, doprint: 1}, 2)
site_labels(1).name = 'HAARP'
site_labels(1).lon  = -145.31944
site_labels(1).lat  = 62.30194


;---Read the SDI data:
;nc_path = 'D:\users\SDI3000\Data\'
nc_path = 'C:\Users\Conde\Main\Poker_SDI\data\'
nc_list = nc_path + ['Poker\PKR 2010_024_Poker_630nm_Red_Sky_Date_01_24.nc', $
                     'Haarp\HRP_2010_024_Elvey_630nm_Red_Sky_Date_01_24.nc']

jslo = -9d99
jshi =  9d99


for j=0, n_elements(nc_list)-1 do begin
    sdi3k_read_netcdf_data,  nc_list(j), metadata=mm, spekfits=spekfits, winds=winds, /preprocess
    timez = (spekfits.start_time + spekfits.end_time)/2.
    if min(timez) gt jslo then jslo = min(timez)
    if max(timez) lt jshi then jshi = max(timez)
    if j eq 0 then begin
       mmarr = mm
       spekfarr = spekfits
       windarr = winds
       nidxarr = intarr(mm.maxrec) + j
       longitarr = winds(0).zone_longitudes
       latitarr = winds(0).zone_latitudes
    endif else begin
       mmarr = [mmarr, mm]
       spekfarr = [spekfarr, spekfits]
       windarr = [windarr, winds]
       nidxarr = [nidxarr, intarr(mm.maxrec) + j]
       longitarr = [longitarr, winds(0).zone_longitudes]
       latitarr  = [latitarr,  winds(0).zone_latitudes]
           endelse
endfor

;---Setup the directory to hold the stills:
dirz =   alldisk_findfiles('\idl_mpegs')
if dirz.nfound eq 0 then return
drivenames = tag_names(dirz)
if dirz.nfound eq 1 then didx = 0 else begin
   drivenames = drivenames(0:dirz.nfound-1)
   mcchoice, 'Drive for Stills?', drivenames, choice
   didx = choice.index
endelse
drive = drivenames(didx)

prefix = ':\idl_mpegs\stills_'
lamstring = strcompress(string(fix(10*mm.wavelength_nm)), /remove_all)
still_folder = drive + prefix + lamstring + '_' + dt_tm_mk(js2jd(0d)+1, jslo, format='0d$_N$_Y$')
if !version.release ne '5.2' then file_mkdir, still_folder else spawn, 'mkdir ' + still_folder

mcchoice, 'Clear Stills?', ['No - Do NOT clear existing stills', 'Yes - DO clear existing stills'], choice, help='WARNING - Answering YES will delete any existing stills in c:\idl_mpegs\stills'
if choice.index then begin
   mcchoice, 'Confirm Clear Stills', ['Cancel - NO! Do NOT clear existing stills', 'Confirm - YES! Really DO clear existing stills'], choice, help='WARNING - Answering YES will delete any existing stills in c:\idl_mpegs\stills'
   if choice.index then begin
      victims  = findfile(still_folder + '\sdi_*.png')
      if victims(0) ne '' then begin
         for j=0,n_elements(victims)-1 do file_delete, victims(j), /verbose
      endif
   endif
endif
tnow = jslo

doing_red = 1
latreach = 8.
lonreach = 16.5
if abs(mmarr(0).wavelength_nm - 557.7) lt 1. then begin
   latreach = 4.5
   lonreach = 9.
   doing_red = 0
endif

latmin = min(mm.latitude)  - latreach
lonmin = min(mm.longitude) - lonreach
latmax = max(mm.latitude)  + latreach
lonmax = max(mm.longitude) + lonreach
cenlat = (latmin + latmax)/2.
cenlon = (lonmin + lonmax)/2.

;--Generate a 1D vector of spacetime coordinates of the observations:
   a_coord = {lats: 0., lons:0., js: 0d, idx: 0}
   obs_coords = replicate(a_coord, n_elements(windarr.meridional_wind))
   obs = replicate(a_coord, n_elements(windarr.meridional_wind))
   lats  = reform(windarr.zone_latitudes,  n_elements(windarr.zone_latitudes))
   lons  = reform(windarr.zone_longitudes, n_elements(windarr.zone_longitudes))
   jsi   = (windarr.start_time + windarr.end_time)/2.
   js    = replicate({js: jsi}, mmarr(0).nzones)
   js    = transpose(js.js)
   js    = reform(js, n_elements(js))
   idx   = replicate({idx: nidxarr}, mmarr(0).nzones)
   idx   = transpose(idx.idx)
   idx   = reform(idx, n_elements(idx))
   obs_coords.lats = lats
   obs_coords.lons = lons
   obs_coords.js   = js
   obs_coords.idx  = idx

;--Generate a 1D vector of observation values:
   loswind = reform(spekfarr.velocity, n_elements(spekfarr.velocity))
   siglosw = reform(spekfarr.sigma_velocity, n_elements(spekfarr.sigma_velocity))
   zen     = reform(windarr.zeniths, n_elements(windarr.zeniths))
   azi     = reform(windarr.azimuths, n_elements(windarr.azimuths))
   an_obs = {loswind: 0., siglosw: 0.}
   obsarr = replicate(an_obs, n_elements(spekfarr.velocity))
   obsarr.loswind = loswind
   obsarr.siglosw = siglosw

;--Generate a 1D vector of spacetime coordinates of the wind output samples:
    nx     = 17
    ny     = 17
    dy     = 1.1*(max(winds(0).zone_latitudes) - min(winds(0).zone_latitudes))/float(nx)
    dx     = dy
    xx     = float(transpose(lindgen(ny,nx)/ny) - nx/2)
    yy     = float(lindgen(nx,ny)/nx - ny/2)
    xx     = cenlon + dx*xx/cos(!dtor*cenlat)
    yy     = cenlat + dy*yy

    keep   = lonarr(n_elements(xx))
    for j=0,n_elements(keep)-1 do begin
        degdist = sqrt((yy(j) - latitarr)^2 + ((xx(j) - longitarr)*cos(!dtor*cenlat))^2)
        nearby  = sort(degdist)
        nearby  = nearby(0:10)
        nearrad = sqrt((cenlat - latitarr(nearby))^2 + ((cenlon - longitarr(nearby))*cos(!dtor*cenlat))^2)
        thisrad = sqrt((cenlat - yy(j))^2 + ((cenlon - xx(j))*cos(!dtor*cenlat))^2)
;        if j gt 64 then stop
;        print, xx(j), yy(j), thisrad, nearrad
        if thisrad le max(nearrad) then keep(j) = 1
    endfor
    kp     = where(keep gt 0)


    xx     = xx(kp)
    yy     = yy(kp)

    out_coords = replicate(a_coord, n_elements(xx))
    out_coords.lons  = xx
    out_coords.lats  = yy
    out_coords.js    = min(obs_coords.js)

;--Get list of asc file names to use:
ascz    = alldisk_findfiles('\poker_asc_data')
ascdisk = tag_names(ascz)
ascdisk = ascdisk(0)
addasc  = ascdisk ne 'NFOUND'
if addasc then begin
   ascdir   = dialog_pickfile(path=ascdisk + ':\poker_asc_data', /dir)
   flis     = findfile(ascdir + '*.fits')
   if flis(0) eq '' then begin
      flis     = findfile(ascdir + '*.jpg')
      asc_gain = 1.
   endif
   flis     = flis(sort(flis))
   if flis(0) eq '' then stop
   sdi3k_asc_refpoints, flis, refpoints, culz, ascrot=mmarr(0).oval_angle
endif

;--Make flare locations:
nxf = 5
nyf = 9
xx   = float(transpose(lindgen(nyf,nxf)/nyf) - nxf/2)
yy   = float(lindgen(nxf, nyf)/nxf           - nyf/2)
fkp  = where(sqrt(xx*xx + 0.3*yy*yy) le 1.2*abs(xx(0,0)))
scl = 0.6
if doing_red then scl = 0.75
flarelons = cenlon + scl*xx*(lonmax - lonmin)/(2*max(xx))
flarelats = cenlat  + scl*yy*(latmax - latmin)/(2*max(yy))
flares    = {lon: flarelons(fkp), lat: flarelats(fkp), use: intarr(n_elements(fkp))+1, js: 0D, rdist: -9999.}
flarebase = flares

scale = 400.

Plot_start:
while !d.window gt 0 do wdelete, !d.window
window, xsize=xpix, ysize=ypix
load_pal, culz
pastel_palette, factor=0.4, dir='dark'
erase, color=culz.bground
map_set, 90., cenlon, 0., /stereo, color=culz.wheat, continents=1, clip=1, $
         limit=[latmin, lonmin, latmax, lonmax], mlinethick=coast_thick, /noborder, /noerase, /hires

MAP_CONTINENTS,/FILL_CONTINENTS,color=culz.chocolate, /horizon, /hires
load_pal, culz
MAP_CONTINENTS,color=culz.wheat, /horizon, /hires

while tnow le jshi do begin
    if !d.x_size ne xpix or !d.y_size ne ypix then begin
       while !d.window gt 0 do wdelete, !d.window
       window, xsize=xpix, ysize=ypix
    endif

    for j=0, n_elements(nc_list)-1 do begin
        thisidx = where(obs_coords.idx eq j)
        deltime = abs(obs_coords(thisidx).js - tnow)
        thistime = where(deltime eq min(deltime))
        thistime = thisidx(thistime)
        nidx     = intarr(n_elements(thistime)) + j
        if j eq 0 then begin
           these = thistime
           thidx = nidx
        endif else begin
           these = [these, thistime]
           thidx = [thidx, nidx]
        endelse
    endfor

    outarr = replicate(obsarr(0), n_elements(these))
    for j=0, n_elements(nc_list)-1 do begin
        oc = obs_coords(these(where(thidx eq j)))
        oc.js = tnow
        sdi3k_spacetime_interpol, obs_coords, oc, obsarr, thout, interp_settings
        outarr(where(thidx eq j)) = thout
    endfor

;    sdi3k_spacetime_interpol, obs_coords, oc, obsarr, outarr, interp_settings


;#####################################################################
; Cut and paste from UCL_VEX:
        print, 'Fitting winds for ', + dt_tm_mak(js2jd(0d)+1, tnow, format='Y$-n$-0d$ h$:m$')
;-------Note that ultimately I need to express x and y in km, not degrees.

;goto, simple_fit

     for j=0,n_elements(out_coords)-1 do begin
              degdist = 0.5 + (sqrt((obs_coords(these).lats - out_coords(j).lats)^2 + ((obs_coords(these).lons - out_coords(j).lons)*cos(!dtor*out_coords(j).lats))^2))
              sdi3k_polywind, outarr.loswind, degdist*outarr.siglosw, $
                              obs_coords(these).lons - out_coords(j).lons, obs_coords(these).lats - out_coords(j).lats, $
                              azi(these), zen(these), order, fitpars, $
                              zonal,   meridional, vertical, sigzon,  sigmer,   sigver,   quality, /horizontal_only, used=used

              if j eq 0 then begin
                 sdi3k_get_poly_wind, out_coords.lons - out_coords(0).lons, out_coords.lats - out_coords(0).lats, order, fitpars, fitwinds
              endif else begin
                 sdi3k_get_poly_wind, 0., 0., order, fitpars, fitw
                 fitwinds.zonal(j) = fitw.zonal(0)
                 fitwinds.meridional(j) = fitw.meridional(0)
              endelse
     endfor

goto, fit_done
simple_fit:
     sdi3k_polywind, outarr.loswind, outarr.siglosw, $
                     obs_coords(these).lons - cenlon, obs_coords(these).lats - cenlat, $
                     azi(these), zen(these), order, fitpars, $
                     zonal,   meridional, vertical, sigzon,  sigmer,   sigver,   quality, /horizontal_only, used=used
     sdi3k_get_poly_wind, out_coords.lons - cenlon, out_coords.lats - cenlat, order, fitpars, fitwinds
fit_done:

    sdi3k_load_asc_fits,      flis, tnow, image, mmarr(0).oval_angle

;------Propagate flare trail:
       for j=n_elements(flares)-1, 0, -1 do begin
           for k=0,n_elements(flares(0).lon)-1 do begin
               if ~flares(j).use(k) then goto, do_not_use
               a_coord.lats = flares(j).lat(k)
               a_coord.lons = flares(j).lon(k)
               a_coord.js   = tnow
;               sdi3k_spacetime_interpol, obs_coords, a_coord, obsarr, outflare, interp_settings
              degdist = 0.5 + (sqrt((obs_coords(these).lats - out_coords(j).lats)^2 + ((obs_coords(these).lons - out_coords(j).lons)*cos(!dtor*out_coords(j).lats))^2))
              sdi3k_polywind, outarr.loswind, degdist*outarr.siglosw, $
                              obs_coords(these).lons - flares(j).lon(k), obs_coords(these).lats - flares(j).lat(k), $
                              azi(these), zen(these), order, flarepars, $
                              zonal,   meridional, vertical, sigzon,  sigmer,   sigver,   quality, /horizontal_only, used=used
              sdi3k_get_poly_wind, 0., 0., order, flarepars, outflare

               rdist = sqrt((obs_coords(these).lats - flares(j).lat(k))^2 + (cos(!dtor*flares(j).lat(k))*(obs_coords(these).lons - flares(j).lon(k)))^2)
               if min(rdist) lt 0.10*(latmax - latmin) then begin
                  zon  = outflare.zonal(0)
                  mer  = outflare.meridional(0)
                  flares(j).lat(k) = flares(j).lat(k) + 0.001*mer*tstep/111.
                  flares(j).lon(k) = flares(j).lon(k) + 0.001*zon*tstep/(111.*cos(!dtor*flares(j).lat(k)))
;                  flares(j).js  = flares(j).js + tstep
               endif else flares(j:*).use(k) = 0
do_not_use:
           endfor
           if total(flares(j).use) eq 0 then veceldel, flares, j
       endfor
       young_enough = where(tnow - flares.js lt 14400, nnn)
       if nnn gt 10 then flares = flares(young_enough)

     erase, color=culz.black
     map_set, 90., cenlon, 0., /stereo, color=culz.wheat, continents=1, clip=1, $
              limit=[latmin, lonmin, latmax, lonmax], mlinethick=coast_thick, /noborder, /noerase, /hires
     load_pal, culz

;------Plot flare trails:
       if tnow - max(flares.js) gt 120 then begin
          flares = [{lon: flarelons(fkp), lat: flarelats(fkp), use: intarr(n_elements(fkp))+1, js: 0D, rdist: -9999.}, flares]
          flares(0).js = tnow
       endif
       for k=0,n_elements(flares(0).lon)-1 do begin
           plots, flares(0).lon(k), flares(0).lat(k), /data
           for j=1,n_elements(flares.lon(0))-1 do begin
               fcul = culz.imgmin + 40+ (culz.imgmax - culz.imgmin - 40 - 2)*(9.0*float(k)/n_elements(fkp) mod 1)
               if flares(j).use(k) then plots, flares(j).lon(k), flares(j).lat(k), /data, /cont, color=fcul, thick=1+6*(cos(0.95*flares(j).js))^2
           endfor
           spot = convert_coord(flarebase(0).lon(k), flarebase(0).lat(k), /data, /to_device)
           Tvcircle, 2, spot(0), spot(1), fcul, thick=4, color=fcul
       endfor
       flare_img = tvrd(/true)

       erase, color=culz.white
       loadct, 0
       for k=0,n_elements(flares(0).lon)-1 do begin
           plots, flares(0).lon(k), flares(0).lat(k), /data
           for j=1,n_elements(flares.lon(0))-1 do begin
               fage = (tnow - flares(j).js)/14400. < 1
               fage = 255.*(fage + fage^2)/2.
               if flares(j).use(k) then plots, flares(j).lon(k), flares(j).lat(k), /data, /cont, color=fage, thick=1+6*(cos(0.95*flares(j).js))^2
           endfor
           spot = convert_coord(flarebase(0).lon(k), flarebase(0).lat(k), /data, /to_device)
           Tvcircle, 2, spot(0), spot(1), 0, color=0, thick=4
           Tvcircle, 1, spot(0), spot(1), 255, thick=2, color=255
       endfor
       flare_age = tvrd(/true)
     load_pal, culz
     pastel_palette, factor=0.4, dir='dark'
     erase, color=culz.bground
     MAP_CONTINENTS,/FILL_CONTINENTS,color=culz.chocolate, /horizon, /hires
     load_pal, culz
     MAP_CONTINENTS,color=culz.wheat, /horizon, /hires

     ;------Display the ASC and temperature images:
     if addasc then begin
       screen_img = tvrd(/true)
       red        = reform(screen_img(0, *, *))
       green      = reform(screen_img(1, *, *))
       blue       = reform(screen_img(2, *, *))
       green(mapvec) = image(ascvec)*asc_gain < 255

       red(mapvec)   = 0
       blue(mapvec)  = 0

       screen_img(0, *, *) = red
       screen_img(1, *, *) = green
       screen_img(2, *, *) = blue
       tv, screen_img, /true
     endif

     for ss=0,n_elements(site_labels)-1 do begin
         if site_labels(ss).doprint then begin
            plots,  site_labels(ss).lon, site_labels(ss).lat, $
                    color=site_labels(ss).color, symsize=site_labels(ss).symsize, $
                    thick=site_labels(ss).thick, psym=site_labels(ss).psym
            dy = 0.065*(max(winds(0).zone_latitudes) - min(winds(0).zone_latitudes))
            xyouts, site_labels(ss).lon, site_labels(ss).lat - dy, site_labels(ss).name, $
                    color=site_labels(ss).color, charsize=site_labels(ss).charsize, $
                    charthick=site_labels(ss).charthick, align=0.5
          endif
     endfor


;     sdi3k_map_windvex, obs_coords(these).lons, obs_coords(these).lats, sin(!dtor*azi(these))*sin(!dtor*zen(these))*outarr.loswind, $
;                                                                        cos(!dtor*azi(these))*sin(!dtor*zen(these))*outarr.loswind, culz.blue, culz.orange, thick=1
     sdi3k_map_windvex, out_coords.lons, out_coords.lats, fitwinds.zonal, fitwinds.meridional, culz.green, culz.white, thick=2
     wait, 0.01
;-------Add the annotation in each of the four corners:
       xyouts, /normal, .04*xshrink, .94,  dt_tm_mk(js2jd(0d)+1, tnow, format='Y$-0n$-0d$'), color=culz.white, charsize=3, charthick=3
       xyouts, /normal, .04*xshrink, .89,  dt_tm_mk(js2jd(0d)+1, tnow, format='h$:m$:s$'),   color=culz.white, charsize=3, charthick=3, align=0
;       xyouts, /normal, .04*xshrink, .03, 'Poker Flat Alaska',                               color=culz.white, charsize=2.2, charthick=2
;       xyouts, /normal, .96*xshrink, .03,  attribution,                                      color=culz.white, charsize=1.5, charthick=2,  align=1
if doing_red then $
       xyouts, /normal, .6, .035, 'Poker Flat/Conde',                               color=culz.white, charsize=3, charthick=3, align=0.5 else $
       xyouts, /normal, .96*xshrink, .035, 'Poker Flat/Conde',                               color=culz.white, charsize=3, charthick=3, alig=1
;       xyouts, /normal, .96*xshrink, .03,  attribution,                                      color=culz.white, charsize=1.5, charthick=2,  align=1

       lat = latmin + 0.97*(latmax - latmin)
       if doing_red then lat = latmin + 0.999*(latmax - latmin)
       lon = 0.515*(lonmin+lonmax)
       lon = 0.485*(lonmin+lonmax)
                hh  = convert_coord(lon, lat, /data, /to_device)
                dlat = convert_coord(lon, lat + 1., /data, /to_device)
                dlon = convert_coord(lon + 1., lat, /data, /to_device)
                lathat = [dlat(0) - hh(0), dlat(1) - hh(1)]
                lonhat = [dlon(0) - hh(0), dlon(1) - hh(1)]/cos(!dtor*lat)
                plots, lon, lat, psym=1, symsize=0.25, thick=3, color=culz.green
                xend = hh(0) + 0.008*100*norm(lonhat)
                yend = hh(1)
                arrow, hh(0), hh(1), xend, yend, color=culz.white, thick=3, hthick=3, hsize=10
                xyouts, /data, lon, lat, '100 m/s ',   color=culz.white,  charsize=2.5, align=1, charthick=3
    clr = '_Green_'
    if doing_red then clr = '_Red_'
    prefix = 'sdi_bistatic'
    frame_name = prefix + clr+ 'DOY' + dt_tm_mk(js2jd(0d)+1, tnow, format='doy$_Y$_(n$_0d$)_at_h$-m$-s$UT') + '.png'
    gif_this, file=still_folder + '\' + frame_name, /png
    tnow = tnow + tstep
    wait, 0.001
endwhile

end
