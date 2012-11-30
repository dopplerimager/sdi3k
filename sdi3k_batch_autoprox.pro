; Example calls:  sdi3k_batch_autoprox, path='d:\users\sdi2000\data\2000_fall\', filter=['*.pf', '*.nc'], calfit='none', skyfit='new', windfit='all', /choose, lookback_seconds=30*86400L
;                 sdi3k_batch_autoprox, path='d:\users\sdi3000\data\spectra', filter=['*.pf', '*.nc', '*.sky', '*.las'], calfit='none', skyfit='none', windfit='all', plotting='all', lookback_seconds=18L*86400L, /choose
;                 sdi3k_batch_autoprox, path='D:\users\SDI3000\Data\HAARP', filter=['*.pf', '*843nm*.nc', '*.sky', '*.las'], calfit='none', skyfit='none', windfit='all', plotting='all', lookback_seconds=18L*86400L, /choose, /ask
;                 sdi3k_batch_autoprox, path='F:\SDIData\Poker\', filter=['pkr_2012*_63*.nc'], calfit='none', skyfit='none', windfit='none', plotting='none', /choose, lookback_seconds=30*86400L

pro sdi3k_get_flat_spec, local_path, color
      color = strupcase(color)
      case color of
           'GREEN': lambda = '5577'
          'SODIUM': lambda = '5890'
             'RED': lambda = '6300'
              'OH': lambda = '8430'
              else: lambda = '6328'
      endcase
      flats    = findfile(local_path + "Wind_flat_field*" + lambda + "*.sav")
      optns = strarr(n_elements(flats))

      menu_items = ['Auto', 'None']
      if flats(0) ne '' then begin
         j=0
         while j lt n_elements(flats) do begin
             wind_flat_field = 'NULL'
             restore, flats(j)
             if size(wind_flat_field, /tname) ne 'STRUCT' then begin
                veceldel, flats, j, /last
                veceldel, optns, j, /last
             endif else begin
                optns(j) = string(10*wind_flat_field.metadata(0).wavelength_nm, format='(i4.4)') + '   '
                optns(j) = optns(j) + dt_tm_mk(js2jd(0d)+1, wind_flat_field.js_average,  format='Y$-N$-0d$') + ' '
                optns(j) = optns(j) + dt_tm_mk(js2jd(0d)+1, wind_flat_field.js_average,  format='(DOY doy$)') + '    [Valid '
                optns(j) = optns(j) + dt_tm_mk(js2jd(0d)+1, wind_flat_field.js_valid(0), format='Y$:doy$') + ' to '
                optns(j) = optns(j) + dt_tm_mk(js2jd(0d)+1, wind_flat_field.js_valid(1), format='Y$:doy$') + ']     Rated '
                optns(j) = optns(j) + strmid(wind_flat_field.stars.name, 0, strpos(wind_flat_field.stars.name, ' ')) + ' by ' + wind_flat_field.operator
                j = j + 1
             endelse
             wait, 0.01
         endwhile
         if n_elements(optns) gt 0 then menu_items = [optns, 'Auto', 'None']
      endif
      mcchoice, 'Wind flat field file for ' + color + '?', menu_items, choice, help='This flat field correction will only be applied to ' + color + ' winds'
      if choice.index lt n_elements(optns) then setenv, 'SDI_' + color + "_ZERO_VELOCITY_FILE=" + flats(choice.index) else $
                                                setenv, 'SDI_' + color + "_ZERO_VELOCITY_FILE=" + choice.name
end

pro sdi3k_batch_autoprox, path=local_path, $
                          filter=filter, $
                          lookback_seconds=lookback_seconds, $
                          chooser=chooser, $
                          skyfit=skyfit, $
                          calfit=calfit, $
                          windfit=windfit, $
                          plotting=plotting, $
                          ascii_export=ascii_export, $
                          plotstage=plotstage, $
                          plot_folder=plot_folder, $
                          ask_flat=ask_flat, $
                          drift_mode=drift_mode, $
                          xy_only=xy_only

if not(keyword_set(local_path))       then local_path       = 'd:\users\sdi3000\data\spectra\'
if not(keyword_set(filter))           then filter           = ['*.pf', '*.nc']
if not(keyword_set(lookback_seconds)) then lookback_seconds = 3*86400L
if not(keyword_set(calfit))           then calfit           = 'all' ; options are: 'all', 'none', or 'new'
if not(keyword_set(skyfit))           then skyfit           = 'all'
if not(keyword_set(windfit))          then windfit          = 'all'
if not(keyword_set(plotting))         then plotting         = 'all'
if not(keyword_set(ascii_export))     then ascii_export     = 'all'
if not(keyword_set(drift_mode))       then drift_mode       = 'data'

calfit   = strupcase(calfit)
skyfit   = strupcase(skyfit)
windfit  = strupcase(windfit)

local_path = file_expand_path(local_path) + '\'

if keyword_set(ask_flat) then begin
   sdi3k_get_flat_spec, local_path, 'GREEN'
   sdi3k_get_flat_spec, local_path, 'RED'
   sdi3k_get_flat_spec, local_path, 'OH'
endif

sdi3k_batch_ncquery, file_desc, path=local_path, filter=filter, /verbose
file_desc = file_desc(where(file_desc.sec_age le lookback_seconds))

skylis = file_desc(where(strupcase(file_desc.metadata.viewtype) eq 'SKY'))
calz   = where(strupcase(file_desc.metadata.viewtype) eq 'CAL', nncal)
if nncal gt 0 then begin
   callis = file_desc(where(strupcase(file_desc.metadata.viewtype) eq 'CAL'))
   skylis = skylis(sort(skylis.name))
   callis = callis(sort(callis.name))
endif

if keyword_set(chooser) then begin
   mcchoice, 'First file to process?', skylis.preferred_name, choice
   lodx = choice.index
   mcchoice, 'Last file to process?',  skylis.preferred_name, choice
   hidx = choice.index
   skylis = skylis(lodx<hidx:hidx>lodx)
endif

for j=0,n_elements(skylis)-1 do begin
       print, 'Processing: ', skylis(j).preferred_name
;       stop
       if nncal gt 0 or getenv('user_specified_insprof') ne '' then begin
          insinf   = mc_fileparse(skylis(j).metadata.path + strmid(skylis(j).insfile, 0, 4) + strmid(skylis(j).insfile, 9, 999))
          insname  = insinf.name_only
          insinf   = mc_fileparse(skylis(j).metadata.path + strmid(callis.insfile,    0, 4) + strmid(callis.insfile,    9, 999))
          insz     = insinf.name_only
          this_ins = where(insz eq insname, nn)
          if nn gt 0 then begin
             use = intarr(nn)
             for nni=0,nn-1 do begin
                 sdi3k_read_netcdf_data, callis(this_ins(nni)).name, metadata=mmi, /close
                 use(nni) = mmi.nzones eq skylis(j).metadata.nzones
             endfor
             best = where(use eq 1)
             best = best(0)
             this_ins = callis(this_ins(best)).name
          endif
          if getenv('user_specified_insprof') ne '' then this_ins = getenv('user_specified_insprof')
          if nn gt 0 or getenv('user_specified_insprof') ne '' then begin
             if strupcase(drift_mode) ne 'DATA' then drift_mode=this_ins
             sdi3k_read_netcdf_data, this_ins, metadata=mm, /close
             doit = size(mm, /tname) eq 'STRUCT'
             if doit then doit = doit and mm.maxrec gt 0
             if doit and skyfit ne 'NONE' then sdi3k_batch_spekfitz, skylis(j).name, this_ins, skip_existing=(skyfit eq 'NEW'), skip_insfit=(calfit eq 'NONE')
          endif
       endif
       sdi3k_read_netcdf_data, skylis(j).name, metadata=mm, /close
       if windfit ne 'NONE' and mm.spekfit_status eq 'Spectra Fitted' then begin

       if (windfit eq 'ALL') or (skylis(j).metadata.windfit_status ne 'Winds Fitted') then sdi3k_batch_windfitz, skylis(j).name, drift_mode=drift_mode
       endif
       if strupcase(strcompress(plotting, /remove)) ne 'NONE' then begin
          sdi3k_batch_plotz, skylis(j).name, skip_existing=(strupcase(strcompress(plotting, /remove)) eq 'NEW'), stage=plotstage, plot_folder=plot_folder, drift_mode=drift_mode, xy_only=xy_only, root_dir='C:\Users\sdi\SDIPlots\', /msis2000, /hwm07
       endif
       if strupcase(strcompress(ascii_export, /remove)) ne 'NONE' then begin
          plot_dir = 'C:\Users\sdi\SDIPlots\'
          lamstring = strcompress(string(fix(10*mm.wavelength_nm)), /remove_all)
          year      = strcompress(string(fix(mm.year)),             /remove_all)
          scode     = strcompress(mm.site_code, /remove_all)
          if strupcase(scode) eq 'PF' then scode = 'PKR'
          md_err = 0
          catch, md_err
          if md_err ne 0 then goto, keep_going
          folder = plot_dir + year + '_' + scode + '_' + lamstring + '\' + 'ASCII_Data' + '\'
          if !version.release ne '5.2' then file_mkdir, folder else spawn, 'mkdir ' + folder
keep_going:
          catch, /cancel
          stp = {export_allsky: 1, $
                 export_skymaps: 1, $
                 export_spectra: 0, $
                 export_wind_gradients: 0, $
                 apply_smoothing: 1, $
                 time_smoothing: 1.1, $
                 space_smoothing: 0.09}

          sdi3k_ascii_export, setup = stp, files = skylis(j).name, outpath = folder, $
                              skip_existing=(strupcase(strcompress(ascii_export, /remove)) eq 'NEW')
       endif

;       if strupcase(strcompress(ascii_export, /remove)) ne 'NONE' then begin
;          plot_dir = 'c:\inetpub\wwwroot\conde\sdiplots\'
;          lamstring = strcompress(string(fix(10*mm.wavelength_nm)), /remove_all)
;          year      = strcompress(string(fix(mm.year)),             /remove_all)
;          scode     = strcompress(mm.site_code, /remove_all)
;          if strupcase(scode) eq 'PF' then scode = 'PKR'
;          md_err = 0
;          catch, md_err
;          if md_err ne 0 then goto, keep_going
;          folder = plot_dir + year + '_' + scode + '_' + lamstring + '\' + 'ASCII_Data' + '\'
;          if !version.release ne '5.2' then file_mkdir, folder else spawn, 'mkdir ' + folder
;keep_going:
;
;          stp = {export_allsky: 1, $
;                 export_skymaps: 1, $
;                 export_spectra: 0, $
;                 apply_smoothing: 1, $
;                 time_smoothing: 1.1, $
;                 space_smoothing: 0.09}
;          sdi3k_ascii_export, setup = stp, files = skylis(j).name, outpath = folder, $
;                              skip_existing=(strupcase(strcompress(ascii_export, /remove)) eq 'NEW')
;       endif
endfor
end
