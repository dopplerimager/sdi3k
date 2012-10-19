pro aero_get_par, aerun, parname, etime, pardat
    pardat  = 'unknown'
    oneline = 'dummy'
    on_ioerror, ioerr

;---Scan file for the chart block for this parameter:
    repeat readf, aerun, oneline until strpos(oneline, parname + 'Chart()') ge 0 or eof(aerun)
ioerr:
    if eof(aerun) then return

;---Get the end time from the plot label strings:
    repeat readf, aerun, oneline until strpos(oneline, 'c.setHorizontalLabels') ge 0 or eof(aerun)
    if eof(aerun) then return
    oneline = strmid(oneline, strpos(oneline, '['), strpos(oneline, ']') - strpos(oneline, '[') + 1)
    sts     = execute('tarr = ' + oneline)
    etime   = tarr(n_elements(tarr) - 1)

;---Get the data array:
    repeat readf, aerun, oneline until strpos(oneline, 'c.add') ge 0 or eof(aerun)
    if eof(aerun) then return
    oneline = strmid(oneline, strpos(oneline, '[') + 1, strpos(oneline, ']') - strpos(oneline, '[') - 1)
    valz = strsplit(oneline, ',', /extract)
    pardat = float(valz)

end

pro read_aerospace_data, fname, aerodat

    openr, aerun, fname, /get_lun

;---Get start date and time from html title block:
    oneline = 'dummy'
    repeat readf, aerun, oneline until strpos(oneline, '<TITLE>') ge 0 or eof(aerun)
    if eof(aerun) then return
    fields = strsplit(oneline, '/', /extract)
    year   = 2000 + fix(strmid(fields(0), strlen(fields(0))-2, 2))
    month  = fix(fields(1))
    day    = fix(strmid(fields(2), 0, 2))
    fields = strsplit(fields(2), ':', /extract)
    shr    = fix(strmid(fields(0), strlen(fields(0))-2, 2))
    smn    = fix(fields(1))
    ssc    = fix(strmid(fields(2), 0, 2))
    startjs = ymds2js(year, month, day, 3600.*shr+60.*smn+ssc)

    aero_get_par, aerun, 'Sunsensors', etime, sunsensors
    aero_get_par, aerun, 'FilterTemp', etime, filtertemp
    aero_get_par, aerun, 'TeleTemp',   etime, teletemp
    aero_get_par, aerun, 'CoolerTemp', etime, coolertemp
    aero_get_par, aerun, 'i4278',      etime, i4278
    aero_get_par, aerun, 'i6300',      etime, i6300
    aero_get_par, aerun, 'i8714',      etime, i8714
    aero_get_par, aerun, 'i8446',      etime, i8446
    aero_get_par, aerun, 'e',          etime, nrg
    aero_get_par, aerun, 'oN2',        etime, O2N2
    aero_get_par, aerun, 'Flux',       etime, flux

    fields  = strsplit(etime, ':', /extract)
    ehr     = fix(fields(0))
    emn     = fix(fields(1))
    endjs   = ymds2js(year, month, day, 3600.*ehr+60.*emn)


    aerec = {js: 0D, sunsensors: 0., filtertemp: 0., teletemp: 0., coolertemp: 0., i4278: 0., i6300: 0., i8714: 0., i8446: 0., nrg: 0., O2N2: 0., flux: 0.}
    aerodat = replicate(aerec, n_elements(i4278))
    aerodat.js = startjs + (endjs - startjs)*findgen(n_elements(i4278))/n_elements(i4278)
    aerodat.sunsensors = sunsensors
    aerodat.filtertemp = filtertemp
    aerodat.teletemp   = teletemp
    aerodat.coolertemp = coolertemp
    aerodat.i4278      = i4278
    aerodat.i6300      = i6300
    aerodat.i8714      = i8714
    aerodat.i8446      = i8446
    aerodat.nrg        = nrg
    aerodat.O2N2       = O2N2
    aerodat.flux       = flux

    close, aerun
    free_lun, aerun
end