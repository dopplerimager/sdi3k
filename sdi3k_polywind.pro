;=================================================================================================
;
;+
; NAME:
;       Polywind
; PURPOSE:
;       This routine fits a horizontal vector wind field
;       to multi-station FPS wind data.  The fit is a 2 or 3-D polynomial
;       field, with polynomial terms up to the order specified by the
;       "order" parameter.
; CALLING:
;       Polywind, los_obs, sigobs,     obs_x,    obs_y,   view_azi,  view_zen,  order, $
;                 zonal,   meridional, vertical, sigzon,  sigmer,    sigver,    quality, $
;                 radians=radz, horizontal_only=hozo
; INPUTS:
;       los_obs:     A 1-D array of line-of-sight horizontal wind observations, positive away.
;       sigobs:      A 1-D array containing the measurement uncertainties in los_obs.
;       obs_x:       A 1-D array containing the x (longitudinal) displacement (in km) of the observed
;                    locations, relative to some convenient origin point.
;       obs_y:       A 1-D array containing the y (latitudinal)  displacement (in km) of the observed
;                    locations, relative to the same origin point.
;       view_azi:    A 1-D array containing the observing azimuth toward each observed
;                    location (eg 0 for a NORTH observation, 90 for an EAST observation, etc).
;                    The units (degrees or radians) depend on the "radians" keyword.
;       view_zen:    A 1-D array containing the observing zenith angle toward each observed
;                    location. The units (degrees or radians) depend on the "radians" keyword.
;       order:       The order of the polynomial fit (1 for a linear fit).
; OUTPUTS:
;       zonal:       A 1-D array of coefficients of the zonal wind.  Array elements are
;                    0 - the uniform wind
;                    1 - the first order variation in the zonal direction
;                    2 - the first order variation in the meridional direction
;                    3 - the second order variation in the zonal direction
;                    4 - the second order variation in the meridional direction
;                    .... etc...
;       meridional:  A 1-D array of coefficients of the meridional wind, with the same structure as
;                    the zonal coefficient array.
;       vertical:    A 1-D array of coefficients of the vertical wind, with the same structure as
;                    the zonal coefficient array.
;       sigzon:      A 1-D array of uncertainties in zonal coefficients.
;       sigmer:      A 1-D array of uncertainties in meridional coefficients.
;       sigver:      A 1-D array of uncertainties in vertical coefficients.
;       quality:     A structure describing the quality of the fit.  Fields in the structure are:
;                    chisq:    The reduced chi-squared parameter of the fit.
;                    df:       The number of statistical degrees of freedom of the fit.
;                    los_fit:  A 1-D array, of the same size as los_obs, containing the fitted
;                              line-of-sight wind estimates.
;                    singular: The number of singular values encountered by SVDFIT (should
;                              be zero for a good fit).
; KEYWORDS:
;       radians:     If set, this indicates that view_azi and view_zen are already in radians,
;                    and do not need converting.
;       horizontal_only: If this keyword is set, the fit will only consider horizontal wind components.
;                    Otherwise, coefficients for the vertical wind are calculated, along with their
;                    corresponding uncertainties.
; PROCEDURE:
;       Chi-squared minimization using singular-value decomposition.
; HISTORY:
;       Written by Mark Conde, Bromley, August 1999.
;-
;
;==============================================================================

function wnd_func, obsidx, nterms
    common wndcom, basis, used
    return, reform(basis(obsidx,used))
end

pro wnd_basis, obs_x, obs_y, obs_azi, obs_zen, nobs, order, vertz
    common wndcom, basis, used

    for i = 0L,nobs-1 do begin
;       Do the zeroth order (uniform wind) separately:
        basis(i, 0) = sin(obs_azi(i))*sin(obs_zen(i))
        basis(i, 1) = cos(obs_azi(i))*sin(obs_zen(i))
        if vertz then basis(i, 2) = cos(obs_zen(i))

        k = vertz+2
;       Now do the spatially varying wind terms:
        for j = 1,order do begin
            for m=j,0,-1 do begin
;               Calculate the basis functions for each power of x and y:
                basis(i, k) = sin(obs_azi(i))*sin(obs_zen(i))*((obs_x(i))^m)*((obs_y(i))^(j-m))
                k = k+1
                basis(i, k) = cos(obs_azi(i))*sin(obs_zen(i))*((obs_x(i))^m)*((obs_y(i))^(j-m))
;                if k eq 3 then basis(i, k) = basis(i, k)*sin(5.*obs_azi(i)) ;#####################
                k = k+1
                if vertz then basis(i, k) = cos(obs_zen(i))*((obs_x(i))^m)*((obs_y(i))^(j-m))
                if vertz then k = k+1
            endfor
        endfor
    endfor
end

pro sdi3k_get_poly_wind, fit_x, fit_y, order, fitpars, result
         if n_elements(fitpars) eq 2*total(indgen(order+1)+1) then vertz=0 else vertz=1

         result = {zonal: findgen(n_elements(fit_x)), $
                   meridional: findgen(n_elements(fit_x)), $
                   vertical: findgen(n_elements(fit_x))}
         nterms  = total(indgen(order+1)+1)
         for n=0,n_elements(fit_x)-1 do begin
             result.zonal(n) = fitpars(0)
             result.meridional(n) = fitpars(1)
             if vertz then result.vertical(n) = fitpars(2)
             k = vertz+2
             for j = 1,order do begin
                 sidx = fix((vertz+2)*total(indgen(j)+1))
                 for m=j,0,-1 do begin
;                     print, j, sidx, m, k
                     result.zonal(n)      = result.zonal(n)      + fitpars(k)*((fit_x(n))^m)*((fit_y(n))^(j-m))
                     k = k + 1
;                     print, j, sidx, m, k
                     result.meridional(n) = result.meridional(n) + fitpars(k)*((fit_x(n))^m)*((fit_y(n))^(j-m))
                     k = k + 1
                     if vertz then result.vertical(n)   = result.vertical(n)   + fitpars(k)*((fit_x(n))^m)*((fit_y(n))^(j-m))
                     if vertz then k = k + 1
                 endfor
             endfor
         endfor
end

pro sdi3k_polywind, los_obs, sigobs,     obs_x,    obs_y,   view_azi, view_zen, order, fitpars, $
              zonal,   meridional, vertical, sigzon,  sigmer,   sigver,   quality, $
              radians=radz, horizontal_only=hozo, used=use_these
    common wndcom, basis,used

    obs_azi = view_azi
    obs_zen = view_zen
    if not(keyword_set(radz)) then obs_azi = !dtor*obs_azi
    if not(keyword_set(radz)) then obs_zen = !dtor*obs_zen
    if keyword_set(hozo) then vertz=0 else vertz=1

    nobs    = n_elements(los_obs)
    nterms  = (vertz + 2)*total(indgen(order+1)+1)
    df      = nobs - nterms
    basis   = fltarr(nobs, nterms)
    if not(keyword_set(use_these)) then used=lindgen(nterms) else used=where(use_these ne 0)

;   Generate the basis functions:
    wnd_basis, obs_x, obs_y, obs_azi, obs_zen, nobs, order, vertz

;   Now fit the parameters:
    vpars = fltarr(n_elements(used))
    fitpars = fltarr(nterms)
    varpars = fltarr(nterms)
    numsing = 0
    los_fit = fltarr(nobs)
    fpars   = svdfit(findgen(nobs), los_obs, n_elements(used), $
                     funct="wnd_func", $
                     weight=fltarr(nobs) + 1./sigobs, $
                     variance=varpars, $
                     singular=numsing, $
                     yfit=los_fit)

     fitpars(used) = fpars
     varpars(used) = vpars

     zonal      = fltarr(nterms/(vertz+2))
     meridional = fltarr(nterms/(vertz+2))
     vertical   = fltarr(nterms/(vertz+2))
     sigzon     = fltarr(nterms/(vertz+2))
     sigmer     = fltarr(nterms/(vertz+2))
     sigver     = fltarr(nterms/(vertz+2))

;    Populate the parameter estimates and uncertainties for uniform wind terms:
     diridx = (vertz+2)*indgen(total(indgen(order+1)+1))
     zonal  = fitpars(diridx)
     meridional = fitpars(diridx + 1)
     if vertz then vertical(0) = fitpars(diridx + 2)

     sigzon    = sqrt(varpars(diridx))
     sigmer(0)    = sqrt(varpars(diridx + 1))
     if vertz then sigver(0)   = sqrt(varpars(diridx + 2))

;    Finally, build the structure to return fit quality info:
     resid     = los_obs - los_fit
     chisq     = total((resid/sigobs)^2)/df
     quality   = {chisq: chisq, df: df, los_fit: los_fit, singular: numsing}
end
