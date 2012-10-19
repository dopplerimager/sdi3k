L1fitmode = 'new'
;L2fitmode = 'new'
L2fitmode = 'all'
plotmode = 'new'

lookback_toolik = 21L
lookback_mawson = 31L
lookback_poker  = 21L
lookback_haarp  = 21L

do_toolik = 1
do_poker  = 1
do_haarp  = 1
do_mawson = 1

file_filter = ['*.pf', '*.nc', '*.sky', '*.las']

setenv, 'SDI_GREEN_ZERO_VELOCITY_FILE=None'
setenv, 'SDI_RED_ZERO_VELOCITY_FILE=None'

if do_toolik then sdi3k_batch_autoprox, path='d:\users\sdi3000\data\Toolik_Lake', $
                  filter=file_filter, $
                  calfit=L1fitmode, skyfit=L1fitmode, windfit=L2fitmode, plotting=plotmode, $
                  lookback_seconds=lookback_toolik*86400L, drift='data'
if do_poker  then sdi3k_batch_autoprox, path='d:\users\sdi3000\data\Poker', $
                  filter=file_filter, $
                  calfit=L1fitmode, skyfit=L1fitmode, windfit=L2fitmode, plotting=plotmode, $
                  lookback_seconds=lookback_poker*86400L, drift='data'
if do_haarp  then sdi3k_batch_autoprox, path='d:\users\sdi3000\data\HAARP', $
                  filter=file_filter, $
                  calfit=L1fitmode, skyfit=L1fitmode, windfit=L2fitmode, plotting=plotmode, $
                  lookback_seconds=lookback_haarp*86400L, drift='data'
if do_mawson then sdi3k_batch_autoprox, path='d:\users\sdi3000\data\Mawson', $
                  filter=file_filter, $
                  calfit=L1fitmode, skyfit=L1fitmode, windfit=L2fitmode, plotting=plotmode, $
                  lookback_seconds=lookback_mawson*86400L, drift='data'
;exit, /no_confirm

end