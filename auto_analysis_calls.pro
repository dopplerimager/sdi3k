sdi3k_batch_autoprox, path='D:\users\SDI3000\Data\Poker\', filter=['*.pf', '*.nc', '*.sky', '*.las'], calfit='all', skyfit='all', windfit='all', plotting='all', lookback_seconds=5L*365L*86400L, /choose, /ask, drift='data'


sdi3k_batch_autoprox, path='D:\users\SDI3000\Data\HAARP\', filter=['*.pf', '*.nc', '*.sky', '*.las'], calfit='all', skyfit='all', windfit='all', plotting='all', lookback_seconds=5L*365L*86400L, /choose, /ask, drift='laser'


sdi3k_batch_autoprox, path='D:\users\SDI3000\Data\Poker\', filter=['*.pf', '*red*.nc','*green*.nc', '*.sky', '*.las'], calfit='none', skyfit='none', windfit='all', plotting='all', $
lookback_seconds=1L*365L*86400L, /choose, /ask, drift='data', plotstage='STAGE_WINDMAPS'


end