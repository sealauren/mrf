pro plotmr
;plot mass vs radius for HIRES 22 paper.

readcol, 'h22_orbital_pars.txt', koi, period, radius, eradius, mass, emass, $
		rho, erho, K, ek, st_rho, est_rho, b, eb, a_rstar, epoch, $
		  delim = ',' ,$
		format = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f' ;
x=strpos(koi,'.')
nkoi = n_elements(koi)
suffix = fltarr(nkoi)

for i=0,nkoi-1 do suffix[i] = strmid(koi[i],x[i],3) 
Tran = where(suffix ne '.10',comp = Non_tran)

; index 'tran' is all of the transiting planets.
; The non-transiting planets will not read in properly since they don't
;	have all of the columns filled.

;radius = radius(tran)
;i = where(radius gt 0.95)
;tran = tran(i)

koi_t = koi[tran]  ; only transiting planets
koi_nt = koi[non_tran] ; NT planets, but they will not index correctly.
period = period(tran)
radius = radius(tran)
eradius = eradius(tran)
mass = mass(tran)
emass = emass(tran)
rho = rho(tran)
erho = erho(tran)
print,'R      M      Rho'
forprint,koi_t,radius,eradius, mass, emass, rho, erho
print,'-------------------------'

plot,radius,rho,ps=8,xtit='Planet Radius (Re)',ytit='Planet Density (g/cc)',yr=[-3,15]

;postscript
!p.charsize=1.8
!p.thick=4
!p.charthick=5
ps_open,'rvsrho'
plot,radius,rho,ps=8,xtit='!6 Planet Radius (Re)',ytit='Planet Density (g/cc)',yr=[-6,15],title='!6Planet Density  vs  Radius'

oploterr,radius,rho,erho

ps_close
spawn,'open rvsrho.ps'

end 
