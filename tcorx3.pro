function tinput,fin,pmin,pmax,dd
	forward_function finput,pdata,sconv,percent
	common qdata
	band=float(band)
	ptmp=pdata(band,pmin,pmax)
	;dcl=percent(band)
	bandx=sconv(band,ptmp[0],ptmp[1],dd)

function finput,fin
	common qdata
	n=long(imax)*jmax
	tmin=long(n*pmin) & tmax=long(n*pmax)
	tmx=sort(tm)
	;print,[tm[tmx[tmin]],tm[tmx[tmax]]]

function pdata,tm,pmin,pmax	
	;dmin = min(tm) + 5
	;dmax = max(tm) - 10 < 210
	;dcl=percent(tm,pmin,pmax,dmin,dmax)
	dcl=percent(tm,pmin,pmax)
	;n=long(imax)*jmax
	;tmin=long(n*pmin) & tmax=long(n*pmax)
	;tmx=sort(tm)
	;print,[tm[tmx[tmin]],tm[tmx[tmax]]]

function sconv,tm,min,max,dmax	

function sdisp,tm,pmin,pmax	

;function percent,img,dmin,dmax,pmin,pmax
;	temp=where((img gt dmin) and (img lt dmax))
function percent,img,pmin,pmax
	common pdata
	;temp=where(vi eq 1)
	temp=where(img le !values.f_infinity,count,complement=temp2)
	tnum=n_elements(temp)
	print,tnum
	timg=img[temp]
	tsort=sort(timg)
	low=long(tnum*pmin) & high=long(tnum*pmax)
	min=timg[tsort[low]] & max=timg[tsort[high]]
	return,[min,max]
end

function quick_image,img,dmin=dmin,dmax=dmax,pmin=pmin,pmax=pmax
	common pdata
	;if n_elements(dmin) eq 0 then dmin = min(img) + 2
	;if n_elements(dmax) eq 0 then dmax = max(img) - 2
	;if n_elements(pmin) eq 0 then pmin=0.03
	;if n_elements(pmax) eq 0 then pmax=0.97
	if n_elements(pmin) eq 0 then pmin=0.001
	if n_elements(pmax) eq 0 then pmax=0.99
	;dscl=percent(img,dmin,dmax,pmin,pmax)
	dscl=percent(img,pmin,pmax)
	return,bytscl(img,min=dscl[0],max=dscl[1])
end


	common qdata
		;print,i,systime(1)-t0
		tmp=where(cls1 eq i)
		cl2i=cls2[tmp] & cl3i=cls3[tmp] & tmxi=tmx[tmp]
	for j=0,d2-1 do begin
		tmp2=where(cl2i eq j,cnt) 
		if(cnt ne 0) then begin
			cl3j=cl3i[tmp2] & tmxj=tmxi[tmp2]
			for k=0,d3-1 do begin
					nper=long(per*cnt)
					tmxx=tmxj[temp]
					stmx=sort(tmxx)
					;ktemp[tmp[tmp2[temp]]]=median(tmxj[temp])
				endif
		endif
	endfor

function aest,tmx,per
	common qdata
		;print,i,systime(1)-t0
		tmp=where(cls1 eq i)
		cl2i=cls2[tmp] & cl3i=cls3[tmp] & tmxi=tmx[tmp] & vii=vi[tmp]
	for j=0,d2-1 do begin
		tmp2=where(cl2i eq j,cnt) 
		if(cnt ne 0) then begin
			cl3j=cl3i[tmp2] & tmxj=tmxi[tmp2] & viij=vii[tmp2]
			for k=0,d3-1 do begin
					nper=long(per*cnt)
					tmxx=tmxj[temp]
					stmx=sort(tmxx)
					;ktemp[tmp[tmp2[temp]]]=median(tmxj[temp])
				endif
		endif
	endfor

function aestx,tmx,per
	common qdata
	d0=max(cls)
		;print,i,systime(1)-t0
		temp=where((cls eq i) and (vi eq 1),cnt)
		if(cnt ne 0) then begin
			nper=long(per*cnt)
			tmxi=tmx[temp]
			stmx=sort(tmxi)
			;ktemp[tmp]=median(tmx[temp])
		endif
	return,ktemp


function best,tm,tma,prad
	temp=where(vi eq 0)
	tmx[temp]=prad
	return,tmx
end

function t_correction,tm,smax,prad,per,dr,flag
	common pdata
	common qdata
	print,"Step 1: Subtraction of assumed path radiance"
	sm=tcrct(tm,dem,prad,hmax)
	print,"Step 2: Estimation of minimum DN for each class"
	wm=aest(sm,per) 
	print,"Step 3: Estimation of path radiance for each pixel"
	tmb=best(tm,wm)                  
	if(flag eq 1) then print,mean(tmb),stdev(tmb)
	tmb=tmb > 0 < smax
	print,"Step 4: Spatial smoothing by median filter"
	print,"       This step takes time !!"
	bmed=median(tmb,dr) 
	if(flag eq 1) then print,mean(bmed),stdev(bmed)
	return,tcrct(tm,dem,bmed,hmax)
end

function t_correction2,tm,prad,per,dr,flag
	common pdata
	common qdata
	print,"Step 1: Subtraction of assumed path radiance"
	sm=tcrct(tm,dem,prad,hmax) > 0.0 
	print,"Step 2: Estimation of minimum DN for each class"
	wm=aest(sm,per) 
	print,"Step 3: Estimation of path radiance for each pixel"
	tmb=best(tm,wm,prad)                  
	if(flag eq 1) then print,mean(tmb),stdev(tmb)
	print,"Step 4: Spatial smoothing by median filter"
	print,"       This step takes time !!"
	bmed=median(tmb,dr) 
	if(flag eq 1) then print,mean(bmed),stdev(bmed)
	;xm=tcrct(tm,dem,bmed,hmax) > 0.0
	;temp=where(vi eq 0)
	;return,xm
	return,tcrct(tm,dem,bmed,hmax)
end

	xmax=n_elements(dem[*,0]) 
	ymax=n_elements(dem[0,*])
	ax=(shift(dem,-1,0)-shift(dem,1,0))/(2*dx)
	ax[0,*]=ax[1,*] & ax[xmax-1,*]=ax[xmax-2,*]
	ay=(shift(dem,0,1)-shift(dem,0,-1))/(2*dy)
	ay[*,0]=ay[*,1] & ay[*,ymax-1]=ay[*,ymax-2]
	tt=!DPI*elv/180.0
	ff=!DPI*azm/180.0
	xa=cos(tt)*sin(ff)
	ya=cos(tt)*cos(ff)
	za=sin(tt)
	return,(-xa*ax-ya*ay+za)/sqrt(1+ax*ax+ay*ay)
end
	forward_function estimate2,iestimate
	common pdata ; 2015/6/8
	common ydata
	rcoef=estimate2(demx,inc,tmx,penvx)
	ref=rcoef[0,*,*]+rcoef[1,*,*]*aerox+rcoef[2,*,*]*aerox^2
	refx=reform(ref,1200,1200)
	if flag then begin ; 2015/6/8
		vi=bytarr(1200,1200)
		temp=where(refx le !values.f_infinity)
		vi[temp]=1
	endif
	wm=aest(ref,per) 
	aeroy=iestimate(wm,rcoef)
	print,mean(refx,/nan),mean(aeroy,/nan)
	print,stddev(refx,/nan),stddev(aeroy,/nan)
end

pro iteratex,demx,inc,tmx,per,flag
	forward_function estimate2,iestimate
	common pdata ; 2015/6/8
	common ydata
	rcoef=estimate2(demx,inc,tmx,penvx)
	ref=rcoef[0,*,*]+rcoef[1,*,*]*aerox+rcoef[2,*,*]*aerox^2
	refx=reform(ref,1200,1200)
	print,mean(refx,/nan),stddev(refx,/nan)
	refx=refx > 0.0 < 1.0
	if flag then begin ; 2015/6/8
		vi=bytarr(1200,1200)
		temp=where(refx le !values.f_infinity)
		vi[temp]=1
	endif
	wm=aestx(ref,per) 
	aeroy=iestimate(wm,rcoef)
	print,mean(aeroy,/nan),stddev(aeroy,/nan)
end

pro mclass
	common pdata
	common qdata
	cls=lonarr(imax,jmax)
	for i=0,imax-1 do for j=0,jmax-1 do $
		cls[i,j]=cls1[i,j]+cls2[i,j]*d2+cls3[i,j]*d2*d3
end

function xmedian,aero,dr
    tempx=where(aero le !values.f_infinity,complement=temp)
    tmean=mean(aero,/nan)
    aero[temp]=tmean
    aero2=smooth(aero,dr)
    aero[temp]=aero2[temp]
    return,median(aero,dr)
end