pro aread,fin
	common xdata
	line=''
	openr,1,fin
	;readf,1,line
	;print,line
	;readf,1,line
	nband=0S & v1=fltarr(2) & v2=fltarr(2) & penv=fltarr(3)
	pquick=fltarr(2)
	dr1=0s & dr2=0s
	while ~eof(1) do begin
		readf,1,line
		case line of
		"sun position:": begin
			readf,1,line
			reads,strmid(line,10),sun_el			readf,1,line
			reads,strmid(line,10),sun_az
			end		"number of bands:": begin
			readf,1,line
			reads,strmid(line,10),nband
			offset = fltarr(nband)
			gain = fltarr(nband)			readf,1,line
			readf,1,line
			reads,strmid(line,10),offset
			readf,1,line
			reads,strmid(line,10),gain
			end		"clear area:":begin			readf,1,line
			reads,strmid(line,12),per			end		"median filter:":begin			readf,1,line
			reads,strmid(line,10),dr1,dr2
			end
		"initial ref:":begin
			readf,1,line
			reads,strmid(line,10),penv
			end
		"initial depth:":begin
			readf,1,line
			reads,strmid(line,10),aero0x
			end
		"vegetation index:":begin
			readf,1,line
			reads,strmid(line,10),v1
			readf,1,line
			reads,strmid(line,10),v2
			end
		"quick look:":begin			readf,1,line
			reads,strmid(line,12),pquick
			end		else: print,line
		endcase
	endwhile	close,1
end

function seppen,fname
	temp1=fltarr(4)
	temp2=fltarr(6)
	temp3=fltarr(7)
	command = 'info.awk '+fname
	spawn,command,result
	reads,result[0],temp1
	reads,result[1],temp2
	reads,result[2],temp3
	return,[temp1,temp2,temp3]
end

function read_height,fname
	data=fltarr(17,20)
	text=''
	temp1=fltarr(4)
	temp2=fltarr(6)
	temp3=fltarr(7)

	openr,1,fname
	for i=0,19 do begin 
		readf,1,text
		readf,1,temp1
		data[0:3,i]=temp1
		readf,1,temp2
		data[4:9,i]=temp2
		readf,1,temp3
		data[10:16,i]=temp3
	endfor
	close,1
	return,data
end

function set_coeff,xx,data
	d_coeff=fltarr(3,17)
	x2=xx^2
	xxx=transpose([[xx],[x2]])
	result=regress(xxx,reform(data[4,*],20),const=c0)
	d_coeff[*,4]=[c0,reform(result)] ; Direct light
	result=regress(xxx,reform(data[5,*],20),const=c0)
	d_coeff[*,5]=[c0,reform(result)] ; Sky light
	result=regress(xxx,reform(data[6,*],20),const=c0)
	d_coeff[*,6]=[c0,reform(result)] ; Environment light
	result=regress(xxx,reform(data[7,*],20),const=c0)
	d_coeff[*,7]=[c0,reform(result)] ; Path radiance
	result=regress(xxx,reform(data[8,*],20),const=c0)
	d_coeff[*,8]=[c0,reform(result)] ; Background radiance
	result=regress(xxx,reform(data[10,*],20),const=c0)
	d_coeff[*,10]=[data[10,0],0,0] ; Gas transmission down
	result=regress(xxx,reform(data[11,*],20),const=c0)
	d_coeff[*,11]=[data[11,0],0,0] ; Gas transmission up
	result=regress(xxx,reform(data[15,*],20),const=c0)
	d_coeff[*,15]=[c0,reform(result)] ; Optical depth 
	result=regress(xxx,reform(data[16,*],20),const=c0)
	d_coeff[*,16]=[c0,reform(result)] ; Spherical albedo
	return,d_coeff
end

function set_coeff2,xx,yy,data
	d_coeff=fltarr(6,17)
	x2=xx^2 & y2=yy^2 & xy=xx*yy
	xxx=transpose([[xx],[yy],[x2],[xy],[y2]])
	result=regress(xxx,reform(data[7,*],10),const=c0)
	d_coeff[*,7]=[c0,reform(result)] ; Path radiance
	result=regress(xxx,reform(data[8,*],10),const=c0)
	d_coeff[*,8]=[c0,reform(result)] ; Background radiance
	result=regress(xxx,reform(data[15,*],10),const=c0)
	d_coeff[*,15]=[c0,reform(result)] ; Optical depth 
	result=regress(xxx,reform(data[16,*],10),const=c0)
	d_coeff[*,16]=[c0,reform(result)] ; Spherical albedo
	return,d_coeff
end


function p_height,data,datax,height,height0
	slope=(data[*,6]-datax)/height0
	return,datax+slope*height
end

pro h_plot,data,datax,n
	xx=findgen(20)*0.5
	plot,xx,data[n,*],yrange=[0.0,0.3]
	oplot,[0.0,3.0],data[n,[0,6]],color=255*256L
	oplot,[0.0,3.0],[datax[n],data[n,6]],color=255*256L
end

function reflectance,data,rad
  path=data[7]
  back=data[8]
  gtrans=data[10]
  edir=data[4]
  esky=data[5]
  eenv=data[6]
  odepth=data[15]
  return,!dpi*(rad-path-back)/gtrans/(edir+esky+eenv)*exp(odepth)
end

function second,x,y
  xx=transpose([[x],[x^2]])
  temp=regress(xx,y,const=c0)
  return,[c0,reform(temp,2)]
end

function e_depth,ref,coeff
  d0=coeff[0]-ref
  d1=coeff[1]
  d2=coeff[2]
  return,-d1/d2/2+sqrt(d1^2/d2^2/4-d0/d2)
end

function adjust,depth,height
	common rdata
	adata=[1.0,depth,depth^2]##transpose(d_coeff)
	;print,reflectance(d_data[*,4],100.0)
	;print,reflectance(adata,100.0)
	return,adata+(-adata+hdata)*height/4.0
end

function adjust2x,depth,height,cosb1,penv1
	common rdata
	hdata2=hdata
	adata=[1.0,depth,depth^2]##transpose(d_coeff)
	adata[4]=adata[4]*cosb1/cosb0 ; direct irradiance
	hdata2[4]=hdata[4]*cosb1/cosb0
	temp=(1-penv0*adata[16])*penv1/(1-penv1*adata[16])/penv0
	temp2=(1-penv0*hdata[16])*penv1/(1-penv1*hdata[16])/penv0
	adata[6]=adata[6]*temp ; env irradiance
	hdata2[6]=hdata[6]*temp2
	adata[8]=adata[8]*temp ; back radiance
	hdata2[8]=hdata[8]*temp2
	;print,reflectance(d_data[*,4],100.0)
	;print,reflectance(adata,100.0)
	return,adata+(-adata+hdata2)*height/4.0
end

function adjust2,depth,height,cosb1,penv1
	common rdata
	hdata2=hdata
	adata=[1.0,depth,depth^2]##transpose(d_coeff)
	adata[4]=adata[4]*cosb1/cosb0 ; direct irradiance
	hdata2[4]=hdata[4]*cosb1/cosb0
	temp=(1-penv0*adata[16])*penv1/(1-penv1*adata[16])/penv0
	temp2=(1-penv0*hdata[16])*penv1/(1-penv1*hdata[16])/penv0
	adata[6]=adata[6]*temp ; env irradiance
	hdata2[6]=hdata[6]*temp2
	adata[8]=adata[8]*temp ; back radiance
	hdata2[8]=hdata[8]*temp2
	return,adata+(-adata+hdata2)*height/4.0
end

function set_ref,height,rad,cosb1,penv1
	common rdata
	ref=fltarr(20)
	for i=0,19 do begin
		adata=adjust2(0.05*i,height,cosb1,penv1)
		ref[i]=reflectance(adata,rad)
	endfor
	;print,ref
	xx=findgen(20)*0.05
	x2=xx^2
	xxx=transpose([[xx],[x2]])
	result=regress(xxx,ref,const=c0)
	return,[c0,reform(result)]
end

function set_ref2,height,rad,cosb1,penv1
	common rdata
	x0=0.2 & x1=0.4 & x2=0.6
	adata=adjust2(x0,height,cosb1,penv1)
	y0=reflectance(adata,rad)
	adata=adjust2(x1,height,cosb1,penv1)
	y1=reflectance(adata,rad)
	adata=adjust2(x2,height,cosb1,penv1)
	y2=reflectance(adata,rad)
	;print,y0,y1,y2
	;c0=(x1*x2*y0-2*x2*x0*y1+x0*x1*y2)/0.08
	c0=3*y0-3*y1+y2
	;c1=(-(x1+x2)*y0+2*(x2+x0)*y1-(x0+x1)*y2)/0.08
	c1=(-y0+1.6*y1-0.6*y2)/0.08
	c2=(y0-2*y1+y2)/0.08
	return,[c0,c1,c2]
end

function set_ref2x,height,rad,cosb1,penv1
	common rdata
	x0=0.2 & x1=0.4
	adata=adjust2(x0,height,cosb1,penv1)
	y0=reflectance(adata,rad)
	adata=adjust2(x1,height,cosb1,penv1)
	y1=reflectance(adata,rad)
	c1=(y1-y0)/(x1-x0)
	c0=(y0*x1-y1*x0)/(x1-x0)
	return,[c0,c1]
end

function estimate,dem,inc,rad,penv1
	rcoef=fltarr(2,1200,1200)
	for j=0,1199 do begin
		;if (j mod 100) eq 0 then print,j
		for i=0,1199 do begin
			temp=set_ref2x(dem[i,j],rad[i,j],inc[i,j],penv1[i,j])
			rcoef[*,i,j]=temp
		endfor
	endfor
	return,rcoef
end

function estimate2,dem,inc,rad,penv1
	rcoef=fltarr(3,1200,1200)
	for j=0,1199 do begin
		;if (j mod 100) eq 0 then print,j
		for i=0,1199 do begin
			temp=set_ref2(dem[i,j],rad[i,j],inc[i,j],penv1[i,j])
			rcoef[*,i,j]=temp
		endfor
	endfor
	return,rcoef
end

function iestimate,wm,rcoef
	aero=fltarr(1200,1200)+!values.f_nan
	a0=reform((rcoef[0,*,*]-wm)/rcoef[2,*,*],1200,1200)
	a1=reform(rcoef[1,*,*]/rcoef[2,*,*],1200,1200)
	for i=0,1199 do begin
	for j=0,1199 do begin
		temp=a1[i,j]^2/4-a0[i,j]
		if temp lt 0 then continue
		x1=-a1[i,j]/2+sqrt(temp)
		x2=-a1[i,j]/2-sqrt(temp)
		if (x1 ge 0) and (x1 le 1.0) then begin
			aero[i,j]=x1
			if (x2 ge 0) and (x2 le 1.0) then aero[i,j]=!values.f_nan
		endif else if (x2 ge 0) and (x2 le 1.0) then aero[i,j]=x2
	endfor		 
	endfor	
	return,aero
end

function re_est,ref_0,wm
	aero2=fltarr(2,1200,1200)
	for j=0,1199 do begin
		if (j mod 100) eq 0 then print,j
		for i=0,1199 do begin
			temp=set_ref2x(dem[i,j],inc[i,j],rad[i,j],aero[i,j])
			aero2[i,j]=temp
		endfor
	endfor
	return,rcoef
end


	

