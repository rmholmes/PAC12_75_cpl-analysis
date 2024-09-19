clear all
close all
%clc
warning off all 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
OGCM_dir 	   = '/media/maillard/HelloWorld/CROCO/SODA/COUPLED_SIMU/';
OGCM_prefix    = 'soda3.4.2_5dy_ocean_reg_npac12_'; %% files pre-processed by process_SODA_outputs.csh (==> DODS format)
OGCM_prefix2   = 'soda3.4.2_5dy_ocean_reg_npac12_withBar_'; %% files names with bartropic U and V used to make CROCRO bry filers
%
Yorig         = 1980;          % reference time for vector time
                               % in croco initial and forcing files
%
Ymin          = 2014;          % first forcing year
Ymax          = 2019;          % last  forcing year
Mmin          = 1;             % first forcing month
Mmax          = 12;             % last  forcing month
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
missval=NaN;
  % Loop on the years and the months
  for Y=Ymin:Ymax
    if Y==Ymin 
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
      mo_max=12;
    end
    for M=mo_min:mo_max
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
	    ' - month ',num2str(M)])
      disp(' ')
      %
      Mm=M-1;Ym=Y;
      if Mm==0
	Mm=12;
	Ym=Y-1;
      end
      Mp=M+1;Yp=Y;
      if Mp==13
	Mp=1;
	Yp=Y+1;
      end
      %

      fname=([OGCM_dir,OGCM_prefix,'Y',num2str(Y),'M',num2str(M),'.nc']);
      fnamebar=([OGCM_dir,OGCM_prefix2,'Y',num2str(Y),'M',num2str(M),'.nc']);

  	disp('    Open SODA dods file')
  	nc2=netcdf(fname);

  	lonT=nc2{'xt_ocean'}(:);
  	latT=nc2{'yt_ocean'}(:);
  	lonU=nc2{'xu_ocean'}(:);
  	latU=nc2{'yu_ocean'}(:);
  	lonV=nc2{'xu_ocean'}(:);
  	latV=nc2{'yu_ocean'}(:);
  	depth=nc2{'st_ocean'}(:);
  	time=nc2{'time'}(:);

  	ssh=nc2{'ssh'}(:);
  	u=nc2{'u'}(:);
  	v=nc2{'v'}(:);
  	tauy=nc2{'taux'}(:);
  	taux=nc2{'tauy'}(:);
  	salt=nc2{'salt'}(:);
  	temp=nc2{'temp'}(:);
  	close(nc2)
	
	u(u==min(min(min(min(u)))))=nan;
	v(v==min(min(min(min(v)))))=nan;
	ssh(ssh==min(min(min(min(ssh)))))=nan;
	temp(temp==min(min(min(min(temp)))))=nan;
	salt(salt==min(min(min(min(salt)))))=nan;
	taux(taux==min(min(min(min(taux)))))=nan;
	tauy(tauy==min(min(min(min(tauy)))))=nan;
	
	min(min(min(min(temp))))

  	disp('    Create SODA file with Ubar and Vbar')

  	missval=NaN;
  	nc=netcdf(fnamebar,'clobber');
  	redef(nc);
  	nc('lonT')=length(lonT);
  	nc('latT')=length(latT);
  	nc('lonU')=length(lonU);
  	nc('latU')=length(latU);
  	nc('lonV')=length(lonV);
  	nc('latV')=length(latV);
  	nc('depth')=length(depth);
  	nc('time')=length(time);
  	nc{'temp'}=ncfloat('time','depth','latT','lonT') ;
  	nc{'temp'}.long_name=ncchar('TEMPERATURE');
  	nc{'temp'}.long_name='TEMPERATURE';
  	nc{'temp'}.units=ncchar('deg. C');
  	nc{'temp'}.units='deg. C';
  	nc{'temp'}.missing_value=missval;
  	nc{'salt'}=ncfloat('time','depth','latT','lonT') ;
  	nc{'salt'}.long_name=ncchar('SALINITY');
  	nc{'salt'}.long_name='SALINITY';
  	nc{'salt'}.units=ncchar('ppt');
  	nc{'salt'}.units='ppt';
  	nc{'salt'}.missing_value=missval;
  	nc{'u'}=ncfloat('time','depth','latU','lonU') ;
  	nc{'u'}.long_name=ncchar('ZONAL VELOCITY');
  	nc{'u'}.long_name='ZONAL VELOCITY';
  	nc{'u'}.units=ncchar('m/sec');
  	nc{'u'}.units='m/sec';
  	nc{'u'}.missing_value=missval;
  	nc{'v'}=ncfloat('time','depth','latV','lonV') ;
  	nc{'v'}.long_name=ncchar('MERIDIONAL VELOCITY');
  	nc{'v'}.long_name='MERIDIONAL VELOCITY';
  	nc{'v'}.units=ncchar('m/sec');
  	nc{'v'}.units='m/sec';
  	nc{'v'}.missing_value=missval;
  	nc{'ubar'}=ncfloat('time','latU','lonU') ;
  	nc{'ubar'}.long_name=ncchar('ZONAL BAROTROPIC VELOCITY');
  	nc{'ubar'}.long_name='ZONAL BAROTROPIC VELOCITY';
  	nc{'ubar'}.units=ncchar('m/sec');
  	nc{'ubar'}.units='m/sec';
  	nc{'ubar'}.missing_value=missval;
  	nc{'vbar'}=ncfloat('time','latV','lonV') ;
  	nc{'vbar'}.long_name=ncchar('MERIDIONAL BAROTROPIC VELOCITY');
  	nc{'vbar'}.long_name='MERIDIONAL BAROTROPIC VELOCITY';
  	nc{'vbar'}.units=ncchar('m/sec');
  	nc{'vbar'}.units='m/sec';
  	nc{'vbar'}.missing_value=missval;
  	nc{'taux'}=ncfloat('time','latU','lonU') ;
  	nc{'taux'}.long_name=ncchar('TAU_X');
  	nc{'taux'}.long_name='TAU_X';
  	nc{'taux'}.units=ncchar('N.m-2');
  	nc{'taux'}.units='N.m-2';
  	nc{'taux'}.missing_value=missval;
  	nc{'tauy'}=ncfloat('time','latV','lonV') ;
  	nc{'tauy'}.long_name=ncchar('TAU_Y');
  	nc{'tauy'}.long_name='TAU_Y';
  	nc{'tauy'}.units=ncchar('N.m-2');
  	nc{'tauy'}.units='N.m-2';
  	nc{'tauy'}.missing_value=missval;
  	nc{'ssh'}=ncfloat('time','latT','lonT') ;
  	nc{'ssh'}.long_name=ncchar('SEA LEVEL HEIGHT');
  	nc{'ssh'}.long_name='SEA LEVEL HEIGHT';
  	nc{'ssh'}.units=ncchar('m');
  	nc{'ssh'}.units='m';
  	nc{'ssh'}.missing_value=missval;
  	nc{'lonT'}=ncdouble('lonT') ;
  	nc{'lonT'}.units=ncchar('degrees_east');
  	nc{'lonT'}.units='degrees_east';
  	nc{'latT'}=ncdouble('latT') ;
  	nc{'latT'}.units=ncchar('degrees_north');
  	nc{'latT'}.units='degrees_north';
  	nc{'lonU'}=ncdouble('lonU') ;
  	nc{'lonU'}.units=ncchar('degrees_east');
  	nc{'lonU'}.units='degrees_east';
  	nc{'latU'}=ncdouble('latU') ;
  	nc{'latU'}.units=ncchar('degrees_north');
  	nc{'latU'}.units='degrees_north';
  	nc{'lonV'}=ncdouble('lonV') ;
  	nc{'lonV'}.units=ncchar('degrees_east');
  	nc{'lonV'}.units='degrees_east';
  	nc{'latV'}=ncdouble('latV') ;
  	nc{'latV'}.units=ncchar('degrees_north');
  	nc{'latV'}.units='degrees_north';
  	nc{'depth'}=ncdouble('depth') ;
  	nc{'depth'}.units=ncchar('meters');
  	nc{'depth'}.units='meters';
  	nc{'time'}=ncdouble('time') ;
  	eval(['nc{''time''}.units = ncchar(''days since 1-Jan-',num2str(Yorig),' 00:00:0.0'');'])
  	eval(['nc{''time''}.units = ''days since 1-Jan-',num2str(Yorig),' 00:00:0.0'';'])
  	endef(nc);
  	%
  	% Fill the file
  	%
  	disp('    Fill the OGCM file')
  	nc{'depth'}(:)=depth;
  	nc{'latT'}(:)=latT;
  	nc{'lonT'}(:)=lonT;
  	nc{'latU'}(:)=latU;
  	nc{'lonU'}(:)=lonU;
  	nc{'latV'}(:)=latV;
  	nc{'lonV'}(:)=lonV;
  	%
  	for tndx=1:length(time)
  	%
  	nc{'time'}(tndx)=time(tndx);
  	%

  	  nc{'ssh'}(tndx,:,:)=squeeze(ssh(tndx,:,:));
  	  nc{'taux'}(tndx,:,:)=squeeze(taux(tndx,:,:));
  	  nc{'tauy'}(tndx,:,:)=squeeze(tauy(tndx,:,:));
  	  u1=squeeze(u(tndx,:,:,:));
  	  v1=squeeze(v(tndx,:,:,:));
  	  nc{'u'}(tndx,:,:,:)=u1;
  	  nc{'v'}(tndx,:,:,:)=v1;
  	  nc{'temp'}(tndx,:,:,:)=squeeze(temp(tndx,:,:,:));
  	  nc{'salt'}(tndx,:,:,:)=squeeze(salt(tndx,:,:,:));

  	%
  	% Compute the barotropic velocities
  	%
  	masku=isfinite(u1);
  	maskv=isfinite(v1);
  	u1(isnan(u1))=0;
  	v1(isnan(v1))=0;
  	dz=gradient(depth);
  	NZ=length(depth);
  	du=0*squeeze(u1(1,:,:));
  	zu=du;
  	dv=0*squeeze(v1(1,:,:));
  	zv=dv;
  	for k=1:NZ
  	  du=du+dz(k)*squeeze(u1(k,:,:));
  	  zu=zu+dz(k)*squeeze(masku(k,:,:));
  	  dv=dv+dz(k)*squeeze(v1(k,:,:));
  	  zv=zv+dz(k)*squeeze(maskv(k,:,:));
  	end
  	du(zu==0)=NaN;
  	dv(zv==0)=NaN;
  	zu(zu==0)=NaN;
  	zv(zv==0)=NaN;
  	ubar=du./zu;
  	vbar=dv./zv;
  	%
  	nc{'ubar'}(tndx,:,:)=ubar;
  	nc{'vbar'}(tndx,:,:)=vbar;
  	%
  	end
  	%
  	close(nc)

	  

end

end



