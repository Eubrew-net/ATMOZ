%% This script computes the Jacobian matrix for Brewer ozone observation
%% uncertainty on measured data irradiance.

clear all;

close all;

clc;

warning off;

disp('--------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code (v9.0)             ')
disp('          Dr. El Gawhary O., Dr. Parra-Rojas F.C. and Redondas A.         ')
disp('                               (2017)                                     ')
disp('--------------------------------------------------------------------------')

% Brewer number
brewer_num = input('Enter the number of the Brewer: ');
brewer_str = num2str(brewer_num);
disp(' ')
switch brewer_str 
	case '157'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '183'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '185'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '17'
		disp('        -You are in Canada-       ')
		disp('Lat:---, Lon=---, Alt:--- m')
	otherwise 
		disp('-----------------------')
		disp('There is no such Brewer')
		disp('-----------------------')
end
disp(' ')

% we introduce the initial date
disp('Input the initial date')
disp('----------------------')
yyyy=input('Enter the year: ');
mm=input('Enter the month: ');
dd=input('Enter the day: ');
period = datenum(yyyy,mm,dd); % matlab time
dates=datestr(period,'yyyy-mm-dd');

disp(' ')

% we introduce the initial date
disp('Input the end date')
disp('------------------')
yyyy_f=input('Enter the year: ');
mm_f=input('Enter the month: ');
dd_f=input('Enter the day: ');
period_f = datenum(yyyy_f,mm_f,dd_f); % matlab time
dates_f=datestr(period_f,'yyyy-mm-dd');

% database urls
url_base='fparra:m23275108M@rbcce.aemet.es/eubrewnet'; % base url
url_func='/data/get/O3L1'; % ozone L1.0 url
url_head='/data/get/ConfigbyDate'; % Config by Date url
url_ds='/data/get/DS'; % Direct Sun url
url_dss='/data/get/DSS'; % Direct Sun Summary url

%%
% ----------------------------O3 Level 1 values------------------------------------------

url_o3l1=['"',url_base,url_func,'?brewerid=',brewer_str,'&date=',dates,'&enddate=',dates_f,'&format=text"']
o3l1 = download(url_o3l1);
%%
Omega = o3l1{12}'; % Calculated Ozone value with Standard algorithm + attenuation filter correction (DU)
t_j = o3l1{7}'; % continuous date index based in Matlab datenum
pre = o3l1{18}'; % Medium Pressure of the Brewer Location (mbar)
temp = o3l1{10}'; % Instrument temperature (ºC)
lon = o3l1{17}'; % Longitude of the Brewer Location (deg)
lat = o3l1{16}'; % Latitude of the Brewer Location (deg)
ms9_db = o3l1{20}'; % MS9, Second double ratio

%%
% this is a temporal artifact to avoid negative values in O3
for i = 1:length(t_j)
	if Omega(i) < 0
		Omega(i) = 0.1;
	end
end

f3 = o3l1{24}'; % corrected measurement of lambda3 (counts/s)


% ----------------------------- Direct Sun ---------------------------------------------------------
url_direct=['"',url_base,url_ds,'?brewerid=',brewer_str,'&date=',dates,'&enddate=',dates_f,'&format=text"']
DS = download(url_direct);

ccl = DS{12}; % number of the slit mask cycles
dark = DS{14}; % raw counts of the dark
rc1 = DS{15}; % raw counts of the incident light at 306.6 nm
rc2 = DS{16}; % raw counts of the incident light at 310.1 nm
rc3 = DS{17}; % raw counts of the incident light at 313.5 nm
rc4 = DS{18}; % raw counts of the incident light at 316.8 nm
rc5 = DS{19}; % raw counts of the incident light at 320.1 nm
nfilpos = DS{7}; % position of the neutral-density filter (0, 64, 128, 192, 256, 320)

% ----------------------------- Direct Sun Summary ---------------------------------------------------------
url_direct=['"',url_base,url_dss,'?brewerid=',brewer_str,'&date=',dates,'&enddate=',dates_f,'&format=text"']
DSS = download(url_direct);

Omega_5_dss=DSS{17}'; % Calculated Ozone value with Standard algorithm + attenuation filter correction from Symmary (DU)

% ratio of the raw counts. It is only a test. It's no real!!
xx_r = -rc2 + 0.5*rc3 + 2.2*rc4 - 1.7*rc5;

save('228_are.mat','Omega','ccl','dark','lat','lon','nfilpos','pre','rc1','rc2','rc3','rc4','rc5','t_j','temp')

return

% --------------------raw to count correction----------------------------------------------------


it = 0.1147; % standard integration time of the Brewer's PMT

% raw to count correction for the dark current
dark_c = 2*dark./(it.*ccl);

% raw to count correction for the incident light for each wavelength
f1_c = 2*rc1./(it.*ccl); % counts/s
f2_c = 2*rc2./(it.*ccl); % counts/s
f3_c = 2*rc3./(it.*ccl); % counts/s
f4_c = 2*rc4./(it.*ccl); % counts/s
f5_c = 2*rc5./(it.*ccl); % counts/s

% uncertainty of the raw to count correction. Expanded (k=2) noise.
u2_1_c = 2*(rc1 + dark)./(it.*ccl);
u2_2_c = 2*(rc2 + dark)./(it.*ccl);
u2_3_c = 2*(rc3 + dark)./(it.*ccl);
u2_4_c = 2*(rc4 + dark)./(it.*ccl);
u2_5_c = 2*(rc5 + dark)./(it.*ccl);

%%
% ----------------- Dark count correction --------------------------------------------------------

% the uncertainty of the dark current is the standard deviation of this
u2_d = (std(dark_c)).^2;

% dark count correction for each wavelenght
f1_d = f1_c - dark_c;
f2_d = f2_c - dark_c;
f3_d = f3_c - dark_c;
f4_d = f4_c - dark_c;
f5_d = f5_c - dark_c;

%%
% filter for small values of the counts
for i=1:length(t_j)
	if f1_d(i) < 2.0
		f1_d(i) = 2.0;
    end
	if f2_d(i) < 2.0
		f2_d(i) = 2.0;
    end
	if f3_d(i) < 2.0
		f3_d(i) = 2.0;
    end
	if f4_d(i) < 2.0
		f4_d(i) = 2.0;
    end
	if f5_d(i) < 2.0
		f5_d(i) = 2.0;
	end
end

%% 

% expanded uncertainty of the measurements including the dark corrention
u2_1_d = u2_1_c + u2_d;
u2_2_d = u2_2_c + u2_d;
u2_3_d = u2_3_c + u2_d;
u2_4_d = u2_4_c + u2_d;
u2_5_d = u2_5_c + u2_d;


% -------------------------dead time correction---------------------------------------
%%
f1_ct=[];
f2_ct=[];
f3_ct=[];
f4_ct=[];
f5_ct=[];

% uncertainty of the deadtime of the Brewer's PMT (Fountoulakis I. et al. 2016)
udt = 1e-9; % ns

% Deadtime correction for each wavelength (9 iterations)
for i = 1:length(t_j)
	f1_ct(1) = f1_d(i);
	f2_ct(1) = f2_d(i);
	f3_ct(1) = f3_d(i);
	f4_ct(1) = f4_d(i);
	f5_ct(1) = f5_d(i);
	for j=2:10
		f1_ct(j)=f1_d(i).*exp(f1_ct(j-1)*dt);
		f2_ct(j)=f2_d(i).*exp(f2_ct(j-1)*dt);
		f3_ct(j)=f3_d(i).*exp(f3_ct(j-1)*dt);
		f4_ct(j)=f4_d(i).*exp(f4_ct(j-1)*dt);
		f5_ct(j)=f5_d(i).*exp(f5_ct(j-1)*dt);
		
		f1ct(i) = f1_ct(j);
		f2ct(i) = f2_ct(j);
		f3ct(i) = f3_ct(j);
		f4ct(i) = f4_ct(j);
		f5ct(i) = f5_ct(j);
			
	end

% expanded uncertainty of the measures including the deadtime correction
	u2_f1_3(i) = ((exp(dt*f1ct(i))./(1-dt*f1ct(i))).^2).*u2_1_d(i) + (f1ct(i).^4)*udt^2;
	u2_f2_3(i) = ((exp(dt*f2ct(i))./(1-dt*f2ct(i))).^2).*u2_2_d(i) + (f2ct(i).^4)*udt^2;
	u2_f3_3(i) = ((exp(dt*f3ct(i))./(1-dt*f3ct(i))).^2).*u2_3_d(i) + (f3ct(i).^4)*udt^2;
	u2_f4_3(i) = ((exp(dt*f4ct(i))./(1-dt*f4ct(i))).^2).*u2_4_d(i) + (f4ct(i).^4)*udt^2;
	u2_f5_3(i) = ((exp(dt*f5ct(i))./(1-dt*f5ct(i))).^2).*u2_5_d(i) + (f5ct(i).^4)*udt^2;

% change to logarithmic space of the measurements	 
	f1_ct_1(i) = 10000*log10(f1ct(i));
	f2_ct_1(i) = 10000*log10(f2ct(i));
	f3_ct_1(i) = 10000*log10(f3ct(i));
	f4_ct_1(i) = 10000*log10(f4ct(i));
	f5_ct_1(i) = 10000*log10(f5ct(i));
	
end

% uncertainty in log space
u2_f1log = (10000*(sqrt(u2_f1_3)./(f1ct*log(10)))).^2;
u2_f2log = (10000*(sqrt(u2_f2_3)./(f2ct*log(10)))).^2;
u2_f3log = (10000*(sqrt(u2_f3_3)./(f3ct*log(10)))).^2;
u2_f4log = (10000*(sqrt(u2_f4_3)./(f4ct*log(10)))).^2;
u2_f5log = (10000*(sqrt(u2_f5_3)./(f5ct*log(10)))).^2;

% provisional uncertainty
u_rel_dt = sqrt((sqrt(u2_f2log)./f2_ct_1).^2 + ((0.5)*sqrt(u2_f3log)./f3_ct_1).^2 + ((2.2)*sqrt(u2_f4log)./f4_ct_1).^2 + ((1.7)*sqrt(u2_f5log)./f5_ct_1).^2);

%%
% -----------------------------Temperature correction--------------------------------------------

% uncertainty of the PMT temperature (the resolution of the PMT is 1ºC)
u_temp = 1/sqrt(3); 

% temperature coefficients and its uncertainties of the IZO2016 campaign
[tc,v,brw] = tblread('temp_izo16.csv',',');
brw_str = str2num(brewer_str);
brw_num = str2num(brw);

% temperature coefficients for each wavelenght
ind = find(brw_num==brw_str);
TC1 = tc(ind,1); % temperature coefficient of lambda1
TC2 = tc(ind,2); % temperature coefficient of lambda2
TC3 = tc(ind,3); % temperature coefficient of lambda3
TC4 = tc(ind,4); % temperature coefficient of lambda4
TC5 = tc(ind,5); % temperature coefficient of lambda5

% uncertainty for each temperature coefficient
utc1 = tc(ind,6); % uncertainty of the temperature coefficient of lambda1
utc2 = tc(ind,7); % uncertainty of the temperature coefficient of lambda2
utc3 = tc(ind,8); % uncertainty of the temperature coefficient of lambda3
utc4 = tc(ind,9); % uncertainty of the temperature coefficient of lambda4
utc5 = tc(ind,10); % uncertainty of the temperature coefficient of lambda5

% PMT temperature correction
f1_temp = f1_ct_1 + temp.*TC1;
f2_temp = f2_ct_1 + temp.*TC2;
f3_temp = f3_ct_1 + temp.*TC3;
f4_temp = f4_ct_1 + temp.*TC4;
f5_temp = f5_ct_1 + temp.*TC5;

% expanded uncertainty of the measurements including the temperature correction
u2_f1t = u2_f1log + (TC1*u_temp).^2 + (temp*utc1).^2;
u2_f2t = u2_f2log + (TC2*u_temp).^2 + (temp*utc2).^2;
u2_f3t = u2_f3log + (TC3*u_temp).^2 + (temp*utc3).^2;
u2_f4t = u2_f4log + (TC4*u_temp).^2 + (temp*utc4).^2;
u2_f5t = u2_f5log + (TC5*u_temp).^2 + (temp*utc5).^2;

% provisional uncertainty
u_rel_temp = sqrt((sqrt(u2_f2t)./f2_temp).^2 + ((0.5)*sqrt(u2_f3t)./f3_temp).^2 + ((2.2)*sqrt(u2_f4t)./f4_temp).^2 + ((1.7)*sqrt(u2_f5t)./f5_temp).^2);

%%
% ------------------------------ Neutral Filters correction -------------------------------------

% neutral filter corrections depending of the filter position
for i = 1:length(t_j)
	if nfilpos(i)./64 == 0
		f1_nd(i) = f1_temp(i) + nfil0;
		f2_nd(i) = f2_temp(i) + nfil0;
		f3_nd(i) = f3_temp(i) + nfil0;
		f4_nd(i) = f4_temp(i) + nfil0;
		f5_nd(i) = f5_temp(i) + nfil0;
	elseif nfilpos(i)./64 == 1
		f1_nd(i) = f1_temp(i) + nfil1;
		f2_nd(i) = f2_temp(i) + nfil1;
		f3_nd(i) = f3_temp(i) + nfil1;
		f4_nd(i) = f4_temp(i) + nfil1;
		f5_nd(i) = f5_temp(i) + nfil1;
	elseif nfilpos(i)./64 == 2
		f1_nd(i) = f1_temp(i) + nfil2;
		f2_nd(i) = f2_temp(i) + nfil2;
		f3_nd(i) = f3_temp(i) + nfil2;
		f4_nd(i) = f4_temp(i) + nfil2;
		f5_nd(i) = f5_temp(i) + nfil2;
	elseif nfilpos(i)./64 == 3
		f1_nd(i) = f1_temp(i) + nfil3;
		f2_nd(i) = f2_temp(i) + nfil3;
		f3_nd(i) = f3_temp(i) + nfil3;
		f4_nd(i) = f4_temp(i) + nfil3;
		f5_nd(i) = f5_temp(i) + nfil3;
	elseif nfilpos(i)./64 == 4
		f1_nd(i) = f1_temp(i) + nfil4;
		f2_nd(i) = f2_temp(i) + nfil4;
		f3_nd(i) = f3_temp(i) + nfil4;
		f4_nd(i) = f4_temp(i) + nfil4;
		f5_nd(i) = f5_temp(i) + nfil4;
	elseif nfilpos(i)./64 == 5
		f1_nd(i) = f1_temp(i) + nfil5;
		f2_nd(i) = f2_temp(i) + nfil5;
		f3_nd(i) = f3_temp(i) + nfil5;
		f4_nd(i) = f4_temp(i) + nfil5;
		f5_nd(i) = f5_temp(i) + nfil5;
	end
end

% double ratio
xx_nf = -f2_nd + 0.5*f3_nd + 2.2*f4_nd - 1.7*f5_nd; 

% uncertainty of the measurements. We consider that the neutral filters uncertainty don't affect 
u_rel_ms9 = sqrt((sqrt(u2_f2t)./f2_nd).^2 + ((0.5)*sqrt(u2_f3t)./f3_nd).^2 + ((2.2)*sqrt(u2_f4t)./f4_nd).^2 + ((1.7)*sqrt(u2_f5t)./f5_nd).^2);

u2_ms9_r = (u_rel_ms9.*xx_nf).^2;

save('17.mat','Omega','ccl','dark','lat','lon','nfilpos','pre','rc1','rc2','rc3','rc4','rc5','t_j','temp')

%%
% ------------------ Definitions of some parameters ------------------------------------------------

Pstan = 1013.25; % standard pressure (mbars)

R = 6371.229e3; % Earth's radius (m)

r = 2370.0; % heigth of Izaña station (m)


%-------------------------------- SZA uncertaiinty---------------------------------------------------
%                             through astronomical formulas
% Spencer, J.W. (1971) Fourier Series Representation of the Position of the Sun. Search, 2, 162-172.
%----------------------------------------------------------------------------------------------------
%%
% from deg to rad
p0 = pi/180;

% time of the measurements (min)
date_vec=datevec(t_j);
days=fix(t_j-datenum(1965,1,1));
t0 = date_vec(:,4)*60.0 + date_vec(:,5) + date_vec(:,6)/60.0;

% time uncertainty (min)
u_t0 = 1/60;

% eccentricity correction factor and unceretainty
t = (days+1)/365.2422;
ElipLong_1965 = 279.4574;
I = (ElipLong_1965 + 360*t + (t0'/1460.97))*p0;
u2_I = ((p0/1460.97)*u_t0)^2;

% equation of time in seconds and the uncertainty
et = 4.2*sin(3*I)-2*cos(2*I)+596.5*sin(2*I)-12.8*sin(4*I)+19.3*cos(3*I)-(102.5+0.142*t).*sin(I)+(0.033*t-429.8).*cos(I);
u2_et = u2_I.*(3*4.2*cos(3*I)+4*sin(2*I)+2*596.5*cos(2*I)-4*12.8*cos(4*I)-3*19.3*sin(3*I)-(102.5+0.142*t).*cos(I)-(0.033*t-429.8).*sin(I)).^2;

% the hour angle and uncertainty
ha = (t0'+(et/60)-720-4*lon)*p0/4;
u2_ha = ((p0/4)*u_t0)^2 +u2_et*(p0/240)^2;

% solar declination and uncertainty
dec = atan(0.4336*sin(I-(p0*et)/240));
u2_dec = (0.4336^2)*u2_I.*(cos(I-(p0*et)/240)./(1+(0.4336*sin(I-(p0*et)/240)).^2)).^2 + (((p0/240)*0.4336)^2)*u2_et.*(cos(I-(p0*et)/240)./(1+(0.4336*sin(I-(p0*et)/240)).^2)).^2;

% solar zenital angle and uncertainty
sza = acos(cos(ha).*cos(dec).*cos(lat*p0)+sin(lat*p0).*sin(dec))/p0;
u2_sza = ((sin(ha).*cos(dec).*cos(lat*p0)).^2).*u2_ha./((p0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*p0)).^2)).^2)+((cos(ha).*sin(dec).*cos(lat*p0)+sin(lat*p0).*cos(dec)).^2).*u2_dec./((p0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*p0)).^2)).^2);

% relative uncertainty of the SZA
urel_sza = 100*sqrt(u2_sza)./sza;

% Figure 1 - SZA relative uncertainty
figure(1)
plot(sza,urel_sza)
title('SZA Relative uncertainty');
xlabel('SZA, grad')
ylabel('Relative uncertainty, %')
grid on
%%

%----------------------------Rayleigh air mass uncertainty----------------------------------------------
heff_r = 5e3; % Rayleigh effective altitude (m)
u_heff_r = 2e2; % uncertainty of the Rayleigh effective altitude (m)

% Rayleigh air mass and uncertainty
sx_r = (R/(R+heff_r))*sin(sza*p0);
m = sec(asin(sx_r)); % air mass

u2_m = ((sx_r.^2)/((1-sx_r.^2).^3)).*(((R^2)/((R+heff_r)^4)).*((sin(sza*p0)).^2)*(u_heff_r^2) + ((R^2)/((R+heff_r)^2)).*((cos(sza*p0)).^2).*(u2_sza));
urel_m = 100*sqrt(u2_m)./m;

% Figure 2 - Rayleigh air mass relative uncertainty 
figure(2)
plot(m,urel_m)
title('Rayleigh air mass relative uncertainty');
xlabel('Rayleigh air mass')
ylabel('Relative uncertainty, %')
grid on

% ---------------------------ozone airmass uncertainty----------------------------------------------
heff_o = 22e3; % Ozone effective altitude (m)
u_heff_o3 = 2e3; % uncertainty of the ozone effective altitude (m)

u_heff_clim = 1.5e3; % uncertainty of the ozone effective altitude through climatology

% combined uncertainty
u_heff_o = sqrt(u_heff_o3^2 + u_heff_clim^2);

% Ozone air mass and uncertainty
sx_o = (R/(R+heff_o))*sin(sza*p0);
mu = sec(asin(sx_o)); % air mass

u2_mu=((sx_o.^2)/((1-sx_o.^2).^3)).*(((R^2)/((R+heff_o)^4)).*((sin(sza*p0)).^2)*(u_heff_o^2) + ((R^2)/((R+heff_o)^2)).*((cos(sza*p0)).^2).*(u2_sza));
urel_mu = 100*sqrt(u2_mu)./mu;

% Figure 3 - Ozone air mass relative uncertainty
figure(3)
plot(mu,urel_mu)
title('Ozone air mass relative uncertainty');
xlabel('Ozone air mass')
ylabel('Relative uncertainty, %')
grid on

% ----------------------------- wavelength error associated to the measurement -------------------------------
% obtained from dispresion

[wv,c,b_wv] = tblread('ms9_wv_test.csv',',');
brw_num_wv = str2num(b_wv);

ind_wv = find(brw_num_wv==brw_str);
k_wv = wv(ind,1); % Wavelength variation (percent)
u_k = wv(ind,2); % uncertainty in wavelength (percent)

% wavelength correction and uncertainty
MS9_wv = xx_nf*(1+k_wv/100);
u2_wv = ((1+k_wv/100)^2).*u2_ms9_r + ((xx_nf/100).^2).*u_k^2;

%%
% ----------------------------- mean of five measurements -----------------------------------------------------

media_ms9=[];
standard=[];
u2_xx_med=[];
u_xx_sta=[];

sza_med=[];

mu_med=[];
u2_mu_med=[];

m_med=[];
u2_m_med=[];

t_j_med=[];

Omega_med=[];

for i=0:1:(size(xx_nf,2)/5)-1
    ms9_m=mean(MS9_wv(1,1+5*i:5+5*i)); % mean of MS9 (each 5 measurements)
	stan_m=std(MS9_wv(1,1+5*i:5+5*i)); % uncertainty
	u2_xx_m = (1/5)*mean(u2_wv(1,1+5*i:5+5*i)); % uncertainty of the mean of MS9
	u_xx_std = std(u2_wv(1,1+5*i:5+5*i));
	media_ms9=[media_ms9;ms9_m];
	standard=[standard;stan_m];	
	u2_xx_med=[u2_xx_med;u2_xx_m];
	u_xx_sta=[u_xx_sta;u_xx_std];
	
	Omega_m=mean(Omega(1,1+5*i:5+5*i)); % mean of ozone each 5 measurements
	Omega_med=[Omega_med;Omega_m];
	
	sza_m = mean(sza(1,1+5*i:5+5*i)); % mean of SZA (each 5 measurements)
	
	mu_m = median(mu(1,1+5*i:5+5*i)); % median of ozone optical mass
	mu_med=[mu_med;mu_m]; 
	u2_mu_m = median(u2_mu(1,1+5*i:5+5*i));
	u2_mu_med = [u2_mu_med;u2_mu_m];
		
	t_j_m = median(mu(1,1+5*i:5+5*i)); % median of time each five measurements
	t_j_med = [t_j_med;t_j_m];
	
	m_m = median(m(1,1+5*i:5+5*i)); % median of Rayleigh optical mass
	m_med=[m_med;m_m];
	u2_m_m = median(u2_m(1,1+5*i:5+5*i));
    u2_m_med = [u2_m_med;u2_m_m];
	
	sza_med=[sza_med;sza_m];
end
ms9_5 = media_ms9';
u2_5 = (standard').^2;
u2_xx_5 = u2_xx_med';
u_xx_5_std = u_xx_sta';

Omega_5 = Omega_med';

sza_5 = sza_med';

mu_5 = mu_med';
u2_mu_5 = u2_mu_med';

m_5 = m_med';
u2_m_5 = u2_m_med';

t_j_5 = t_j_med';

pre_5 = mean(pre)*ones(1,length(t_j)/5);

% MS9
MS9_r = ms9_5;

%----------------------------cross correlation matrix calculation (Omar)----------------------------------------------

% Measurements substracting the ETC
N_meas = (MS9_r)+1; 

A_nominal = zeros(1,length(sza)/5);
Omega_nominal = zeros(1,length(sza)/5);
mu_nominal = zeros(1,length(sza)/5);
B_nominal = zeros(1,length(sza)/5);
m_nominal = zeros(1,length(sza)/5);
pre_nominal = zeros(1,length(sza)/5);
etc_nominal = zeros(1,length(sza)/5);
%%
for r=1:length(sza)/5
                    %% Nominal values parameters
                    
                    fprintf('status computation = %2.1f  \n',(r/length(sza)/5)*100);

                    A_nominal(r) = A(r);

                    Omega_nominal(r) = Omega_5(r);

                    mu_nominal(r) = mu_5(r);

                    B_nominal(r) = B(r);
					
					m_nominal(r) = m_5(r);
					
					pre_nominal(r) = pre_5(r);
					
					etc_nominal(r) = etc(r);


                    %% Here Levenberg algorithm would start.

                    % The present functions computes the chi2 distribution of

                    % intialization
                    flag  = 0; % variable to be set to 1 when the stop condition is reached

                    count = 0; % variable that counts how many times chi2 did not have changes
                               % higher than 

                    %% Choose a starting value for the parameter l in the algorithm

                    l = 0.0001;

                    % -------------------------------------------------------------------------
                    % sigma is the noise level. It will be assumed to be the same in every
                    % pixel
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                    sigma = 1e-4;
                   
                    % Now we generate the first chi2

                    [chi2_min, N] = generate_chi2(N_meas(r),sigma,A_nominal(r),Omega_nominal(r),mu_nominal(r), B_nominal(r), m_nominal(r),pre_nominal(r),etc_nominal(r));

                    [alpha,beta] = compute_jacobian(N, N_meas(r), A_nominal(r),Omega_nominal(r),mu_nominal(r), B_nominal(r),m_nominal(r),pre_nominal(r),etc_nominal(r),sigma);


                    %% Compute Jacobian and beta vector necessary for the iterative algorithm

                    % Iterations start from here

                    while (count~=2) 

                            % Compute the alpha_prime which is required in the algorithm
                            % basis in 6 dimensions

                            alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4) alpha(1,5) alpha(1,6) alpha(1,7);...
                                            alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4) alpha(2,5) alpha(2,6) alpha(2,7);...
                                            alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4) alpha(3,5) alpha(3,6) alpha(3,7);...
                                            alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l) alpha(4,5) alpha(4,6) alpha(4,7);...
											alpha(5,1) alpha(5,2) alpha(5,3) alpha(5,4) alpha(5,5)*(1+l) alpha(5,6) alpha(5,7);...
											alpha(6,1) alpha(6,2) alpha(6,3) alpha(6,4) alpha(6,5) alpha(6,6)*(1+l) alpha(6,7);
											alpha(7,1) alpha(7,2) alpha(7,3) alpha(7,4) alpha(7,5) alpha(7,6) alpha(7,7)*(1+l)];

                            % Then we compute the parameters variations

                            delta_param = alpha_prime\beta';

                            %Compute the new parameters

                            A_nominal_new(r) = A_nominal(r)  + delta_param(1);

                            Omega_nominal_new(r)  = Omega_nominal(r)  + delta_param(2);

                            mu_nominal_new(r)  = mu_nominal(r)  + delta_param(3);

                            B_nominal_new(r)  = B_nominal(r)  + delta_param(4);
							
							m_nominal_new(r)  = m_nominal(r)  + delta_param(5);
							
							pre_nominal_new(r)  = pre_nominal(r)  + delta_param(6);
							
							etc_nominal_new(r)  = etc_nominal(r)  + delta_param(7);

                            % Compute the new chi2

                            [chi2_new, N] = generate_chi2(N_meas(r),sigma,A_nominal_new(r),Omega_nominal_new(r),mu_nominal_new(r),B_nominal_new(r),m_nominal_new(r),pre_nominal_new(r),etc_nominal_new(r));

                            diff_chi2 = (chi2_new-chi2_min);

                            if (diff_chi2<0)

								l = l/10

                                A_nominal(r)  = A_nominal_new(r); 

                                Omega_nominal(r)  = Omega_nominal_new(r); 

                                mu_nominal(r)  = mu_nominal_new(r); 

                                B_nominal(r)  = B_nominal_new(r); 
								
								m_nominal(r)  = m_nominal_new(r);
								
								pre_nominal(r)  = pre_nominal_new(r);
								
								etc_nominal(r)  = etc_nominal_new(r);

                                [alpha, beta] = compute_jacobian(N,N_meas(r),A_nominal(r),Omega_nominal(r),mu_nominal(r),B_nominal(r),m_nominal(r),pre_nominal(r), etc_nominal(r),sigma);

                                chi2_aux = chi2_min;

                                chi2_min = chi2_new;

                                if (abs(diff_chi2/chi2_aux)<0.00000001)

                                    count = count +1;

                                end

                            else

                                count = 0;

                                l = l*10

                            end

                            if (diff_chi2.*1e10==0)

                                count =2;
                             
							end
                    
						end


                    l = 0;


                    alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4) alpha(1,5) alpha(1,6) alpha(1,7);...
                                            alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4) alpha(2,5) alpha(2,6) alpha(2,7);...
                                            alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4) alpha(3,5) alpha(3,6) alpha(3,7);...
                                            alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l) alpha(4,5) alpha(4,6) alpha(4,7);...
											alpha(5,1) alpha(5,2) alpha(5,3) alpha(5,4) alpha(5,5)*(1+l) alpha(5,6) alpha(5,7);...
											alpha(6,1) alpha(6,2) alpha(6,3) alpha(6,4) alpha(6,5) alpha(6,6)*(1+l) alpha(6,7);...
											alpha(7,1) alpha(7,2) alpha(7,3) alpha(7,4) alpha(7,5) alpha(7,6) alpha(7,7)*(1+l)];
                    
                    
                    alpha_prime_corrected = alpha_prime+0.3.*eye(7,7);
                    
                    C  = inv(alpha_prime_corrected);
                    
                    uncert_A(r)  = sigma.*sqrt(abs(C(1,1)));

                    uncert_Omega(r)  = sigma.*sqrt(abs(C(2,2)));

                    uncert_mu(r)  = sigma.*sqrt(abs(C(3,3)));

                    uncert_B(r)  = sigma.*sqrt(abs(C(4,4)));
					
					uncert_m(r)  = sigma.*sqrt(abs(C(5,5)));
					
					uncert_pre(r)  = sigma.*sqrt(abs(C(6,6)));
					
					uncert_etc(r)  =sigma.*sqrt(abs(C(7,7)));
                     
                     

                    rho(r,1,1) = C(1,1)/(sqrt(abs(C(1,1)*C(1,1))));
                    rho(r,1,2) = C(1,2)/(sqrt(abs(C(1,1)*C(2,2))));
                    rho(r,1,3) = C(1,3)/(sqrt(abs(C(1,1)*C(3,3))));
                    rho(r,1,4) = C(1,4)/(sqrt(abs(C(1,1)*C(4,4))));
					rho(r,1,5) = C(1,5)/(sqrt(abs(C(1,1)*C(5,5))));
					rho(r,1,6) = C(1,6)/(sqrt(abs(C(1,1)*C(6,6))));
					rho(r,1,7) = C(1,7)/(sqrt(abs(C(1,1)*C(7,7))));
					
                    rho(r,2,1) = C(2,1)/(sqrt(abs(C(2,2)*C(1,1))));
                    rho(r,2,2) = C(2,2)/(sqrt(abs(C(2,2)*C(2,2))));
                    rho(r,2,3) = C(2,3)/(sqrt(abs(C(2,2)*C(3,3))));
                    rho(r,2,4) = C(2,4)/(sqrt(abs(C(2,2)*C(4,4))));
					rho(r,2,5) = C(2,5)/(sqrt(abs(C(2,2)*C(5,5))));
					rho(r,2,6) = C(2,6)/(sqrt(abs(C(2,2)*C(6,6))));
					rho(r,2,7) = C(2,7)/(sqrt(abs(C(2,2)*C(7,7))));
                    
					rho(r,3,1) = C(3,1)/(sqrt(abs(C(3,3)*C(1,1))));
                    rho(r,3,2) = C(3,2)/(sqrt(abs(C(3,3)*C(2,2))));
                    rho(r,3,3) = C(3,3)/(sqrt(abs(C(3,3)*C(3,3))));
                    rho(r,3,4) = C(3,4)/(sqrt(abs(C(3,3)*C(4,4))));
					rho(r,3,5) = C(3,5)/(sqrt(abs(C(3,3)*C(5,5))));
					rho(r,3,6) = C(3,6)/(sqrt(abs(C(3,3)*C(6,6))));
					rho(r,3,7) = C(3,7)/(sqrt(abs(C(3,3)*C(7,7))));
					
					rho(r,4,1) = C(4,1)/(sqrt(abs(C(4,4)*C(1,1))));
                    rho(r,4,2) = C(4,2)/(sqrt(abs(C(4,4)*C(2,2))));
                    rho(r,4,3) = C(4,3)/(sqrt(abs(C(4,4)*C(3,3))));
                    rho(r,4,4) = C(4,4)/(sqrt(abs(C(4,4)*C(4,4))));
					rho(r,4,5) = C(4,5)/(sqrt(abs(C(4,4)*C(5,5))));
					rho(r,4,6) = C(4,6)/(sqrt(abs(C(4,4)*C(6,6))));
					rho(r,4,7) = C(4,7)/(sqrt(abs(C(4,4)*C(7,7))));
					
					rho(r,5,1) = C(5,1)/(sqrt(abs(C(5,5)*C(1,1))));
					rho(r,5,2) = C(5,2)/(sqrt(abs(C(5,5)*C(2,2))));
					rho(r,5,3) = C(5,3)/(sqrt(abs(C(5,5)*C(3,3))));
					rho(r,5,4) = C(5,4)/(sqrt(abs(C(5,5)*C(4,4))));
					rho(r,5,5) = C(5,5)/(sqrt(abs(C(5,5)*C(5,5))));
					rho(r,5,6) = C(5,6)/(sqrt(abs(C(5,5)*C(6,6))));
					rho(r,5,7) = C(5,7)/(sqrt(abs(C(5,5)*C(7,7))));
					
					rho(r,6,1) = C(6,1)/(sqrt(abs(C(6,6)*C(1,1))));
					rho(r,6,2) = C(6,2)/(sqrt(abs(C(6,6)*C(2,2))));
					rho(r,6,3) = C(6,3)/(sqrt(abs(C(6,6)*C(3,3))));
					rho(r,6,4) = C(6,4)/(sqrt(abs(C(6,6)*C(4,4))));
					rho(r,6,5) = C(6,5)/(sqrt(abs(C(6,6)*C(5,5))));
					rho(r,6,6) = C(6,6)/(sqrt(abs(C(6,6)*C(6,6))));
					rho(r,6,7) = C(6,7)/(sqrt(abs(C(6,6)*C(7,7))));
					
					rho(r,7,1) = C(7,1)/(sqrt(abs(C(7,7)*C(1,1))));
					rho(r,7,2) = C(7,2)/(sqrt(abs(C(7,7)*C(2,2))));
					rho(r,7,3) = C(7,3)/(sqrt(abs(C(7,7)*C(3,3))));
					rho(r,7,4) = C(7,4)/(sqrt(abs(C(7,7)*C(4,4))));
					rho(r,7,5) = C(7,5)/(sqrt(abs(C(7,7)*C(5,5))));
					rho(r,7,6) = C(7,6)/(sqrt(abs(C(7,7)*C(6,6))));
					rho(r,7,7) = C(7,7)/(sqrt(abs(C(7,7)*C(7,7))));

end

% Figure 4 - Cross Correlation between A and Omega
figure(4)
plot(t_j_5,rho(:,1,2));
datetick('x',dateFormat)
title('Correlation A-\Omega');
xlabel('Time')
ylabel('Correlation')

% Figure 5 - Cross Correlation between mu and Omega
figure(5)
plot(t_j_5,rho(:,2,3));
datetick('x',dateFormat)
title('Correlation Omega-\mu');
xlabel('Time')
ylabel('Correlation')

% Figure 6 - Cross Correlation between A and mu
figure(6)
plot(t_j_5,rho(:,1,3));
datetick('x',dateFormat)
title('Correlation \mu-A');
xlabel('Time')
ylabel('Correlation')

% Figure 7 - Cross Correlation between A and B
figure(7)
plot(t_j_5,rho(:,1,4));
datetick('x',dateFormat)
title('Correlation A-B');
xlabel('Time')
ylabel('Correlation')

% Figure 8 - Cross Correlation between B and Omega
figure(8)
plot(t_j_5,rho(:,2,4));
datetick('x',dateFormat)
title('Correlation B-\Omega');
xlabel('Time')
ylabel('Correlation')

% Figure 9 - Cross Correlation between B and mu
figure(9)
plot(t_j_5,rho(:,3,4));
datetick('x',dateFormat)
title('Correlation B-\mu');
xlabel('Time')
ylabel('Correlation')

% Figure 10 - Cross Correlation between A and m
figure(10)
plot(t_j_5,rho(:,1,5));
datetick('x',dateFormat)
title('Correlation A-m');
xlabel('Time')
ylabel('Correlation')

% Figure 11 - Cross Correlation between A and p
figure(11)
plot(t_j_5,rho(:,1,6));
datetick('x',dateFormat)
title('Correlation A-p');
xlabel('Time')
ylabel('Correlation')

% Figure 12 - Cross Correlation between A and ETC
figure(12)
plot(t_j_5,rho(:,1,7));
datetick('x',dateFormat)
title('Correlation A-ETC');
xlabel('Time')
ylabel('Correlation')

% MS9 with Rayleigh correction
meas = N_meas;
o3_meas = (meas-m_5.*pre_5.*B/Pstan)./(A.*mu_5);

% O3 without Rayleigh correction
o3_meas_r = meas./(A.*mu_5);
%%

% uncertainty of the pressure
u_pre = 15.0; % mbar

% uncertainty of the measurements
u_meas = sqrt(u2_xx_5);

% Uncertainty of the ETC (counts/s)
u_etc_f = 5*ones(1,length(t_j_5));

% Uncertainty of the ozone absorption coefficient
%-------------------------------------------------

% from climatology
unc_teff = 5.0;

% percentage gradient for the BP and SG cross sections (Redondas A. et al., 2014)
grad_bp = 9.3601e-2;
grad_sg = 9.6391e-3;

[u_bp,u_sg]=sigma2(brewer_str); % fuction for the cross secion uncertainty
u_boot = cross(brewer_str); % function for the Bootstraping method uncertainty

%%
% variation of the ozone absorption coefficient with the wavelength obtained from dispersion (absolute)
if brewer_num == 17
	Astep = 0.0012;
else
	Astep = 0.00109;
end

% uncertainty of the ozone absorption coefficient for each cross section
u_A_bp = sqrt(u_bp^2 + u_boot^2 + (unc_teff*grad_bp*mean(A)/100)^2 + (Astep)^2);
u_A_sg = sqrt(u_sg^2 + u_boot^2 + (unc_teff*grad_sg*mean(A)/100)^2 + (Astep)^2);

% Fix relative uncertainty (percent/100) and uncertainty of the Rayleigh coefficient
urel_B = 0.01;
u_B = urel_B*B;

%%
% Ozone uncertainty with fixed uncertainties for BP cross section
a = (mu_5.*A).^2;

b = 2*mu_5.*A.*(Omega_5.*A.*sqrt(u2_mu_5).*rho(:,2,3)' + mu_5.*Omega_5.*u_A_bp.*rho(:,1,2)' + (pre_5/Pstan).*B.*sqrt(u2_m_5).*rho(:,1,5)' +...  
    (m_5.*B/Pstan)*u_pre.*rho(:,2,6)' + m_5.*(pre_5/Pstan).*u_B.*rho(:,2,4)' + u_etc_f.*rho(:,2,7)');

c = u2_mu_5.*(Omega_5.*A).^2 + (mu_5.*Omega_5.*u_A_bp).^2 + ((B.*pre_5/Pstan).^2).*u2_m_5 + (m_5.*(pre_5/Pstan).*u_B).^2 +...
    ((m_5.*B./Pstan).^2).*u_pre.^2 + u_etc_f.^2 + 2*(Omega_5.^2).*A.*mu_5.*sqrt(u2_mu_5).*u_A_bp.*rho(:,1,3)' +...
	2*Omega_5.*A.*m_5.*(pre_5/Pstan).*sqrt(u2_mu_5).*u_B.*rho(:,3,4)' + 2*mu_5.*Omega_5.*m_5.*(pre_5/Pstan).*u_A_bp.*u_B.*rho(:,1,4)'+...
    2*(Omega_5.*A.*pre_5.*B/Pstan).*sqrt(u2_mu_5).*sqrt(u2_m_5).*rho(:,3,5)' + 2*(Omega_5.*A.*m_5.*B/Pstan).*sqrt(u2_mu_5).*u_pre.*rho(:,3,6)'+...
    2*(mu_5.*Omega_5.*pre_5.*B./Pstan).*u_A_bp.*sqrt(u2_m_5).*rho(:,1,5)' + 2*(mu_5.*Omega_5.*m_5.*B./Pstan).*u_A_bp.*u_pre.*rho(:,1,6)' +...
    2*(pre_5.*B.*B.*m_5./(Pstan^2)).*sqrt(u2_m_5).*u_pre.*rho(:,5,6)' + 2*(B.*m_5.*pre_5.*pre_5/(Pstan^2)).*sqrt(u2_m_5).*u_B.*rho(:,4,5)' +...
    2*(B.*pre_5.*m_5.*m_5./(Pstan^2)).*u_pre.*u_B.*rho(:,4,6)' + 2*Omega_5.*A.*u_etc_f.*sqrt(u2_mu_5).*rho(:,3,7)'+...
    2*mu_5.*Omega_5.*u_etc_f.*u_A_bp.*rho(:,1,7)' + 2*B.*(pre_5/Pstan).*u_etc_f.*sqrt(u2_m_5).*rho(:,5,7)'+...
    2*(m_5.*B/Pstan).*u_etc_f.*u_pre.*rho(:,6,7)' + 2*(m_5.*pre_5/Pstan).*u_etc_f.*u_B.*rho(:,4,7)' + u2_xx_5 ;

uo3_pos = abs((-b+sqrt(b.^2 - 4*a.*c))./(2*a));
uo3_neg = abs((-b-sqrt(b.^2 - 4*a.*c))./(2*a));

% we select the lower uncertainty 
if uo3_pos > uo3_neg
	uo3 = uo3_neg;
elseif uo3_pos < uo3_neg 
	uo3 = uo3_pos;
else 
	uo3 = uo3_pos;
end

% Figure 13 - Ozone Absolute uncertainty
figure(13)
plot(Omega_5.*mu_5,uo3,'+')
title('Ozone absolute uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
grid on

% Figure 14 - Ozone Relative uncertainty
figure(14)
plot(sza_5,100.*uo3./Omega_5,'+')
title('Ozone relative uncertainty');
xlabel('Solar Zenith Angle, deg')
ylabel('Relative Uncertainty, %')
grid on

%%
% Ozone uncertainty with fixed uncertainties for SG cross section
a = (mu_5.*A).^2;

b_sg = 2*mu_5.*A.*(Omega_5.*A.*sqrt(u2_mu_5).*rho(:,2,3)' + mu_5.*Omega_5.*u_A_sg.*rho(:,1,2)' + (pre_5/Pstan).*B.*sqrt(u2_m_5).*rho(:,1,5)' +...  
    (m_5.*B/Pstan)*u_pre.*rho(:,2,6)' + m_5.*(pre_5/Pstan).*u_B.*rho(:,2,4)' + u_etc_f.*rho(:,2,7)');

c_sg = u2_mu_5.*(Omega_5.*A).^2 + (mu_5.*Omega_5.*u_A_sg).^2 + ((B.*pre_5/Pstan).^2).*u2_m_5 + (m_5.*(pre_5/Pstan).*u_B).^2 +...
    ((m_5.*B./Pstan).^2).*u_pre.^2 + u_etc_f.^2 + 2*(Omega_5.^2).*A.*mu_5.*sqrt(u2_mu_5).*u_A_sg.*rho(:,1,3)' +...
	2*Omega_5.*A.*m_5.*(pre_5/Pstan).*sqrt(u2_mu_5).*u_B.*rho(:,3,4)' + 2*mu_5.*Omega_5.*m_5.*(pre_5/Pstan).*u_A_sg.*u_B.*rho(:,1,4)'+...
    2*(Omega_5.*A.*pre_5.*B/Pstan).*sqrt(u2_mu_5).*sqrt(u2_m_5).*rho(:,3,5)' + 2*(Omega_5.*A.*m_5.*B/Pstan).*sqrt(u2_mu_5).*u_pre.*rho(:,3,6)'+...
    2*(mu_5.*Omega_5.*pre_5.*B./Pstan).*u_A_sg.*sqrt(u2_m_5).*rho(:,1,5)' + 2*(mu_5.*Omega_5.*m_5.*B./Pstan).*u_A_sg.*u_pre.*rho(:,1,6)' +...
    2*(pre_5.*B.*B.*m_5./(Pstan^2)).*sqrt(u2_m_5).*u_pre.*rho(:,5,6)' + 2*(B.*m_5.*pre_5.*pre_5/(Pstan^2)).*sqrt(u2_m_5).*u_B.*rho(:,4,5)' +...
    2*(B.*pre_5.*m_5.*m_5./(Pstan^2)).*u_pre.*u_B.*rho(:,4,6)' + 2*Omega_5.*A.*u_etc_f.*sqrt(u2_mu_5).*rho(:,3,7)'+...
    2*mu_5.*Omega_5.*u_etc_f.*u_A_sg.*rho(:,1,7)' + 2*B.*(pre_5/Pstan).*u_etc_f.*sqrt(u2_m_5).*rho(:,5,7)'+...
    2*(m_5.*B/Pstan).*u_etc_f.*u_pre.*rho(:,6,7)' + 2*(m_5.*pre_5/Pstan).*u_etc_f.*u_B.*rho(:,4,7)' + u2_xx_5 ;

uo3_pos_sg = abs((-b_sg+sqrt(b_sg.^2 - 4*a.*c_sg))./(2*a));
uo3_neg_sg = abs((-b_sg-sqrt(b_sg.^2 - 4*a.*c_sg))./(2*a));

% we select the lower uncertainty 
if uo3_pos_sg > uo3_neg_sg
	uo3_sg = uo3_neg_sg;
elseif uo3_pos_sg < uo3_neg_sg 
	uo3_sg = uo3_pos_sg;
else 
	uo3_sg = uo3_pos_sg;
end

% Figure 15 - Ozone Relative uncertainty versus BP and SG cross sections
figure(15)
plot(sza_5,100.*uo3./Omega_5,'+',sza_5,100.*uo3_sg./Omega_5,'+')
title('Ozone relative uncertainty');
xlabel('Solar Zenith Angle, deg')
ylabel('Relative Uncertainty, %')
legend('BP cross section','SG cross section')
grid on

% Figura 16 - Nominal Ozone
figure(16)
plot(sza_5,Omega_5,'+',sza_5,Omega_nominal,'+')
title('Ozone vs. Nominal Ozone');
xlabel('Solar Zenith Angle, deg')
ylabel('TOC, DU')
legend('Omega-med','Omega-nominal')
grid on


max(100*abs(Omega_5 - Omega_nominal)./Omega_5)
