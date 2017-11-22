%% This script computes the Jacobian matrix for Brewer ozone observation
%% uncertainty on measured data irradiance.

clear all;

close all;

clc;

warning off;

disp('-------------------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code for IZO16 Campaign             ')
disp('                     Triad Brewers and dates: 12/09/2016 -- 25/09/2016                     ')
disp('              Dr. Parra-Rojas F.C., Dr. El Gawhary O., Redondas A. and Egli L.       ')
disp('                                       (2017)                                     ')
disp('-------------------------------------------------------------------------------------')

disp('Press any key to continue')
pause;
disp(' ')

% Select the Brewer number
brewer_num = input('Enter the number of the Brewer: ');
brewer_str = num2str(brewer_num);
disp(' ')

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
	otherwise 
		disp('-----------------------')
		disp('There is no such Brewer')
		disp('-----------------------')
		break
end
disp(' ')


% Select the data base (1 -> Bass and Paur, 2 -> Serdyuchenko and Gorshelev)
DS_num = input('Enter the number of the data set (1 or 2): ');
disp(' ')

% load the data of the Brewer 185 only to calculate ETC transfer
ref = load('ref16.mat');

% load the input of the Brewer number
if brewer_str == '157'
	cal = load('br157.mat');
elseif brewer_str == '183'
	cal = load('br183.mat');
elseif brewer_str == '185'
	cal = load('br185.mat');
else
	disp('Wrong brewer number')
	break
end

% this is a temporal artifact to avoid negative values in O3
for i = 1:length(cal.t_j)
	if cal.Omega(i) < 0
		cal.Omega(i) = 0.1;
	end
end

% Format of the date to plot
dateFormat=16;

% ----------------------------------Config by Date values-----------------------------------

% for the Brewers 157, 183, 185 respectively

A1_b = [0.3395 0.341 0.342]*10; % ozone absorption coefficient for the triad (atm cm)^-1
etc_b = [1615 1630 1575]; % ETC for the triad (counts s^-1)
dt_b = [2.6 2.3 2.9]*1e-8; % dead time the triad (nsec)

if brewer_str == '157'
	A1 = A1_b(1);
	etc1 = etc_b(1);
	dt = dt_b(1);
elseif brewer_str == '183'
	A1 = A1_b(2);
	etc1 = etc_b(2);
	dt = dt_b(2);
elseif brewer_str == '185'
	A1 = A1_b(3);
	etc1 = etc_b(3);
	dt = dt_b(3);
end

% ozone absorption coefficients in an array of the data size
A = A1*ones(1,length(cal.t_j)/5);


% Extraterrestrial constant ratio in an array of the data size
etc = etc1*ones(1,length(cal.t_j)/5);

% Brewer weighting coefficients at the different wavelengths
w = [0.0 -1.0 0.5 2.2 -1.7];

% Rayleigh scattering coefficients at the different wavelenghts (atm^-1)

if DS_num == 1
	BE = [4870 4620 4410 4220 4040];
elseif DS_num == 2
	disp('still not included')
	break
end

% Rayleigh scattering coefficients in an array of the data size (=1)
B = sum(w.*BE)*ones(1,length(cal.t_j)/5);

% New variables
Omega = cal.Omega;
t = cal.t_j;
pre = cal.pre;
temp = cal.temp;
lon = cal.lon;
lat = cal.lat;
ccl = cal.ccl;
dark = cal.dark;
rc2 = cal.rc2;
rc3 = cal.rc3;
rc4 = cal.rc4;
rc5 = cal.rc5;
nfilpos = cal.nfilpos;

% --------------------raw to count correction----------------------------------------------------

it = 0.1147; % standard integration time of the Brewer's PMT

% raw to count correction for the dark current
dark_c = 2*dark./(it.*ccl);

% raw to count correction for the incident light for each wavelength
f2_c = 2*rc2./(it.*ccl); % counts/s
f3_c = 2*rc3./(it.*ccl); % counts/s
f4_c = 2*rc4./(it.*ccl); % counts/s
f5_c = 2*rc5./(it.*ccl); % counts/s

% uncertainty of the raw to count correction. Expanded (k=1) noise.
u2_2_c = (f2_c + dark_c)./(it.*ccl);
u2_3_c = (f3_c + dark_c)./(it.*ccl);
u2_4_c = (f4_c + dark_c)./(it.*ccl);
u2_5_c = (f5_c + dark_c)./(it.*ccl);

% ----------------- Dark count correction --------------------------------------------------------

% the uncertainty of the dark current is the standard deviation of this
u2_d = (std(dark_c)).^2;

% dark count correction for each wavelenght
f2_d = f2_c - dark_c;
f3_d = f3_c - dark_c;
f4_d = f4_c - dark_c;
f5_d = f5_c - dark_c;

% filter for small values of the counts
for i=1:length(t)
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

% expanded uncertainty of the measurements including the dark corrention
u2_2_d = u2_2_c + u2_d;
u2_3_d = u2_3_c + u2_d;
u2_4_d = u2_4_c + u2_d;
u2_5_d = u2_5_c + u2_d;

% -------------------------dead time correction---------------------------------------
%%
f2_ct=[];
f3_ct=[];
f4_ct=[];
f5_ct=[];

% uncertainty of the deadtime of the Brewer's PMT (Fountoulakis I. et al. 2016)
udt = 1e-9; % ns

% Deadtime correction for each wavelength (9 iterations)
for i = 1:length(t)
	f2_ct(1) = f2_d(i);
	f3_ct(1) = f3_d(i);
	f4_ct(1) = f4_d(i);
	f5_ct(1) = f5_d(i);
	for j=2:10
		f2_ct(j)=f2_d(i).*exp(f2_ct(j-1)*dt);
		f3_ct(j)=f3_d(i).*exp(f3_ct(j-1)*dt);
		f4_ct(j)=f4_d(i).*exp(f4_ct(j-1)*dt);
		f5_ct(j)=f5_d(i).*exp(f5_ct(j-1)*dt);
		
		f2ct(i) = f2_ct(j);
		f3ct(i) = f3_ct(j);
		f4ct(i) = f4_ct(j);
		f5ct(i) = f5_ct(j);
			
	end

% expanded uncertainty of the measures including the deadtime correction
	u2_f2_3(i) = ((exp(dt*f2ct(i))./(1-dt*f2ct(i))).^2).*u2_2_d(i) + (f2ct(i).^4)*udt^2;
	u2_f3_3(i) = ((exp(dt*f3ct(i))./(1-dt*f3ct(i))).^2).*u2_3_d(i) + (f3ct(i).^4)*udt^2;
	u2_f4_3(i) = ((exp(dt*f4ct(i))./(1-dt*f4ct(i))).^2).*u2_4_d(i) + (f4ct(i).^4)*udt^2;
	u2_f5_3(i) = ((exp(dt*f5ct(i))./(1-dt*f5ct(i))).^2).*u2_5_d(i) + (f5ct(i).^4)*udt^2;

% change to logarithmic space of the measurements	 
	f2_ct_1(i) = 10000*log10(f2ct(i));
	f3_ct_1(i) = 10000*log10(f3ct(i));
	f4_ct_1(i) = 10000*log10(f4ct(i));
	f5_ct_1(i) = 10000*log10(f5ct(i));
	
end

u2_f2log = (10000*(sqrt(u2_f2_3)./(f2ct*log(10)))).^2;
u2_f3log = (10000*(sqrt(u2_f3_3)./(f3ct*log(10)))).^2;
u2_f4log = (10000*(sqrt(u2_f4_3)./(f4ct*log(10)))).^2;
u2_f5log = (10000*(sqrt(u2_f5_3)./(f5ct*log(10)))).^2;

% -----------------------------Temperature correction--------------------------------------------

% uncertainty of the PMT temperature (the resolution of the PMT is 1ºC)
u_temp = 1/sqrt(3); 

% temperature coefficients and its uncertainties of the IZO2016 campaign
[tc,v,brw] = tblread('temp_izo16.csv',',');
brw_str = str2num(brewer_str);
brw_num = str2num(brw);

% temperature coefficients for each wavelenght
ind = find(brw_num==brw_str);
TC2 = tc(ind,2);
TC3 = tc(ind,3);
TC4 = tc(ind,4);
TC5 = tc(ind,5);

% uncertainty for each temperature coefficient
utc2 = tc(ind,7);
utc3 = tc(ind,8);
utc4 = tc(ind,9);
utc5 = tc(ind,10);

% PMT temperature correction
f2_temp = f2_ct_1 + temp.*TC2;
f3_temp = f3_ct_1 + temp.*TC3;
f4_temp = f4_ct_1 + temp.*TC4;
f5_temp = f5_ct_1 + temp.*TC5;

% expanded uncertainty of the measurements including the temperature correction
u2_f2t = u2_f2log + (TC2*u_temp).^2 + (temp*utc2).^2;
u2_f3t = u2_f3log + (TC3*u_temp).^2 + (temp*utc3).^2;
u2_f4t = u2_f4log + (TC4*u_temp).^2 + (temp*utc4).^2;
u2_f5t = u2_f5log + (TC5*u_temp).^2 + (temp*utc5).^2;

MS9_temp = -f2_temp + 0.5*f3_temp + 2.2*f4_temp - 1.7*f5_temp;

u_rel_ms9 = sqrt((sqrt(u2_f2t)./f2_temp).^2 + ((0.5)*sqrt(u2_f3t)./f3_temp).^2 + ((2.2)*sqrt(u2_f4t)./f4_temp).^2 + ((1.7)*sqrt(u2_f5t)./f5_temp).^2);

u2_ms9_temp = (u_rel_ms9.*MS9_temp).^2;

% ------------------------------ Neutral Filters correction -------------------------------------

% neutral filter corrections depending of the filter position

% Neutral filter coefficients and its uncertainties of the IZO2016 campaign
[nd,v,brw] = tblread('nd_izo16.csv',',');
brw_str = str2num(brewer_str);
brw_num = str2num(brw);

% neutral filter
ind = find(brw_num==brw_str);
NF0 = 0;
NF1 = nd(ind,6);
NF2 = nd(ind,7);
NF3 = nd(ind,8);
NF4 = nd(ind,9);
NF5 = nd(ind,10);

% uncertainty for each neutral filter
unf0 = 0;
unf1 = (nd(ind,11) + nd(ind,16))/2;
unf2 = (nd(ind,12) + nd(ind,17))/2;
unf3 = (nd(ind,13) + nd(ind,18))/2;
unf4 = (nd(ind,14) + nd(ind,19))/2;
unf5 = (nd(ind,15) + nd(ind,20))/2;

for i = 1:length(t)
	if nfilpos(i)./64 == 0
		MS9_nf(i) = MS9_temp(i) + NF0;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf0^2;
	elseif nfilpos(i)./64 == 1
		MS9_nf(i) = MS9_temp(i) + NF1;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf1^2;
	elseif nfilpos(i)./64 == 2
		MS9_nf(i) = MS9_temp(i) + NF2;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf2^2;
	elseif nfilpos(i)./64 == 3
		MS9_nf(i) = MS9_temp(i) + NF3;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf3^2;
	elseif nfilpos(i)./64 == 4
		MS9_nf(i) = MS9_temp(i) + NF4;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf4^2;
	elseif nfilpos(i)./64 == 5
		MS9_nf(i) = MS9_temp(i) + NF5;
		u2_ms9_nf(i) = u2_ms9_temp(i) + unf5^2;
	end
end

% MS9
xx_nf = MS9_nf;

% uncertainty of the measurements. The neutral filters uncertainty is not included
u2_ms9_r = u2_ms9_nf;

disp('Measurement uncertainty: OK')

% ------------------ Definitions of some parameters ------------------------------------------------

Pstan = 1013.25; % standard pressure (mbars)

R = 6371.229e3; % Earth's radius (m)

%-------------------------------- SZA uncertaiinty---------------------------------------------------
%                             through astronomical formulas
% Spencer, J.W. (1971) Fourier Series Representation of the Position of the Sun. Search, 2, 162-172.
%----------------------------------------------------------------------------------------------------
%%
% from deg to rad
p0 = pi/180;

% time of the measurements (min)
date_vec=datevec(t);
days=fix(t-datenum(1965,1,1));
t0 = date_vec(:,4)*60.0 + date_vec(:,5) + date_vec(:,6)/60.0;

% time uncertainty (min)
u_t0 = 1/60;

% eccentricity correction factor and unceretainty
t_e = (days+1)/365.2422;
ElipLong_1965 = 279.4574;
I = (ElipLong_1965 + 360*t_e + (t0'/1460.97))*p0;
u2_I = ((p0/1460.97)*u_t0)^2;

% equation of time in seconds and the uncertainty
et = 4.2*sin(3*I)-2*cos(2*I)+596.5*sin(2*I)-12.8*sin(4*I)+19.3*cos(3*I)-(102.5+0.142*t_e).*sin(I)+(0.033*t_e-429.8).*cos(I);
u2_et = u2_I.*(3*4.2*cos(3*I)+4*sin(2*I)+2*596.5*cos(2*I)-4*12.8*cos(4*I)-3*19.3*sin(3*I)-(102.5+0.142*t_e).*cos(I)-(0.033*t_e-429.8).*sin(I)).^2;

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

disp('Solar Zenith Angle uncertainty: OK')

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
if DS_num == 1
	u_heff_r = 2.0e2; % uncertainty of the Rayleigh effective altitude (m)
elseif DS_num == 2
	u_heff_r = 0.005*heff_r; % uncertainty of the Rayleigh effective altitude (m)
end

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

disp('Rayleigh air mass uncertainty: OK')

% ---------------------------ozone airmass uncertainty----------------------------------------------
heff_o = 22e3; % Ozone effective altitude (m)
%u_heff_o3 = 5e2; % uncertainty of the ozone effective altitude (m)

%u_heff_clim = 1e3; % uncertainty of the ozone effective altitude through climatology

if DS_num == 1
	u_heff_o = 2.0e3; % uncertainty of the Rayleigh effective altitude (m)
elseif DS_num == 2
	u_heff_o = 1.0e3; % uncertainty of the Rayleigh effective altitude (m)
end

% combined uncertainty
%u_heff_o = sqrt(u_heff_o3^2 + u_heff_clim^2);

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

disp('Ozone air mass uncertainty: OK')

% ----------------------------- wavelength error associated to the measurement -------------------------------
% obtained from dispresion

[wv,c,b_wv] = tblread('ms9_wv_test.csv',',');
brw_num_wv = str2num(b_wv);

ind_wv = find(brw_num_wv==brw_str);
k_wv = wv(ind,1);
u_k = wv(ind,2);

MS9_wv = xx_nf*(1+k_wv/100);
u2_wv = ((1+k_wv/100)^2).*u2_ms9_r + ((xx_nf/100).^2).*u_k^2;

disp('Wavelength error associated to the measurement: OK')

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
    ms9_m=mean(MS9_wv(1,1+5*i:5+5*i));
	stan_m=std(MS9_wv(1,1+5*i:5+5*i));
	u2_xx_m = (1/5)*mean(u2_wv(1,1+5*i:5+5*i));
	u_xx_std = std(u2_wv(1,1+5*i:5+5*i));
	media_ms9=[media_ms9;ms9_m];
	standard=[standard;stan_m];	
	u2_xx_med=[u2_xx_med;u2_xx_m];
	u_xx_sta=[u_xx_sta;u_xx_std];
	
	Omega_m=mean(Omega(1,1+5*i:5+5*i));
	Omega_med=[Omega_med;Omega_m];
	
	sza_m = mean(sza(1,1+5*i:5+5*i));
	
	mu_m = median(mu(1,1+5*i:5+5*i));
	mu_med=[mu_med;mu_m];
	u2_mu_m = median(u2_mu(1,1+5*i:5+5*i));
	u2_mu_med = [u2_mu_med;u2_mu_m];
		
	t_j_m = median(t(1,1+5*i:5+5*i));
	t_j_med = [t_j_med;t_j_m];
	
	m_m = median(m(1,1+5*i:5+5*i));
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

pre_5 = mean(pre)*ones(1,length(t)/5);

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

								l = l/10;

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

                                l = l*10;

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

disp('Cross correlations: OK')

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


% O3 with Rayleigh correction
meas = N_meas;
o3_meas = (meas-etc-m_5.*pre_5.*B/1013.25)./(A.*mu_5);

% O3 without Rayleigh correction
o3_meas_r = (meas-etc)./(A.*mu_5);
%%

% uncertainty of the pressure (mbar)
if DS_num == 1
	u_pre = 15.0; 
elseif DS_num == 2
	u_pre = 1.3;
end

% uncertainty of the measurements
u_meas = sqrt(u2_xx_5);

% ------------------------ Uncertainty of the ozone absorption coefficient -------------------------------

% Climatology (kelvin)
if DS_num == 1
	unc_teff = 0.0;
elseif DS_num == 2
	unc_teff = 1.0;
end

% percentage gradient for the BP and SG cross sections (Redondas A. et al., 2014)
grad_bp = 9.3601e-2;
grad_sg = 9.6391e-3;

% cross section uncertainty
[u_bp,u_sg]=sigma2_test(brewer_str);
disp('Cross section uncertainty: OK')

% Bootstraping method uncertainty
u_boot = cross(brewer_str);
disp('Integration method uncertainty: OK')

% variation of the ozone absorption coefficient with the wavelength obtained from dispersion (absolute)
Astep = 0.00109;

% uncertainty of the ozone absorption coefficient for each cross section
if DS_num == 1
	u_A = sqrt(u_bp^2 + u_boot^2 + (unc_teff*grad_bp*mean(A)/100)^2 + (Astep)^2);
	u_A_bp = u_A;
elseif DS_num == 2
	u_A = sqrt(u_sg^2 + u_boot^2 + (unc_teff*grad_sg*mean(A)/100)^2 + (Astep)^2);
	u_A_sg = u_A;
end
disp('Ozone absorption uncertainty: OK')

% ----------------------------------- Uncertainty of the ETC (counts/s) --------------------------------------
if brewer_num == 185
	u_etc = 10;
else
	% load the trensfer instrument data (t, O3, mu, SZA, MS9, sigma2MS9)
	a = [t_j_5' Omega_5' mu_5' sza_5' ms9_5' u2_xx_5'];

	% Load the reference Brewer data (t, O3, mu, SZA, MS9, Sigma2MS9, sigmaO3_BP, sigmaO3_SG)
	bc = load('ref16.mat');
	b = [bc.t_j_5' bc.Omega_5' bc.mu_5' bc.sza_5' bc.ms9_5' bc.u2_xx_5' bc.uo3' bc.uo3_sg'];

	% 5 minutes syncronization
	% dat = t(ins),t(inst)-t(ref),O3(inst),mu(inst),SZA(inst),MS9(inst),SigmasMS9(inst),O3(ref),mu(ref),SZA(ref),MS9(ref)
	% Sigma2MS9(ref),sigmaO3_BP(ref),sigmaO3_SG(ref)
	[x,r,rp,ra,dat,ox,osc_smooth_ini]=ratio_min_ozone(a,b,5,brewer_str,'185',0);

	t_etc = dat(:,2);

	diff_time = diff(dat(:,1));
	for i = 1:length(diff_time)
		if diff_time(i) == 0
			diff_time(i) = 0.01;
		end
	end

	diff_o3 = diff(dat(:,8))./diff_time; % O3 difference
	diff_o3 = [diff_o3(1),diff_o3'];

	% Figure 18 - ETC fit calculation
	figure(18)
	plot(dat(:,4),dat(:,6),'+')
	[h,Lr]=rline;
	etc_cal = Lr(2);

	% expanded uncertainties of ETC transfer
	if DS_num == 1
		u2_etc_bp = (A(1)*dat(:,4).*(dat(:,13)+t_etc.*diff_o3')).^2 + (dat(:,8).*dat(:,4)*u_A_bp).^2 + mean(u2_mu_5)*(dat(:,8).*A(1)).^2 + dat(:,7);
		u_etc = median(sqrt(u2_etc_bp))./sqrt(length(dat(:,4)));
	elseif DS_num == 2
		u2_etc_sg = (A(1)*dat(:,4).*(dat(:,14)+t_etc.*diff_o3')).^2 + (dat(:,8).*dat(:,4)*u_A_sg).^2 + mean(u2_mu_5)*(dat(:,8).*A(1)).^2 + dat(:,7);
		u_etc = mean(sqrt(u2_etc_sg));
	end
end
disp('ETC uncertainty: OK')

% --------------------------- Rayleigh coefficients uncertainty ----------------------------------------------
if DS_num == 1
	u_B = 0.01*B;
elseif DS_num == 2
	disp('still not included')
	break
end
disp('Rayleigh coefficient uncertainty: OK')

% -------------------------------- OVERALL UNCERTAINTY --------------------------------------------------------
% Ozone uncertainty with fixed uncertainties for BP cross section
a = (mu_5.*A).^2;

b = 2*mu_5.*A.*(Omega_5.*A.*sqrt(u2_mu_5).*rho(:,2,3)' + mu_5.*Omega_5.*u_A.*rho(:,1,2)' + (pre_5/Pstan).*B.*sqrt(u2_m_5).*rho(:,1,5)' +...  
    (m_5.*B/Pstan)*u_pre.*rho(:,2,6)' + m_5.*(pre_5/Pstan).*u_B.*rho(:,2,4)' - u_etc.*rho(:,2,7)');

c = u2_mu_5.*(Omega_5.*A).^2 + (mu_5.*Omega_5.*u_A).^2 + ((B.*pre_5/Pstan).^2).*u2_m_5 + (m_5.*(pre_5/Pstan).*u_B).^2 +...
    ((m_5.*B./Pstan).^2).*u_pre.^2 + u_etc.^2 + 2*(Omega_5.^2).*A.*mu_5.*sqrt(u2_mu_5).*u_A.*rho(:,1,3)' +...
	2*Omega_5.*A.*m_5.*(pre_5/Pstan).*sqrt(u2_mu_5).*u_B.*rho(:,3,4)' + 2*mu_5.*Omega_5.*m_5.*(pre_5/Pstan).*u_A.*u_B.*rho(:,1,4)'+...
    2*(Omega_5.*A.*pre_5.*B/Pstan).*sqrt(u2_mu_5).*sqrt(u2_m_5).*rho(:,3,5)' + 2*(Omega_5.*A.*m_5.*B/Pstan).*sqrt(u2_mu_5).*u_pre.*rho(:,3,6)'+...
    2*(mu_5.*Omega_5.*pre_5.*B./Pstan).*u_A.*sqrt(u2_m_5).*rho(:,1,5)' + 2*(mu_5.*Omega_5.*m_5.*B./Pstan).*u_A.*u_pre.*rho(:,1,6)' +...
    2*(pre_5.*B.*B.*m_5./(Pstan^2)).*sqrt(u2_m_5).*u_pre.*rho(:,5,6)' + 2*(B.*m_5.*pre_5.*pre_5/(Pstan^2)).*sqrt(u2_m_5).*u_B.*rho(:,4,5)' +...
    2*(B.*pre_5.*m_5.*m_5./(Pstan^2)).*u_pre.*u_B.*rho(:,4,6)' - 2*Omega_5.*A.*u_etc.*sqrt(u2_mu_5).*rho(:,3,7)'-...
    2*mu_5.*Omega_5.*u_etc.*u_A.*rho(:,1,7)' - 2*B.*(pre_5/Pstan).*u_etc.*sqrt(u2_m_5).*rho(:,5,7)'-...
    2*(m_5.*B/Pstan).*u_etc.*u_pre.*rho(:,6,7)' - 2*(m_5.*pre_5/Pstan).*u_etc.*u_B.*rho(:,4,7)' - u2_xx_5 ;

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

% Figure 19 - Ozone Absolute uncertainty
figure(19)
plot(Omega_5.*mu_5,uo3,'+')
set(gca,'Ylim',[0,20])
title('Ozone absolute uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
grid on

% Figure 20 - Ozone Relative uncertainty
figure(20)
plot(t_j_5,100.*uo3./Omega_5,'+')
datetick('x',dateFormat)
set(gca,'Ylim',[0,10])
title('Ozone relative uncertainty');
xlabel('Time')
ylabel('Relative Uncertainty, %')
grid on

% Figure 21 - Ozone Relative uncertainty versus BP and SG cross sections
figure(21)
plot(t_j_5,100.*2*uo3./Omega_5,'+')
datetick('x',dateFormat)
set(gca,'Ylim',[0,20])
title('Ozone relative uncertainty (2\sigma)');
xlabel('Time')
ylabel('Relative Uncertainty, %')
grid on

% Figura 22 - Nominal Ozone
figure(22)
plot(sza_5,Omega_5,'+',sza_5,Omega_nominal,'+')
title('Ozone vs. Nominal Ozone');
xlabel('Solar Zenith Angle, deg')
ylabel('TOC, DU')
legend('Omega-med','Omega-nominal')
grid on

%t_save = t_j_5';
%uo3_save = uo3';
%t0_save = t0';
%o3_save = Omega_5';
%sza_save = sza_5';
%ums9_save = (sqrt(u2_xx_5))';
%ms9_save = MS9_r';
%save('izo16_185.mat','t_save','t0_save','uo3_save','o3_save','sza_save','ms9_save','ums9_save')

disp(' ')
disp('OVERALL OZONE UNCERTAINTY: OK!!')
disp(' ')
disp('Thanks!!')