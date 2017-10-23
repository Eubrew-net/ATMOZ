%% This script computes the Jacobian matrix for Brewer ozone observation
%% uncertainty on measured data irradiance.

clear all;

close all;

clc;

warning off;

disp('--------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code (v6.0)             ')
disp('          Dr. El Gawhary O., Dr. Parra-Rojas F.C. and Redondas A.         ')
disp('                               (2017)                                     ')
disp('--------------------------------------------------------------------------')

% Brewer number
brewer_num = input('Enter the number of the Brewer: ');
brewer_str = num2str(brewer_num);
disp(' ')
switch brewer_str 
	case '214'
		disp('      -You are in Sodankyla-       ')
		disp('Lat:67.3675, Lon=-26.633, Alt:100.0 m')
	case '37'
		disp('      -You are in Sodankyla-       ')
		disp('Lat:67.3675, Lon=-26.633, Alt:100.0 m')
	case '6'
		disp('      -You are in Vindeln-       ')
		disp('Lat:64.244, Lon=-19.767, Alt:225.0 m')
	case '107'
		disp('      -You are in Jokioinen-       ')
		disp('Lat:60.814, Lon=-23.499, Alt:106.0 m')
	case '202'
		disp('      -You are in Sondrestrom-       ')
		disp('Lat:66.996, Lon=50.6214, Alt:150.0 m')
	case '53'
		disp('      -You are in Sondrestrom-       ')
		disp('Lat:66.996, Lon=50.6214, Alt:150.0 m')
	case '171'
		disp('      -You are in Fairbanks-       ')
		disp('Lat:64.819, Lon=147.869, Alt:138.0 m')
	case '128'
		disp('     -You are in Norrkoping-       ')
		disp('Lat:58.58, Lon=-16.15, Alt:43.0 m')
	case '228'
		disp('      -You are in Copenhagen-       ')
		disp('Lat:55.7185, Lon=-12.569, Alt:90.0 m')
	case '82'
		disp('      -You are in Copenhagen-       ')
		disp('Lat:55.7185, Lon=-12.569, Alt:90.0 m')
	case '44'
		disp('      -You are in Obninsk-       ')
		disp('Lat:55.099, Lon=-36.6066, Alt:100.0 m')
	case '43'
		disp('      -You are in Kislovodsk-       ')
		disp('Lat:43.733, Lon=-42.661, Alt:2070.0 m')
	case '188'
		disp('      -You are in Ankara-       ')
		disp('Lat:39.95, Lon=-32.88, Alt:913.0 m')
	case '86'
		disp('      -You are in Thessaloniki-      ')
		disp('Lat:40.634, Lon=-22.956, Alt:60.0 m')
	case '5'
		disp('      -You are in Thessaloniki-      ')
		disp('Lat:40.634, Lon=-22.956, Alt:60.0 m')
	case '207'
		disp('      -You are in Warsaw-       ')
		disp('Lat:52.246, Lon=-20.94, Alt:120.0 m')
	case '225'
		disp('   -You are in Poprad-Ganovce-       ')
		disp('Lat:49.03, Lon=-20.32, Alt:706.0 m')
	case '97'
		disp('   -You are in Poprad-Ganovce-       ')
		disp('Lat:49.03, Lon=-20.32, Alt:706.0 m')
	case '152'
		disp('  -You are in Budapest (Lorinc)-   ')
		disp('Lat:47.43, Lon=-19.18, Alt:139.0 m')
	case '184'
		disp('      -You are in Hradec Kralove-       ')
		disp('Lat:50.1772, Lon=-15.8386, Alt:285.0 m')
	case '98'
		disp('      -You are in Hradec Kralove-       ')
		disp('Lat:50.1772, Lon=-15.8386, Alt:285.0 m')
	case '118'
		disp('      -You are in Lidenberg-       ')
		disp('Lat:52.21, Lon=-14.12, Alt:112.0 m')
	case '78'
		disp('      -You are in Lidenberg-       ')
		disp('Lat:52.21, Lon=-14.12, Alt:112.0 m')
	case '30'
		disp('      -You are in Lidenberg-       ')
		disp('Lat:52.21, Lon=-14.12, Alt:112.0 m')
	case '67'
		disp('      -You are in Rome-       ')
		disp('Lat:41.9, Lon=-12.5, Alt:75.0 m')
	case '10'
		disp('  -You are in Hohenpeissenberg-       ')
		disp('Lat:47.8, Lon=-11.01, Alt:985.0 m')
	case '156'
		disp('       -You are in Arosa-       ')
		disp('Lat:46.78, Lon=-9.67, Alt:1840.0 m')
	case '163'
		disp('       -You are in Davos-       ')
		disp('Lat:46.8, Lon=-9.83, Alt:1560.0 m')
	case '72'
		disp('       -You are in Davos-       ')
		disp('Lat:46.8, Lon=-9.83, Alt:1560.0 m')
	case '66'
		disp('        -You are in Aosta-       ')
		disp('Lat:45.7422, Lon=-7.357, Alt:569.0 m')
	case '201'
		disp('      -You are in Tamanrasset-       ')
		disp('Lat:22.79, Lon=-5.529, Alt:1320.0 m')
	case '178'
		disp('        -You are in Uccle-       ')
		disp('Lat:50.799, Lon=-4.357, Alt:100.0 m')
	case '16'
		disp('        -You are in Uccle-       ')
		disp('Lat:50.799, Lon=-4.357, Alt:100.0 m')
	case '117'
		disp('      -You are in Murcia-       ')
		disp('Lat:38.0, Lon=1.17, Alt:69.0 m')
	case '166'
		disp('      -You are in Zaragoza-       ')
		disp('Lat:41.63, Lon=0.914, Alt:250.0 m')
	case '75'
		disp('      -You are in Reading-       ')
		disp('Lat:51.44, Lon=0.94, Alt:61.0 m')
	case '126'
		disp('   -You are in Manchester-       ')
		disp('Lat:53.47, Lon=2.23, Alt:76.0 m')
	case '172'
		disp('   -You are in Manchester-       ')
		disp('Lat:53.47, Lon=2.23, Alt:76.0 m')
	case '186'
		disp('  -You are in Madrid (Barajas)-   ')
		disp('Lat:40.45, Lon=3.72, Alt:680.0 m')
	case '70'
		disp('  -You are in Madrid (Barajas)-   ')
		disp('Lat:40.45, Lon=3.72, Alt:680.0 m')
	case '150'
		disp('  -You are in El Arenosillo-    ')
		disp('Lat:37.1, Lon=6.73, Alt:41.0 m')
	case '151'
		disp('     -You are in La Coruna-       ')
		disp('Lat:43.33, Lon=8.47, Alt:62.0 m')
	case '88'
		disp('   -You are in Valentia-       ')
		disp('Lat:51.938, Lon=10.25, Alt:14.0 m')
	case '227'
		disp('   -You are in Valentia-       ')
		disp('Lat:51.938, Lon=10.25, Alt:14.0 m')
	case '33'
		disp('   -You are in Santa Cruz de Tenerife-   ')
		disp('    Lat:28.47, Lon=16.25, Alt:52.0 m')
	case '157'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '183'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '185'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '48'
		disp('        -You are in Funchal-       ')
		disp('Lat:32.644, Lon=16.887, Alt:0.0 m')
	case '155'
		disp('        -You are in Carrasco International Airport-       ')
		disp('            Lat:-34.86, Lon=56.0, Alt:90.0 m')
	case '229'
		disp('        -You are in Rio Gallegos-       ')
		disp('     Lat:-51.6, Lon=69.32, Alt:5.0 m')
	case '180'
		disp('     -You are in Punta Arenas-       ')
		disp('Lat:-53.137, Lon=70.88, Alt:22.0 m')
	case '218'
		disp('       -You are in Izobamba-       ')
		disp('Lat:-0.366, Lon=78.55, Alt:3058.0 m')
	case '179'
		disp('        -You are in Hobart-       ')
		disp('Lat:-42.904, Lon=-147.327, Alt:20.0 m')
	case '232'
		disp('        -You are in Singapore-       ')
		disp('Lat:1.368, Lon=-103.982, Alt:14.0 m')
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
period = (datenum(yyyy,mm,dd)); % matlab time
dates=datestr(period,'yyyy-mm-dd');

disp(' ')

% we introduce the initial date
disp('Input the end date')
disp('------------------')
yyyy_f=input('Enter the year: ');
mm_f=input('Enter the month: ');
dd_f=input('Enter the day: ');
period_f = (datenum(yyyy_f,mm_f,dd_f)); % matlab time
dates_f=datestr(period_f,'yyyy-mm-dd');

% database urls
url_base='fparra:m23275108M@rbcce.aemet.es/eubrewnet'; % base url
url_func='/data/get/O3L1'; % ozone L1.0 url
url_head='/data/get/ConfigbyDate'; % Config by Date url
url_ds='/data/get/DS'; % Direct Sun url

% ----------------------------O3 Level 1 values------------------------------------------
for i=1:length(period)
	dates(i,:)
	url_str=[url_base,url_func,'?brewerid=',brewer_str,'&date=',dates(i,:),'&enddate=',dates_f(i,:)]
	curl_str=['curl -s --connect-timeout 120 "',url_str,'"'];
	[status,data_json]=system(curl_str);
	O3L1=loadjson(data_json);
end

for j = 2:length(O3L1)
	Omega(j-1) = cell2num(O3L1{1,j}(10)); % Calculated Ozone value with Standard algorithm + attenuation filter correction (DU)
	t_j(j-1) = cell2num(O3L1{1,j}(5)); % continuous date index based in Matlab datenum
	pre(j-1) = cell2num(O3L1{1,j}(16)); % Medium Pressure of the Brewer Location (mbar)
	temp(j-1) = cell2num(O3L1{1,j}(8)); % Instrument temperature (ºC)
	lon(j-1) = cell2num(O3L1{1,j}(15)); % Longitude of the Brewer Location (deg)
	lat(j-1) = cell2num(O3L1{1,j}(14)); % Latitude of the Brewer Location (deg)
	ms9_db(j-1) = cell2num(O3L1{1,j}(18)); % MS9, Second double ratio
	
	f3(j-1) = cell2num(O3L1{1,j}(23)); % corrected measurement of lambda3 (counts/s)
end

% ----------------------------------Config by Date values-----------------------------------
dateFormat=16;

for i=1:length(period)
	dates(i,:)
	url_str=[url_base,url_head,'?brewerid=',brewer_str,'&date=',dates(i,:),'&enddate=',dates_f(i,:)]
	curl_str=['curl -s --connect-timeout 120 "',url_str,'"'];
	[status,data_json]=system(curl_str);
	CbD=loadjson(data_json);
end

for k = 2:length(CbD)
	A1(k) = cell2num(CbD{1,k}(12))*10; % Ozone absorption coefficient (atm cm)^-1
	etc1(k) = cell2num(CbD{1,k}(15)); % extraterrestrial constant obtained by Langley extrapolations (185) or transfer (others)
	t_h(k) = cell2num(CbD{1,k}(3)); % date of the configuration
	dt(k) = cell2num(CbD{1,k}(17)); % deadtime of the PMT (ns)
	nfil0(k) = cell2num(CbD{1,k}(21)); % the attenuation value of the neutral-density filter (with no filter) 
	nfil1(k) = cell2num(CbD{1,k}(22)); % the attenuation value of the neutral-density filter 1
	nfil2(k) = cell2num(CbD{1,k}(23)); % the attenuation value of the neutral-density filter 2
	nfil3(k) = cell2num(CbD{1,k}(24)); % the attenuation value of the neutral-density filter 3
	nfil4(k) = cell2num(CbD{1,k}(25)); % the attenuation value of the neutral-density filter 4
	nfil5(k) = cell2num(CbD{1,k}(26)); % the attenuation value of the neutral-density filter 5
end

% ozone absorption coefficients in an array of the data size (atm cm)^-1
A = A1(1,2)*ones(1,length(t_j));

% Extraterrestrial constant ratio in an array of the data size (counts/s)
etc = etc1(1,2)*ones(1,length(t_j));

% Brewer weighting coefficients at the different wavelengths
w = [0.0 -1.0 0.5 2.2 -1.7];

% Rayleigh scattering coefficients at the different wavelenghts (atm^-1)
BE = [4870 4620 4410 4220 4040];

% Rayleigh scattering coefficients in an array of the data size (=1) (atm)^-1
B = sum(w.*BE)*ones(1,length(t_j));

% ----------------------------- Direct Sun ---------------------------------------------------------

for i=1:length(period)
	dates(i,:)
	url_str=[url_base,url_ds,'?brewerid=',brewer_str,'&date=',dates(i,:),'&enddate=',dates_f(i,:)]
	curl_str=['curl -s --connect-timeout 120 "',url_str,'"'];
	[status,data_json]=system(curl_str);
	DS=loadjson(data_json);
end

for k = 2:length(DS)
	ccl(k-1) = cell2num(DS{1,k}(9)); % number of the slit mask cycles
	dark(k-1) = cell2num(DS{1,k}(11)); % raw counts of the dark
	rc1(k-1) = cell2num(DS{1,k}(12)); % raw counts of the incident light at 306.6 nm
	rc2(k-1) = cell2num(DS{1,k}(13)); % raw counts of the incident light at 310.1 nm
	rc3(k-1) = cell2num(DS{1,k}(14)); % raw counts of the incident light at 313.5 nm
	rc4(k-1) = cell2num(DS{1,k}(15)); % raw counts of the incident light at 316.8 nm
	rc5(k-1) = cell2num(DS{1,k}(16)); % raw counts of the incident light at 320.1 nm
	nfilpos(k-1) = cell2num(DS{1,k}(4)); % position of the neutral-density filter (0, 64, 128, 192, 256, 320)
end

% ratio of the raw counts. It is only a test. It's no real!!
xx_r = -rc2 + 0.5*rc3 + 2.2*rc4 - 1.7*rc5;

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

% ratio of the counts. It is only a test. It's no real!!
xx_c = -f2_c + 0.5*f3_c + 2.2*f4_c - 1.7*f5_c;

% uncertainty of the ratio of the counts. It is only a test. It's no real!!
u2_xx_c = u2_2_c + (0.5^2)*u2_3_c + (2.2^2)*u2_4_c + (1.7^2)*u2_5_c;
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

% filter for small values of the counts
for i=1:length(t_j)
	if f1_d(i) < 2.0
		f1_d(i) = 2.0;
	elseif f2_d(i) < 2.0
		f2_d(i) = 2.0;
	elseif f3_d(i) < 2.0
		f3_d(i) = 2.0;
	elseif f4_d(i) < 2.0
		f4_d(i) = 2.0;
	elseif f5_d(i) < 2.0
		f5_d(i) = 2.0;
	end
end

% expanded uncertainty of the measurements including the dark corrention
u2_1_d = u2_1_c + u2_d;
u2_2_d = u2_2_c + u2_d;
u2_3_d = u2_3_c + u2_d;
u2_4_d = u2_4_c + u2_d;
u2_5_d = u2_5_c + u2_d;

% ratio of the counts. It is only a test. It's no real!!
xx_d = -f2_d + 0.5*f3_d + 2.2*f4_d - 1.7*f5_d;

% uncertainty of the ratio of the counts. It is only a test. It's no real!!
u2_xx_d = u2_2_d + (0.5^2)*u2_3_d + (2.2^2)*u2_4_d + (1.7^2)*u2_5_d;

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
		f1_ct(j)=f1_d(i).*exp(f1_ct(j-1)*dt(2));
		f2_ct(j)=f2_d(i).*exp(f2_ct(j-1)*dt(2));
		f3_ct(j)=f3_d(i).*exp(f3_ct(j-1)*dt(2));
		f4_ct(j)=f4_d(i).*exp(f4_ct(j-1)*dt(2));
		f5_ct(j)=f5_d(i).*exp(f5_ct(j-1)*dt(2));
		
		f1ct = f1_ct(j);
		f2ct = f2_ct(j);
		f3ct = f3_ct(j);
		f4ct = f4_ct(j);
		f5ct = f5_ct(j);
			
	end

% expanded uncertainty of the measures including the deadtime correction
	u2_f1_3(i) = ((exp(dt(2)*f1ct)./(1-dt(2)*f1ct)).^2).*u2_1_d(i) + (f1ct.^4)*udt^2;
	u2_f2_3(i) = ((exp(dt(2)*f2ct)./(1-dt(2)*f2ct)).^2).*u2_2_d(i) + (f2ct.^4)*udt^2;
	u2_f3_3(i) = ((exp(dt(2)*f3ct)./(1-dt(2)*f3ct)).^2).*u2_3_d(i) + (f3ct.^4)*udt^2;
	u2_f4_3(i) = ((exp(dt(2)*f4ct)./(1-dt(2)*f4ct)).^2).*u2_4_d(i) + (f4ct.^4)*udt^2;
	u2_f5_3(i) = ((exp(dt(2)*f5ct)./(1-dt(2)*f5ct)).^2).*u2_5_d(i) + (f5ct.^4)*udt^2;

% change to logarithmic space of the measurements	 
	f1_ct_1(i) = 10000*log10(f1ct);
	f2_ct_1(i) = 10000*log10(f2ct);
	f3_ct_1(i) = 10000*log10(f3ct);
	f4_ct_1(i) = 10000*log10(f4ct);
	f5_ct_1(i) = 10000*log10(f5ct);

% change to logarithmic space of the uncertainties (wrong)
	u2_f1dt(i) = 10000*log10(u2_f1_3(i));
	u2_f2dt(i) = 10000*log10(u2_f2_3(i));
	u2_f3dt(i) = 10000*log10(u2_f3_3(i));
	u2_f4dt(i) = 10000*log10(u2_f4_3(i));
	u2_f5dt(i) = 10000*log10(u2_f5_3(i));
	
end

% ratio of the counts. It is only a test. It's no real!!
xx_dt = -f2_ct_1 + 0.5*f3_ct_1 + 2.2*f4_ct_1 - 1.7*f5_ct_1;

% uncertainty of the ratio of the counts. It is only a test. It's no real!!
u2_xx_dt = u2_f2dt + (0.5^2)*u2_f3dt + (2.2^2)*u2_f4dt + (1.7^2)*u2_f5dt;


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
u2_f1t = u2_f1dt + (TC1*u_temp).^2 + (temp*utc1).^2;
u2_f2t = u2_f2dt + (TC2*u_temp).^2 + (temp*utc2).^2;
u2_f3t = u2_f3dt + (TC3*u_temp).^2 + (temp*utc3).^2;
u2_f4t = u2_f4dt + (TC4*u_temp).^2 + (temp*utc4).^2;
u2_f5t = u2_f5dt + (TC5*u_temp).^2 + (temp*utc5).^2;

% ratio of the counts. It is only a test. It's no real!!
xx_temp = -f2_temp + 0.5*f3_temp + 2.2*f4_temp - 1.7*f5_temp;

% uncertainty of the ratio of the counts. It is only a test. It's no real!!
u2_xx_temp = u2_f2t + (0.5^2)*u2_f3t + (2.2^2)*u2_f4t + (1.7^2)*u2_f5t;
%%
% ------------------------------ Neutral Filters correction -------------------------------------

% neutral filter corrections depending of the filter position
for i = 1:length(t_j)
	if nfilpos(i)./64 == 0
		f1_nd(i) = f1_temp(i) + nfil0(2);
		f2_nd(i) = f2_temp(i) + nfil0(2);
		f3_nd(i) = f3_temp(i) + nfil0(2);
		f4_nd(i) = f4_temp(i) + nfil0(2);
		f5_nd(i) = f5_temp(i) + nfil0(2);
	elseif nfilpos(i)./64 == 1
		f1_nd(i) = f1_temp(i) + nfil1(2);
		f2_nd(i) = f2_temp(i) + nfil1(2);
		f3_nd(i) = f3_temp(i) + nfil1(2);
		f4_nd(i) = f4_temp(i) + nfil1(2);
		f5_nd(i) = f5_temp(i) + nfil1(2);
	elseif nfilpos(i)./64 == 2
		f1_nd(i) = f1_temp(i) + nfil2(2);
		f2_nd(i) = f2_temp(i) + nfil2(2);
		f3_nd(i) = f3_temp(i) + nfil2(2);
		f4_nd(i) = f4_temp(i) + nfil2(2);
		f5_nd(i) = f5_temp(i) + nfil2(2);
	elseif nfilpos(i)./64 == 3
		f1_nd(i) = f1_temp(i) + nfil3(2);
		f2_nd(i) = f2_temp(i) + nfil3(2);
		f3_nd(i) = f3_temp(i) + nfil3(2);
		f4_nd(i) = f4_temp(i) + nfil3(2);
		f5_nd(i) = f5_temp(i) + nfil3(2);
	elseif nfilpos(i)./64 == 4
		f1_nd(i) = f1_temp(i) + nfil4(2);
		f2_nd(i) = f2_temp(i) + nfil4(2);
		f3_nd(i) = f3_temp(i) + nfil4(2);
		f4_nd(i) = f4_temp(i) + nfil4(2);
		f5_nd(i) = f5_temp(i) + nfil4(2);
	elseif nfilpos(i)./64 == 5
		f1_nd(i) = f1_temp(i) + nfil5(2);
		f2_nd(i) = f2_temp(i) + nfil5(2);
		f3_nd(i) = f3_temp(i) + nfil5(2);
		f4_nd(i) = f4_temp(i) + nfil5(2);
		f5_nd(i) = f5_temp(i) + nfil5(2);
	end
end

% ratio of the counts.
xx_nf = -f2_nd + 0.5*f3_nd + 2.2*f4_nd - 1.7*f5_nd;

% uncertainty of the measurements. We consider that the neutral filters uncertainty don't affect 
%u2_ms9_r = u2_xx_temp;
u2_xx_temp1 = -u2_f2t + (0.5)*u2_f3t + (2.2)*u2_f4t - (1.7)*u2_f5t; % it's not correct
u2_ms9_r = u2_xx_temp1;
%%
% ------------------ Definitions of some parameters ------------------------------------------------

P0 = 1013.25; % standard pressure (mbars)

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
m = sec(asin(sx_r));

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
u_heff_o = 2e3; % uncertainty of the ozone effective altitude (m)

% Ozone air mass and uncertainty
sx_o = (R/(R+heff_o))*sin(sza*p0);
mu = sec(asin(sx_o));

u2_mu=((sx_o.^2)/((1-sx_o.^2).^3)).*(((R^2)/((R+heff_o)^4)).*((sin(sza*p0)).^2)*(u_heff_o^2) + ((R^2)/((R+heff_o)^2)).*((cos(sza*p0)).^2).*(u2_sza));
urel_mu = 100*sqrt(u2_mu)./mu;

% Figure 3 - Ozone air mass relative uncertainty
figure(3)
plot(mu,urel_mu)
title('Ozone air mass relative uncertainty');
xlabel('Ozone air mass')
ylabel('Relative uncertainty, %')
grid on

%%
% ----------------------------- mean of five measurements -----------------------------------------------------

media=[];
standard=[];
sza_med=[];
mu_med=[];
m_med=[];
u2_xx_med=[];

for i=0:1:(size(xx_nf,2)/5)-1
    ser=mean(xx_nf(1,1+5*i:5+5*i)); % mean of MS9 (each 5 measurements)
	stan_m=std(xx_nf(1,1+5*i:5+5*i)); % uncertainty

	sza_m = mean(sza(1,1+5*i:5+5*i)); % mean of SZA (each 5 measurements)
	mu_m = median(mu(1,1+5*i:5+5*i)); % median of ozone optical mass
	m_m = median(m(1,1+5*i:5+5*i)); % median of Rayleigh optical mass 
	
	u2_xx_m = (1/5)*mean(u2_ms9_r(1,1+5*i:5+5*i)); % uncertainty of the mean of MS9
    
	media=[media;ser];
	standard=[standard;stan_m];
	sza_med=[sza_med;sza_m];
	mu_med=[mu_med;mu_m];
	m_med=[m_med;m_m];
	
	u2_xx_med=[u2_xx_med;u2_xx_m];
end

o3_5 = media';
u2_5 = (standard').^2;
sza_5 = sza_med';
mu_5 = mu_med';
m_5 = m_med';
u2_xx_5 = u2_xx_med';

MS9_r = xx_nf;


%----------------------------cross correlation matrix calculation (Omar)----------------------------------------------

% Measurements substracting the ETC
N_meas = (MS9_r-etc)+1; 


A_nominal = zeros(1,length(sza));
Omega_nominal = zeros(1,length(sza));
mu_nominal = zeros(1,length(sza));
B_nominal = zeros(1,length(sza));
%%
for p=1:length(sza)
                    %% Nominal values parameters
                    
                    fprintf('status computation = %2.1f  \n',(p/length(sza))*100);

                    A_nominal(p) = A(p);

                    Omega_nominal(p) = Omega(p);

                    mu_nominal(p) = mu(p);

                    B_nominal(p) = B(p);


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

                    [chi2_min, N] = generate_chi2(N_meas(p),sigma,A_nominal(p),Omega_nominal(p),mu_nominal(p), B_nominal(p));

                    [alpha,beta] = compute_jacobian(N, N_meas(p), A_nominal(p),Omega_nominal(p),mu_nominal(p), B_nominal(p),sigma);


                    %% Compute Jacobian and beta vector necessary for the iterative algorithm

                    % Iterations start from here

                    while (count~=2) 

                            % Compute the alpha_prime which is required in the algorithm
                            % basis in 5 dimensions

                            alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4);...
                                            alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4);...
                                            alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4);...
                                            alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l)];

                            % Then we compute the parameters variations

                            delta_param = alpha_prime\beta';

                            %Compute the new parameters

                            A_nominal_new(p) = A_nominal(p)  + delta_param(1);

                            Omega_nominal_new(p)  = Omega_nominal(p)  + delta_param(2);

                            mu_nominal_new(p)  = mu_nominal(p)  + delta_param(3);

                            B_nominal_new(p)  = B_nominal(p)  + delta_param(4);

                            % Compute the new chi2

                            [chi2_new, N] = generate_chi2(N_meas(p),sigma,A_nominal_new(p) ,Omega_nominal_new(p) ,mu_nominal_new(p) , B_nominal_new(p));

                            diff_chi2 = (chi2_new-chi2_min);

                            if (diff_chi2<0)

								l = l/10

                                A_nominal(p)  = A_nominal_new(p); 

                                Omega_nominal(p)  = Omega_nominal_new(p); 

                                mu_nominal(p)  = mu_nominal_new(p); 

                                B_nominal(p)  = B_nominal_new(p); 

                                [alpha, beta] = compute_jacobian(N, N_meas(p), A_nominal(p) ,Omega_nominal(p) ,mu_nominal(p) , B_nominal(p) , sigma);

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


                    alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4);...
                                            alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4);...
                                            alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4);...
                                            alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l)];
                    
                    
                    alpha_prime_corrected = alpha_prime+0.3.*eye(4,4);
                    
                    C  = inv(alpha_prime_corrected);
                    
                    uncert_A(p)  = sigma.*sqrt(abs(C(1,1)));

                    uncert_Omega(p)  = sigma.*sqrt(abs(C(2,2)));

                    uncert_mu(p)  = sigma.*sqrt(abs(C(3,3)));

                    uncert_B(p)  = sigma.*sqrt(abs(C(4,4)));
                     
                     

                    rho(p,1,1) = C(1,1)/(sqrt(abs(C(1,1)*C(1,1))));
                    rho(p,1,2) = C(1,2)/(sqrt(abs(C(1,1)*C(2,2))));
                    rho(p,1,3) = C(1,3)/(sqrt(abs(C(1,1)*C(3,3))));
                    rho(p,1,4) = C(1,4)/(sqrt(abs(C(1,1)*C(4,4))));
                    rho(p,2,1) = C(2,1)/(sqrt(abs(C(2,2)*C(1,1))));
                    rho(p,2,2) = C(2,2)/(sqrt(abs(C(2,2)*C(2,2))));
                    rho(p,2,3) = C(2,3)/(sqrt(abs(C(2,2)*C(3,3))));
                    rho(p,2,4) = C(1,1)/(sqrt(abs(C(1,1)*C(1,1))));
                    rho(p,3,1) = C(3,1)/(sqrt(abs(C(3,3)*C(1,1))));
                    rho(p,3,2) = C(3,2)/(sqrt(abs(C(3,3)*C(2,2))));
                    rho(p,3,3) = C(3,3)/(sqrt(abs(C(3,3)*C(3,3))));
                    rho(p,3,4) = C(3,4)/(sqrt(abs(C(3,3)*C(4,4))));
                    rho(p,4,1) = C(4,1)/(sqrt(abs(C(4,4)*C(1,1))));
                    rho(p,4,2) = C(4,2)/(sqrt(abs(C(4,4)*C(2,2))));
                    rho(p,4,3) = C(4,3)/(sqrt(abs(C(4,4)*C(3,3))));
                    rho(p,4,4) = C(4,4)/(sqrt(abs(C(4,4)*C(4,4))));

end

% Figure 4 - Cross Correlation between A and Omega
figure(4)
plot(t_j,rho(:,1,2));
datetick('x',dateFormat)
title('Correlation A-\Omega');
xlabel('Time')
ylabel('Correlation')

% Figure 5 - Cross Correlation between mu and Omega
figure(5)
plot(t_j,rho(:,2,3));
datetick('x',dateFormat)
title('Correlation Omega-\mu');
xlabel('Time')
ylabel('Correlation')

% Figure 6 - Cross Correlation between A and mu
figure(6)
plot(t_j,rho(:,1,3));
datetick('x',dateFormat)
title('Correlation \mu-A');
xlabel('Time')
ylabel('Correlation')

% Figure 7 - Cross Correlation between A and B
figure(7)
plot(t_j,rho(:,1,4));
datetick('x',dateFormat)
title('Correlation A-B');
xlabel('Time')
ylabel('Correlation')

% Figure 8 - Cross Correlation between B and Omega
figure(8)
plot(t_j,rho(:,2,4));
datetick('x',dateFormat)
title('Correlation B-\Omega');
xlabel('Time')
ylabel('Correlation')

% Figure 9 - Cross Correlation between B and mu
figure(9)
plot(t_j,rho(:,3,4));
datetick('x',dateFormat)
title('Correlation B-\mu');
xlabel('Time')
ylabel('Correlation')

%%
% MS9 with Rayleigh correction
meas = N_meas; % counts/s
o3_meas = (meas-m.*pre.*B/1013.15)./(A.*mu); % counts/s

% O3 without Rayleigh correction
o3_meas_r = meas./(A.*mu); % counts/s
%%

% uncertainty of the measurements
u_meas = sqrt(u2_ms9_r); % counts/s

% Uncertainty of the ETC (counts/s)
u_etc_f = 5*ones(1,length(t_j));

% Uncertainty of the ozone absorption coefficient (Bootstraping Method)

u_A = cross(brewer_str); % (atm cm)^-1

% Fix relative uncertainty (percent/100) and uncertainty of the Rayleigh coefficient
urel_B = 0.01;
u_B = urel_B*B; % (atm)^-1

%%
% Uncertainty of ozone with fixed uncertainties
a = (mu.*A).^2;
b = 2*mu.*A.*(Omega.*A.*sqrt(u2_mu).*rho(:,2,3)' + mu.*Omega.*u_A.*rho(:,1,2)' + m.*(pre/1013).*u_B.*rho(:,2,4)');
c = u2_mu.*(Omega.*A).^2 + (mu.*Omega.*u_A).^2 + (m.*(pre/1013).*u_B).^2 + 2*(Omega.^2).*A.*mu.*sqrt(u2_mu).*u_A.*rho(:,1,3)' +...
2*Omega.*A.*m.*(pre/1013).*sqrt(u2_mu).*u_B.*rho(:,3,4)' + 2*mu.*Omega.*m.*(pre/1013).*u_A.*u_B.*rho(:,1,4)' - u2_ms9_r - u_etc_f.^2;

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

% Figure 10 - Absolute ozone uncertainty
figure(10)
plot(Omega.*mu,uo3)
title('Ozone absolute uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
grid on

% ---------------- ozone uncertainty contributions ------------------------------------
% -------------------------------------------------------------------------------------
% ----------------- measurements ------------------------------------------------------

u_F_rel = [0.01:0.005:0.04];

for l = 1:7
    for i = 1:length(t_j)
        u_F_abs(l,i) = MS9_r(i)*u_F_rel(l);
		a_f(i) = (mu(i).*A(i)).^2;
		b_f(i) = 2*mu(i).*A(i).*(Omega(i).*A(i).*sqrt(u2_mu(i)).*rho(i,2,3)' + mu(i).*Omega(i).*u_A.*rho(i,1,2)' + m(i).*(pre(i)/1013).*u_B(i).*rho(i,2,4)');
		c_f(l,i) = u2_mu(i).*(Omega(i).*A(i)).^2 + (mu(i).*Omega(i).*u_A).^2 + (m(i).*(pre(i)/1013).*u_B(i)).^2 + 2*(Omega(i).^2).*A(i).*mu(i).*sqrt(u2_mu(i)).*u_A.*rho(i,1,3)' +...
		2*Omega(i).*A(i).*m(i).*(pre(i)/1013).*sqrt(u2_mu(i)).*u_B(i).*rho(i,3,4)' + 2*mu(i).*Omega(i).*m(i).*(pre(i)/1013).*u_A.*u_B(i).*rho(i,1,4)' - u_F_abs(l,i).^2 - u_etc_f(i).^2;
    
		uo3_f(l,i) = abs((-b_f(i)+sqrt(b_f(i).^2 - 4*a_f(i).*c_f(l,i)))./(2*a_f(i)));
	end
end

% Figure 11 - Ozone absolute uncertainty vs. measurement uncertainty
figure(11)
plot(Omega.*mu,uo3_f(1:7,:))
title('Ozone absolute uncertainty vs. measurement uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on

% ------------------------------- ozone absorption --------------------------------------------

u_A_rel = [0.01:0.005:0.04];

for l = 1:7
    for i = 1:length(t_j)
        u_A_abs(l,i) = A(i)*u_A_rel(l);
		a_a(i) = (mu(i).*A(i)).^2;
		b_a(l,i) = 2*mu(i).*A(i).*(Omega(i).*A(i).*sqrt(u2_mu(i)).*rho(i,2,3)' + mu(i).*Omega(i).*u_A_abs(l,i).*rho(i,1,2)' + m(i).*(pre(i)/1013).*u_B(i).*rho(i,2,4)');
		c_a(l,i) = u2_mu(i).*(Omega(i).*A(i)).^2 + (mu(i).*Omega(i).*u_A_abs(l,i)).^2 + (m(i).*(pre(i)/1013).*u_B(i)).^2 + 2*(Omega(i).^2).*A(i).*mu(i).*sqrt(u2_mu(i)).*u_A_abs(l,i).*rho(i,1,3)' +...
		2*Omega(i).*A(i).*m(i).*(pre(i)/1013).*sqrt(u2_mu(i)).*u_B(i).*rho(i,3,4)' + 2*mu(i).*Omega(i).*m(i).*(pre(i)/1013).*u_A_abs(l,i).*u_B(i).*rho(i,1,4)' - u2_ms9_r(i) - u_etc_f(i).^2;
    
		uo3_a(l,i) = abs((-b_a(l,i)+sqrt(b_a(l,i).^2 - 4*a_a(i).*c_a(l,i)))./(2*a_a(i)));
	end
end

% Figure 12 - Ozone absolute uncertainty vs. absorption uncertainty
figure(12)
plot(Omega.*mu,uo3_a(1:7,:))
title('Ozone absolute uncertainty vs. absorption uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on

% ------------------------------- Rayleigh coefficients --------------------------------------------

u_B_rel = [0.01:0.005:0.04];

for l = 1:7
    for i = 1:length(t_j)
        u_B_abs(l,i) = B(i)*u_B_rel(l);
		a_b(i) = (mu(i).*A(i)).^2;
		b_b(l,i) = 2*mu(i).*A(i).*(Omega(i).*A(i).*sqrt(u2_mu(i)).*rho(i,2,3)' + mu(i).*Omega(i).*u_A.*rho(i,1,2)' + m(i).*(pre(i)/1013).*u_B_abs(l,i).*rho(i,2,4)');
		c_b(l,i) = u2_mu(i).*(Omega(i).*A(i)).^2 + (mu(i).*Omega(i).*u_A).^2 + (m(i).*(pre(i)/1013).*u_B_abs(l,i)).^2 + 2*(Omega(i).^2).*A(i).*mu(i).*sqrt(u2_mu(i)).*u_A.*rho(i,1,3)' +...
		2*Omega(i).*A(i).*m(i).*(pre(i)/1013).*sqrt(u2_mu(i)).*u_B_abs(l,i).*rho(i,3,4)' + 2*mu(i).*Omega(i).*m(i).*(pre(i)/1013).*u_A.*u_B_abs(l,i).*rho(i,1,4)' - u2_ms9_r(i) - u_etc_f(i).^2;
    
		uo3_b(l,i) = abs((-b_b(l,i)+sqrt(b_b(l,i).^2 - 4*a_b(i).*c_b(l,i)))./(2*a_b(i)));
	end
end

% Figure 13 - Ozone absolute uncertainty vs. Rayleigh uncertainty
figure(13)
plot(Omega.*mu,uo3_b(1:7,:))
title('Ozone absolute uncertainty vs. Rayleigh uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on

% ------------------------------- ETC --------------------------------------------------

u_ETC_rel = [0.01:0.005:0.04];

for l = 1:7
    for i = 1:length(t_j)
        u_ETC_abs(l,i) = etc(i)*u_ETC_rel(l);
		a_e(i) = (mu(i).*A(i)).^2;
		b_e(i) = 2*mu(i).*A(i).*(Omega(i).*A(i).*sqrt(u2_mu(i)).*rho(i,2,3)' + mu(i).*Omega(i).*u_A.*rho(i,1,2)' + m(i).*(pre(i)/1013).*u_B(i).*rho(i,2,4)');
		c_e(l,i) = u2_mu(i).*(Omega(i).*A(i)).^2 + (mu(i).*Omega(i).*u_A).^2 + (m(i).*(pre(i)/1013).*u_B(i)).^2 + 2*(Omega(i).^2).*A(i).*mu(i).*sqrt(u2_mu(i)).*u_A.*rho(i,1,3)' +...
		2*Omega(i).*A(i).*m(i).*(pre(i)/1013).*sqrt(u2_mu(i)).*u_B(i).*rho(i,3,4)' + 2*mu(i).*Omega(i).*m(i).*(pre(i)/1013).*u_A.*u_B(i).*rho(i,1,4)' - u2_ms9_r(i)-u_ETC_abs(l,i).^2;
    
		uo3_e(l,i) = abs((-b_e(i)+sqrt(b_e(i).^2 - 4*a_e(i).*c_e(l,i)))./(2*a_e(i)));
	end
end

% Figure 14 - Ozone absolute uncertainty vs. ETC uncertainty
figure(14)
plot(Omega.*mu,uo3_e(1:7,:))
title('Ozone absolute uncertainty vs. ETC uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on

% ------------------------------- ozone air mass --------------------------------------------------

u_mu_rel = [0.01:0.005:0.04];

for l = 1:7
    for i = 1:length(t_j)
        u_mu_abs(l,i) = mu(i)*u_mu_rel(l);
		a_mu(i) = (mu(i).*A(i)).^2;
		b_mu(l,i) = 2*mu(i).*A(i).*(Omega(i).*A(i).*u_mu_abs(l,i).*rho(i,2,3)' + mu(i).*Omega(i).*u_A.*rho(i,1,2)' + m(i).*(pre(i)/1013).*u_B(i).*rho(i,2,4)');
		c_mu(l,i) = (u_mu_abs(l,i).^2).*(Omega(i).*A(i)).^2 + (mu(i).*Omega(i).*u_A).^2 + (m(i).*(pre(i)/1013).*u_B(i)).^2 + 2*(Omega(i).^2).*A(i).*mu(i).*u_mu_abs(l,i).*u_A.*rho(i,1,3)' +...
		2*Omega(i).*A(i).*m(i).*(pre(i)/1013).*u_mu_abs(l,i).*u_B(i).*rho(i,3,4)' + 2*mu(i).*Omega(i).*m(i).*(pre(i)/1013).*u_A.*u_B(i).*rho(i,1,4)' - u2_ms9_r(i) - u_etc_f(i).^2;
    
		uo3_mu(l,i) = abs((-b_mu(l,i)+sqrt(b_mu(l,i).^2 - 4*a_mu(i).*c_mu(l,i)))./(2*a_mu(i)));
	end
end

% Figure 15 - Ozone absolute uncertainty vs. ozone air mass uncertainty
figure(15)
plot(Omega.*mu,uo3_mu(1:7,:))
title('Ozone absolute uncertainty vs. ozone air mass uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on