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
