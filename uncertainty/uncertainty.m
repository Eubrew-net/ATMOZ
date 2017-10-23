%% This script computes the Jacobian matrix for Brewer ozone observation
%% uncertainty on measured data irradiance.

clear all;

close all;

clc;

warning off;

disp('--------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code (v2.0)             ')
disp('          Dr. El Gawhary O., Dr. Parra-Rojas F.C. and Redondas A.         ')
disp('                               (2017)                                     ')
disp('--------------------------------------------------------------------------')

% Brewer number and date of measurement
brewer_str = input('Enter the number of the Brewer (between quotes): ');
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
		disp('   Lat:---, Lon=---, Alt:--- m')
	otherwise 
		disp('-----------------------')
		disp('There is no such Brewer')
		disp('-----------------------')
end

% input the date to study
yyyy=input('Enter the year: ');
mm=input('Enter the month: ');
dd=input('Enter the day: ');
period = (datenum(yyyy,mm,dd)); % matlab time
dates=datestr(period,'yyyy-mm-dd');

% database urls
url_base='fparra:m23275108M@rbcce.aemet.es/eubrewnet'; % base url
url_func='/data/get/O3L1'; % O3 Level 1.0 url
url_head='/data/get/BHeader'; % B Header url

% ----------------------------O3 Level 1 values------------------------------------------
for i=1:length(period)
	dates(i,:)
	url_str=[url_base,url_func,'?brewerid=',brewer_str,'&date=',dates(i,:)]
	curl_str=['curl -s --connect-timeout 120 "',url_str,'"'];
	[status,data_json]=system(curl_str);
	O3L1=loadjson(data_json);
end

for j = 2:length(O3L1)
	Z_deg(j-1) = cell2num(O3L1{1,j}(6)); % Solar Zenith Angle (deg)
	mu(j-1) = cell2num(O3L1{1,j}(7)); % ozone optical air mass
	Omega(j-1) = cell2num(O3L1{1,j}(10)); % ozone Level 1.0 (DU)
	t_j(j-1) = cell2num(O3L1{1,j}(5)); % % continuous date index based in Matlab datenum
	pre(j-1) = cell2num(O3L1{1,j}(16)); % pressure of the station (mbar)
	ms9(j-1) = cell2num(O3L1{1,j}(18)); % double ratio MS9 (counts/s)
	temp(j-1) = cell2num(O3L1{1,j}(8)); % temperarure of the station (Celsius)
	lon(j-1) = cell2num(O3L1{1,j}(15)); % Longitude of the observatory (deg)
	lat(j-1) = cell2num(O3L1{1,j}(14));	% Latitude of the observatory (deg)
end

% ----------------------------------B Header values-----------------------------------
dateFormat=16;

for i=1:length(period)
	dates(i,:)
	url_str=[url_base,url_head,'?brewerid=',brewer_str,'&date=',dates(i,:)]
	curl_str=['curl -s --connect-timeout 120 "',url_str,'"'];
	[status,data_json]=system(curl_str);
	BH=loadjson(data_json);
end

for k = 2:length(BH)
	A1(k) = cell2num(BH{1,k}(14))*10; % Ozone absorption coefficient (atm cm)^-1
	etc1(k) = cell2num(BH{1,k}(17)); % Extraterrestrial constant (counts/s)
	t_h(k) = cell2num(BH{1,k}(2)); % continuous date index based in Matlab datenum
end

% ozone absorption coefficients
A = A1(1,2)*ones(1,length(Z_deg)); % Ozone absorption coefficient for all the times (atm cm)^-1

% Extraterrestrial constant ratio
etc = etc1(1,2)*ones(1,length(Z_deg)); % Extraterrestrial constant for all the times (counts/s)

w = [0.0 -1.0 0.5 2.2 -1.7]; % statistical weights for earch wavelength
BE = [4807 4620 4410 4220 4040]; % Rayleight coefficients for each wavelength (atm^-1)
B = sum(w.*BE)*ones(1,length(Z_deg)); % weighted Rayleight coefficient for all the times (=1) (atm^-1) 

% ------------------ Definitions of some parameters ------------------------------------------------

P0 = 1013; % standard pressure (mbars)

R = 6371.229e3; % Earth's radius (m)

h = 5e3; % altitude of the Rayleigh disperssion layer (m)

r = 0; % altitude of the station in Huelva (m)

Z_rad = Z_deg*pi/180; % Solar Zenith Angle (rad)

m = sec(asin((R/(R+h))*sin(Z_rad))); % Rayleight optical air mass 

%-------------------------------- SZA uncertaiinty---------------------------------------------------
%                             through astronomical formulas
% Spencer, J.W. (1971) Fourier Series Representation of the Position of the Sun. Search, 2, 162-172.
%----------------------------------------------------------------------------------------------------

% from deg to rad
p0 = pi/180;

date_vec=datevec(t_j);
days=fix(t_j-datenum(1965,1,1));
t0 = date_vec(:,4)*60.0 + date_vec(:,5) + date_vec(:,6)/60.0;

% time uncertainty (min)
u_t0 = 1/60;

% eccentricity correction factor
t = (days+1)/365.2422;

ElipLong_1965 = 279.4574;

I = (ElipLong_1965 + 360*t + (t0'/1460.97))*p0;
u_I = (p0/1460.97)*u_t0;

% equation of time in seconds and the uncertainty
et = 4.2*sin(3*I)-2*cos(2*I)+596.5*sin(2*I)-12.8*sin(4*I)+19.3*cos(3*I)-(102.5+0.142*t).*sin(I)+(0.033*t-429.8).*cos(I);
u_et = u_I.*(3*4.2*cos(3*I)+4*sin(2*I)+2*596.5*cos(2*I)-4*12.8*cos(4*I)-3*19.3*sin(3*I)-(102.5+0.142*t).*cos(I)-(0.033*t-429.8).*sin(I));

% the hour angle and the uncertainty
ha = (t0'+(et/60)-720-4*lon)*p0/4;
u_ha = (p0/4)*u_t0+(p0/240)*u_et;

% solar declination and uncertainty
dec = atan(0.4336*sin(I-(p0*et)/240));
u_dec = 0.4336*u_I.*cos(I-(p0*et)/240)./(1+(0.4336*sin(I-(p0*et)/240)).^2)-(p0/240)*0.4336*u_et.*cos(I-(p0*et)/240)./(1+(0.4336*sin(I-(p0*et)/240)).^2);

% solar zenital angle and uncertainty
sza = acos(cos(ha).*cos(dec).*cos(lat*p0)+sin(lat*p0).*sin(dec))/p0;
u_sza = -sin(ha).*cos(dec).*cos(lat*p0).*u_ha./(p0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*p0)).^2))-(cos(ha).*sin(dec).*cos(lat*p0)+sin(lat*p0).*cos(dec)).*u_dec./(p0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*p0)).^2));
urel_sza = 100*u_sza./sza;

% Figure 1 - SZA relative uncertainty
figure(1)
plot(sza,urel_sza)
title('SZA Relative uncertainty');
xlabel('SZA, grad')
ylabel('Relative uncertainty, %')
grid on

%----------------------------Rayleigh air mass uncertainty----------------------------------------------
heff_r = 5e3; % Rayleigh effective altitude (m)
u_heff_r = 2e2; % uncertainty of the Rayleigh effective altitude (m)

% Rayleigh air mass and uncertainty
air_m = sec(asin((R/(R+h))*sin(Z_rad))); % air mass
u_am_r=(((R^2)*((sin(sza*p0)))*(R+heff_r))./((R^2)+(heff_r^2)+(R^2)*(sin(sza*p0)).^2).^(3/2)).*(u_sza+(((sin(sza*p0)).^2)/(R+heff_r)^2)*u_heff_r);
urel_am_r = 100*u_am_r./air_m;

% Figure 2 - Rayleigh air mass relative uncertainty 
figure(2)
plot(air_m,urel_am_r)
title('Rayleigh air mass relative uncertainty');
xlabel('Rayleigh air mass')
ylabel('Relative uncertainty, %')
grid on

% ---------------------------ozone airmass uncertainty----------------------------------------------
heff_o = 22e3; % Ozone effective altitude (m)
u_heff_o = 2e3; % Ozone effective altitude uncertainty (m)

% Ozone air mass and uncertainty
air_m_o = sec(asin((R/(R+heff_o))*sin(Z_rad))); % air mass
u_am_o=(((R^2)*((sin(sza*p0)))*(R+heff_o))./((R^2)+(heff_o^2)+(R^2)*(sin(sza*p0)).^2).^(3/2)).*sqrt((cos(sza*p0).*u_sza).^2+(((sin(sza*p0)/(R+heff_o))*u_heff_o).^2));
urel_am_o = 100*u_am_o./air_m_o;

% Figure 3 - Ozone air mass relative uncertainty
figure(3)
plot(air_m_o,urel_am_o)
title('Ozone air mass relative uncertainty');
xlabel('Ozone air mass')
ylabel('Relative uncertainty, %')
grid on

%----------------------------cross correlation matrix calculation----------------------------------------------
N_meas = (ms9-etc)+1;

% derivative of mu with respest to Z

dev_mu = (R+h)*(R+r)^2.*sin(Z_rad).*cos(Z_rad)./((R+h)^2-(R+r)^2.*(sin(Z_rad)).^2).^(3/2);
dev_m = sin(Z_rad)./cos(Z_rad).^2.*[(-1+0.0018167)+2*0.002875.*(sec(Z_rad)-1)+3*0.0008083.*(sec(Z_rad)-1).^2];

A_nominal = zeros(1,length(Z_rad));
Omega_nominal = zeros(1,length(Z_rad));
mu_nominal = zeros(1,length(Z_rad));
B_nominal = zeros(1,length(Z_rad));

for p=1:length(Z_rad)
                    %% Nominal values parameters
                    
                    fprintf('status computation = %2.1f  \n',(p/length(Z_deg))*100);

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

% Relative uncertainty of the measurement (percent/100)
urel_f = [0.01:0.005:0.04];

% Relative uncertainty of the ozone absorption (percent/100)
urel_A = [0.003:0.001:0.006];

% Fix relative uncertainty of the ozone absorption (percert/100)
urel_ff = 0.01;

% MS9 with no Rayleigh correction and uncertainty
meas = ms9-air_m.*B.*pre/1013; % counts/s
u_meas = urel_ff*meas; % counts/s

% Fix relative uncertainty of the ETC (counts/s)
u_etc_f = 5*ones(1,length(Z_deg));

% Relative uncertainty of the ETC (counts/s)
u_etc = [5:1:10];

% Fix relative uncertainty (percent/100) and uncertainty of the ozone absorption coefficient
urel_A_f = 0.003;
u_A = urel_A_f*A; % (atm cm)^-1

% Fix relative uncertainty (percent/100) and uncertainty of the Rayleigh coefficient
urel_B = 0.01;
u_B = urel_B*B; % (atm^-1)

%%
% Uncertainty of ozone with fixed uncertainties
y_omega = - mu.*A;
y_mu = -A.*Omega;
y_B = -air_m.*pre/1013;
y_A = -mu.*Omega;

a = y_omega.^2;
b = 2*y_omega.*(y_mu.*rho(:,2,3)'.*u_am_o + y_B.*rho(:,2,4)'.*u_B + y_A.*rho(:,1,2)'.*u_A);
c = ((y_mu.^2).*(u_am_o.^2) + (y_A.^2).*(u_A.^2) + (y_B.^2).*(u_B.^2) + 2*y_mu.*y_A.*rho(:,1,3)'.*u_am_o.*u_A +...
2*y_mu.*y_B.*rho(:,3,4)'.*u_am_o.*u_B + 2*y_B.*y_A.*rho(:,1,4)'.*u_B.*u_A - u_meas.^2 - u_etc_f.^2);

uo3 = (-b+sqrt(b.^2 - 4*a.*c))./(2*a);


% Figure 10 - Absolute ozone uncertainty
figure(10)
plot(Omega.*mu,uo3)
title('Ozone absolute uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
grid on

%%
% -----------------------uncertainty contriburtion------------------------------------

% Figure 11 - Contribution of each uncertainty to the total relative uncertainty 
figure(11)
plot(Omega.*mu,a,Omega.*mu,b,Omega.*mu,c,Omega.*mu,uo3)
title('Contribution of Ozone absolute uncertainty');
xlabel('Ozone, Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('u-a', 'u-b', 'u-c', 'u-tot')
grid on
%%

% ------- measurement uncertainty variation -------------------------------------------
%% 
for l = 1:length(urel_f)
    for i = 1:length(Z_deg)
        u_msf(l,i) = urel_f(l)*meas(i);
		u2_o3_f_1(l,i) = (-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A(i)) +...
		sqrt((-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A(i))).^2 -...
		4*(y_omega(i).^2).*((y_mu(i).^2).*(u_am_o(i).^2) + (y_A(i).^2).*(u_A(i).^2) + (y_B(i).^2).*(u_B(i).^2) + 2*y_mu(i).*y_A(i).*rho(i,1,3)'.*u_am_o(i).*u_A(i) +...
		2*y_mu(i).*y_B(i).*rho(i,3,4)'.*u_am_o(i).*u_B(i) + 2*y_B(i).*y_A(i).*rho(i,1,4)'.*u_B(i).*u_A(i) - u_msf(l,i).^2 - u_etc_f(i).^2)))./(2*y_omega(i).^2);
    end
end

% Figure 12 - Ozone absolute uncertainty vs. measurement uncertainty
figure(12)
plot(Omega.*mu,u2_o3_f_1(1:7,:))
title('Ozone absolute uncertainty vs. measurement uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('1%', '1.5%', '2%', '2.5%', '3%', '3.5%', '4%')
grid on
%%

% ------- etc uncertainty variation -------------------------------------------
%% 
for l = 1:length(u_etc)
    for i = 1:length(Z_deg)
        u2_o3_f_2(l,i) = (-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A(i)) +...
		sqrt((-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A(i))).^2 -...
		4*(y_omega(i).^2).*((y_mu(i).^2).*(u_am_o(i).^2) + (y_A(i).^2).*(u_A(i).^2) + (y_B(i).^2).*(u_B(i).^2) + 2*y_mu(i).*y_A(i).*rho(i,1,3)'.*u_am_o(i).*u_A(i) +...
		2*y_mu(i).*y_B(i).*rho(i,3,4)'.*u_am_o(i).*u_B(i) + 2*y_B(i).*y_A(i).*rho(i,1,4)'.*u_B(i).*u_A(i) - u_meas(i).^2 - u_etc(l).^2)))./(2*y_omega(i).^2);
    end
end

% Figure 13 - Ozone absolute uncertainty vs. ETC uncertainty
figure(13)
plot(Omega.*mu,u2_o3_f_2(1:6,:))
title('Ozone absolute uncertainty vs. etc uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('5u', '6u', '7u', '8u', '9u', '10u')
grid on
%%

% ------- ozone absorption uncertainty variation -------------------------------------------
%% 
for l = 1:length(urel_A)
    for i = 1:length(Z_deg)
		u_A_f(l,i) = urel_A(l)*A(i);
		u2_o3_f_3(l,i) = (-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A_f(l,i)) +...
		sqrt((-2*y_omega(i).*(y_mu(i).*rho(i,2,3)'.*u_am_o(i) + y_B(i).*rho(i,2,4)'.*u_B(i) + y_A(i).*rho(i,1,2)'.*u_A_f(l,i))).^2 -...
		4*(y_omega(i).^2).*((y_mu(i).^2).*(u_am_o(i).^2) + (y_A(i).^2).*(u_A_f(l,i).^2) + (y_B(i).^2).*(u_B(i).^2) + 2*y_mu(i).*y_A(i).*rho(i,1,3)'.*u_am_o(i).*u_A_f(l,i) +...
		2*y_mu(i).*y_B(i).*rho(i,3,4)'.*u_am_o(i).*u_B(i) + 2*y_B(i).*y_A(i).*rho(i,1,4)'.*u_B(i).*u_A_f(l,i) - u_meas(i).^2 - u_etc_f(i).^2)))./(2*y_omega(i).^2);
    end
end

% Figure 14 - Ozone absolute uncertainty vs. ozone absorption uncertainty
figure(14)
plot(Omega.*mu,u2_o3_f_3(1:4,:))
title('Ozone absolute uncertainty vs. ozone absorption uncertainty');
xlabel('Ozone Slant Column, DU')
ylabel('Absolute Uncertainty, DU')
legend('0.3%', '0.4%', '0.5%', '0.6%')
grid on
%%