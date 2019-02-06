%% This script computes the Jacobian matrix for Brewer ozone observation
%% uncertainty on measured data irradiance.

clear all;

close all;

clc;

warning off;

disp('-------------------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code for IZO16 Campaign             ')
disp('                     Dobsons and dates: 12/09/2016 -- 25/09/2016                     ')
disp('              Dr. Parra-Rojas F.C., Dr. El Gawhary O., Redondas A. and Egli L.       ')
disp('                                       (2019)                                     ')
disp('-------------------------------------------------------------------------------------')

disp('Press any key to continue')
pause;
disp(' ')

DS = input('Enter the number of the data set (1 or 2): ');
disp(' ')

% Select the Brewer number
dobson_num = input('Enter the number of the Dobson: ');
dobson_str = num2str(dobson_num);
disp(' ')

read_dobson;
 d=Dobson2012;
 
 % modifica el tiempo
 d(:,2)=24*(d(:,1)-fix(d(:,1))); % GMT hour
 td=round(d(:,1)*24*60/5);
 
 % datos AOD
 load('aod.mat');
 taod=round(aod(:,1)*24*60/5);
 dobson=scan_join([td,d],[taod,aod(:,[1,8])]);
 dobson(1:find(~isnan(dobson(:,2)),1,'first')-1,:)=[];

 dobson=scan_join([td,d],[td,d]);
 
t = dobson(:,2);
m = dobson(:,23);
mu = dobson(:,22);
Omega = dobson(:,24);

% Izana
lat = 28.308;
lon = 16.499;

if DS == 1
	u_teff = 3.5;
	grad = 1.3285e-1;
	u_rel_beta = 0.01;
	u_p = 15;
	heff_o = 26-(lat/10);
	u_heff_o = 2.24e3;
elseif DS == 2
	u_teff = 1.0;
	grad = 1.0384e-1;
	u_rel_beta = 0.005;
	u_p = 1.3;
	heff_o = 24.24e3;
	u_heff_o = 0.5e3;
end

heff_r = 5e3;
u_heff_r = 2e2;
	

%% Initialization of variables

%% differential cross sections as given by Dobson manual
delta_alpha2 = 1.806;
delta_alpha1 = 0.374;


%% Rayleigh differential cross sections as given by Dobson manual

delta_beta2 = 0.114;
delta_beta1 = 0.104;


%% Definition of influence parameters.

R = 6371.229e3; % Earth's radius in m

h = 22e3; % height Ozone layer in m

r = 2360; % height station in m Izana

Z_rad= asin((R+h).*(sqrt(mu.^2-1))./((R+r).*mu));

Z_deg = Z_rad*180/pi;

% Z_deg = zenith_040_1; % zenith in degrees
% 
% Z_rad = Z_deg*pi/180; % zenith in rad

%Omega = ozone_040_1; % in atm-cm

%mu = air_m3_040_1; % Ozone airmass

P0 = 1013.25; % atmospheric pressure in mbars

%P = 7.6988; %0.830*1013.25; % station pressure in mbars
P = 770; %0.830*1013.25; % station pressure in mbars

%-------------------------------- SZA uncertaiinty---------------------------------------------------
%                             through astronomical formulas
% Spencer, J.W. (1971) Fourier Series Representation of the Position of the Sun. Search, 2, 162-172.
%----------------------------------------------------------------------------------------------------
%%
% from deg to rad
pp0 = pi/180;

% time of the measurements (min)
date_vec=datevec(t);
days=fix(t-datenum(1965,1,1));
t0 = date_vec(:,4)*60.0 + date_vec(:,5) + date_vec(:,6)/60.0;

% time uncertainty (min)
u_t0 = 1/60;

% eccentricity correction factor and unceretainty
t_e = (days+1)/365.2422;
ElipLong_1965 = 279.4574;
I = (ElipLong_1965 + 360*t_e + (t0./1460.97))*pp0;
u2_I = ((pp0/1460.97)*u_t0)^2;

% equation of time in seconds and the uncertainty
et = 4.2*sin(3*I)-2*cos(2*I)+596.5*sin(2*I)-12.8*sin(4*I)+19.3*cos(3*I)-(102.5+0.142*t_e).*sin(I)+(0.033*t_e-429.8).*cos(I);
u2_et = u2_I.*(3*4.2*cos(3*I)+4*sin(2*I)+2*596.5*cos(2*I)-4*12.8*cos(4*I)-3*19.3*sin(3*I)-(102.5+0.142*t_e).*cos(I)-(0.033*t_e-429.8).*sin(I)).^2;

% the hour angle and uncertainty
ha = (t0+(et/60)-720-4*lon)*pp0/4;
u2_ha = ((pp0/4)*u_t0)^2 +u2_et*(pp0/240)^2;

% solar declination and uncertainty
dec = atan(0.4336*sin(I-(pp0*et)/240));
u2_dec = (0.4336^2)*u2_I.*(cos(I-(pp0*et)/240)./(1+(0.4336*sin(I-(pp0*et)/240)).^2)).^2 + (((pp0/240)*0.4336)^2)*u2_et.*(cos(I-(pp0*et)/240)./(1+(0.4336*sin(I-(pp0*et)/240)).^2)).^2;

% solar zenital angle and uncertainty
sza = acos(cos(ha).*cos(dec).*cos(lat*pp0)+sin(lat*pp0).*sin(dec))/pp0;
u2_sza = ((sin(ha).*cos(dec).*cos(lat*pp0)).^2).*u2_ha./((pp0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*pp0)).^2)).^2)+((cos(ha).*sin(dec).*cos(lat*pp0)+sin(lat*pp0).*cos(dec)).^2).*u2_dec./((pp0*sqrt(1-(cos(ha).*cos(dec).*cos(lat*pp0)).^2)).^2);

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

sx_r = (R/(R+heff_r))*sin(sza*pp0);
u2_m = ((sx_r.^2)./((1-sx_r.^2).^3)).*(((R^2)/((R+heff_r)^4)).*((sin(sza*pp0)).^2).*(u_heff_r^2) + ((R^2)/((R+heff_r)^2)).*((cos(sza*pp0)).^2).*(u2_sza));
urel_m = 100*sqrt(u2_m)./m;

sx_o = (R/(R+heff_o))*sin(sza*pp0);
u2_mu=((sx_o.^2)./((1-sx_o.^2).^3)).*(((R^2)/((R+heff_o)^4)).*((sin(sza*pp0)).^2)*(u_heff_o^2) + ((R^2)/((R+heff_o)^2)).*((cos(sza*pp0)).^2).*(u2_sza));
urel_mu = 100*sqrt(u2_mu)./mu;


%m = air_mass_040_1; % airmass

delta_delta = 0;

% derivative of mu with respect to Z

dev_mu = -(R+h)*(R+r)^2.*sin(Z_rad).*cos(Z_rad)./((R+h)^2-(R+r)^2.*(sin(Z_rad)).^2).^(3/2);

dev_m = -sin(Z_rad)./(cos(Z_rad)).^2-0.0018167.*(-sin(Z_rad)./(cos(Z_rad)).^2)-0.002875.*(2.*(sec(Z_rad)-1).*(-sin(Z_rad)./(cos(Z_rad)).^2))-...
    0.0008083.*(3.*(sec(Z_rad)-1).^2).*(sec(Z_rad)-1).^2.*(-sin(Z_rad)./(cos(Z_rad)).^2);

dev_sec = sin(Z_rad)/(cos(Z_rad)).^2;
	
 

alpha_matrix = zeros(length(Z_rad),10,10);

%% Jacobian matrix. It deals with a 8x8 matrix

J11 = (mu.*Omega).^2;
J12 = -(mu.*Omega).^2;
J13 = mu.^2.*Omega.*(delta_alpha2-delta_alpha1);
J14 = (Omega.^2).*(delta_alpha2-delta_alpha1).*mu.*dev_mu;
J15 = mu.*Omega.*m.*P./P0;
J16 = -(mu.*Omega).*m.*P./P0;
J17 = mu.*Omega.*(delta_beta2-delta_beta1).*dev_m.*P./P0;
J18 = mu.*Omega.*(delta_beta2-delta_beta1).*m./P0;
J19 = mu.*Omega.*sec(Z_rad);
J110 = mu.*Omega.*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%----------------------------------------


J21 = J12;
J22 = (mu.*Omega).^2;
J23 = -mu.^2.*Omega.*(delta_alpha2-delta_alpha1);
J24 = -(Omega.^2).*(delta_alpha2-delta_alpha1).*mu.*dev_mu;
J25 = -mu.*Omega.*m.*P./P0;
J26 = mu.*Omega.*m.*P./P0;
J27 = -mu.*Omega.*(delta_beta2-delta_beta1).*dev_m.*P./P0;
J28 = -mu.*Omega.*(delta_beta2-delta_beta1).*m./P0;
J29 = -mu.*Omega.*sec(Z_rad);
J210 =  -mu.*Omega.*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%---------------------------------------------

% J31 = J13;
% J32 = J23;
% J33 = ((delta_alpha2-delta_alpha1).*mu).^2;
% J34 = -(m.*P./P0).^2;
% J35 = (m.*P./P0).^2;
% J36 = -m.*P./P0.*(delta_beta2-delta_beta1)./P0;
% J37 = -m.*P./P0;
% J38 = -m.*P./P0.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

J31 = J13;
J32 = J23;
J33 = ((delta_alpha2-delta_alpha1).*mu).^2;
J34 = ((delta_alpha2-delta_alpha1).^2).*mu.*Omega.*dev_mu;
J35 = (delta_alpha2-delta_alpha1).*mu.*m.*P/P0;
J36 = -(delta_alpha2-delta_alpha1).*mu.*m.*P/P0;
J37 = (delta_alpha2-delta_alpha1).*mu.*dev_m.*(delta_beta2-delta_beta1).*P/P0;
J38 = (delta_alpha2-delta_alpha1).*mu.*(delta_beta2-delta_beta1).*m./P0;
J39 = (delta_alpha2-delta_alpha1).*mu.*sec(Z_rad);
J310 = (delta_alpha2-delta_alpha1).*mu.*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%-----------------------------------------------

%-----------------------------------------------

J41 = J14;
J42 = J24;
J43 = J34;
J44 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).^2;
J45 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*(m.*P./P0);
J46 = -(Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*(m.*P./P0);
J47 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*(delta_beta2-delta_beta1).*dev_m.*P/P0;
J48 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*(delta_beta2-delta_beta1).*m./P0;
J49 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*sec(Z_rad);
J410 = (Omega.*(delta_alpha2-delta_alpha1).*dev_mu).*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%------------------------------------------------

J51 = J15;
J52 = J25;
J53 = J35;
J54 = J45;
J55 = (m.*P./P0).^2;
J56 = -(m.*P./P0).^2;
J57 = m.*P./P0.*(delta_beta2-delta_beta1).*dev_m.*P/P0;
J58 = m.*P./P0.*(delta_beta2-delta_beta1).*m./P0;
J59 = m.*P./P0.*sec(Z_rad);
J510 = m.*P./P0.*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%--------------------------------------------------


J61 = J16;
J62 = J26;
J63 = J36;
J64 = J46;
J65 = J56;
J66 = -(m.*P/P0).^2;
J67 = -m.*((P/P0).^2).*(delta_beta2-delta_beta1).*dev_m;
J68 = -m.*(P/P0).*(delta_beta2-delta_beta1).*m./P0;
J69 = -m.*(P/P0).*sec(Z_rad);
J610 = -m.*(P/P0).*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%------------------------------------------------

J71 = J17;
J72 = J27;
J73 = J37;
J74 = J47;
J75 = J57;
J76 = J67;
J77 = ((delta_beta2-delta_beta1).*dev_m.*P/P0).^2;
J78 = ((delta_beta2-delta_beta1).^2).*m.*dev_m.*P/(P0.^2);
J79 = ((delta_beta2-delta_beta1).*dev_m.*P/P0).*sec(Z_rad);
J710 = ((delta_beta2-delta_beta1).*dev_m.*P/P0).*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

%-----------------------------------------------

J81 = J18;
J82 = J28;
J83 = J38;
J84 = J48;
J85 = J58;
J86 = J68;
J87 = J78;
J88 = ((delta_beta2-delta_beta1).*m./P0).^2;
J89 = ((delta_beta2-delta_beta1).*m./P0).*sec(Z_rad);
J810 = ((delta_beta2-delta_beta1).*m./P0).*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

J91 = J19;
J92 = J29;
J93 = J39;
J94 = J49;
J95 = J59;
J96 = J69;
J97 = J79;
J98 = J89;
J99 = (sec(Z_rad)).^2;
J910 = (sec(Z_rad)).*delta_delta.*(sin(Z_rad)./cos(Z_rad).^2);

J101 = J110;
J102 = J210;
J103 = J310;
J104 = J410;
J105 = J510;
J106 = J610;
J107 = J710;
J108 = J810;
J109 = J910
J1010 = (delta_delta.*(sin(Z_rad)./cos(Z_rad))).^2;



%------------------------------------------------

for i=1:length(Z_rad)
    
    
    
    fprintf('status computation = %2.1f  \n',(i/length(Z_rad))*100);

    alpha_matrix= [squeeze(J11(i)) squeeze(J12(i)) squeeze(J13(i)) squeeze(J14(i)) squeeze(J15(i)) squeeze(J16(i)) squeeze(J17(i)) squeeze(J18(i)) squeeze(J19(i)) squeeze(J110(i));
            squeeze(J21(i)) squeeze(J22(i)) squeeze(J23(i)) squeeze(J24(i)) squeeze(J25(i)) squeeze(J26(i)) squeeze(J27(i)) squeeze(J28(i)) squeeze(J29(i)) squeeze(J210(i));
            squeeze(J31(i)) squeeze(J32(i)) squeeze(J33(i)) squeeze(J34(i)) squeeze(J35(i)) squeeze(J36(i)) squeeze(J37(i)) squeeze(J38(i)) squeeze(J39(i)) squeeze(J310(i));
            squeeze(J41(i)) squeeze(J42(i)) squeeze(J43(i)) squeeze(J44(i)) squeeze(J45(i)) squeeze(J46(i)) squeeze(J47(i)) squeeze(J48(i)) squeeze(J49(i)) squeeze(J410(i));
            squeeze(J51(i)) squeeze(J52(i)) squeeze(J53(i)) squeeze(J54(i)) squeeze(J55(i)) squeeze(J56(i)) squeeze(J57(i)) squeeze(J58(i)) squeeze(J59(i)) squeeze(J510(i));
            squeeze(J61(i)) squeeze(J62(i)) squeeze(J63(i)) squeeze(J64(i)) squeeze(J65(i)) squeeze(J66(i)) squeeze(J67(i)) squeeze(J68(i)) squeeze(J69(i)) squeeze(J610(i));
            squeeze(J71(i)) squeeze(J72(i)) squeeze(J73(i)) squeeze(J74(i)) squeeze(J75(i)) squeeze(J76(i)) squeeze(J77(i)) squeeze(J78(i)) squeeze(J79(i)) squeeze(J710(i));
            squeeze(J81(i)) squeeze(J82(i)) squeeze(J83(i)) squeeze(J84(i)) squeeze(J85(i)) squeeze(J86(i)) squeeze(J87(i)) squeeze(J88(i)) squeeze(J89(i)) squeeze(J810(i));
			squeeze(J91(i)) squeeze(J92(i)) squeeze(J93(i)) squeeze(J94(i)) squeeze(J95(i)) squeeze(J96(i)) squeeze(J97(i)) squeeze(J98(i)) squeeze(J99(i)) squeeze(J910(i));
			squeeze(J101(i)) squeeze(J102(i)) squeeze(J103(i)) squeeze(J104(i)) squeeze(J105(i)) squeeze(J106(i)) squeeze(J107(i)) squeeze(J108(i)) squeeze(J109(i)) squeeze(J1010(i))];
        
        
        

%--------------------------------------------------------

             %% Determination of the correlation matrix
             
             
            alpha_prime_corrected = alpha_matrix+0.2*eye(10,10);
                    
            C  = inv(alpha_prime_corrected);
            

            for j=1:10

                for k=1:10

                    rho(i,j,k) = C(j,k)/(sqrt(C(j,j)*C(k,k)));


                end
            end


        
              
end  


    %% Now we compute the total uncertainty based on the input matrix provided by Ulf Kohler
    %% Relative uncertainties

	
    u_rel_wv = 0.01;
    u_rel_xs = 0.01;
    u_rel_temp = 0.01;

	u_rel_alpha1 = sqrt((u_teff.*grad.*delta_alpha1./100).^2 + u_rel_wv.^2 + u_rel_xs.^2 + u_rel_temp.^2);
	u_rel_alpha2 = sqrt((u_teff.*grad.*delta_alpha2./100).^2 + u_rel_wv.^2 + u_rel_xs.^2 + u_rel_temp.^2);

    u_rel_betest = 0.0001;
	
	u_rel_beta1 = sqrt(u_rel_beta.^2 + u_rel_betest.^2);
	u_rel_beta2 = u_rel_beta1;
    
    u_rel_P = u_p./P;
    
    u_rel_delta = 0.0002;
	
	u_rel_sec = (sin(Z_rad)./((cos(Z_rad)).^2)).*sqrt(u2_sza)./sec(Z_rad);

%    u_rel_etc = sqrt(0.01^2+0.005^2+0.005^2);

    u_rel_mu = urel_mu./100;
	
	u_rel_m = urel_m./100;

    u_rel_wedge = 0.005;

    u_rel_straylight = sqrt(0.01^2+0.001^2);
    
    u_rel_y = u_rel_wedge;
    
    
    u_rel_interference_instruments = 0.001;
    
    u_rel_operat_impre = 0.01;
    
	
    dy_da1 = mu.*Omega;
    
    dy_da2 = -mu.*Omega;
    
    dy_da3 = ((delta_alpha2-delta_alpha1).*mu);
	
	dy_da4 = Omega.*(delta_alpha2-delta_alpha1);
    
    dy_da5 = m.*P./P0;
    
    dy_da6 = -m.*P./P0;
	
	dy_da7 = (delta_beta2-delta_beta1).*P/P0;
    
    dy_da8 = (delta_beta2-delta_beta1).*m./P0;
    
    dy_da9 = sec(Z_rad);
    
    dy_da10 = delta_delta;
    
    a = dy_da3.^2;

    b = 2.*(dy_da1.*dy_da3.*u_rel_alpha2.*rho(:,1,3)+dy_da2.*dy_da3.*u_rel_alpha1.*rho(:,2,3)+dy_da4.*dy_da3.*u_rel_mu.*rho(:,4,3)+...
			dy_da5.*dy_da3.*u_rel_beta2.*rho(:,5,3)+dy_da6.*dy_da3.*u_rel_beta1.*rho(:,6,3)+dy_da7.*dy_da3.*u_rel_m.*rho(:,7,3)+...
			dy_da8.*dy_da3.*u_rel_P.*rho(:,8,3)+dy_da9.*dy_da3.*u_rel_delta.*rho(:,9,3)+dy_da10.*dy_da3.*u_rel_sec.*rho(:,10,3));

	c = ((dy_da1.*u_rel_alpha2).^2)+((dy_da2.*u_rel_alpha1).^2) + ((dy_da4.*u_rel_mu).^2) + ((dy_da5.*u_rel_beta2).^2) + ((dy_da6.*u_rel_beta1).^2)+...
		((dy_da7.*u_rel_m).^2) + ((dy_da8.*u_rel_P).^2) + ((dy_da9.*u_rel_delta).^2) + ((dy_da10.*u_rel_sec).^2)+...
		2.*(dy_da1.*dy_da2.*u_rel_alpha2.*u_rel_alpha1.*rho(:,1,2) + dy_da1.*dy_da4.*u_rel_alpha2.*u_rel_mu.*rho(:,1,4)+...
		dy_da1.*dy_da5.*u_rel_alpha2.*u_rel_beta2.*rho(:,1,5) + dy_da1.*dy_da6.*u_rel_alpha2.*u_rel_beta1.*rho(:,1,6)+...
		dy_da1.*dy_da7.*u_rel_alpha2.*u_rel_m.*rho(:,1,7) + dy_da1.*dy_da8.*u_rel_alpha2.*u_rel_P.*rho(:,1,8)+...
		dy_da1.*dy_da9.*u_rel_alpha2.*u_rel_delta.*rho(:,1,9) + dy_da1.*dy_da10.*u_rel_alpha2.*u_rel_sec.*rho(:,1,10)+...
		dy_da2.*dy_da4.*u_rel_alpha1.*u_rel_mu.*rho(:,2,4) + dy_da2.*dy_da5.*u_rel_alpha1.*u_rel_beta2.*rho(:,2,5)+...
		dy_da2.*dy_da6.*u_rel_alpha1.*u_rel_beta1.*rho(:,2,6) + dy_da2.*dy_da7.*u_rel_alpha1.*u_rel_m.*rho(:,2,7)+...
		dy_da2.*dy_da8.*u_rel_alpha1.*u_rel_P.*rho(:,2,8) + dy_da2.*dy_da9.*u_rel_alpha1.*u_rel_delta.*rho(:,2,9)+...
		dy_da2.*dy_da10.*u_rel_alpha1.*u_rel_sec.*rho(:,2,10)+... 
		dy_da4.*dy_da5.*u_rel_mu.*u_rel_beta2.*rho(:,4,5) + dy_da4.*dy_da6.*u_rel_mu.*u_rel_beta1.*rho(:,4,6)+...
		dy_da4.*dy_da7.*u_rel_mu.*u_rel_m.*rho(:,4,7) + dy_da4.*dy_da8.*u_rel_mu.*u_rel_P.*rho(:,4,8)+...
		dy_da4.*dy_da9.*u_rel_mu.*u_rel_delta.*rho(:,4,9) + dy_da4.*dy_da10.*u_rel_mu.*u_rel_sec.*rho(:,4,10)+...
		dy_da5.*dy_da6.*u_rel_beta2.*u_rel_beta1.*rho(:,5,6) + dy_da5.*dy_da7.*u_rel_beta2.*u_rel_m.*rho(:,5,7)+...
		dy_da5.*dy_da8.*u_rel_beta2.*u_rel_P.*rho(:,5,8) + dy_da5.*dy_da9.*u_rel_beta2.*u_rel_delta.*rho(:,5,9)+...
		dy_da5.*dy_da10.*u_rel_beta2.*u_rel_sec.*rho(:,5,10)+...
		dy_da6.*dy_da7.*u_rel_beta1.*u_rel_m.*rho(:,6,7) + dy_da6.*dy_da8.*u_rel_beta1.*u_rel_P.*rho(:,6,8)+...
		dy_da6.*dy_da9.*u_rel_beta1.*u_rel_delta.*rho(:,6,9) + dy_da6.*dy_da10.*u_rel_beta1.*u_rel_sec.*rho(:,6,10)+...
		dy_da7.*dy_da8.*u_rel_m.*u_rel_P.*rho(:,7,8) + dy_da7.*dy_da9.*u_rel_m.*u_rel_delta.*rho(:,7,9)+...
		dy_da7.*dy_da10.*u_rel_m.*u_rel_sec.*rho(:,7,10)+...
		dy_da8.*dy_da9.*u_rel_P.*u_rel_delta.*rho(:,8,9) + dy_da8.*dy_da10.*u_rel_P.*u_rel_sec.*rho(:,8,10)+...
		dy_da9.*dy_da10.*u_rel_delta.*u_rel_sec.*rho(:,9,10)) - 2.*(u_rel_y.^2); 
	
	
	
        
   u_toc_final_1 = abs((-b+sqrt(b.^2-4.*a.*c))./(2.*a));
   u_toc_final_2 = abs((-b-sqrt(b.^2-4.*a.*c))./(2.*a));
   
   u_toc_final__combined = sqrt(u_toc_final_2.^2+u_rel_interference_instruments.^2+u_rel_operat_impre.^2)  ; 
 

%% outputs


figure(20)
plot(Z_deg,100.*u_toc_final__combined./Omega,'+');
ylabel('relative TOC uncertainty');
xlabel('Zenith angle')

figure(21)
plot(Z_deg,Omega);
hold on;
plot(Z_deg,u_toc_final__combined,'r');
xlabel('Zenith angle');




%% Generate output plots

%%Degrees of correlation
% figure(1)
% plot(Z_deg,rho(:,1,2))
% 
% figure(2)
% plot(Z_deg,rho(:,1,3))
% 
% figure(3)
% plot(Z_deg,rho(:,1,4))
% 
% figure(4)
% plot(Z_deg,rho(:,1,5))
% 
% figure(5)
% plot(Z_deg,rho(:,1,6))
% 
% figure(6)
% plot(Z_deg,rho(:,1,7))
% 
% figure(7)
% plot(Z_deg,rho(:,1,8))
% 
% figure(8)
% plot(Z_deg,rho(:,2,3))
% 
% figure(9)
% plot(Z_deg,rho(:,2,4))
% 
% figure(10)
% plot(Z_deg,rho(:,5,8))