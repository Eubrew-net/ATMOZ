clear all;

close all;

clc;

warning off;

disp('-------------------------------------------------------------------------------------')
disp('                 Welcome to the ATMOZ uncertainty code for IZO16 Campaign            ')
disp('                     Dobsons and dates: 12/09/2016 -- 25/09/2016                     ')
disp('                 Parra-Rojas F.C., El Gawhary O., Redondas A., Kohler U.,            ')
disp('                              Egli L. and Grobner J.                                 ')
disp('                                       (2019)                                        ')
disp('-------------------------------------------------------------------------------------')

disp('Press any key to continue')
pause;
disp(' ')
 
 % lee los datos Dobson
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

 %% try to reproduce dobson calculations
PAIR={'C','D','A','CD','AD'};

 
 %    C      D     A     CD    AD
 A1=[0.833,0.374,1.806,0.459,1.432];
 B1=[0.109,0.104,0.114,0.005,0.010];
 
 % calculamos ozono
 i=0:2;
 N=dobson(:,7+i*5);N=[N,N(:,1)-N(:,2),N(:,3)-N(:,2)]; 
 mu_dob=dobson(:,4+i*5);mu_dob=[mu_dob,dobson(:,[19,22])];
 m_dob=dobson(:,5+i*5);m_dob=[m_dob,dobson(:,[20,23])];
 xf=dobson(:,8+5*i);xf=[xf,dobson(:,[21,24])]; 
 
 RC=matmul(B1,m_dob)*770/1013;
 
 o3d= matmul(A1,mu_dob);
 x=(N/100-RC)./o3d;
 
 figure(1);
 plot(dobson(:,2),100*abs(xf(:,4)-x(:,4))/xf(:,4));
 title('Relative uncertainty of Xcd');
 
 figure(2);
 plot(dobson(:,2),100*abs(xf(:,5)-x(:,5))/xf(:,5));
 title('Relative uncertainty of Xad');
 
 DS_num = input('Enter the number of the data set (1 or 2): ');
 disp(' ');
 
 dobson_num = input('Enter the number of the Dobson: ');
 dobson_str = num2str(dobson_num);
 disp(' ');
 
 disp(' ');
 
 switch dobson_str
	case '64'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '74'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	case '83'
		disp('        -You are in Izana-       ')
		disp('Lat:28.308, Lon=16.499, Alt:2370.0 m')
	otherwise 
		disp('-----------------------')
		disp('There is no such Dobson')
		disp('-----------------------')
		return
end

disp(' ')

% load the input of the Dobson number
if dobson_str == '64'
%	cal = load('br157.mat');
	disp('Dobson D64')
elseif dobson_str == '74'
%	cal = load('br183.mat');
	disp('Dobson D74')
elseif dobson_str == '83'
%	cal = load('br183.mat');
	disp('Dobson D83')
else
	disp('Wrong Dobson number')
	return
end

disp(' ');



% -------------------------------- incertidumbre de la medida --------------------------------------
 
u2_Nc = (0.005*N(:,1)).^2;
u2_Nd = (0.005*N(:,2)).^2;
u2_Na = (0.005*N(:,3)).^2;

u2_Ncd = u2_Nc + u2_Nd;
u2_Nad = u2_Na + u2_Nd;


 
 % ------------------ Definitions of some parameters ------------------------------------------------

Pstan = 1013.25; % standard pressure (mbars)

R = 6371.229e3; % Earth's radius (m)

t=dobson(:,2);

lon = 28.308;

lat = 16.499;

%-------------------------------- SZA uncertainty---------------------------------------------------
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

% eccentricity correction factor and uncertainty
t_e = (days+1)/365.2422;
ElipLong_1965 = 279.4574;
I = (ElipLong_1965 + 360*t_e + (t0/1460.97))*p0;
u2_I = ((p0/1460.97)*u_t0)^2;

% equation of time in seconds and the uncertainty
et = 4.2*sin(3*I)-2*cos(2*I)+596.5*sin(2*I)-12.8*sin(4*I)+19.3*cos(3*I)-(102.5+0.142*t_e).*sin(I)+(0.033*t_e-429.8).*cos(I);
u2_et = u2_I.*(3*4.2*cos(3*I)+4*sin(2*I)+2*596.5*cos(2*I)-4*12.8*cos(4*I)-3*19.3*sin(3*I)-(102.5+0.142*t_e).*cos(I)-(0.033*t_e-429.8).*sin(I)).^2;

% the hour angle and uncertainty
ha = (t0+(et/60)-720-4*lon)*p0/4;
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
figure(3)
plot(sza,urel_sza)
title('SZA Relative uncertainty');
xlabel('SZA, grad')
ylabel('Relative uncertainty, %')
grid on

%----------------------------Rayleigh air mass uncertainty----------------------------------------------
heff_r = 5e3; % Rayleigh effective altitude (m)
u_heff_r = 2e2;

% Rayleigh air mass and uncertainty
sx_r = (R/(R+heff_r))*sin(sza*p0);
m = sec(asin(sx_r)); % air mass

u2_m = ((sx_r.^2)./((1-sx_r.^2).^3)).*(((R^2)/((R+heff_r)^4)).*((sin(sza*p0)).^2)*(u_heff_r^2) + ((R^2)/((R+heff_r)^2)).*((cos(sza*p0)).^2).*(u2_sza));
urel_m = 100*sqrt(u2_m)./m;

urel_m_A = urel_m;
urel_m_C = urel_m; 
urel_m_D = urel_m;

disp('Rayleigh air mass uncertainty: OK')

% Figure 2 - Rayleigh air mass relative uncertainty 
figure(4)
plot(m,urel_m)
title('Rayleigh air mass relative uncertainty');
xlabel('Rayleigh air mass')
ylabel('Relative uncertainty, %')
grid on


% ---------------------------ozone airmass uncertainty----------------------------------------------
 % Ozone effective altitude (m)
if DS_num == 1
	heff_o = 26-lat/10;
	u_heff_o3 = 2.24e3;
	u_heff_clim = 0; % uncertainty of the ozone effective altitude (m)
elseif DS_num == 2
	heff_o = 24.24e3;
	u_heff_o3 = 0.5e3;
	u_heff_clim = 0; % uncertainty of the ozone effective altitude (m)
end

% combined uncertainty
u_heff_o = sqrt(u_heff_o3^2 + u_heff_clim^2);

% Ozone air mass and uncertainty
sx_o = (R/(R+heff_o))*sin(sza*p0);
mu = sec(asin(sx_o)); % air mass

u2_mu=((sx_o.^2)./((1-sx_o.^2).^3)).*(((R^2)/((R+heff_o)^4)).*((sin(sza*p0)).^2)*(u_heff_o^2) + ((R^2)/((R+heff_o)^2)).*((cos(sza*p0)).^2).*(u2_sza));
urel_mu = 100*sqrt(u2_mu)./mu;

u2_mu_A = u2_mu;
u2_mu_C = u2_mu;
u2_mu_D = u2_mu;

% Figure 3 - Ozone air mass relative uncertainty
figure(5)
plot(mu,urel_mu)
title('Ozone air mass relative uncertainty');
xlabel('Ozone air mass')
ylabel('Relative uncertainty, %')
grid on

disp('Ozone air mass uncertainty: OK')

AC = A1(1)*ones(1,length(sza));
BC = B1(1)*ones(1,length(sza));
pre = 770*ones(1,length(sza));

NC_j = N(:,1);
AC_j = zeros(1,length(sza));
OmegaC_j = zeros(1,length(sza));
muC_j = zeros(1,length(sza));
BC_j = zeros(1,length(sza));
mC_j = zeros(1,length(sza));
pre_j = zeros(1,length(sza));

for r=1:length(sza)
	fprintf('status computation = %2.1f  \n',(r/length(sza))*100);
	
	AC_j(r) = AC(r);

    OmegaC_j(r) = xf(r,1);

    muC_j(r) = mu_dob(r,1);

    BC_j(r) = BC(r);
					
	mC_j(r) = m_dob(r,1);
					
	pre_j(r) = pre(r);
	
	
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
	
	[chi2_min, N] = generate_chi2(NC_j(r),sigma,AC_j(r),OmegaC_j(r),muC_j(r), BC_j(r), mC_j(r),pre_j(r));

    [alpha,beta] = compute_jacobian(N, NC_j(r), AC_j(r),OmegaC_j(r),muC_j(r), BC_j(r),mC_j(r),pre_j(r),sigma);

	while (count~=2) 

    % Compute the alpha_prime which is required in the algorithm
    % basis in 6 dimensions

		alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4) alpha(1,5) alpha(1,6);...
					   alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4) alpha(2,5) alpha(2,6);...
                       alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4) alpha(3,5) alpha(3,6);...
                       alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l) alpha(4,5) alpha(4,6);...
				       alpha(5,1) alpha(5,2) alpha(5,3) alpha(5,4) alpha(5,5)*(1+l) alpha(5,6);...
				       alpha(6,1) alpha(6,2) alpha(6,3) alpha(6,4) alpha(6,5) alpha(6,6)*(1+l)];
	
		% Then we compute the parameters variations

        delta_param = alpha_prime\beta';

        %Compute the new parameters

        AC_j_new(r) = AC_j(r)  + delta_param(1);

        OmegaC_j_new(r)  = OmegaC_j(r)  + delta_param(2);

        muC_j_new(r)  = muC_j(r)  + delta_param(3);

        BC_j_new(r)  = BC_j(r)  + delta_param(4);
			
		mC_j_new(r)  = mC_j(r)  + delta_param(5);
							
		pre_j_new(r)  = pre_j(r)  + delta_param(6);
		
		% Compute the new chi2
		
		[chi2_new, N] = generate_chi2(NC_j(r),sigma,AC_j_new(r),OmegaC_j_new(r),muC_j_new(r),BC_j_new(r),mC_j_new(r),pre_j_new(r));

        diff_chi2 = (chi2_new-chi2_min);
		
		if (diff_chi2<0)

			l = l/10;

			AC_j(r)  = AC_j_new(r); 

			OmegaC_j(r)  = OmegaC_j_new(r); 

			muC_j(r)  = muC_j_new(r); 

			BC_j(r)  = BC_j_new(r); 
								
			mC_j(r)  = mC_j_new(r);
								
			pre_j(r)  = pre_j_new(r);
								

			[alpha, beta] = compute_jacobian(N,NC_j(r),AC_j(r),OmegaC_j(r),muC_j(r),BC_j(r),mC_j(r),pre_j(r),sigma);

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
	
	alpha_prime = [alpha(1,1)*(1+l) alpha(1,2) alpha(1,3) alpha(1,4) alpha(1,5) alpha(1,6);...
                   alpha(2,1) alpha(2,2)*(1+l) alpha(2,3) alpha(2,4) alpha(2,5) alpha(2,6);...
                   alpha(3,1) alpha(3,2) alpha(3,3)*(1+l) alpha(3,4) alpha(3,5) alpha(3,6);...
                   alpha(4,1) alpha(4,2) alpha(4,3) alpha(4,4)*(1+l) alpha(4,5) alpha(4,6);...
				   alpha(5,1) alpha(5,2) alpha(5,3) alpha(5,4) alpha(5,5)*(1+l) alpha(5,6);...
				   alpha(6,1) alpha(6,2) alpha(6,3) alpha(6,4) alpha(6,5) alpha(6,6)*(1+l)];

	
	alpha_prime_corrected = alpha_prime+0.3.*eye(7,7);
                    
    C  = inv(alpha_prime_corrected);
                    
    uncert_AC(r)  = sigma.*sqrt(abs(C(1,1)));

    uncert_OmegaC(r)  = sigma.*sqrt(abs(C(2,2)));

    uncert_muC(r)  = sigma.*sqrt(abs(C(3,3)));

    uncert_BC(r)  = sigma.*sqrt(abs(C(4,4)));
					
	uncert_mC(r)  = sigma.*sqrt(abs(C(5,5)));
					
	uncert_pre(r)  = sigma.*sqrt(abs(C(6,6)));
					
	
	rho(r,1,1) = C(1,1)/(sqrt(abs(C(1,1)*C(1,1))));
    rho(r,1,2) = C(1,2)/(sqrt(abs(C(1,1)*C(2,2))));
    rho(r,1,3) = C(1,3)/(sqrt(abs(C(1,1)*C(3,3))));
    rho(r,1,4) = C(1,4)/(sqrt(abs(C(1,1)*C(4,4))));
	rho(r,1,5) = C(1,5)/(sqrt(abs(C(1,1)*C(5,5))));
	rho(r,1,6) = C(1,6)/(sqrt(abs(C(1,1)*C(6,6))));
					
    rho(r,2,1) = C(2,1)/(sqrt(abs(C(2,2)*C(1,1))));
    rho(r,2,2) = C(2,2)/(sqrt(abs(C(2,2)*C(2,2))));
    rho(r,2,3) = C(2,3)/(sqrt(abs(C(2,2)*C(3,3))));
    rho(r,2,4) = C(2,4)/(sqrt(abs(C(2,2)*C(4,4))));
	rho(r,2,5) = C(2,5)/(sqrt(abs(C(2,2)*C(5,5))));
	rho(r,2,6) = C(2,6)/(sqrt(abs(C(2,2)*C(6,6))));
	                
	rho(r,3,1) = C(3,1)/(sqrt(abs(C(3,3)*C(1,1))));
    rho(r,3,2) = C(3,2)/(sqrt(abs(C(3,3)*C(2,2))));
    rho(r,3,3) = C(3,3)/(sqrt(abs(C(3,3)*C(3,3))));
    rho(r,3,4) = C(3,4)/(sqrt(abs(C(3,3)*C(4,4))));
	rho(r,3,5) = C(3,5)/(sqrt(abs(C(3,3)*C(5,5))));
	rho(r,3,6) = C(3,6)/(sqrt(abs(C(3,3)*C(6,6))));
					
	rho(r,4,1) = C(4,1)/(sqrt(abs(C(4,4)*C(1,1))));
    rho(r,4,2) = C(4,2)/(sqrt(abs(C(4,4)*C(2,2))));
    rho(r,4,3) = C(4,3)/(sqrt(abs(C(4,4)*C(3,3))));
    rho(r,4,4) = C(4,4)/(sqrt(abs(C(4,4)*C(4,4))));
	rho(r,4,5) = C(4,5)/(sqrt(abs(C(4,4)*C(5,5))));
	rho(r,4,6) = C(4,6)/(sqrt(abs(C(4,4)*C(6,6))));
					
	rho(r,5,1) = C(5,1)/(sqrt(abs(C(5,5)*C(1,1))));
	rho(r,5,2) = C(5,2)/(sqrt(abs(C(5,5)*C(2,2))));
	rho(r,5,3) = C(5,3)/(sqrt(abs(C(5,5)*C(3,3))));
	rho(r,5,4) = C(5,4)/(sqrt(abs(C(5,5)*C(4,4))));
	rho(r,5,5) = C(5,5)/(sqrt(abs(C(5,5)*C(5,5))));
	rho(r,5,6) = C(5,6)/(sqrt(abs(C(5,5)*C(6,6))));
					
	rho(r,6,1) = C(6,1)/(sqrt(abs(C(6,6)*C(1,1))));
	rho(r,6,2) = C(6,2)/(sqrt(abs(C(6,6)*C(2,2))));
	rho(r,6,3) = C(6,3)/(sqrt(abs(C(6,6)*C(3,3))));
	rho(r,6,4) = C(6,4)/(sqrt(abs(C(6,6)*C(4,4))));
	rho(r,6,5) = C(6,5)/(sqrt(abs(C(6,6)*C(5,5))));
	rho(r,6,6) = C(6,6)/(sqrt(abs(C(6,6)*C(6,6))));
					
end

disp('Cross correlations: OK')

% Figure 4 - Cross Correlation between A and Omega
figure(6)
plot(dobson(:,2),rho(:,1,2));
datetick('x',dateFormat)
title('Correlation A-\Omega');
xlabel('Time')
ylabel('Correlation')
