%% This script computes the Jacobian matrix for Lambert-Beer law for O3 determination
%% Then, it computes the correlation matrix from the model to estimate the degrees of correlation
%% Finally, it computes the final uncertainty by including all the uncertainty sources for Dobson instruments


clear all;

clc;
warning off;


%% Loading instrument data
% load Omega.mat; % Ozone values as given by Dobson measurements
% load mu.mat; % mu vaues as given by Dobson measurements
% load m.mat; % air mass values as given by Dobson meas.

% data_dobson = importdata('data_1.txt');
% Omega = data_dobson.data(:,3);
% mu = data_dobson.data(:,1);
% m = data_dobson.data(:,2);
rdobson;
Omega = dobson_ad(:,4);
mu = dobson_ad(:,2);
m = dobson_ad(:,4);



%% Initialization of variables

%% differential cross sections as given by Dobson manual couple AD! Eventually adjust!
delta_alpha2 = 1.8;
delta_alpha1 = 0.3636;

%% Rayleigh differential cross sectino as given by Dobson manual couple AD! Eventually adjust!

delta_beta2 = 0.114;
delta_beta1 = 0.104;

%% Definition of influence parameters.

R = 6371.229e3; % Earth's radius in m

h = 22e3; % heigth Ozone layer in m

r = 2360; % heigth station in m Izana

Z_rad= asin((R+h).*(sqrt(mu.^2-1))./((R+r).*mu));

Z_deg = Z_rad*180/pi;

% Z_deg = zenith_040_1; % zenith in degrees
% 
% Z_rad = Z_deg*pi/180; % zenith in rad

%Omega = ozone_040_1; % in atm-cm

%mu = air_m3_040_1; % Ozone airmass

P0 = 1013; % atmospehric pressure in mbars

P = 769.88; %0.830*1013.25; % staion pressure in mbars


%m = air_mass_040_1; % airmass

delta_delta = 0;

% derivative of mu with respest to Z

dev_mu = -(R+h)*(R+r)^2.*sin(Z_rad).*cos(Z_rad)./((R+h)^2-(R+r)^2.*(sin(Z_rad)).^2).^(3/2);

dev_m = -sin(Z_rad)./(cos(Z_rad)).^2-0.0018167.*(-sin(Z_rad)./(cos(Z_rad)).^2)-0.002875.*(2.*(sec(Z_rad)-1).*(-sin(Z_rad)./(cos(Z_rad)).^2))-...
    0.0008083.*(3.*(sec(Z_rad)-1).^2).*(sec(Z_rad)-1).^2.*(-sin(Z_rad)./(cos(Z_rad)).^2);

%% 

alpha_matrix = zeros(length(Z_rad),8,8);

%% Jacobian matrix. It deals with a 8x8 matrix

J11 = (mu.*Omega).^2;
J12 = -(mu.*Omega).^2;
J13 = mu.^2.*Omega.*(delta_alpha2-delta_alpha1);
J14 = mu.*Omega.*m.*P./P0;
J15 = -(mu.*Omega).*m.*P./P0;
J16 = mu.*Omega.*(delta_beta2-delta_beta1).*m./P0;
J17 = mu.*Omega.*sec(Z_rad);
J18 = mu.*Omega.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%----------------------------------------


J21 = J12;
J22 = (mu.*Omega).^2;
J23 = -mu.^2.*Omega.*(delta_alpha2-delta_alpha1);
J24 = -mu.*Omega.*m.*P./P0;
J25 = mu.*Omega.*m.*P./P0;
J26 = -mu.*Omega.*(delta_beta2-delta_beta1).*m./P0;
J27 = -mu.*Omega.*sec(Z_rad);
J28 =  -mu.*Omega.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%---------------------------------------------


J31 = J13;
J32 = J23;
J33 = ((delta_alpha2-delta_alpha1).*mu).^2;
J34 = (delta_alpha2-delta_alpha1).*mu.*m.*P/P0;
J35 = -(delta_alpha2-delta_alpha1).*mu.*m.*P/P0;
J36 = (delta_alpha2-delta_alpha1).*mu.*m.*(delta_beta2-delta_beta1)/P0;
J37 = (delta_alpha2-delta_alpha1).*mu.*sec(Z_rad);
J38 = (delta_alpha2-delta_alpha1).*mu.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));



%-----------------------------------------------

J41 = J14;
J42 = J24;
J43 = J34;
J44 = (m.*P./P0).^2;
J45 = -(m.*P./P0).^2;
J46 = m.*P./P0.*(delta_beta2-delta_beta1).*m./P0;
J47 = m.*P./P0.*sec(Z_rad);
J48 = m.*P./P0.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%------------------------------------------------

J51 = J15;
J52 = J25;
J53 = J35;
J54 = J45;
J55 = (m.*P./P0).^2;
J56 = -m.*P./P0.*(delta_beta2-delta_beta1).*m./P0;
J57 = -m.*P./P0.*sec(Z_rad);
J58 = -m.*P./P0.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%--------------------------------------------------


J61 = J16;
J62 = J26;
J63 = J36;
J64 = J46;
J65 = J56;
J66 = ones(1,length(Z_rad))'.*((delta_beta2-delta_beta1).*m./P0).^2;
J67 = (delta_beta2-delta_beta1)./P0.*sec(Z_rad);
J68 = (delta_beta2-delta_beta1)./P0.*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%------------------------------------------------

J71 = J17;
J72 = J27;
J73 = J37;
J74 = J47;
J75 = J57;
J76 = J67;
J77 = (sec(Z_rad)).^2;
J78 = sec(Z_rad).*((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2));

%-----------------------------------------------

J81 = J18;
J82 = J28;
J83 = J38;
J84 = J48;
J85 = J58;
J86 = J68;
J87 = J78;
J88 = (((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2))).^2;


%------------------------------------------------

for m=1:length(Z_rad)
    
    
    
    fprintf('status computation = %2.1f  \n',(m/length(Z_rad))*100);

    alpha_matrix= [squeeze(J11(m)) squeeze(J12(m)) squeeze(J13(m)) squeeze(J14(m)) squeeze(J15(m)) squeeze(J16(m)) squeeze(J17(m)) squeeze(J18(m));
            squeeze(J21(m)) squeeze(J22(m)) squeeze(J23(m)) squeeze(J24(m)) squeeze(J25(m)) squeeze(J26(m)) squeeze(J27(m)) squeeze(J28(m));
            squeeze(J31(m)) squeeze(J32(m)) squeeze(J33(m)) squeeze(J34(m)) squeeze(J35(m)) squeeze(J36(m)) squeeze(J37(m)) squeeze(J38(m));
            squeeze(J41(m)) squeeze(J42(m)) squeeze(J43(m)) squeeze(J44(m)) squeeze(J45(m)) squeeze(J46(m)) squeeze(J47(m)) squeeze(J48(m));
            squeeze(J51(m)) squeeze(J52(m)) squeeze(J53(m)) squeeze(J54(m)) squeeze(J55(m)) squeeze(J56(m)) squeeze(J57(m)) squeeze(J58(m));
            squeeze(J61(m)) squeeze(J62(m)) squeeze(J63(m)) squeeze(J64(m)) squeeze(J65(m)) squeeze(J66(m)) squeeze(J67(m)) squeeze(J68(m));
            squeeze(J71(m)) squeeze(J72(m)) squeeze(J73(m)) squeeze(J74(m)) squeeze(J75(m)) squeeze(J76(m)) squeeze(J77(m)) squeeze(J78(m));
            squeeze(J81(m)) squeeze(J82(m)) squeeze(J83(m)) squeeze(J84(m)) squeeze(J85(m)) squeeze(J86(m)) squeeze(J87(m)) squeeze(J88(m))];
        
        
        

%--------------------------------------------------------

             %% Determination of the correlation matrix
             
             
            alpha_prime_corrected = alpha_matrix+0.2*eye(8,8);
                    
            C  = inv(alpha_prime_corrected);
            

            for j=1:8

                for k=1:8

                    rho(m,j,k) = C(j,k)/(sqrt(C(j,j)*C(k,k)));


                end
            end


        
              
end  


    %% Now we compute the total uncertainty based on the input matrix provided by Ulf Kohler
    %% Relative uncertainties
    




    u_rel_alpha = 0.01; % wavelength adjustment and  bandwidth effect
    u_rel_o3 = 0.01; % ozone absorption
    u_rel_temp = 0.01; % temperature ozone layer
    
    u_rel_combined_abs = sqrt(u_rel_alpha^2+u_rel_o3.^2+u_rel_temp^2);
%---------------------------------------------------------
    

    u_rel_beta = 0.001; %atmospheric scattering

%----------------------------------------------------------
 

    u_rel_P = 0.0016; % relative Atmospheric relative uncertainty
    
%----------------------------------------------------------
    
    u_rel_aerosol = 0.002; %particle scattering
 
%----------------------------------------------------------

    u_rel_etc = sqrt(0.01^2+0.005^2+0.005^2+0.01^2); %direct intercomparison, Langley method for abs. cal. of reference, spectral stability and
                                              %relative optical path through the ozone layer

%-------------------------------------------------------


    u_rel_wedge = 0.005; %wedge calibration
    
%--------------------------------------------------------    


    u_rel_straylight = sqrt(0.01^2+0.001^2); %stray light (internal and external)
    
    
%-------------------------------------------------------
    
    u_rel_y = sqrt(u_rel_etc.^2+u_rel_wedge.^2); % uncertainty readings
    
%-------------------------------------------------------
    
    u_rel_Z = 0.00005; % uncertainty Zenith angle
    
%------------------------------------------------------
    
    u_rel_interference_atmopsheric_cost = 0.001; %interference with other atmospheric constituents
    
%------------------------------------------------------

    
    u_rel_operat_impre = 0.01; %operational precision

    
    %%-----------------------------
    
%% Here we compute the expression of the TOC uncertainty, as function of all the other uncertainties contributions and degrees of correlations.
%% In order to do this, we need to compute all the sensitivity coefficients. 

% Sensitivity coefficients
    
    dy_da1 = mu.*Omega;
    
    dy_da2 = -mu.*Omega;
    
    dy_da3 = ((delta_alpha2-delta_alpha1).*mu);
    
    dy_da4 = m.*P./P0;
    
    dy_da5 = -m.*P./P0;
    
    dy_da6 = (delta_beta2-delta_beta1).*m./P0;
    
    dy_da7 = sec(Z_rad);
    
    dy_da8 = (((delta_alpha2-delta_alpha1).*Omega.*dev_mu+(delta_beta2-delta_beta1).*P./P0.*dev_m+delta_delta.*(-sin(Z_rad)./cos(Z_rad).^2)));
    
 %%------------------------------------------------------------------------------------------------------------------------------------------
 
 %% Coefficients of the second order equation to solve
    
    A1 = dy_da3.^2;

    A2 = 2.*((u_rel_combined_abs.*dy_da1.*dy_da3.*rho(:,1,3))+u_rel_combined_abs.*dy_da2.*dy_da3.*rho(:,2,3)+...
            u_rel_beta.*dy_da4.*dy_da3.*rho(:,4,3)+u_rel_beta.*dy_da5.*dy_da3.*rho(:,5,3)+u_rel_P.*dy_da6.*dy_da3.*rho(:,6,3)+...
            u_rel_aerosol.*dy_da7.*dy_da3.*rho(:,7,3)+u_rel_Z.*dy_da8.*dy_da3.*rho(:,8,3));
        
    A3 = 2.*(u_rel_combined_abs.*u_rel_combined_abs.*dy_da1.*dy_da2.*rho(:,1,2)+u_rel_combined_abs.*u_rel_beta.*dy_da1.*dy_da4.*rho(:,1,4)+...
            u_rel_combined_abs.*u_rel_beta.*dy_da1.*dy_da5.*rho(:,1,5)+u_rel_combined_abs.*u_rel_P.*dy_da1.*dy_da6.*rho(:,1,6)+...
            u_rel_combined_abs.*u_rel_aerosol.*dy_da1.*dy_da7.*rho(:,1,7)+u_rel_combined_abs.*u_rel_Z.*dy_da1.*dy_da8.*rho(:,1,8)+...
            u_rel_combined_abs.*u_rel_combined_abs.*dy_da2.*dy_da1.*rho(:,2,1)+u_rel_combined_abs.*u_rel_beta.*dy_da2.*dy_da4.*rho(:,2,4)+...
            u_rel_combined_abs.*u_rel_beta.*dy_da2.*dy_da5.*rho(:,2,5)+u_rel_combined_abs.*u_rel_P.*dy_da2.*dy_da6.*rho(:,2,6)+...
            u_rel_combined_abs.*u_rel_aerosol.*dy_da2.*dy_da7.*rho(:,2,7)+u_rel_combined_abs.*u_rel_Z.*dy_da2.*dy_da8.*rho(:,2,8)+...
            u_rel_beta.*u_rel_combined_abs.*dy_da4.*dy_da1.*rho(:,4,1)+u_rel_beta.*u_rel_combined_abs.*dy_da4.*dy_da2.*rho(:,4,2)+...
            u_rel_beta.*u_rel_beta.*dy_da4.*dy_da5.*rho(:,4,5)+u_rel_beta.*u_rel_P.*dy_da4.*dy_da6.*rho(:,4,6)+...
            u_rel_beta.*u_rel_aerosol.*dy_da4.*dy_da7.*rho(:,4,7)+u_rel_beta.*u_rel_Z.*dy_da4.*dy_da8.*rho(:,4,8)+...
            u_rel_beta.*u_rel_combined_abs.*dy_da5.*dy_da1.*rho(:,5,1)+u_rel_beta.*u_rel_combined_abs.*dy_da5.*dy_da2.*rho(:,5,2)+...
            u_rel_beta.*u_rel_beta.*dy_da5.*dy_da4.*rho(:,5,4)+u_rel_beta.*u_rel_P.*dy_da5.*dy_da6.*rho(:,5,6)+...
            u_rel_beta.*u_rel_aerosol.*dy_da5.*dy_da7.*rho(:,5,7)+u_rel_beta.*u_rel_Z.*dy_da5.*dy_da8.*rho(:,5,8)+...
            u_rel_P.*u_rel_combined_abs.*dy_da6.*dy_da1.*rho(:,6,1)+u_rel_P.*u_rel_combined_abs.*dy_da6.*dy_da2.*rho(:,6,2)+...
            u_rel_P.*u_rel_beta.*dy_da6.*dy_da4.*rho(:,6,4)+u_rel_P.*u_rel_beta.*dy_da6.*dy_da5.*rho(:,6,5)+...
            u_rel_P.*u_rel_aerosol.*dy_da6.*dy_da7.*rho(:,6,7)+u_rel_P.*u_rel_Z.*dy_da6.*dy_da8.*rho(:,6,8)+...
            u_rel_aerosol.*u_rel_combined_abs.*dy_da7.*dy_da1.*rho(:,7,1)+u_rel_aerosol.*u_rel_combined_abs.*dy_da7.*dy_da2.*rho(:,7,2)+...
            u_rel_aerosol.*u_rel_beta.*dy_da7.*dy_da4.*rho(:,7,4)+u_rel_aerosol.*u_rel_beta.*dy_da7.*dy_da5.*rho(:,7,5)+...
            u_rel_aerosol.*u_rel_P.*dy_da7.*dy_da6.*rho(:,7,6)+u_rel_aerosol.*u_rel_Z.*dy_da7.*dy_da8.*rho(:,7,8)+...
            u_rel_Z.*u_rel_combined_abs.*dy_da8.*dy_da1.*rho(:,8,1)+u_rel_Z.*u_rel_combined_abs.*dy_da8.*dy_da2.*rho(:,8,2)+...
            u_rel_Z.*u_rel_beta.*dy_da8.*dy_da4.*rho(:,8,4)+u_rel_Z.*u_rel_beta.*dy_da8.*dy_da5.*rho(:,8,5)+...
            u_rel_Z.*u_rel_P.*dy_da8.*dy_da6.*rho(:,8,6)+u_rel_Z.*u_rel_aerosol.*dy_da8.*dy_da7.*rho(:,8,7))-u_rel_y.^2;
        
%----------------------------------------------------------------------------------------------------------------------------------

%% Two solutions
        
   u_toc_final_1 = abs((-A2+sqrt(A2.^2-4.*A1.*A3))./(2.*A1));
   u_toc_final_2 = abs((-A2-sqrt(A2.^2-4.*A1.*A3))./(2.*A1));
   
%% Final result   
   
   u_toc_final__combined = sqrt(u_toc_final_2.^2+u_rel_interference_atmopsheric_cost.^2+u_rel_operat_impre.^2); 
 

%% outputs


figure(20)
plot(Z_deg,u_toc_final__combined./Omega.*100,'.');
ylabel('relative TOC uncertainty %');
xlabel('Zenith angle')

figure(21)
plot(Z_deg,Omega.*1000,'ob');
hold on;
plot(Z_deg,u_toc_final__combined.*1000,'xr');
xlabel('Zenith angle');

figure(22)
plot(Z_deg,u_toc_final__combined.*1000,'.r');
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