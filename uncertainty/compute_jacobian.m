% this function computes the uncertainties tespect to 3 parameters:
% MIDCD, height and SWA


function [alpha, beta] = compute_jacobian(N, N_meas, A_nominal,Omega_nominal,mu_nominal, B_nominal, sigma);



    



% parameters variations used in sensitivity analysis

delta_A = 1e-8;

delta_Omega = 0.00001;

delta_mu = 1e-10;

delta_B = 1e-10;





% new values for the parameters

A_max = A_nominal+delta_A ;

A_min = A_nominal-delta_A ;


Omega_max = Omega_nominal+delta_Omega ; % swa in degrees

Omega_min = Omega_nominal-delta_Omega ;

mu_max = mu_nominal+delta_mu; % swa in degrees

mu_min = mu_nominal-delta_mu ;



B_max = B_nominal+delta_B;

B_min = B_nominal-delta_B;


%--------------------------The main loop starts here
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%  1) Variation of A



N_forwards = A_max*mu_nominal*Omega_nominal+B_nominal;

N_backwards = A_min*mu_nominal*Omega_nominal+B_nominal;

  
derivative_A = (N_forwards-N_backwards)./(2*delta_A);
  
%  2) Variation of Omega



N_forwards = A_nominal*mu_nominal*Omega_max+B_nominal;

N_backwards = A_nominal*mu_nominal*Omega_min+B_nominal;

  
derivative_Omega = (N_forwards-N_backwards)./(2*delta_Omega);


%  3) Variation of mu



N_forwards = A_nominal*mu_max*Omega_nominal+B_nominal;

N_backwards = A_nominal*mu_min*Omega_nominal+B_nominal;

  
derivative_mu = (N_forwards-N_backwards)./(2*delta_mu);

%  3) Variation of B



N_forwards = A_nominal*mu_nominal*Omega_nominal+B_max;

N_backwards = A_nominal*mu_nominal*Omega_nominal+B_min;

  
derivative_B = (N_forwards-N_backwards)./(2*delta_B);
%--------------------------------------------------------------------------



% compute the jacobian
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



% compute the alpha and beta matrices

alpha11 = derivative_A*derivative_A;

alpha12 = derivative_A*derivative_Omega;

alpha13 = derivative_A*derivative_mu;

alpha14 = derivative_A*derivative_B;



alpha21 = derivative_Omega*derivative_A;

alpha22 = derivative_Omega*derivative_Omega;

alpha23 = derivative_Omega*derivative_mu;

alpha24 = derivative_Omega*derivative_B;



alpha31 = derivative_mu*derivative_A;

alpha32 = derivative_mu*derivative_Omega;

alpha33 = derivative_mu*derivative_mu;

alpha34 = derivative_mu*derivative_B;



alpha41 = derivative_B*derivative_A;

alpha42 = derivative_B*derivative_Omega;

alpha43 = derivative_B*derivative_mu;

alpha44 = derivative_B*derivative_B;









alpha = [alpha11 alpha12 alpha13 alpha14; alpha21 alpha22 alpha23 alpha24;...
    alpha31 alpha32 alpha33 alpha34; alpha41 alpha42 alpha43 alpha44]./(sigma.^2);

beta1 = (N_meas-N)*derivative_A./(sigma.^2);
beta2 = (N_meas-N)*derivative_Omega./(sigma.^2);
beta3 = (N_meas-N)*derivative_mu./(sigma.^2);
beta4 = (N_meas-N)*derivative_B./(sigma.^2);



beta = [beta1 beta2 beta3 beta4];



end


