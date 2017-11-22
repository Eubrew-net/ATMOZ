% this function computes the uncertainties tespect to 7 parameters:
% MIDCD, height and SWA


function [alpha, beta] = compute_jacobian(N, N_meas, A_nominal,Omega_nominal,mu_nominal, B_nominal,m_nominal,p_nominal,etc_nominal, sigma);



% parameters variations used in sensitivity analysis

delta_A = 1e-8;

delta_Omega = 1e-5;

delta_mu = 1e-10;

delta_B = 1e-10;

delta_m = 1e-10;

delta_p = 1e-8;

delta_etc = 1e-10;


% new values for the parameters

A_max = A_nominal+delta_A ;

A_min = A_nominal-delta_A ;


Omega_max = Omega_nominal+delta_Omega ; % swa in degrees

Omega_min = Omega_nominal-delta_Omega ;


mu_max = mu_nominal+delta_mu; % swa in degrees

mu_min = mu_nominal-delta_mu ;



B_max = B_nominal+delta_B;

B_min = B_nominal-delta_B;


m_max = m_nominal+delta_m;

m_min = m_nominal-delta_m;


p_max = p_nominal+delta_p;

p_min = p_nominal-delta_p;


etc_max = etc_nominal+delta_etc;

etc_min = etc_nominal-delta_etc;

%--------------------------The main loop starts here
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%  1) Variation of A



N_forwards = etc_nominal - 10*A_max*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_nominal/1013.25;

N_backwards = etc_nominal - 10*A_min*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_nominal/1013.25;

  
derivative_A = (N_forwards-N_backwards)./(2*delta_A);
  
%  2) Variation of Omega



N_forwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_max - B_nominal*m_nominal*p_nominal/1013.25;

N_backwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_min - B_nominal*m_nominal*p_nominal/1013.25;

  
derivative_Omega = (N_forwards-N_backwards)./(2*delta_Omega);


%  3) Variation of mu



N_forwards = etc_nominal - 10*A_nominal*mu_max*Omega_nominal - B_nominal*m_nominal*p_nominal/1013.25;

N_backwards = etc_nominal - 10*A_nominal*mu_min*Omega_nominal - B_nominal*m_nominal*p_nominal/1013.25;

  
derivative_mu = (N_forwards-N_backwards)./(2*delta_mu);

%  4) Variation of B



N_forwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_max*m_nominal*p_nominal/1013.25;

N_backwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_min*m_nominal*p_nominal/1013.25;

  
derivative_B = (N_forwards-N_backwards)./(2*delta_B);

%  5) Variation of m

N_forwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_max*p_nominal/1013.25;

N_backwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_min*p_nominal/1013.25;


derivative_m = (N_forwards-N_backwards)./(2*delta_m);

%  6) Variation of p

N_forwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_max/1013.25;

N_backwards = etc_nominal - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_min/1013.25;


derivative_p = (N_forwards-N_backwards)./(2*delta_p);

%  7) Variation of ETC

N_forwards = etc_max - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_nominal/1013;

N_backwards = etc_min - 10*A_nominal*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_nominal/1013;

derivative_etc = (N_forwards-N_backwards)./(2*delta_etc); 

%--------------------------------------------------------------------------



% compute the jacobian
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------



% compute the alpha and beta matrices

alpha11 = derivative_A*derivative_A;

alpha12 = derivative_A*derivative_Omega;

alpha13 = derivative_A*derivative_mu;

alpha14 = derivative_A*derivative_B;

alpha15 = derivative_A*derivative_m;

alpha16 = derivative_A*derivative_p;

alpha17 = derivative_A*derivative_etc;


alpha21 = derivative_Omega*derivative_A;

alpha22 = derivative_Omega*derivative_Omega;

alpha23 = derivative_Omega*derivative_mu;

alpha24 = derivative_Omega*derivative_B;

alpha25 = derivative_Omega*derivative_m;

alpha26 = derivative_Omega*derivative_p;

alpha27 = derivative_Omega*derivative_etc;



alpha31 = derivative_mu*derivative_A;

alpha32 = derivative_mu*derivative_Omega;

alpha33 = derivative_mu*derivative_mu;

alpha34 = derivative_mu*derivative_B;

alpha35 = derivative_mu*derivative_m;

alpha36 = derivative_mu*derivative_p;

alpha37 = derivative_mu*derivative_etc;



alpha41 = derivative_B*derivative_A;

alpha42 = derivative_B*derivative_Omega;

alpha43 = derivative_B*derivative_mu;

alpha44 = derivative_B*derivative_B;

alpha45 = derivative_B*derivative_m;

alpha46 = derivative_B*derivative_p;

alpha47 = derivative_B*derivative_etc;



alpha51 = derivative_m*derivative_A;

alpha52 = derivative_m*derivative_Omega;

alpha53 = derivative_m*derivative_mu;

alpha54 = derivative_m*derivative_B;

alpha55 = derivative_m*derivative_m;

alpha56 = derivative_m*derivative_p;

alpha57 = derivative_m*derivative_etc;



alpha61 = derivative_p*derivative_A;

alpha62 = derivative_p*derivative_Omega;

alpha63 = derivative_p*derivative_mu;

alpha64 = derivative_p*derivative_B;

alpha65 = derivative_p*derivative_m;

alpha66 = derivative_p*derivative_p;

alpha67 = derivative_p*derivative_etc;


alpha71 = derivative_etc*derivative_A;

alpha72 = derivative_etc*derivative_Omega;

alpha73 = derivative_etc*derivative_mu;

alpha74 = derivative_etc*derivative_B;

alpha75 = derivative_etc*derivative_m;

alpha76 = derivative_etc*derivative_p;

alpha77 = derivative_etc*derivative_etc;



alpha = [alpha11 alpha12 alpha13 alpha14 alpha15 alpha16 alpha17; alpha21 alpha22 alpha23 alpha24 alpha25 alpha26 alpha27;...
    alpha31 alpha32 alpha33 alpha34 alpha35 alpha36 alpha37; alpha41 alpha42 alpha43 alpha44 alpha45 alpha46 alpha47;...
	alpha51 alpha52 alpha53 alpha54 alpha55 alpha56 alpha57; alpha61 alpha62 alpha63 alpha64 alpha65 alpha66 alpha67;...
	alpha71 alpha72 alpha73 alpha74 alpha75 alpha76 alpha77]./(sigma.^2);

beta1 = (N_meas-N)*derivative_A./(sigma.^2);
beta2 = (N_meas-N)*derivative_Omega./(sigma.^2);
beta3 = (N_meas-N)*derivative_mu./(sigma.^2);
beta4 = (N_meas-N)*derivative_B./(sigma.^2);
beta5 = (N_meas-N)*derivative_m./(sigma.^2);
beta6 = (N_meas-N)*derivative_p./(sigma.^2);
beta7 = (N_meas-N)*derivative_etc./(sigma.^2);



beta = [beta1 beta2 beta3 beta4 beta5 beta6 beta7];



end


