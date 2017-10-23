% This function generates a chi2 given a set of shape parameters


function [chi2_min, N] = generate_chi2(N_meas,sigma,A_nominal,Omega_nominal,mu_nominal,B_nominal,m_nominal,p_nominal,etc_nominal);

    
    

     N = etc_nominal - A_nominal*mu_nominal*Omega_nominal - B_nominal*m_nominal*p_nominal/1013.25;


 
     
     % Then we compute the final chi2
     
    
     
     chi2_min = [((N_meas-N)./(sigma)).^2];

    
     
     
     
     
end