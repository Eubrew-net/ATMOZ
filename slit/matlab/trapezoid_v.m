function y=trapezoid_v(x,lamda,a1,a3)
%input  
% x wavelength in nm
% lamda central wavelenths of the slit
% a1 base of slit a o in nm
% a1 top  of slit at 1 in nm

 s=sign(x-lamda); 
 y=(a1-2*s.*(x-lamda))./(a1-a3);
  
  
 y(x>lamda+a1/2 | x<lamda-a1/2)=0; 
 y(x>lamda-a3/2 & x<lamda+a3/2)=1;
