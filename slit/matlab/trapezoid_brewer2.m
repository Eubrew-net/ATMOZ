function y=trapezoid_brewer(x,lamda,F,H)
%input  
% x wavelength in nm
% lamda central wavelenths of the slit
% H= height default =0.8
% a1 base of slit a o in nm
% a1 top  of slit at 1 in nm
if nargin==3
    H=0.87;
end
 s=-sign(x-lamda); 
 y=(s.*(x-lamda)/(2*F))+1;
  
  
 y(x>lamda+2*F | x<lamda-2*F)=0; 
 y(x>lamda-(1-H)*2*F & x<lamda+ (1-H)*F)=H;
 y=y/max(y);
