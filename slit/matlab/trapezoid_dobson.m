function y=trapezoid_dobson(x,lamda,a1,a2,a3)
%input  
% x wavelength in nm
% lamda central wavelenths of the slit
% a1 base of slit a o in nm
% a2   of slit at 0.5
% a3 top of slit at 1

 s=sign(x-lamda); 
 y=NaN*x;
 %base
 x2=abs((lamda-x))>=a1/2;
 y(x2)=0;
 % first slope
 x2=abs((x-lamda))<a1/2  & abs((x-lamda))>=a2/2 ;
 y(x2)=0.5*(a1-2*s(x2).*(x(x2)-lamda))./(a1-a2);
 %second slope
 x2=abs((x-lamda))>a3/2  & abs((x-lamda))<a2/2 ;
 y(x2)=0.5+0.5*((a2-2*s(x2).*(x(x2)-lamda))./(a2-a3));
 %top
 y(x>=lamda-a3/2 & x<=lamda+a3/2)=1;
 
  
  