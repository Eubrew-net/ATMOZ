function [wc,fwhm]=slit_interp(slit,prec,plot)
% 
x=min(slit(:,1)):prec:max(slit(:,1));
h=interp1(slit(:,1),slit(:,2:end),x');
dh=[h(1,:);diff(h)];

h1=h.*sign(dh); 

[a0,c0]=max(abs(h1));
[a2,b1]=min(abs(h1-0.5));
[a3,b2]=min(abs(h1+0.5));

fwhm=x(b2)-x(b1);
wc=x(b1)+fwhm/2;
%mv=x(c0);