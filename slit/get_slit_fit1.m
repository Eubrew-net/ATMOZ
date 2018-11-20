function slit = get_slit_fit1( fwhm, central_wl, base,wavelengthsF )
%   get_slit_fit: function used to produce a fit of the slit
%   fwhm: full width at half maximum
%   central_wl: central wavelength
%   base: the level of the base (e.g. 10^-5)
%   width2: width of the gaussian describing stray light
%   wavelengthsF: wl range and step of the output datapoints (in nm - e.g. 300:0.01:360) 

wavelengths=-120:0.01:120; % wavelengths for initial calculation 
b_param=-log10(base); 


% get the fit for the central part
width1=(fwhm+0.0816)./1.93;
gin=(fwhm./(2.*width1));
slope1=-0.3665./log(gin);
dx=wavelengths;%-central_wl;

slit=1*(1+(1.*((exp(-(abs(dx./width1)).^slope1))-1)));

slit=interp1(wavelengths,slit,wavelengthsF-central_wl); % interpolate to the desired wavelengths
slit(slit<=base)=base;

%%
% get the fit for the stray light wavelengths 

% if b_param<=5.5 % use different parameterization for single and double Brewers
%     width2=(0.5793*(b_param^2))+(-3.589*b_param)+5.228;
% elseif b_param>5.5
%        a =   5.497e-07;
%        b =       1.701;
%        c =  -1.569e+07;
%        d =      -3.766;
%     width2= a.*exp(b.*b_param) + c.*exp(d.*b_param);
% end
% 
% aparam=b_param+1.3;
% slope2=(abs(width2)./83.38).^(1/3.824);
% y_strl=10.^(1.*(1+aparam.*(exp(-(abs(dx/width2)).^slope2)-1)));
% 
% % if the fit for the base of the stray light fit differs significantly from
% % the real base, force it to the real base. 
% mnbase=nanmean(y_strl((dx>-40 & dx<-15)|(dx>15 & dx<40)));
% 
% lgdif=log10(base)-log10(mnbase);
% if abs(lgdif)>0.25
%     aparam=aparam-lgdif;
%     y_strl=10.^(1.*(1+aparam.*(exp(-(abs(dx/width2)).^slope2)-1)));
% end
% 
% % get limits seperating the two functions
% ylim=0.01;
% 
% % find the xdatapoints where the y value is closest to the ylim
% dify1=abs(y_centr(dx<0)-ylim);
% zz1=find(dify1==nanmin(dify1));
% dxa=dx(dx<0);
% xlim1=dxa(zz1);
% dify2=abs(y_centr(dx>=0)-ylim);
% zz2=find(dify2==nanmin(dify2));
% dxb=dx(dx>=0);
% xlim2=dxb(zz2);
% 
% % get final slit
% slit=[y_strl(dx<xlim1),y_centr(dx>=xlim1 & dx<=xlim2),y_strl(dx>xlim2)];
% 
% slit=interp1(wavelengths,slit,wavelengthsF-central_wl); % interpolate to the desired wavelengths
end





