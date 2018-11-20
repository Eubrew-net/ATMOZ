function [ o3xsec,o3xj ] = o3xsec_int( o3x,slit )
%function [ o3xsec ] = o3xsec_int( bp,temp )
% integrate the spectra to the slit resolution

%xspec=scan_join(o3x,slit);
%interpolate the slit to spectra resolution  

% remove NaN
slit=slit(~isnan(slit(:,2)),:);
%
yslit=interp1(slit(:,1),slit(:,2),o3x(:,1),'linear');
%%
o3xj=scan_join(o3x,[o3x(:,1),yslit]);
% NaN to zero
o3xj(isnan(o3xj))=0;
% resolution of the o3x
o3xj(o3xj(:,2)==0,:)=[];

%% integrate
o3xsec=trapz(o3xj(:,1),o3xj(:,2).*o3xj(:,3))/trapz(o3xj(:,1),o3xj(:,3));
