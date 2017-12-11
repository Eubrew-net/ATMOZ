function [o3x,x,y]=xcs_cal_nm(peak,o3_set)
% calculate the five crosssection study
% wavelengths in Amstrong !!
%load(fullfile('..','..','o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3x=[];
for k=1:size(o3_set,2)-1  %o3abs

  ka=~isnan(o3_set(:,k+1));  % clean
  x=o3_set(ka,1); %nm   
  cs=o3_set(ka,k+1);  % o3xs 
    
  if size(peak,1)==1 % wc and fwhm  
    %syntetic 
    y=trapezoid_brewer(x,peak(1),peak(2)/2,.87);
    o3x(k)=trapz(x,cs.*y)./trapz(x,y);
  else     %% measured
             % slit interpolated to cross section resolution
      try
        y1=interp1(peak(:,1),peak(:,2),x);
        c=[x,y1,cs];
        ki=all(~isnan(c)')';
        c=c(ki,:);
        o3x(k)=trapz(c(:,1),c(:,2).*c(:,3))./trapz(c(:,1),c(:,2));
      catch
        o3x(k)=NaN;
      end
  end
end
%semilogy(x,cs.*y,'-')