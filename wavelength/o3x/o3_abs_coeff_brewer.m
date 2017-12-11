clear all
%D:\MEDIA\Biblioteca\dropbox\Dropbox\Dobson\brewer_o3x
if ismac
  addpath(genpath('/Volumes/IMAC/CODE/are2011/matlab'))
  addpath(genpath('../matlab')) 
else
   addpath(genpath('..\matlab')) 
   addpath(genpath('d:\code\are2011\matlab'));
end

addpath(genpath(fullfile('..','..','matlab')));
load(fullfile('..','..','o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
[dsp_sum_are,cab]=liesfile('dsp_sum_are2011.txt',1,23);
dsp_summ=dsp_sum_are(:,1:15);
bp_r=dsp_sum_are(:,[1,2,16:end]);
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

o3_set(:,2:end)=o3_set(:,2:end)/log(10);

%% Brewer slits from dispersion
n_samples=size(dsp_summ,1);
o3_set(o3_set(:,1)>400,:)=[];
n_wv=size(o3_set,1);
slits=zeros(n_wv(1),6+1,n_samples);
wv=[];fwhm=[];
figure;hold all;
for i=1:n_samples;
 brw_t=[];
  for j=1:6  % brewer six sltis
   x=o3_set(:,1);%dsp_summ(in_,3+j)+[-10:.01:10];
   y=trapezoid_brewer(x,dsp_summ(i,3+j)/10,dsp_summ(i,9+j)/10,.87);
   brw_t=scan_join(brw_t,[x,y]);;
   fwhm(i,j)=dsp_summ(i,9+j)/10;
   wv(i,j)=dsp_summ(i,3+j)/10;
  end
  slits(:,:,i)=brw_t; 
end
ploty(brw_t,'-.')

%% procedure,slits sintetic
%
% simple aproximation (brewer method)
o3x_sum=NaN*zeros(6,4,n_samples);
o3_abs=NaN*zeros(n_samples,4);

for i=1:n_samples
    figure(1);hold all;
    brw_t=slits(:,:,i);
    c=scan_join_v2(brw_t,o3_set);
    %c=scan_join(brw_t,o3_set);
    
    semilogy(c(:,1),c(:,2:end),'.'); hold on;
    set(gca,'Ylim',[10^-3,10]);
    set(gca,'Xlim',[298,325]);
    semilogyy(brw_t);
    title(sprintf(' %s %d',datestr(dsp_summ(i,1)),dsp_summ(i,2)));
    %DBWV= [305.5000  311.5000  317.6000  325.4000  332.4000  339.8000];
    %% c(1) wv, 2:7 slits  8:11 b&b
    o3x=NaN*zeros(6,4);
    o3x2=o3x;
    figure(2);hold all;
    for j=1:6  %slits
    
        subplot(2,3,j);hold all;
        for k=1:4  %o3abs
            
         r=~isnan(c(:,7+k)); % resolution of o3abs         
         s=c(r,j+1); %slit
         l=c(r,1); %wv
         o=c(r,7+k); %o3 cross
         s(isnan(s))=0;
         o3x(j,k)=trapz(l,s.*o)./trapz(l,s);
         o3x2(j,k)=o3xsec_int([l,o],[l,s]);
         semilogy(l,s.*o,'.')
    end
    if(j==1) legend('Brewer b&p','IGACO BP','D&M','IUP'); end
    title(num2str(o3x(j,:)));
    try
    set(gca,'XLim',wv(i,j)+[-1,1])
    catch
    end    
    end
    o3x_sum(:,:,i)=o3x;
    o3x(isnan(o3x))=0;
    o3abs(:,i)=-(o3x')*O3W';
    o3abs(o3abs(:,i)==0,i)=NaN;
    suptitle([sprintf(' %s Brw %d o3abs=',datestr(dsp_summ(i,1)),dsp_summ(i,2)),num2str(o3abs(:,i)')]);
    %printfiles_append(f1,f2,'o3abs_cal')
end

save o3x_sum o3x_sum;
save o3abs o3abs;

o3abs_c=bp_r(:,5:end)*-O3W(1:end)';
o3x_bp=squeeze(o3x_sum(:,1,:));

%% calculation check
figure;
 plot(o3x_bp'-bp_r(:,5:end))
%%
figure
gscatter(1:length(o3abs_c),[o3abs_c-o3abs(1,:)'],dsp_summ(:,2))
 
%%
figure
plot(100*(o3abs_c-o3abs(1,:)')./o3abs(1,:)')

%%
figure; 
 gscatter(1:length(o3abs_c),100*(o3abs_c-o3abs(1,:)')./o3abs(1,:)',dsp_summ(:,2));
title(' difference % operational process, this work process')
%% ratio for slit
f=figure;hold all;
r_slits=NaN*zeros(6,4,n_samples);o3x_r=squeeze(o3x_sum(1:6,1,:));;for i=1:4,o3x_1=squeeze(o3x_sum(1:6,i,:));r_slits(:,i,:)=o3x_1./o3x_r;end
for j=1:6,
    subplot(2,3,j);
    for i=2:4,
         figure(f+i-2),
         subplot(2,3,j);
        aux=sortrows([dsp_summ(:,9),squeeze(r_slits(j,i,:))],1);
        ploty(aux);
        hold all;
    end;
end
%%
 plot(dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')'),'.')
rline
          
legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')


%%

 m=grpstats([dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')')'],fix(dsp_summ(:,9)*25));
[m,s]=grpstats([dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')')'],fix(dsp_summ(:,9)*25));
legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')

%%
% %d_slit_meas=d_slit;
% figure;
% gplotmatrix(wv(:,2:end),fwhm(:,2:end),dsp_summ(:,2),'brg','+o*.xsd^v<>ph',10)
% %%
% 
% for i=1:6,
%     figure; 
%     scatterhist(wv(:,i),fwhm(:,i));
%     ylabel('fwhm nm');
%     xlabel('wavelength nm');
%     title(sprintf(' Brewer wavelength and fwhm distribution (Set 64 test) Slit #%d',i))
% end
% 
% %%
% 
% 
% for i=1:6,
%     figure; 
%      boxplot(fwhm(:,i),dsp_summ(:,2));
%     ylabel('fwhm nm');
%     xlabel('brewer serial #');
%     title(sprintf(' Brewer fwhm distribution (Set 20 brewer) Slit #%d',i))
% end
% %%
% %%
% 
% 
% for i=1:6,
%     figure; 
%     boxplot(wv(:,i),dsp_summ(:,2));
%     ylabel('wavelength nm');
%     xlabel('brewer serial #');
%     title(sprintf(' Brewer wavelength distribution (Set 20 brewer) Slit #%d',i))
% end

%% results
hist(matdiv(o3abs(2:end,:),o3abs(1,:))',100);
title(' Ozone cross Section ratio to Brewer B&P ');
legend('IGACO BP','D&M','IUP');

%%
figure;
plot(matdiv(o3abs(2:end,:),o3abs(1,:))');
 legend('IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')
 xlabel('# Dispersion tests')
 ylabel('Brewer Ozone Absortion coeffient atm cm^-1')
 title('Ratio to Brewer B&P ');
set(gca,'Ylim',[0.97,1.05]);
grid;

%%

figure;
[m,s,n,gname]=grpstats(o3abs',dsp_summ(:,2),0.5);
 legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')
 xlabel('Brewer')
 ylabel('Brewer Ozone Absortion coeffient atm cm^-1')
 title('Ratio to Brewer B&P ');

%%
figure;
[m,s,n,gname]=grpstats(matdiv(o3abs,o3abs(1,:))',dsp_summ(:,2),0.5);
 legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP')
 xlabel('Brewer')
 ylabel('Ratio')
 title('Ratio to Brewer B&P ');
 
%%
figure

 
%%
legend_o3abs={'Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP'};
rp=round(10000*matdiv(matadd(o3abs,-o3abs(1,:)),o3abs(1,:)))'/100;
displaytable([mean(rp);std(rp)],legend_o3abs,15,'.3f');

%%
legend_o3abs={'Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP'};
ratio=round(10000*matdiv(o3abs,o3abs(1,:)))'/10000;
displaytable([mean(ratio);std(ratio);max(ratio);min(ratio);range(ratio)],legend_o3abs,20,'.4f');

%%
boxplot(ratio(:,2:end),'labels',legend_o3abs(2:end))
title('Ratio to Brewer B&P');
% title('Dobson Slits and Ozone cross Section')
% legend('Brewer b&p','IGACO BP','D&M','IUP')
% plot(DBWV',o3x,'p')
% vline(DBWV');
% load  benhard.dat
% k=~isnan(benhard(:,2));
% plot(benhard(k,1),benhard(k,5),'o');
% %% comparison
% 
% [a,b]=ismember(DBWV,benhard(:,1));
% comp_k=benhard(b,[1,5]); % Khomyr
% 



% %%plot(x,y)
% %%
% b_p=[];d_m=[];iup=[];
% figure;hold all;
% for sl=1:6
%   slit=d_slit(:,[1:2]+(sl-1)*2);
%   slit(isnan(slit(:,1)),:)=[];
%   y1=interp1(slit(:,1),slit(:,2),o3bp_45(:,1));
%   y1=y1/max(y1);
%   y1(isnan(y1))=0;
%   b_p(sl)=trapz(o3bp_45(:,1),y1.*o3bp_45(:,2))/trapz(o3bp_45(:,1),y1);
%   %daumont
%   y2=interp1(slit(:,1),slit(:,2),daumont_xs(:,1));
%   y2=y2/max(y2);
%   y2(isnan(y2))=0;
%   d_m(sl)=trapz(daumont_xs(:,1),y2.*daumont_xs(:,2))/trapz(daumont_xs(:,1),y2);
%   %% uip
%   y3=interp1(slit(:,1),slit(:,2),o3abs_46(:,1));
%   y3=y3/max(y3);
%   y3(isnan(y3))=0;
%   o3abs_46(isnan(o3abs_46(:,2)),2)=0;
%   iup(sl)=trapz(o3abs_46(:,1),y3.*o3abs_46(:,2))/trapz(o3abs_46(:,1),y3);
%   
%   
%   plot(daumont_xs(:,1),[y2*5,daumont_xs(:,2)]);
%   plot(o3abs_46(:,1),[y3*5,o3abs_46(:,2)]);
%   plot(o3bp_45(:,1),[y1*5,o3bp_45(:,2)]);
%   set(gca,'Xlim',[300,350])
%   set(gca,'Ylim',[0,6])
% end
% 
% 
% %%
%  %DBWV=sort([305.5,325.4,311.5,332.4 ,317.6,339.8]) 
%  DBWV= [305.5000  311.5000  317.6000  325.4000  332.4000  339.8000];
%  disp(b_p)
%  disp(d_m)
%  disp(iup)
% %% % A 305.5/325.4; C: 311.5/332.4;  D1): 317.6/339.8
%  
%  A_bp=b_p(1)-b_p(4);
%  A_dm=d_m(1)-d_m(4);
%  A_iup=iup(1)-iup(4);
%  disp('A pair daumont / b&p')
%  disp(A_dm./A_bp)
%  disp('A pair iup / b&p')
%  disp(A_iup./A_bp)
%  %%
%  C_bp=b_p(2)-b_p(5);
%  C_dm=d_m(2)-d_m(5);
%  C_iup=iup(2)-iup(5);
%  
%  disp('C pair daumont / b&p')
%  disp(C_dm./C_bp)
%  disp('C pair iup / b&p')
%  disp(C_iup./C_bp)
%  
%  D_bp=b_p(3)-b_p(6);
%  D_dm=d_m(3)-d_m(6);
%  D_iup=iup(3)-iup(6);
%  
%  disp('D pair daumont / b&p')
%  disp(D_dm./D_bp)
%  disp('D pair iup / b&p')
%  disp(D_iup./D_bp)
%  
%  
%  % AD
%  AD_bp=A_bp-D_bp;
%  AD_dm=A_dm-D_dm;
%  AD_iup=A_iup-D_iup;
%  
%  disp('AD pair daumont / b&p')
%  disp(AD_dm./AD_bp)
%  disp('AD pair iup / b&p')
%  disp(AD_iup./AD_bp)
%  
%  % CD
%  CD_bp=C_bp-D_bp;
%  CD_dm=C_dm-D_dm;
%  CD_iup=C_iup-D_iup;
%  
%  disp('CD pair daumont / b&p')
%  disp(CD_dm/CD_bp)
%  
%  disp('CD pair IUP / b&p')
%  disp(CD_iup/CD_bp)
%  