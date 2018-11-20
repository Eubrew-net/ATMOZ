addpath(genpath('matlab')) 

%addpath(genpath(fullfile('..','..','matlab')));
load(fullfile('.','o3x_atmoz','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45; %i45
% o3_set=scan_join(o3_set,o3x_45); % Brewer
% o3_set=scan_join(o3_set,o3bp_45); %igaco
% o3_set=scan_join(o3_set,daumont_45); % Daumont
% o3_set=scan_join(o3_set,o3sgw7_45);   %uip
% o3_set=scan_join(o3_set,sdt.t_45);   %uip quad
% o3_set=scan_join(o3_set,BPB_46); % Bernhard 
% o3_set=scan_join(o3_set,o3gw7_45);   %uip 2016
xs_legend={'Brw','B&P','DMB','UIP','UIQ','B05','GW7'};

%[dsp_sum_are,cab]=liesfile('dsp_sum_are2011.txt',1,23);
[dsp_sum_are,cab]=liesfile(fullfile('.','dsp_sum.txt'),1,23);
%%
[C,ia,ib]=unique(dsp_sum_are(:,1:2),'rows','stable');
dsp_summ=dsp_sum_are(ia,1:15);
%dsp_summ=dsp_sum_are(:,1:15);
bp_r=dsp_sum_are(ia,[1,2,16:end]);


dsp_m=nanmean(dsp_summ);
% 017 nan values replaced by mean brewer=1;
dsp_summ(2,:)=dsp_m;
dsp_summ(2,2)=1;
dsp_summ(2,1)=dsp_summ(1,1);
% values not available on 14
dsp_summ(1,4)=dsp_summ(2,4);
dsp_summ(1,10)=dsp_summ(2,10);




O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

o3_set(:,2:end)=o3_set(:,2:end)/log(10);
%% Brewer slits with stray-light
slit_=load('slit_are2015.dat');
brw_=slit_(1,2:end);
slit_n=slit_(2,2:end);
%fwhm
slit_fwhm=slit_(3,2:end);
%stray light rejection
slit_slr=10.^slit_(4,2:end);

% we dont have the original wavelenghs -> nominal as first approx
n_samples=length(unique(brw_));

% added one column for 0 slit
 fw=reshape(slit_fwhm,5,[])';
 fw=[fw(:,1),fw];
 slr=reshape(slit_slr,5,[])';
 slr=[slr(:,1),slr];
 %sl=reshape(slit_(5:end,2:end),5,[])';

 sl=reshape(slit_(5:end,2:end),[],7*5)';
 wsl=slit_(5:end,1);
% xsetction
o3_set(o3_set(:,1)>400 | o3_set(:,1)<280,:)=[];
n_wv=size(o3_set,1);
slits=zeros(n_wv(1),6+1,n_samples);

% Stray Light

o3_sum_sl=NaN*zeros(6,7,n_samples);
o3_abs_sl=NaN*zeros(n_samples,5);

o3x_gaus=NaN*zeros(6,7,n_samples);
o3_gaus=NaN*zeros(n_samples,5);


wv_sl(n_samples,6)=0;
fwhm_sl(n_samples,6)=0;

for i=1:n_samples
    figure(1);hold all;
    o3x=NaN*zeros(6,4);
    %o3x2=o3x;
    o3x_sl=o3x;
    
    %figure
    for j=1:6  %slit  -
        %figure(j);subplot(2,3,j);hold off
         for k=1:7  %o3abs
             ka=~isnan(o3_set(:,k+1));
             x=o3_set(ka,1);
             cs=o3_set(ka,k+1);
             %y0=trapezoid_brewer(x,BWN(j),fw(i,j)/10/2,.87);
             y=(get_slit_fit1(fw(i,j)/10,BWN(j),1E-10,x));
             y1=(get_slit_fit2(fw(i,j)/10,BWN(j),slr(i,j),x));
            % if isreal(y1)==0
            %     y1=real(y1);
            % end
            % if j==1
            % y2=interp1(wsl/10,sl((i-1)*5+j,:),x,'linear''extrap')
            %  else
            % y2=interp1(wsl/19,sl((i-1)*5+j,:),x,'linear''extrap');
            % end
             wv(i,j)=BWN(j);%dsp_summ(i,3+j)/10;
             fwhm(i,j)=fw(i,j)/10/2;%dsp_summ(i,9+j)/10/2;
             o3x(j,k)=trapz(x,cs.*y)./trapz(x,y);
             %disp(i)
             o3x_sl(j,k)=trapz(x,cs.*y1)./trapz(x,y1);
             %o3x2(j,k)=o3xsec_int([l,o],[l,s]);
             figure(1)
             semilogy(x,cs.*y,'x');hold on
             semilogy(x,cs.*y1,'-');
             title('');
    end
    %if(j==1) legend('Brewer b&p','IGACO BP','D&M','IUP','orient','horizontal'); end
    %title(num2str(o3x(j,:)));
    %try
    %set(gca,'XLim',wv(i,j)+[-1,1])
    %catch
    %end    
    end
    
    
    
    o3x_gaus(:,:,i)=o3x;
    o3x(isnan(o3x))=0;
    
    o3_gaus(:,i)=-(o3x')*O3W';
    o3_gaus(o3_gaus(:,i)==0,i)=NaN;
    
    o3x_sum_sl(:,:,i)=o3x_sl;
    o3x_sl(isnan(o3x_sl))=0;
    
    o3_abs_sl(:,i)=-(o3x_sl')*O3W'
    o3_abs_sl(o3_abs_sl(:,i)==0,i)=NaN;
    
    %suptitle([sprintf(' %s Brw %d o3abs=',datestr(dsp_summ(i,1)),dsp_summ(i,2)),num2str(o3abs(:,i)')]);
    %printfiles_append(f1,f2,'o3abs_cal')
end


save o3_sum_sl o3_sum_sl;
save o3_abs_sl o3_abs_sl;


save o3x_gaus o3x_gaus;
save o3_gaus o3_gaus;

save wk;


o3abs_c=bp_r(:,5:end)*-O3W(1:end)';
o3x_bp=squeeze(o3x_sum(:,1,:));

%%
figure;plot(100*(o3_abs_sl-o3_gaus)./o3_gaus,'.-');
legend(num2str(unique(brw_)'));
grid
set(gca,'xticklabel',xs_legend)
%%
figure;plot(100*((o3_abs_sl./o3_gaus)-1)','.-')
set(gca,'xtick',1:7,'xticklabel',num2str(unique(brw_)'))
legend(xs_legend)
grid


%% Brewer slits from dispersion
n_samples=size(dsp_summ,1);
o3_set(o3_set(:,1)>340 | o3_set(:,1)<280,:)=[];
n_wv=size(o3_set,1);
slits=zeros(n_wv(1),6+1,n_samples);

%% procedure,slits sintetic
%
% simple aproximation (brewer method)
o3x_sum=NaN*zeros(6,7,n_samples);
o3_abs=NaN*zeros(n_samples,5);
wv(n_samples,6)=0;
fwhm(n_samples,6)=0;

for i=1:n_samples
    %figure(1);hold all;
    o3x=NaN*zeros(6,4);
    o3x2=o3x;
    %figure(2);hold all;
    for j=1:6  %slit
        subplot(2,3,j);hold all;
         for k=1:7  %o3abs
             ka=~isnan(o3_set(:,k+1));
             x=o3_set(ka,1);
             cs=o3_set(ka,k+1);
             y=trapezoid_brewer2(x,dsp_summ(i,3+j)/10,dsp_summ(i,9+j)/10/2,.87);
             wv(i,j)=dsp_summ(i,3+j)/10;fwhm(i,j)=dsp_summ(i,9+j)/10/2;
             o3x(j,k)=trapz(x,cs.*y)./trapz(x,y);
            %o3x2(j,k)=o3xsec_int([l,o],[l,s]);
            %semilogy(x,cs.*y,'-')
    end
    %if(j==1) legend('Brewer b&p','IGACO BP','D&M','IUP','orient','horizontal'); end
    %title(num2str(o3x(j,:)));
    %try
    %set(gca,'XLim',wv(i,j)+[-1,1])
    %catch
    %end    
    end
    o3x_sum(:,:,i)=o3x;
    o3x(isnan(o3x))=0;
    o3abs(:,i)=-(o3x')*O3W';
    o3abs(o3abs(:,i)==0,i)=NaN;
    %suptitle([sprintf(' %s Brw %d o3abs=',datestr(dsp_summ(i,1)),dsp_summ(i,2)),num2str(o3abs(:,i)')]);
    %printfiles_append(f1,f2,'o3abs_cal')
end

save o3x_sum o3x_sum;
save o3abs o3abs;
save wk;

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
% %% ratio for slit
% f=figure;hold all;
% r_slits=NaN*zeros(6,4,n_samples);o3x_r=squeeze(o3x_sum(1:6,1,:));;for i=1:4,o3x_1=squeeze(o3x_sum(1:6,i,:));r_slits(:,i,:)=o3x_1./o3x_r;end
% for j=1:6,
%     subplot(2,3,j);
%     for i=2:4,
%          figure(f+i-2),
%          subplot(2,3,j);
%         aux=sortrows([dsp_summ(:,6+j),squeeze(r_slits(j,i,:))],1);
%         ploty(aux);
%         hold all;
%     end;
% end
%%
% plot(dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')'),'.')
%rline
%          
% legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')

%% AROSA

x=[040,072,156];
ja=ismember(dsp_summ(:,2),x);
arosa_b=[dsp_summ(ja,2)';o3abs(:,ja)];
printmatrix(grpstats(arosa_b',arosa_b(1,:)')',4)
%%

% m=grpstats([dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')')'],fix(dsp_summ(:,9)*25));
[m,s]=grpstats([dsp_summ(:,9),matadd(matdiv(o3abs,o3abs(1,:)),- nanmean(matdiv(o3abs,o3abs(1,:))')')'],fix(dsp_summ(:,9)*25));
legend('Brewer Bass & Paur','IGACO B&P','Daumont & Malicet','IUP','IUPQ','orientation','horizontal')
hold on
%% nominal
figure;
scatterhist(0.1*(dsp_summ(:,9)-median(dsp_summ(:,9))),o3abs_c-nanmedian(o3abs_c))

ylabel('O_3 abs coeff - O_3 abs coeff_{(nominal)} (amt  cm)^{-1}')
%title({'BOp. diferences to nominal Brewer vs wavelength diferences to nominal'})
xlabel('wavelength \lambda_{6} -\lambda_{6 nominal} (nm)')

%robust_line_
b=regress((o3abs_c-nanmedian(o3abs_c)),0.1*(dsp_summ(:,9)-median(dsp_summ(:,9))))
refline(b);
set(gca,'Xlim',0.1*[-0.5,0.5]);
set(gca,'Xtick',0.1*(-0.5:0.1:0.5));
set(gca,'Ylim',[-5.5,5.5]*10^-3);

set(gcf,'Tag','Brw_BP_nominal_wv');
printfiles_publication(gcf,'.','LineWidth','auto','Width',18,'Height',12)



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
figure
hist(matdiv(o3abs(2:end,:),o3abs(1,:))',100);
title(' Ozone cross Section ratio to Brewer B&P ');
legend(xs_legend(2:end));





%% The change in ozone is correlated in waveleng 4 and in less in wv 3
%  c= a + b(l4-lo)
%  
figure
a=[];b={}; tab=[];
j=(abs(dsp_summ(:,9)-mean(dsp_summ(:,9)))<1);
for ii=1:size(o3abs,1)
    [a(ii,:),b{ii}]=robustfit(.1*(dsp_summ(j,9)-mean(dsp_summ(j,9))),...
        (o3abs(ii,j)./o3abs(1,j)));
    
    rs{ii}=regstats(o3abs(ii,j)./o3abs(1,j),.1*(dsp_summ(j,9)-mean(dsp_summ(j,9))),'linear');
    %plot(dsp_summ(j,9)-mean(dsp_summ(j,9)),...
    %    (o3abs(ii,j)./o3abs(1,j))-mean(o3abs(ii,j)./o3abs(1,j)),'.');
    
    h(ii)=plot(.1*(dsp_summ(j,9)-mean(dsp_summ(j,9))),...
        (o3abs(ii,j)./o3abs(1,j)),'.');

    hold all
    refline(a(ii,end:-1:1));
    tab(ii,:)=[rs{ii}.beta',rs{ii}.rsquare];
end
set(h(3),'MarkerSize',5,'Marker','p')
set(h(4),'MarkerSize',5,'Marker','s')
set(h(5),'Marker','x')
set(h(2),'Marker','+') 
%title({'Brewer BOp DXS ratios to different cross sections  vs  wavelength diferences to nominal'})
legend(h(2:end),'IGQ4','DBM','IUP','IUP Q','orientation','horizontal',2);
ylabel('A_{i,xs} / A_{i,op} ');
%(amt cm)^{-1}
%xlabel('wavelength S4 - wavelength So4_{nominal} (A)')
set(gca,'Xlim',.1*[-0.5,0.5]);
set(gca,'ylim',[0.98,1.04]);
%ylabel('O_3 abs coeff - O_3 abs coeff_{(nominal)}')
xlabel('wavelength \lambda_{6} -\lambda_{6 nominal} (nm)')
grid
set(gcf,'Tag','Brw_ratio_to_BP_nominal_wv');
printfiles_publication(gcf,'.','LineWidth','auto','Width',17,'Height',13)
disp('Figure 5')
   
%% 
% figure
% ii=3
% plot(dsp_summ(j,9)-mean(dsp_summ(j,9)),...
%         (o3abs(ii,j)./o3abs(1,j))-mean(o3abs(ii,j)./o3abs(1,j)),'.');
% set(gca,'Xlim',[-0.5,0.5]);
% rl1=refline(a(ii,2));   
% 
% title('{ B&P Brewer operative B.oP ratio to DBM  difference to nominal
% to diferences to nominal vs wavelength diferences to nominal')
% robust_line_
% set(gca,'Xlim',[-0.5,0.5]);
% 



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
rotateticklabel(gca)
%legend('BrwBP','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')
legend_o3abs={'BOp','IGQ4','DBM','IUP','UIPQ'};
legend(legend_o3abs,'Orientation','Horizontal');

 l1=xlabel('')
 ylabel(' A_{xs} {(atm cm)^{ - 1}}')
 %title('Brewer ozone absorption coeffient RBCC-E dataset ');
 title('');
set(gcf,'Tag','Brw_DXS');
printfiles_publication(gcf,'.','LineWidth','auto','Width',16,'Height',12)

%$$c{cm^{ - 1}}$$

%%
figure;
[m,s,n,gname]=grpstats(matdiv(o3abs,o3abs(1,:))',dsp_summ(:,2),0.5);
rotateticklabel(gca)
%legend('BrwBP','IGACO B&P','Daumont & Malicet','IUP','orientation','horizontal')
legend_o3abs={'B.op','IGAQ','DBM','IUP','UIPQ'};
legend(legend_o3abs,'Orientation','Horizontal');
% rotate
 xlabel(' ')
 ylabel('Ratio')
 title('Ratio to Brewer B&P ');
 
%%


 
%%
legend_o3abs={'B.op','IGAQ','DBM','IUP','UIPQ'};
rp=round(10000*matdiv(matadd(o3abs,-o3abs(1,:)),o3abs(1,:)))'/100;
displaytable([mean(rp);std(rp)],legend_o3abs,15,'.3f');

%
ratio=round(10000*matdiv(o3abs,o3abs(1,:)))'/10000;
%f=fopen('table_brewer_ratio_o3abs.csv','wt');
brewer_stats=[nanmean(ratio);nanstd(ratio);max(ratio);min(ratio);range(ratio)];
printmatrix(round(brewer_stats*1000)/1000)
displaytable([nanmean(ratio);nanstd(ratio);max(ratio);min(ratio);range(ratio)],...
    legend_o3abs,20,'.4f',{'mean','std','max','min','range'},...
    1,'\t');
%fclose(f)
ratios_from_excel=[
  %  	185	157	EC set	Scarnatto	40	72	156	RBCC-E
0.9695	0.969	0.97	0.979324604	0.938547486	0.91689008	0.932795699	
0.978479197	 NaN NaN		0.975961538	1.087378641	0.947368421	1.00872093	
]

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