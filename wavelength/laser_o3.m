%laser ozone mode
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
addpath(genpath('matlab'))
addpath(genpath('/Users/aredondas/Dropbox/ATMOZ/PTB/ptb/matlab'));
%startup
%
[j1]=readb_j3dk('bdata185/B02616.185','j1');
[j2]=readb_j3dk('bdata185/B02616.185','j2');
[j3]=readb_j3dk('bdata185/B02616.185','j3');
[j4]=readb_j3dk('bdata185/B02616.185','j4');
[j5]=readb_j3dk('bdata185/B02616.185','j5');
[j6]=readb_j3dk('bdata185/B02616.185','j6');
j3=j3(:,111:end); % elimianmos las pruebas

m=cat(3,j1,j2,j3,j4,j5,j6);
% 6:12  raw counts for the 7 slits
% 15 waveleng registered by lasser
% 17:23 counts/sec dark corrected
%
%% data ozone mode
%w_time,  set_wavelength, Monitor_signal, Brewer, Brewer_dark, Brewer_dark_corr, Nonlin_corr, Brewer_nonlin_corr,Brewer_nonlin_corr/Mon. 
% date, time
% fecha=datenum(2015,1,26)+m(2,:)/60/24;
% wv=m(15,:);
% mon=wv*NaN;
% brw=m(6:12,:);
% brw_cs=m(17:23,:);
% for i=1:7
% %brwl=[brw(:,1)';j3(20,:)./monitor_brewer(j3(20,:))]';
% brw_lc(:,i)=m(16+i,:)./monitor_brewer(m(16+i,:));
% end
% 
% 
%sli=[0,2:6];
%
% for j=1:6
% figure
%  ha = tight_subplot(3,2,[.03 .1],[.2 .01],[.1 .1])
%     for i=1:6
%     jn=m(:,:,j);
%     axes(ha(i))
%     %subplot(3,2,i); 
%     semilogy(jn(15,:),jn([sli(i)+6,7],:),'-x');
%     set(ha(i),'Ylim',[0,1E6]);
%     grid;
%     set(ha(1:4),'XTickLabel','');
%     %set(ha(3:2:6),'YTickLabel','')
% end
% suptitle(sprintf('Brewer raw counts Exp #%d',j));
% end
%%


figure
semilogy(m(15,:),m(6:12,:),'-')
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer Raw counts')
title('Brewer ozone mode, (Raw counts)')



% Brewer counts second
figure
semilogy(m(15,:),m(17:end,:),'-')
legend(cellstr(num2str([0,2:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode (counts/sec)');

%% Brewer counts second corrected
figure
r=m(17:end,:)./monitor_brewer_m(m(17:end,:));

plot(m(15,:),r./max(r')','-');
legend(cellstr(num2str([0,2:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  Linearity Corrected (counts/sec)')


%% dark count ?
figure
slit=0
semilogy(m(15,:),m(6+slit,:),'r-o')
hold on
semilogy(m(15,:),m(18,:),'k-')
semilogy(m(15,:),m(17+slit,:),'-.','linewidth',2)

dk=m(18,:)./monitor_brewer(m(18,:));
dk_c=(m(17+slit,:)./monitor_brewer(m(17+slit,:)))-dk;

semilogy(m(15,:),dk_c,'b')


xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  Slit #0')
legend('raw counts','dark','counts/second','nl corrected')
grid

set(gcf,'Tag','Laser_Brewer_dark')
printfiles_report(gcf,'.');

%%
%% dark count ?
figure
slit=2
semilogy(m(15,:),m(6+slit,:),'r-o')
hold on
semilogy(m(15,:),m(18,:),'k-')
semilogy(m(15,:),m(17+slit,:),'-.','linewidth',2)

dk=m(18,:)./monitor_brewer(m(18,:));
dk_c=(m(17+slit,:)./monitor_brewer(m(17+slit,:)))-dk;

semilogy(m(15,:),dk_c,'b')


xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  Slit #')
legend('raw counts','dark','counts/second','nl corrected')
grid

%set(gcf,'Tag','Laser_Brewer_dark')
%printfiles_report(gcf,'.');

%%
slit=4
figure
semilogy(m(15,:),m(6+slit,:)/max(m(6+slit,:)),'r-o')
hold on
semilogy(m(15,:),m(7,:)/max(m(7,:)),'k-')
semilogy(m(15,:),m(17+slit,:)/max(m(17+slit,:)),'-.','linewidth',2)

dk=m(18,:)./monitor_brewer(m(18,:));
c_dk=(m(17+slit,:)./monitor_brewer(m(17+slit,:)))-dk;
semilogy(m(15,:),c_dk/max(c_dk),'b')
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  Slit #')
grid on
legend('raw counts','dark','counts/second','nl corrected')
%% linearity correction
% Saulius Data
%% data ozone mode
%w_time,  set_wavelength, Monitor_signal, Brewer, Brewer_dark, Brewer_dark_corr, Nonlin_corr, Brewer_nonlin_corr,Brewer_nonlin_corr/Mon. 
% date, time
load laser_ptb/LaserScan_310nm.txt;
opo310=LaserScan_310nm;
brw=[j3(15,:)/10;j3(14,:)]';
%brwn=[brw(:,1),brw(:,2)./j3(5,:)'];
brwc=[brw(:,1),j3(20,:)'];
brwl=[brw(:,1)';j3(20,:)./monitor_brewer(j3(20,:))]';
f=figure
opo=opo310;
%semilogy(brw(:,1),[brwc(:,2)/max(brwc(:,2)),brwl(:,2)/max(brwl(:,2))],'-o')
plot(brw(:,1),[brwc(:,2)/max(brwc(:,2)),brwl(:,2)/max(brwl(:,2))],'-o')

hold all
%semilogy(opo(:,1),[opo(:,2)/max(opo(:,2)),opo(:,3)/max(opo(:,3))],':p')
plot(opo(:,1),[opo(:,2)/max(opo(:,2)),opo(:,3)/max(opo(:,3))],':p')
hold on
legend('brw c/s','brw nl correcred','ptb c/cy','ptb nl')
grid;
title('Laser scan 310 nm');

% wavelength it looks pretty simmilar  % CENTER OF MASS
%  4 set saulious/alberto uncorrected/corrected
%  
f=double(f)
cwv=[];fwhm=[];
[cwv(1),fwhm(1)]=slit_fit([opo(:,1),opo(:,2)],4,1,f);
[cwv(2),fwhm(2)]=slit_fit([opo(:,1),opo(:,3)],4,1,f);
[cwv(3),fwhm(3)]=slit_fit([brw(:,1),brwc(:,2)],4,1,f);
[cwv(4),fwhm(4)]=slit_fit([brw(:,1),brwl(:,2)],4,1,f);

printmatrix([cwv,cwv-cwv(2)],4)
printmatrix([fwhm,fwhm-fwhm(2)],4)
% there is small differences on central wavelengs 0.001 nm and bigger on
% fwhm ~.02
%%
brw=[j3(15,:)/10;j3(14,:)]';
brwc=[brw(:,1),j3(20,:)'];
brwl=[brw(:,1)';j3(20,:)./monitor_brewer(j3(20,:))]';



% wavelength it looks pretty simmilar  % CENTER OF MASS
%  4 set saulious/alberto uncorrected/corrected
%  
f=double(f)
cwv=[];fwhm=[];
for i=1:4
[cwv(i,1),fwhm(i,1)]=slit_fit([brw(:,1),brwc(:,2)],i,0,f);
[cwv(i,2),fwhm(i,2)]=slit_fit([brw(:,1),brwl(:,2)],i,0,f);
end

printmatrix([cwv,cwv-cwv(4)],4)
printmatrix([fwhm,fwhm-fwhm(4)],4)

f=figure
%semilogy(brw(:,1),[brwc(:,2)/max(brwc(:,2)),brwl(:,2)/max(brwl(:,2))],'-o')
h=plot(brw(:,1),[brwc(:,2)/max(brwc(:,2)),brwl(:,2)/max(brwl(:,2))],'-o')
title('Laser scan 310 nm');
legend(h,'brw c/s','brw nl correcred')
grid;
xlabel('wavelength (nm)')
ylabel('normalized counts/second')
set(f,'Tag','Corrected vs uncorrected')
printfiles_report(f,'figures');




%% Figure 
figure
y=trapezoid_brewer(opo(:,1),cwv(2),fwhm(1,2)/2);
plot(opo(:,1),opo(:,3)/max(opo(:,3)),opo(:,1),y)
legend('Laser Measurment','Brewer Slit');



%%
% figure
% semilogy(m(15,:),m(17,:)/max(m(17,:)),'-.')
% hold on
% semilogy(m(15,:),m(6,:)/max(m(6,:)),'r-o')
% %semilogy(m(15,:),m(7,:),'k-')
% %legend('counts/second','raw counts','dark')
%% slit measurements

a=squeeze(m(15,:,:)); x=a(:);
d=squeeze(m(7,:,:));  dk=d(:);% dark
z=dk./monitor_brewer(dk)';
f=figure;
f=double(f);
semilogy(x,dk/max(dk),'.'); hold on;
semilogy(x,z/max(z),'.r'); 
legend('dark','dark lc')
%%



s=[];s1=[];s0=[];wc=[];fw=[];
d=squeeze(m(7,:,:)); 
dk=d(:);
cy=squeeze(m(5,:,:));  cy=cy(:); 
a=squeeze(m(15,:,:)); x=a(:);
for i=[1,3:7]
  
   b=squeeze(m(16+i,:,:)); 
   y=b(:);% counts/second 
   c=squeeze(m(5+i,:,:)); 
   c=(c(:)-dk) ./cy; % counts/cy
   
   z=y./monitor_brewer(y)';
   z1=c./monitor_brewer_c(c)';
   
   s0=scan_join(s0,[x,c/max(c)]); % raw
   
   s=scan_join(s,[x,z/max(z)]);   % lc 1
   
   s1=scan_join(s1,[x,z1/max(z1)]);  % lc 
   
   [wc(i,1),fw(i,1)]=slit_fit([x,z],4,1,f);
   [wc(i,2),fw(i,2)]=slit_fit([x,z1/max(z1)],4);
   
end

legend('dark','c/s','raw','nl c/s','nl raw')
%% ozone cross section calculation


%%
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SGQ','SG1','BPB'}

O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];
color_s={'rx-','bo','k:','g.'}
% o3_set=[];
% o3_set=scan_join(o3_set,o3x_45); % Brewer
% o3_set=scan_join(o3_set,o3bp_45); %igaco
% o3_set=scan_join(o3_set,daumont_45); % Daumont
% o3_set=scan_join(o3_set,o3sgw_45);   %uip
% o3_set=scan_join(o3_set,sdt.t_45);   %uip quad
% o3_set=scan_join(o3_set,uip.t_45);   %uip v1
% o3_set=scan_join(o3_set,BPB_46); % Bernhard 
% Brewer standard , parametrizated slit
%  6 slits  
%  4 cross-section
%  2 samples 
%      1 Brewer standard slit
%      2 Brewer measurement
%      3 Brewer PTB(?) 
n_samples=2;

% simple aproximation (brewer method)
o3x_sum=NaN*zeros(6,4,n_samples);
o3_abs=NaN*zeros(n_samples,4);
wv(n_samples,6)=0;
fwhm(n_samples,6)=0;

for i=1  %:n_samples
    fx=figure;
    %hold all;
    %ha = tight_subplot(2,3,[.2 .2],[.1 .1],[.1 .1])
    o3x=NaN*zeros(6,4);
    o3x2=o3x;
    o3x1=o3x;
    %figure(2);
    hold all;
    for j=1:6  %slit
        sn=[1,3:7]; % slit number
         %figure
         subplot(2,3,j);hold all;
         %axes(ha(j));
         l=[];
         for k=1:4  %o3abs
             ka=~isnan(o3_set(:,k+1));
             x=o3_set(ka,1); %wv nm
             cs=o3_set(ka,k+1);  % o3xs
             %sintetic 
             y=trapezoid_brewer(x,wc(sn(j),i)/10,fw(sn(j),i)/10/2,.87);
             y0=trapezoid_brewer(x,wc(sn(j),i)/10,fw(sn(j),i)/10/2,1);
             %wv(i,j)=dsp_summ(i,3+j)/10;
             %fwhm(i,j)=dsp_summ(i,9+j)/10;
             o3x(j,k)=trapz(x,cs.*y)./trapz(x,y);
             o3x1(j,k)=trapz(x,cs.*y0)./trapz(x,y0);
             %% measured
             % slit interpolated to cross section resolution
             y1=interp1(s(:,1),s(:,j+1),x*10);
             c=[x,y1,cs];
             ki=all(~isnan(c)')';
             c=c(ki,:);
             o3x2(j,k)=trapz(c(:,1),c(:,2).*c(:,3))./trapz(c(:,1),c(:,2));
             
             yyaxis left
             plot(x,y,'b-x',x,y0,'k:.',c(:,1),c(:,2),'r:',s(:,1)/10,s(:,j+1),'ro')
             set(gca,'XLim',wc(sn(j),1)/10+[-1,1]);
             hold all
             if mod(j,3)==1 ylabel('normalized'); end
             yyaxis right
             hold all
             l(k)=plot(x,cs,color_s{k});
             if mod(j,3)==0 ylabel('(atm cm)^-1'); end
             if j>3 xlabel(' wavelength (nm)'); end
             box on;
    end
    if j==1
        [lx,xl]=legend(l,'Brewer B&P','IGACO B&P','DMB','SGW','Orientation','Horizontal');
    end
    %title(num2str(o3x(j,:)));
    100*(o3x(j,:)-o3x2(j,:))./o3x(j,:)
    try
    %set(gca,'XLim',wc(i,j)+[-1,1])
    catch
    end    
    end
    o3x_sum(:,:,1)=o3x;
    o3x_sum(:,:,2)=o3x2;
    
    o3x(isnan(o3x))=0;
    o3x1(isnan(o3x1))=0;
    o3x2(isnan(o3x2))=0;
    
    o3abs(:,i)=-(o3x')*O3W'/log(10);
    o3abs1(:,i)=-(o3x1')*O3W'/log(10);
    o3abs2(:,i)=-(o3x2')*O3W'/log(10);
    o3abs(o3abs(:,i)==0,i)=NaN;
    %suptitle([sprintf(' %s Brw %d o3abs=','PTB '),num2str(o3abs(:,i)')]);
    %printfiles_append(f1,f2,'o3abs_cal')
end

set(lx,'Position',[0.22232      0.94643      0.57768     0.040476]);
set(fx,'Tag','Laser_Brewer_ozone_mode')
printfiles_report(fx,'figures','Width',20,'Height',15,'LineWidth',1)

%figure
%plot(100*(o3abs-o3abs2)./o3abs)


table_abs=round([o3abs,o3abs1,o3abs2],4);