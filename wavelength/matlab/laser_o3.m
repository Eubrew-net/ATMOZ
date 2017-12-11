%laser ozone mode


[j1]=readb_j3('bdata185/B02616.185','j1');
[j2]=readb_j3('bdata185/B02616.185','j2');
[j3]=readb_j3('bdata185/B02616.185','j3');
[j4]=readb_j3('bdata185/B02616.185','j4');
[j5]=readb_j3('bdata185/B02616.185','j5');
[j6]=readb_j3('bdata185/B02616.185','j6');
j3=j3(:,111:end); % elimianmos las pruebas

m=cat(3,j1,j2,j3,j4,j5,j6);
%% 6:12  raw counts for the 7 slits
%% 15 waveleng registered by lasser
%% 17:23 counts/sec dark corrected


%%
figure
semilogy(m(15,:),m(6:12,:))
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer Raw counts')
title('Brewer ozone mode')

%% Brewer counts second
figure
semilogy(m(15,:),m(17:end,:))
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode Dark corrected(counts/sec)')

%% dark count ?
figure
semilogy(m(15,:),m(17,:),'-.')
hold on
semilogy(m(15,:),m(6,:),'r-o')
semilogy(m(15,:),m(7,:),'k-')
legend('counts/second','raw counts','dark')
%% linearity correction
% Saulius Data
load LaserScan_310nm.txt;
opo310=LaserScan_310nm;
brw=[j3(15,:)/10;j3(14,:)]';
brwn=[brw(:,1),brw(:,2)./j3(5,:)'];
brwc=[brw(:,1),j3(20,:)'];
brwl=[brw(:,1)';j3(20,:)./monitor_brewer(j3(20,:))]';
f=figure
opo=opo310;
semilogy(brw(:,1),[brwc(:,2)/max(brwc(:,2)),brwl(:,2)/max(brwl(:,2))],'-o')
hold all
semilogy(opo(:,1),[opo(:,2)/max(opo(:,2)),opo(:,3)/max(opo(:,3))],':p')
hold on
legend('brw c/s','brw nl correcred','s c/cy','s nl')
grid;
title('Laser scan 310 nm');

%% wavelength it looks pretty simmilar
%  4 set saulious/alberto uncorrected/corrected
%  
f=double(f)
cwv=[];fwhm=[];
[cwv(1),fwhm(1)]=slitfit([opo(:,1),opo(:,2)],4,1,f);
[cwv(2),fwhm(2)]=slitfit([opo(:,1),opo(:,3)],4,1,f);
[cwv(3),fwhm(3)]=slitfit([brw(:,1),brwc(:,2)],4,1,f);
[cwv(4),fwhm(4)]=slitfit([brw(:,1),brwl(:,2)],4,1,f);

cwv-cwv(2)
fwhm-fwhm(2)
% there is small differences on central wavelengs 0.001 nm and bigger on
% fwhm ~.02
%% Figure 
figure
y=trapezoid_brewer(opo(:,1),cwv(2)/2,fwhm(2));
plot(opo(:,1),opo(:,3)/max(opo(:,3)),opo(:,1),y)
legend('Laser Measurment','Brewer Slit');



%%
figure
semilogy(m(15,:),m(17,:)/max(m(17,:)),'-.')
hold on
semilogy(m(15,:),m(6,:)/max(m(6,:)),'r-o')
%semilogy(m(15,:),m(7,:),'k-')
%legend('counts/second','raw counts','dark')
%% slit measurements

a=squeeze(m(15,:,:)); x=a(:);
b=squeeze(m(7,:,:));  y=b(:);% dark
z=y./monitor_brewer(y)';
f=figure;f=double(f);
semilogy(x,y); hold on;
%semilogy(x,z./max(z),'.r'); iguales
s=[];s1=[];s0=[];wc=[];fw=[];
for i=1:7
   a=squeeze(m(15,:,:)); x=a(:);
   b=squeeze(m(16+i,:,:));  y=b(:);% dark corrected  
   c=squeeze(m(5+i,:,:)); c=c(:); % raw counts
   z=y./monitor_brewer(y)';
   z1=c./monitor_brewer(c)';
   s0=scan_join(s,[x,c/max(c)]);
   s=scan_join(s,[x,z/max(z)]);
   s1=scan_join(s1,[x,z1/max(z1)]);
   
   [wc(i,1),fw(i,1)]=slitfit([x,z],4,1,f);
   [wc(i,2),fw(i,2)]=slitfit([x,z1],0);
   
end
%plot(x,y);
%plot(x,c)
%plot(x,z)
%plot(x,z1)

legend('dark','c/s','raw','nl c/s','nl raw')
%% ozone cross section calculation
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SBQ','SG1','BPB'}
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
n_samples=2;

% simple aproximation (brewer method)
o3x_sum=NaN*zeros(6,4,n_samples);
o3_abs=NaN*zeros(n_samples,4);
wv(n_samples,6)=0;
fwhm(n_samples,6)=0;

for i=1:n_samples
    figure(1);
    hold all;
    o3x=NaN*zeros(6,4);
    o3x2=o3x;
    figure(2);
    hold all;
    for j=1:6  %slit
        subplot(2,3,j);hold all;
         for k=1:4  %o3abs
             ka=~isnan(o3_set(:,k+1));
             x=o3_set(ka,1); %wv nm
             cs=o3_set(ka,k+1);  % o3xs
             %sintetic 
             y=trapezoid_brewer(x,dsp_summ(i,3+j)/10,dsp_summ(i,9+j)/10,.87);
             wv(i,j)=dsp_summ(i,3+j)/10;fwhm(i,j)=dsp_summ(i,9+j)/10;
             o3x(j,k)=trapz(x,cs.*y)./trapz(x,y);
            %o3x2(j,k)=o3xsec_int([l,o],[l,s]);
            semilogy(x,cs.*y,'-')
    end
    if(j==1) legend('Brewer b&p','IGACO BP','D&M','IUP','orient','horizontal'); end
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