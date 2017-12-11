%laser ozone mode99999u
J={};
addpath(genpath('~/CODE/rbcce.aemet.es/iberonesia/matlab')) 
addpath(genpath('matlab')) 
 path(fullfile(pwd,'matlab'),path)
 set(0,'defaultfigurecolor',[1 1 1])
 
%% 
for i=1:10
   J{i}=readb_j3dk('bdata185/B02716.185',sprintf('y%d',i));
   if i==10
   J{i}=readb_j3dk('bdata185/B02716.185','ya');
   end
   J{i}(2,:)=J{i}(2,:)/24/60+datenum(2016,0,27);
   if i>8
       J{i}(2,:)=J{i}(2,:)+1;
   end
end

for i=11:24    
    ch=char(87+i);
    sprintf('y%c',ch)
    J{i}=readb_j3dk('bdata185/B02816.185',sprintf('y%c',ch));
    J{i}(2,:)=J{i}(2,:)/24/60+datenum(2016,0,28);
end

% ozone
n=1;
for i=25:30   
    
    J{i}=readb_j3dk('bdata185/B02616.185',sprintf('j%d',n));
    J{i}(2,:)=J{i}(2,:)/24/60+datenum(2016,0,26);
    n=n+1;
    J{i}(16,:)=1020;  %calc_step
end
J{27}=J{27}(:,end-76:end); 

% ozone 2
n=1;
for i=31:35   
    
    J{i}=readb_j3dk('bdata185/B02716.185',sprintf('h%d',n));
    J{i}(2,:)=J{i}(2,:)/24/60+datenum(2016,0,27);
    n=n+1;
    %J{i}(16,:)=1020;  %calc_step
end
J{31}=J{31}(:,77:end); 

% revisar bad mic
%J{i+1}=readb_j3('bdata185/B02716.185',sprintf('j%d',3));
%J{i+1}(2,:)=J{i+1}(2,:)/24/60+datenum(2016,0,27);


%% measurement check
s=cellfun(@(x) size(x,2)',J)
%%


% j2 repetido
aux1=J{2}(:,1:77);
aux2=J{2}(:,78:end);
J{2}=aux1;
J{end+1}=aux2;
% elimianmos las pruebas
J{4}=J{4}(:,1+93-77:end);
%% medidas repetidas
for x=21:24
x2=J{x};
x3=x2(:,1:77);
x4=x2(:,end-76:end);
J{x}=x3;
J{end+1}=x4;
end
s=cellfun(@(x) size(x,2)',J)
m=cat(2,J{:});
% 23 columnas, 77 wavelengns, 29 measurements

sl=reshape(m,23,77,[]);
timex=squeeze(mean(sl(2,:,:)));
timex(8)
times=squeeze(std(sl(2,:,:)));
%% examples
f=figure
h={};

ha = tight_subplot(2,3,[.05 .05],[.1 .2],[.1 .1])
            %for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
            %set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
j_o3=find(sl(16,1,:)==1020)
for i=1:6
%figure
%subplot(2,3,i)
axes(ha(i));
h{i}=semilogyy(sl([15,6:12],:,j_o3(i))')
%legend(mmcellstr(sprintf('slit #%d |',[0:6]')));
grid on
if i>3 xlabel('wavelength (A)'); end
if mod(i,3)==1 ylabel('Brewer raw counts'); end
set(h{i}(2),'LineWidth',2,'Marker','x');
set(ha(i),'Ytick',[10         100        1000       10000]);
 axis('tight')
end
%set(ha(1:4),'XTickLabel','');
%set(ha(2:2:6),'YTickLabel','')
l1=legend(h{1},mmcellstr(sprintf('slit #%d |',[0:6]')),'Orientation','Horizontal');
set(l1,'Position',[0.077608      0.84228      0.83661     0.040429]);
suptitle(sprintf('Brewer ozone mode raw counts  mic=%d',1020))

set(f,'tag','laser_log');
printfiles_report(gcf,'figures','Width',18,'Height',18);

%% sort by time
[i,j]=sort(timex);
datestr(timex(j))
Js=J(j);
m=cat(2,Js{:});
sl=reshape(m,23,77,[]);
timex=squeeze(mean(sl(2,:,:)));
%% micrometer set
mic_=unique(sl(16,:,:))

%% Raw counts=
figure
semilogy(m(15,:),m(6:12,:),'.')
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer Raw counts')
title('Brewer ozone mode')

%% & Brewer counts second
figure
%norf=max(max(m(17:end,:)));
norf=max(m(17:end,:)');
semilogy(m(15,:),matdiv(m(17:end,:),norf'),'.-')
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  (counts/sec)')
yyaxis right
plot(m(15,:),m(7,:),'k-')

hold on
%% no linearity correction
nl=m(17:end,:);
for i=1:7
  nl(i,:)=monitor_brewer(nl(i,:));
end
dk=m(18,:)./monitor_brewer(m(18,:));

lc=m(17:end,:)./nl;
lc=matadd(lc,-dk);
%norf=max(max(lc));
norf=max(lc');
figure
r=sortrows([m(15,:)',matdiv(lc',norf)],1);
ly=semilogy(m(15,:)',matdiv(lc',norf),'.');
set(gca,'Ylim',[1E-8,1])
 set(ly(2),'Marker','x')
grid
xlabel('wavelength (set)')
ylabel('Brewer normalized  counts/sec')
title('Brewer ozone mode corrected  (counts/sec)')
l=legend(cellstr(num2str([0:6]')))
yyaxis right
p=plot(m(15,:),m(7,:),'k:.','MarkerSize',10);
ylabel('Dark raw counts')
l.String{end}='Dark'

% linearity correction
%% for counts/sec >12000
 for i=1:7
     for j=1:size(sl,3)
       %dk=sl(7,:,j)./monitor_brewer(sl(7,:,j));  
       sl(16+i,:,j)= sl(16+i,:,j)./monitor_brewer(sl(16+i,:,j))/norf(i);
     end
 end

%%
%% ozone cross section calculation
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SGQ','SG1','BPB'}

O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

%% examples


%% wv, fwhm per mic position
wv=[]; fwhm=[]; l=[];o3xs_m=[];o3xs_p=[];
for i=1:length(mic_)
    fg(i)=figure;
    n=1
   for j=[1,3:7]
    subplot(2,3,n);
    hold on
    x1=sl([15,16+j],sl(16,:,:)==mic_(i))';
    t1=sl(2,sl(16,:,:)==mic_(i))';
    %ploty(x1);
    % only peak 
    [ix,ij]=max(x1(:,2));
    tmax(i,j)=t1(ij);
    x2=x1(ij-34:ij+34,:);  % some errors on the file
    plot(x2(:,1),x2(:,2)/max(x2(:,2)));
    [wv(i,j),fwhm(i,j)]=slit_fit(x2,0,1,fg(i));
    
    [wv_(i,j),fwhm_(i,j)]=slit_fit(x2,1,0,fg(i));
    o3xs_m(i,j,:)=xcs_cal(x2,o3_set);
    o3xs_p(i,j,:)=xcs_cal([wv(i,j),fwhm(i,j)],o3_set);
    n=n+1;
    hold on
    %title(datestr(tmax(i,j)));
    title([datestr(tmax(i,j),'HH:MM'),sprintf(' wv  %.1f  fwhm %.1f',wv(i,j),fwhm(i,j))]);
    line=[tmax(i,j),mic_(i),j-1,wv(i,j),fwhm(i,j),o3xs_m(i,j,4),o3xs_p(i,j,4)]
    l=[l;line];
    box on;
    grid on;
    xlabel('wavelength')
    axis('tight')
    set(gca,'Ylim',[0,1.2]);
   end
   %suptitle(num2str(mic_(i)));
end
%    

set(fg(1),'Tag','Laser_scan_dsp')
printfiles_report(fg(1),'figures','Width',20,'Height',12)



%%
wv_=[]; fwhm_=[]; mic=[];slit=[];o3xsm=[];o3xsp=[];
k=0;
sl(18,:,:)=0;
for i=1:size(sl,3)
    nr=max(max(sl([17:end],:,i)));
    for j=[1,3:7]    
     if max(sl([16+j],:,i))/nr>0.3  % normalized
          if mod(k,6)==0
           f=figure;
           subplot(2,3,1);
           hold on;
           k=1;
          else
            k=k+1;
            subplot(2,3,k);
          end
         
          [wv_(i,j),fwhm_(i,j)]=slit_fit(sl([15,16+j],:,i)',0,1,f);
          mic(i,j)=unique(sl(16,:,i));
          slit(i,j)=j;
          o3xsm(i,j,:)=xcs_cal(sl([15,16+j],:,i)',o3_set);
          o3xsp(i,j,:)=xcs_cal([wv_(i,j),fwhm_(i,j)],o3_set);
          title([datestr(timex(i),'HH:MM'),sprintf('  %d %d %.1f  %.1f',slit(i,j),mic(i,j),wv_(i,j),fwhm_(i,j))]);
          axis('tight')
     else
          wv_(i,j)=0;
          fwhm_(i,j)=0;
     end  
 
    end
end
fwhm_(isnan(fwhm_))=0;

%%
laserscan=[timex,sum(mic,2),sum(slit,2)-1,sum(wv_,2),sum(fwhm_,2),sum(o3xsm(:,:,4),2),sum(o3xsp(:,:,4),2)];
save laserscan laserscan
o3_laser=sortrows(laserscan(laserscan(:,2)==1020,:),3);
save o3_laser o3_laser
%% PTB comparison
ptb_scan =[
       3100.6          5.7
       3135.1          6.1
       3168.1            6
       3200.1          5.9
       3031.9          6.1
         3063          6.1
       3100.5          5.8
       3133.4            6
       3163.5          5.9
       3199.8          5.6
       3233.1          5.9
       3264.9          5.7
       3295.7          5.5
       3109.3            6
       3139.7          5.8
       3139.7            6
       3176.4          5.7
       3209.9            6
       3242.1          5.8
       3273.2          5.7
       3255.8          5.9
       3284.7          5.7
       3319.5          5.4
       3351.2          5.7
       3381.6          5.5
       3410.9          5.5
         3393          5.5
       3420.5          5.5
       3453.4          5.3
       3483.5          5.5
         3512          5.3
       3539.5          5.2
       3490.7          5.4
       3517.1          5.4
       3548.7          5.1
       3577.6          5.3
       3604.9          5.1
       3631.1            5
       3548.7          5.2
       3577.6          5.3
       3604.9          5.1
       3631.1          5.1];
   
   %% Falta una medida ??
   ptb_scan(7,:)=[];
   ptb_scan(12,:)=[];
o3x_ptb=[];   
for i=1:size(ptb_scan,1)
    o3x_ptb(i,1:7)=xcs_cal([ptb_scan(i,1),ptb_scan(i,2)],o3_set);
    ptb_scan(i,3)=o3x_ptb(i,4);
end
   
   % round to 0.1 Amstrong
%laserscan=[timex,round(sum(wv_,2),1),round(sum(fwhm_,2),1)];

j=abs(laserscan(:,4)-ptb_scan(:,1))<0.2; %outliers
figure
plot(laserscan(j,5)-ptb_scan(j,2))
figure
plot(laserscan(j,4)-ptb_scan(j,1))

figure
plot(laserscan(:,4),laserscan(:,6)./ptb_scan(:,3),'o')
grid
title('Efective Serduchenko ozone cross section difference RBCC-E vs PTB')

%%
lscan=[laserscan,ptb_scan];
save lscan lscan


%% Ozone calculation
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];

ozone_slits=sortrows([laserscan(1:6,:),ptb_scan(1:6,:)],3)
o3abs=-ozone_slits(:,[6:7,end])'*O3W'/log(10)

printmatrix(ozone_slits,5)
printmatrix(o3abs,5)
printmatrix(o3abs/o3abs(2))



