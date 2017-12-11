%  Columns 1 through 11
%    'w_time_s'    'WL_set'    'Mon'    'Brewer'    'BrewerDark'    'Brewer_dk'    'NLK'    'Brewer_dk_NL'    'BrewNL_Mon'    'Brewermicrom'    'Peak_steps'
%  Columns 12 through 14
%    'FWHM_steps'    'slit'    'mwl'
% w_time_s (in seconds), WL_set (laser wavelength that was set), Mon (Monitor detector signal), Brewer (recorded Brewer counts), 
% BrewerDark (Brewer dark counts), Brewer_dk (Brewer-BrewerDark)/4 - divide by 4 for the compatibility of nonlineartity corrections),
% NLK (Nonlinearity correction), Brewer_dk_NL        (Brewer_dk*NLK), BrewNL_Mon (Brewer_dk_NL/Mon),       
% Brewermicrom (Brewer grating position in steps), Peak_steps (Determined peak position of the bandpass function in steps),   
% FWHM_steps  (Determined FWHM in steps),   Slit (Brewer used slit)  wavelength (as was assigned by the brewer with the current calibration)
 filename='laser_ptb/BrewerScan_20160129.txt'
 
 addpath(genpath('~/CODE/rbcce.aemet.es/iberonesia/matlab')) 
 path(fullfile(pwd,'matlab'),path)
 addpath(genpath(fullfile('matlab')));
 %filename='laser_ptb/BrewerScan_20160127a.txt'
 
 [a,b]=read_ptb_ab(filename);
 filename=strrep(filename,'_',' ');
 [C,IA,IC] = unique(a(:,[2,13]),'rows','legacy');
 [wv,wa,wb]=unique(C(:,1));
 
 %% NL correction
 % Dark corrected a convert to c/s 0.1144 is the integration time
 IT=0.1147;
 CY=1;  % cycles on DSP
 F_0=a(:,6);  % CYCLES ALREADY CORRECTED
 DT=2.9E-8;
 %Fcs= 2*matdiv(matadd(F.raw,-F.dark),CY)/IT;
 Fcs= 2*F_0/CY/IT;
 Fcs(Fcs<0)=2;
 Fcs(Fcs>1E07)=1E07;  
 % dead time correction
  F0=Fcs;
  for j=1:9
    Fcs=F0.*exp(Fcs*DT);   
  end
 
 FNLC=Fcs./monitor_brewer(Fcs)';
 a=[a,Fcs,FNLC];
 b=[b,'Brw c/s','Brw c/s (corr)'];
 b=strrep(b,'_',' ')
 %% check
 lnes=10;
 s=a(IC==lnes,:);
 sc=Fcs(IC==lnes,:);
 snl=FNLC(IC==lnes,:);
 ln=C(lnes,1);
 slits=C(lnes,2);
 x=1:length(s(:,6));
 %
 figure
 plot(s(:,6),1./monitor_brewer(s(:,6))','r-x');
 %legend('PTB','RBCC-E');
  xlabel(b(6));
  %ylabel(b(7));
 %title('RBCC-E correction');
  %legend('PTB','RBCC-E');
  title({filename, [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
  grid
  
%% choose one scan
figure
mmplotyy(x,s(:,6),'o',sc,'-x')
legend([b(6),'Brewer counts/second']);
ylabel(b(6));
mmplotyy('RBCC-E counts/secsetond');
title({[filename], [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
%%
figure
semilogy(x,s(:,6)./max(s(:,6)),'-x',x,sc./max(sc),'o')

%%
figure
mmplotyy(s(:,10),matdiv(s(:,8:9),max(s(:,8:9))),'o',s(:,6)/max(s(:,6)),'-x')
suptitle({[filename], [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
legend([b(8),b(9),b(6)]);
title('PTB analysis');
%%
figure
mmplotyy(x,matdiv(s(:,8:9),max(s(:,8:9))),'o',matdiv([sc.*s(:,7),snl],max([sc.*s(:,7),snl])),'-x')
legend([strrep(b(8:9),'_','-'),'Brewer c/s * NLK','IZO NLK']);
ylabel(b(8));
mmplotyy('RBCC-E counts/second');
suptitle({[filename], [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
%%
figure
mmplotyy(s(:,10)',matdiv([s(:,6),sc],max([s(:,6),sc])),'o',...
    matdiv([sc.*s(:,7),snl],max([sc.*s(:,7),snl])),'-x')
legend([b(6),'Brewer c/s ','Brewer c/s * NLK','Brewer NLK']);
ylabel([b(6),'  & Brewer c/s ']);
xlabel(b(10));
mmplotyy('RBCC-E counts/second');
suptitle({[filename], [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
title('PTB vs RBCC NL corr');
%%
figure
mmplotyy(s(:,10)',matdiv(s(:,8:9),max(s(:,8:9))),'o',matdiv([sc.*s(:,7),snl],max([sc.*s(:,7),snl])),'-x')
legend([b(8:9),'Brewer c/s * NLK','RBCC_E NLK']);
ylabel([b{8},' ',b{9}]);
mmplotyy('RBCC-E counts/second');
suptitle({[filename], [num2str(ln),'nm  slit=',num2str(slits),' line ',num2str(lnes)]}); 
title('PTB vs RBCC NL corr');



%% analisys

maxlines=length(wv)
dsp_arrU=zeros(maxlines,6);   % contains data for Up scan
dsp_arrD=zeros(maxlines,6);   % for DOWN scan
DSP=zeros(maxlines,6);        % mean up/dw scan
wl=zeros(maxlines,1);         % wavelength
DSPstd=zeros(maxlines,6);     % std up/dw scan
fwhm=zeros(maxlines,6);       % mean up/dw scab
fwhmstd=zeros(maxlines,6);    % std up/dw scan
backlash=zeros(maxlines,6);   % difference Up/Dw scan
resup=zeros(maxlines,6);      % resolution
scans={};
scans{maxlines,6}=[];
daycnt=1;

K=1
for i=1:length(C)
     ln=C(i,1);
     slits=C(i,2);
     lnes=wb(i);
     
     figure
     s=a(IC==i,:);
     mic=s(:,10);
     [u,iu]=unique(mic,'legacy');
     f=histc(mic,u);
     f=ismember(mic,u(f>1));
     % only up & down values
     scans{lnes,slits+1}=s(f,:);
     scan_analyzed=s(f,([10,16]));
     
     mmplotyy(s(f,10),s(f,16),'o',s(f,15),'x')
     hold on
    try
     [dsp_arrU(lnes,slits+1),dsp_arrD(lnes,slits+1),tempfw,resup(lnes,slits+1)]=dsp_brewer(scan_analyzed(:,:),ln,0.2,0.8,1);
     if isnan(tempfw)
         %ploty(scans{lnes,slits+1}(2:end-1,[10,8]),'-o')
     end
    catch
        %ploty(scans{lnes,slits+1}(2:end-1,[10,8]),'-x')
    end
     title({filename, [num2str(ln),' nm  slit=',num2str(slits),' line=',num2str(lnes)]}); 

                temp(daycnt,:)=[dsp_arrU(lnes,slits+1),dsp_arrD(lnes,slits+1)];
                tempb(daycnt,:)=tempfw;
                tempc(daycnt)=dsp_arrU(lnes,slits+1)-dsp_arrD(lnes,slits+1); % backlash

       
            %names{lnes,slits+1}=file_name;
            fwhm(lnes,slits+1)=mean(tempb(:));
            fwhmstd(lnes,slits+1)=std(tempb(:));
            DSP(lnes,slits+1)=mean(temp(:,1));
            DSPstd(lnes,slits+1)=std(temp(:));
            backlash(lnes,slits+1)=mean(tempc(:));
end

%%
IZO=struct('fwhm',fwhm,'fwhmstd',fwhmstd,'dsp',DSP,'wv',wv,'backlash',backlash,'scans',scans);
wl=wv*10;
%% standard
[fwl,fstps,pwl,pstps]=normaldsp(wl,DSP);  % quad for every slit



dsp=DSP;
dspstd=DSPstd;
FWHM=fwhm;
FWHM(FWHM==0)=NaN;
DSP(DSP==0)=NaN;
%
slits_l= mmcellstr(sprintf('slit#%d \n',1:6));
year=16;
comment='PTB';
brewnb=185;
day=026;
fpath='.';
fname=sprintf('dspnorm_%03d%02d_%s.%03d',day(1),year,comment,brewnb);
ozone_pos=1020;
mic_zero=1733;



%

fnamealldsp=sprintf('ozonedsp_%03d%02d_%s.%03d',day(1),year,comment,brewnb);
save(fullfile(fpath,fnamealldsp),'wl','dsp','dspstd','fwhm','fwhmstd','backlash','fwl','fstps','pwl','pstps');

data=pwl(:,end:-1:1)';
mydata=data(:);
savefmt(fullfile(fpath,fname),[mydata(4:end);mydata(1:3)],'','%.7e');
DSP_QUAD=[mydata(4:end);mydata(1:3)];

disp('Use polyval(pwl(2,:),wl) for calculating normal wavelengths')

% now calculate ozonecoeffs
outfname_cuadratic=sprintf('opos%03d%02d_%s.%03d',day(1),year,comment,brewnb);
fprintf('Saving ozonecoeffs to %s\n',fullfile(fpath,outfname_cuadratic));

CALC_STEP=1020; ozonepos=1020;%fnameN='O3f30012.185';
[res,detail,salida]=ozonecoeff3(fullfile(fpath,fnamealldsp),[ozonepos+mic_zero,ozonepos],fullfile(fpath,fname),fullfile(fpath,outfname_cuadratic));

jcal=find((res(:,1)==CALC_STEP),1);
jumk=find((res(:,1)==CALC_STEP),1,'last');

% ozone wv and fwhm
detail(1:2,:,:)


%%
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SGQ','SG1','BPB'}

O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];



%% QUAD
mic_range=[1020
        2087
        2423
        4174
        6213
        7711];
    sitn=[0,2:6]
bs=[];    
for i=1:length(mic_range)
  for j=1:6 
    o3xs_q=xcs_cal(detail(1:2,j,i)',o3_set);
    l=[mic_range(i),sitn(j),detail(1,j,i),detail(2,j,i),o3xs_q];
    bs=[bs;l];
  end
end
%% CUBIC

fname7=sprintf('dsp_%03d%02d_%s.%03d',day(1),year,comment,brewnb);
fprintf('saving powfiu7 to %s\n',fullfile(fpath,fname7));
polynb=3;
%dsp(~dsp)=NaN;
dsp(isnan(dsp))=0;
nslit=min(sum(~isnan(dsp)));
if nslit<=polynb
   polynb=nslit-1;
end
pos0=powfiu7;
if nanmean(fwhm(fwhm(:,1)~=0,1))<65,pos0=powfiu7(1);end

[fstd,slitpos,rmsf]=powfiu7(pos0,wl,dsp,polynb,[],[],[],fullfile(fpath,fname7));

f2=figure;
plot(wl,fstd(:,:,end),'o-'); set(gca,'FontSize',11,'FontWeight','bold');
xlabel('wavelength ')
ylabel('residual (A)')
title('Brewer Scan laser lines, cubic');
grid
set(f2,'Tag','brewer_scan_cubic_residual')
 printfiles_report(f2,'figures','Width',12,'Height',12)
disp('Use brstps2 to calculate steps and wavelengths');


stps17=brstps2(wl,1,[],fullfile(fpath,fname7));

outfname_cubic=sprintf('opos_pow7_%03d%02d_%s.%03d',day(1),year,comment,brewnb);

[res2,detail2,salida2]=ozonecoeff3(fullfile(fpath,fnamealldsp),[ozonepos+mic_zero,ozonepos],...
    fullfile(fpath,fname7),fullfile(fpath,outfname_cubic));

mic_range=[1020
        2087
        2423
        4174
        6213
        7711];
    sitn=[0,2:6];
bc=[];    
for i=1:length(mic_range)
  for j=1:6 
    o3xs_c=xcs_cal(detail2(1:2,j,i)',o3_set);
    l=[mic_range(i),sitn(j),detail2(1,j,i),detail2(2,j,i),o3xs_c];
    bc=[bc;l];
  end
end
%%
brw_scan_cubic=bc;
brw_scan_quad=bs;

save brw_scan_cubic brw_scan_cubic
save brw_scan_quad brw_scan_quad



%%
% figure;
% 
% plot(tst(:,4),tst(:,4)-tst(:,12))
% hold on
% plot(tst(:,4),tst(:,4)-tst(:,8))
% legend('cuadratic','cubic')
% grid
% title('Central Wavelength Laser Scan vs Brewer Scan')