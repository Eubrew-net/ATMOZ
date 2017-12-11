%use local files
 addpath(genpath('~/CODE/rbcce.aemet.es/iberonesia/matlab')) 
 path(fullfile(pwd,'matlab'),path)

%% all:dsp

day=020:28;
year=16;
brew=185;
lines=1:20;
minslit=0;
maxslit=6;
dspp='DSP/185_16_020/';
usewl=[];

% lamparas Hg Cd
[wl,DSP,DSPstd,fwhm,fwhmstd,backlash,resup,filenames]=alldsp_dev(day,year,brew,lines,minslit,maxslit,usewl,dspp);

DSP(DSP==0)=NaN
DSPstd(DSPstd==0)=NaN
fwhm(fwhm==0)=NaN
fwhmstd(fwhmstd==0)=NaN

[fwl,fstps,pwl,pstps]=normaldsp(wl,DSP);  % quad for every slit
[fwlc,fstpsc,pwlc,pstpsc]=cubicdsp(wl,DSP);  % quad for every slit

%%
f1=figure
plot(wl,fwlc,'-'); hold on;
plot(wl,fwl,':');
title('Residuals quadratic (doted) & cubic dispersion relations');
grid
xlabel('wavelength [A] ');
ylabel('residual [A]')
set(f1,'Tag','Lamp_residuals');
printfiles_report(f1,'figures','Width',18,'Height',18)

dsp=DSP;dspstd=DSPstd;

FWHM=fwhm;
FWHM(FWHM==0)=NaN;
DSP(DSP==0)=NaN;
%%
slits_l= mmcellstr(sprintf('slit#%d \n',1:6));
year=16;
comment='LAMP';
brewnb=185;
day=020;
fpath='.';
fname=sprintf('dspnorm_%03d%02d_%s.%03d',day(1),year,comment,brewnb)
ozone_pos=1020;
mic_zero=1733;

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
CALC_STEP=1020; ozonepos=1020;
[res,detail,salida]=ozonecoeff3(fullfile(fpath,fnamealldsp),[ozonepos+mic_zero,ozonepos],fullfile(fpath,fname),fullfile(fpath,outfname_cuadratic));

%%
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SGQ','SG1','BPB'}

O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];
mic_range=[1020
        2087
        2423
        4174
        6213
        7711];
    sitn=[0,2:6]
lamp=[];    
for i=1:length(mic_range)
  for j=1:6 
    o3xs_q=xcs_cal(detail(1:2,j,i)',o3_set);
    l=[mic_range(i),sitn(j),detail(1,j,i),detail(2,j,i),o3xs_q]
    lamp=[lamp;l];
  end
end
lamp_q=lamp;
save lamp_q lamp_q
 
%% cubic
fname7=sprintf('dsp_%03d%02d_%s.%03d',day(1),year,comment,brewnb);
fprintf('saving powfiu7 to %s\n',fullfile(fpath,fname7));
polynb=3;
dsp(isnan(dsp))=0;
%dsp(~dsp)=NaN;
nslit=min(sum(~isnan(dsp)));
if nslit<=polynb
   polynb=nslit-1;
end
pos0=powfiu7;
if nanmean(fwhm(fwhm(:,1)~=0,1))<65,pos0=powfiu7(1);end
dsp(isnan(dsp))=0;

[fstd,slitpos,rmsf]=powfiu7(pos0,wl,dsp,polynb,[],[],[],fullfile(fpath,fname7));
%%
%f2=figure;
%plot(wl,fstd(:,:,end),'o-'); 
f2=figure;
plot(wl,fstd(:,:,end),'o-'); set(gca,'FontSize',11,'FontWeight','bold');
xlabel('wavelength ')
ylabel('residual (A)')
title('Brewer Scan laser lines, cubic');
grid
set(f2,'Tag','brewer_scan_cubic_residual')
printfiles_report(f2,'figures','Width',12,'Height',12)

%%




disp('Use brstps2 to calculate steps and wavelengths');


outfname_cubic=sprintf('opos_pow7_%03d%02d_%s.%03d',day(1),year,comment,brewnb);

[res2,detail2,salida2]=ozonecoeff3(fullfile(fpath,fnamealldsp),[ozonepos+mic_zero,ozonepos],...
    fullfile(fpath,fname7),fullfile(fpath,outfname_cubic));

mic_range=[1020
        2087
        2423
        4174
        6213
        7711];
    sitn=[0,2:6]
lamp_c=[];    
for i=1:length(mic_range)
  for j=1:6 
    o3xs_c=xcs_cal(detail2(1:2,j,i)',o3_set);
    l=[mic_range(i),sitn(j),detail2(1,j,i),detail2(2,j,i),o3xs_c];
    lamp_c=[lamp_c;l];
  end
end
save lamp_c lamp_c
%[wl,DSP_,DSPstd_,fwhm_,fwhmstd_,backlash_,resup_,filenames_]=alldsp_dev_nl(day,year,brew,lines,minslit,maxslit,usewl,dspp);
%[wl,DSP,DSPstd,fwhm,fwhmstd,backlash,resup]=alldsp_dev(day,year,brew,lines,minslit,maxslit,usewl,dspp);

%% 
dspp='DSP/185_16_021/';
