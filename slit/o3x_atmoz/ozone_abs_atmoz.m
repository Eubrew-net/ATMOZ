%% Ozone absortion cross section for ATMOZ
% 
% 
%Original files from:
%http://igaco-o3.fmi.fi/ACSO/cross_sections.html  
%
% Three Bass &Paur cross sections:
% Brewer Operative (BOp) : 
%   Is the B&P at 228K without any adjustment.  
%   Not available at IGACO web page but at -Mainz Spectral Database-.
%  
%%
% IGACO B&P IGQ4 
%  quadracit coeffients on the file  on the file "Bp.par" 
%
%  The Individual temperatures files from ACSO web page are not used 
%  This files do not agree with B&P paper and do not include -45 C set. 
%  This dataset appears to be spectra at selected temperatures calculated 
%  from the polynomial fitted to the original data excluding 218K (Weber 2011),
% 
%%     
%Bernhard 2005 (B05),
%  B&P corrected by temperature dependence and extended 
%  with Molina & Molina cross section. 
%  Includes dynamical effects convolving the XS with temperature profile 
%  and ozone profile and top of the atmosphere solar spectrum.
%%  
%Daumont, Brion & Malicet (DBM),
%   Measured at 5 temperatures, including 228K 
%%   
%Serdychenco et all (IUP), 
%  The newly determined data set from Bremen University, 
%  Institute for e Enviromental Physics (IUP) are also available at IGACO 
%  with ten temperatures files and the quadratic fit (IUPQ).
%
%
% The cross section converted to (atm cm)-1
% The wavelengs are converted  to air

% original files

% Bass & Paur Brewer (operative)
%bpb.file
% IGCACO B&P
%bbi.file
% Bernhard 2005 (B&P)
%b05.file
% Daumont Brion Malicet
%dbm.file
% UIP Serdyuchenco Groseleb Weber
%sgw.file




%%
T0=273.15
L=2.69e19;% molecules per cm2 -> cm a (atm cm)-1
%
load slit_dob83_brw  % brewer and dobson slits





%% temperature coefficients  by serdyuchenco
% sigma(T,lambda)=10^(-20)*a_0(lambda)*[1+a_1(lambda)*T+a_2(lambda)*T^2]
% lambda: wavelength in nm
% T: temperature in degree Celsius
% sigma: ozone cross section in cm^2
% vacuum to air converted with va2air (libradtran)
load serdyunchenkoo3temp.air
st=serdyunchenkoo3temp;
o3abs_temp=293:-10:193;
sdt.file='serdyunchenkoo3temp.air';
sdt.lamda=st(:,1);
sdt.tempc=matmul(st(:,2),[ones(size(st(:,1))),st(:,3:4)]);
sdt.temp=[o3abs_temp,228];
o3sdt=polyvac(sdt.tempc(:,3:-1:1)',(sdt.temp-273.15));
sdt.o3x=o3sdt'*L*1E-20;
sdt.t_45=[sdt.lamda,L*1E-20*polyvac(sdt.tempc(:,3:-1:1)',-45)'];
save sdt sdt
%%
%% Serdyuchenco version 0713 (Fourier filtered)
%
load SGW0713.air.dat;
sgw7_air=SGW0713_air;
L=2.69e19;% molecules per cm2 -> cm a (atm cm)-1
sgw7_air(:,2:end)=sgw7_air(:,2:end)*L;
sgw7_temp=293:-10:193;
% interpolate to -45
o3sgw7_i45=[sgw7_air(:,1),interp1(sgw7_temp,sgw7_air(:,2:end)',228.15)'];
% interpolate to -46.13
o3sgw7_i46=[sgw7_air(:,1),interp1(sgw7_temp,sgw7_air(:,2:end)',227)'];


P_sgw7=polyfic(sgw7_temp-T0,sgw7_air(:,2:end)',2);
%o3uip_45=[o3abs_air(:,1),interp1(o3abs_temp,o3abs_air(:,2:end)',228.15)'];
o3sgw7_45=[sgw7_air(:,1),polyvac(P_sgw7,228.15-T0)'];
o3sgw7_46=[sgw7_air(:,1),polyvac(P_sgw7,227-T0)'];




% o3 temperature
sgw7.lamda=sgw7_air(:,1);
sgw7.o3x=[sgw7_air(:,2:end),o3sgw7_45(:,2)];
sgw7.temp=[sgw7_temp,228];
sgw7.t_45=o3sgw7_45;
sgw7.t_46=o3sgw7_46;
sgw7.ti_45=o3sgw7_i45;
sgw7.ti_46=o3sgw7_i46;
sgw7.q=P_sgw7;
sgw7.file='SGW0713_air.dat'

save sgw7 sgw7

%%
%%
%% gorshelev version 0917 (ATMOZ)
%
load gorshelev_201709_v3.air;
gw7_air=gorshelev_201709_v3;
L=2.69e19;% molecules per cm2 -> cm a (atm cm)-1
gw7_air(:,2:end)=gw7_air(:,2:end)*L;
gw7_temp=293:-10:193;
% interpolate to -45
o3gw7_i45=[gw7_air(:,1),interp1(gw7_temp,gw7_air(:,2:end)',228.15)'];
% interpolate to -46.13
o3gw7_i46=[gw7_air(:,1),interp1(gw7_temp,gw7_air(:,2:end)',227)'];


P_gw7=polyfic(gw7_temp-T0,gw7_air(:,2:end)',2);
%o3uip_45=[o3abs_air(:,1),interp1(o3abs_temp,o3abs_air(:,2:end)',228.15)'];
o3gw7_45=[gw7_air(:,1),polyvac(P_gw7,228.15-T0)'];
o3gw7_46=[gw7_air(:,1),polyvac(P_gw7,227-T0)'];




% o3 temperature
gw7.lamda=sgw7_air(:,1);
gw7.o3x=[sgw7_air(:,2:end),o3sgw7_45(:,2)];
gw7.temp=[sgw7_temp,228];
gw7.t_45=o3sgw7_45;
gw7.t_46=o3sgw7_46;
gw7.ti_45=o3sgw7_i45;
gw7.ti_46=o3sgw7_i46;
gw7.q=P_sgw7;
gw7.file='gorshelev_201709_v3.air'

save gw7 gw7

%%



%%
figure
h=plot(gw7.lamda,sgw7.o3x,'-',slit.dob83(:,1),slit.dob83(:,2:end)*3,'.');
legend([cellstr(num2str(gw7.temp'))])
set(gca,'Xlim',[300,350])
hold on
%ploty(o3sgw7_45,'k.-.');
%set(h,'linewidth',2)
title('Bremen O_3 XS');
ylabel('(atm cm)^-1')
xlabel(' wv nm');
insetLoc = [315,340;4,10];
rectLoc = [315,320;0.5 1.5];

inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],...
    'LineWidth',1);

legend([cellstr(num2str(o3abs_temp'));'228'])

 set(gcf,'Tag',' DOBSON_IUP ');
 %printfiles_publication(gcf,'.','LineWidth','auto','Width',16,'Height',12)


%%
figure
h=plot(sdt.lamda,sdt.o3x,'-',slit.dob83(:,1),slit.dob83(:,2:end)*3,':.');
set(gca,'Xlim',[300,350])
legend([cellstr(num2str(gw7.temp'))])
hold on
ylabel('(atm cm)^-1')
xlabel(' wv nm');
insetLoc = [315,340;4,10];
rectLoc = [315,320;0.5 1.5];

inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],...
    'LineWidth',1);

legend([cellstr(num2str(o3abs_temp'));'228'])

 set(gcf,'Tag',' DOBSON_IUP_V3 ');
 %printfiles_publication(gcf,'.','LineWidth','auto','Width',16,'Height',12)


%%
%
f1=figure
h=plot(gw7.lamda,gw7.o3x,'-',...
    slit.brw(:,1),slit.brw(:,2:end)*3,'b-',...
    slit.dob83(:,1),slit.dob83(:,2:end)*3,'r-');

legend(h,cellstr(num2str(o3abs_temp')))
set(gca,'Xlim',[300,350])
hold on
title('Bremem 2017 ');
ylabel('(atm cm)^-1')
xlabel(' wv nm');
insetLoc = [315,340;4,11];
rectLoc = [310,330;0.1 3.1];

inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],...
    'LineWidth',1);

legend([cellstr(num2str(o3abs_temp'));'228'])

set(gcf,'Tag','IUP cross section and Dobson-Brewer Slits');
%printfiles_publication(gcf,'.','LineWidth','auto','Width',16,'Height',12);

%%
figure
orient portrait
[j,k]=ismember(sdt.lamda,gw7.lamda);
h=plot(gw7.lamda(k),100*(gw7.o3x(k,:)-sdt.o3x(j,:))./sdt.o3x(j,:),'-');
interactivelegend(h,cellstr(num2str(o3abs_temp')))
set(gca,'Xlim',[305,340])
set(gca,'Ylim',[-10,5])
title('IUP Residuals of the quadratic aproximation');
ylabel(' % ')
xlabel(' wv nm');
insetLoc = [310,335;-9,-4];
rectLoc = [312,315; -3.5,3.5];
legend([cellstr(num2str(gw7.temp'))])
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);

set(gcf,'Tag','Residuals of IUP quadratic approximation');
 %printfiles_publication(gcf,'.','LineWidth','scaled','Width',16,'Height',12)



%% Brewer B&P
f=load('ozxsec.dat');
L=2.69e19;% molecules per cm2 -> cm a (atm cm)-1
f(f==-9.99)=NaN;
f(:,2:end)=f(:,2:end)*L;
%k=1.38062e-23; % bolzmann
% brewer soft o3x1T=f(:,4)*1.013*1e5/(k*273.1)*1e-6; %in atm /cm
o3x1T=f(:,4);
o3x1wl=f(:,1);
brw_bp=[o3x1wl,o3x1T];
o3x_temp=[ 203.15       218.15       228.15       243.15       273.15       298.15];
% interpolate to -46
o3x_46=[f(:,1),interp1(o3x_temp,f(:,2:end)',227)'];
o3x_45=[f(:,1),o3x1T];
%
 bpb.lamda=f(:,1);
 bpb.temp=o3x_temp;
 bpb.o3x=f(:,2:end);
 bpb.t_45=o3x_45;
 bpb.file='ozxsec.dat';
 save bpb bpb
%%
figure;
hold all;
h=plot(bpb.lamda,bpb.o3x,'-');
interactivelegend(h,cellstr(num2str(bpb.temp')))
set(gca,'Xlim',[300,350])

ploty(o3x_46,'k.-.');
set(h,'linewidth',2)
legend([cellstr(num2str(o3x_temp'));'227'])
title('Brewer B&P');
ylabel('(atm cm)^-1');
xlabel(' wavelength (nm) ');
insetLoc = [315,340;4,10];
rectLoc = [312,315;1.25 2];
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);
box on;

%% B&P IGACO
load bp.par
o3bp_45=[bp(:,1),L*1E-20*(bp(:,2)+bp(:,3)*-45+bp(:,4)*-45*-45)];
o3bp_46=[bp(:,1),L*1E-20*(bp(:,2)+bp(:,3)*-46.3+bp(:,4)*-46.3*-46.3)];

%asumimos sgw7 temperaturas
bpi.file='bp.par';
bpi.lamda=bp(:,1);
bpi.tempc=bp(:,2:4);
bpi.t_45=o3bp_45;
temps=[193   203   213   218   223   228   233   243   253   263   273   283   293   295   298];
bpi.temp=temps; %'IUP temp'
%bpi.temp=sgw7.temp; %'IUP temp'

o3paur=polyvac(bp(:,4:-1:2)',(temps-273.15));
bpi.o3x=o3paur'*L*1E-20;
save bpi bpi

%%
figure;
hold all;
h=plot(bpi.lamda,bpi.o3x','-');
interactivelegend(h,cellstr(num2str(bpi.temp')))
set(gca,'Xlim',[300,350])

ploty(o3bp_45,'k.-.');
set(h,'linewidth',2)
legend([cellstr(num2str(bpi.temp'));'227'])
title('IGACO B&P');
ylabel('(atm cm)^-1')
insetLoc = [315,340;4,10];
rectLoc = [312,315;1.25 2];
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);
box on;
xlabel(' wavelength (nm) ');

%% B&P Bernhard 2005
bernhardBP=xlsread('Germar2005BP.xls',1);
BPB_46=[bernhardBP(:,1),bernhardBP(:,2)*log(10)]; % is already normalized
b05.t_45=BPB_46(2:end,:);
b05.t_46=BPB_46(2:end,:);
b05.o3x=BPB_46(2:end,2);
b05.lamda=BPB_46(2:end,1);
b05.temp=227;
b05.file='Germar2005BP.xls';


save b05 b05


%% DAUMONT
%% % here load daumont
%1e-20  1  1  0  1  5  2  218  3  228  4  243  5  273  6  295  6
%Daumont o3 cross sections for 5 temperatures, originals in omsao_o3_218k_brion.dat etc. for 218, 228, 243, 273 and 295 K
%from Daumont et al., J.Atmos.Chem.15, 145(1992) and Malicet et al., J.Atmos.Chem.21, 263(1995), data are 1e-20cm2 per molecule


[data,head]=liesfile('o3daumont5temp.abs',4,6);
%1e-20cm2 per molecule 
data(:,2:end)=data(:,2:end)*L*1e-20;
data(:,1)=data(:,1);
% 228 K -45
%[k,jk]=min(abs(daumont_xs(:,1)-o3x1wl(1)));
%hold on;plot(o3x1wl,o3x1T./max(o3x1T),'r.');
daumont_temp=[218  228  243  273  295]; 
dbm.lamda=data(:,1);
dbm.o3x=data(:,2:end);
dbm.temp=daumont_temp;
dbm.file='o3daumont5temp.abs';

daumont_45=data(:,[1,3]);
daumont_46=[dbm.lamda,interp1(dbm.temp,dbm.o3x',227)'];
dbm.t_45=daumont_45;
dbm.tempc=polyfic(dbm.temp'-T0,dbm.o3x',2);


save dbm dbm

%%  output files

t_bpb=array2table([bpb.lamda,bpb.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',bpb.temp))']);
writetable(t_bpb,'Brewer_operative_B&P');

t_bpi=array2table([bpi.lamda,bpi.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',bpi.temp))']);
writetable(t_bpi,'IGACO_B&P_qtemp');

t_b05=array2table([b05.lamda,b05.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',b05.temp))']);
writetable(t_b05,'Bernhard_2005_B&P');

t_dbm=array2table([dbm.lamda,dbm.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',dbm.temp))']);
writetable(t_dbm,'DMB_');

t_sgw=array2table([sgw7.lamda,sgw7.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',sgw7.temp))']);
writetable(t_sgw,'SGW_');

t_gw=array2table([gw7.lamda,sgw7.o3x],'VariableNames',['wv_nm',mmcellstr(sprintf('T_%.f|',gw7.temp))']);
writetable(t_sgw,'GW_');

%%
figure;
hold all;
h=plot(dbm.lamda,dbm.o3x','.');
set(gca,'Xlim',[300,350])

ploty(daumont_45,'k.-.');
set(h,'linewidth',2)
legend([cellstr(num2str(dbm.temp'));'228'])
title('Daumont Malicet Brion IGACO ');
ylabel('(atm cm)^-1')
insetLoc = [315,340;4,10];
rectLoc = [312,315;1.25 2];
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);
box on;

%% B&P set to -46
b_p_46=[];
b_p_46=scan_join(o3x_46,o3bp_46);
b_p_46=scan_join(b_p_46,BPB_46);

%%B&P set to -45
b_p_45=[];
b_p_45=scan_join(o3x_45,o3bp_45);
b_p_45=scan_join(b_p_45,BPB_46);


%% Final Set at -45 

o3_set=[];
o3_set=scan_join(o3_set,o3x_45); % Brewer
o3_set=scan_join(o3_set,o3bp_45); %igaco
o3_set=scan_join(o3_set,daumont_45); % Daumont
o3_set=scan_join(o3_set,o3sgw7_45);   %uip
o3_set=scan_join(o3_set,sdt.t_45);   %uip quad
o3_set=scan_join(o3_set,BPB_46); % Bernhard 
o3_set=scan_join(o3_set,o3gw7_45);   %uip quad
o3_set_45=o3_set;
save o3_set_45 o3_set_45;
%% Bremen 2014/2017
%
figure
jc=find(all(~isnan(o3_set(:,[1,5,8])'))');
yyaxis left
plot(o3_set(jc,1),100*matdiv(matadd(o3_set(jc,5),-o3_set(jc,8)),o3_set(jc,8)),'o:')
title('Bremen cross section at -45C Ratio to 2017');

ylabel('ratio %')
xlabel('wavelengh nm');
grid
yyaxis right
h=plot(slit.brw(:,1),.8+slit.brw(:,2:end)/5,'b-',...
    slit.dob83(:,1),.8+slit.dob83(:,2:end)/5,'r-');
%hline(1.007,'r','Bernhard to IGACO mean= 1.007')
set(gca,'Xlim',[290,350])
legend({'2014 vs 2017'})
%% Final Set at -46 

o3_set_46=[];
o3_set_46=scan_join(o3_set_46,o3x_46); % Brewer
o3_set_46=scan_join(o3_set_46,o3bp_46); %igaco
o3_set_46=scan_join(o3_set_46,daumont_46); % Daumont
o3_set_46=scan_join(o3_set_46,o3sgw7_46);   %uip ->updated !!
o3_set_46=scan_join(o3_set_46,BPB_46); % Bernhard ,BPB_46); % Bernhard 
o3_set_46=scan_join(o3_set_46,o3gw7_46); % Bremen 2017; % Bernhard 
save o3_set_46 o3_set_46  

%%



%%
figure;
h=plot(o3_set(:,1),o3_set(:,2:end),'o',...
      slit.brw(:,1),.8+slit.brw(:,2:end)/5,'b-',...
    slit.dob83(:,1),.8+slit.dob83(:,2:end)/5,'r-');
set(h(3),'MarkerSize',3,'Marker','p')
set(h(4),'MarkerSize',3,'Marker','s')
set(h(1),'Marker','x')
set(h(2),'Marker','+')
title('');
set(gca,'Xlim',[298,350])
insetLoc = [309,345;4,12];
rectLoc = [316,319;.8 1];

inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);
legend('B Op.','IGQ4','DBM','IUP','B&P B05','GW_','Orientation','Horizontal')
title('O_3 ozone cross sections T=228K ');
ylabel('(atm cm)^-1')
xlabel('nm');
set(h(7:18),'Visible','off')

set(gcf,'Tag','O3 cross section T228K');
%printfiles_publication(gcf,'.','LineWidth','auto','Width',16,'Height',12)

%%

figure;
hold all;
h=plot(b_p_46(:,1),b_p_46(:,2:end),'.',...
    slit.brw(:,1),slit.brw(:,2:end)*2,'b-',...
    slit.dob83(:,1),slit.dob83(:,2:end)*2,'r-');
interactivelegend(h,{'Brewer','IGACO','Benhard'})
set(gca,'Xlim',[300,350])

legend({'Brewer','IGACO','Benhard'})
title('Bass & Paur cross section at -46.3C?');
ylabel('(atm cm)^-1')
insetLoc = [315,340;4,9];
rectLoc = [312,315;1.25 2];
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);




%% common values
figure
jc=find(all(~isnan(b_p_46(:,2:end)'))');
plot(b_p_46(jc,1),matdiv(b_p_46(jc,[2,4]),b_p_46(jc,3)),'o:')
title('Bass & Paur cross section at -46.3C: Ratio to IGACO');
legend({'Brewer','Benhard'})
ylabel('ratio')
xlabel('wavelengh nm');
hline(1.007,'r','Bernhard to IGACO mean= 1.007')

%% Set to -45
figure;
hold all;
h=plot(b_p_45(:,1),b_p_45(:,2:end));
interactivelegend(h,{'Brewer','IGACO','Benhard'})
set(gca,'Xlim',[300,350])

legend({'Brewer','IGACO','Benhard'})
title('Bass & Paur cross section at -45C  (Bernhard -46.3)');
ylabel('(atm cm)^-1')
insetLoc = [315,340;4,9];
rectLoc = [312,314;1.4 1.9];
inset2DAbsolute(gca,insetLoc,'ZOOMN',rectLoc,'Color',[0 0 0],'LineWidth',1);


%% common values
figure
jc=find(all(~isnan(b_p_45(:,2:end)'))');
plot(b_p_45(jc,1),matdiv(b_p_45(jc,[2,4]),b_p_45(jc,3)),'o:')
title('Bass & Paur cross section at -45C? (Bernhard -46): Ratio to IGACO');
legend({'Brewer','Benhard'})
ylabel('ratio')
xlabel('wavelengh nm')


