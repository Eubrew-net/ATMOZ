startup;
addpath(genpath(fullfile('~','CODE','iberonesia','matlab')));

%[sl_raw,TC,slraw_c]=readb_sl_rawl('bdata185/B*.185');
%sl_raw=sortrows(sl_raw,1);
%save  sl_raw186.mat  sl_raw
load sl_raw185.mat
sl_legend ={'fexha' 'hg' 'idx'	'temp'...
                      'filter1','filter2','min','sli0','slit1','cy'...
                        'raw L0'	'rL1'	'rL2'	'rL3'	'rL4'	'rL5'	'rL6'...  
                        'counts/second TC=0 cL0'	'cL1'	'cL2'	'cL3'	'cL4'	'cL5'	'cL6'...
                        'ms4'	'ms5'	'ms6'	'ms7'	'ms8'	'ms9'	'R1'	'R5'};
 
 T=sl_raw(:,4);                   
 D=sl_raw(:,1);
 F0= sl_raw(:,18:24);
 R=sl_raw(:,25:32);  %calculated using the configuration of Bfile
 
 
 %save sl_analysis
 %% Ratios calculation
% Weight definition for the seven slits
% slit 0 used for hg calibration slit 1-> dark
O3W=[  0.00      0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0.00     -1.00    0.00    0.00    4.20   -3.20];
WN=[302.1,306.3,310.1,313.5,316.8,320.1];
% MS8 SO2 ms9 o3 en el soft del brewer.
% single ratios used in brewer software
rms4=[0  -1  0  0  1  0];
rms5=[0   0 -1  0  1  0];
rms6=[0   0  0 -1  1  0];
rms7=[0   0  0  0 -1  1];
% matriz de ratios
% Ratios=F*W;0.
W=[rms4;rms5;rms6;rms7;SO2W;O3W]';
F0= sl_raw(:,18:24);
F0=F0(:,[1,3:end]); % rm dark
R0=F0*W;
%Fc=round(log10(F0)*10^4); already scaled
load cf
Fc=F0;
Fcc=F0-(T*cf);
RC=Fcc*W;

leg_ratios={'MS4','MS5','MS6','MS7','SO_2','O_3'};
leg_slits={'S0','S2','S3','S4','S5','S6'}

%%
%%
figure
mmplotyy_temp(D,matdiv(F0,F0(1,:)),T,'.',[0,50])
mmplotyy('temperature')
legend(cellstr(num2str(1:6,'%d')'))
datetick('keeplimits','keepticks')
grid
title('Brewer#185 internal lamp normalized intensity')
vline(datenum(2016,1,[18,30]),'-')
vline_v(datenum(2016,1,13),'-','Lamp Replacement')
%vline_v(datenum(2015,5,[25,45]),'k-','Huelva Campaign')
vline_v(datenum(2015,5,[25,40]),'k-')
legend(mmcellstr(sprintf('Slit #%d|',1:6)))

%%
figure
plot(D,T,'x')
title('Brewer#185 internal temperature')
vline_v(datenum(2015,12,[10,20]),'-','Laboratory')
vline_v(datenum(2016,1,[18,30]),'-','Chamber experiment')
vline_v(datenum(2015,5,[25,40]),'k-','Huelva Campaign')
datetick
ylabel('Temperature')
grid
%%
figure
time_periods=[datenum(2015,4,15),...
    datenum(2015,6,4),...
    datenum(2015,6,14),...
    datenum(2016,1,14)];
gr=group_time(D,time_periods);


gscatter(T,F0(:,2)/F0(1,2),gr,'','o')
title('Normalized ratio Slit #1')
grid








%% PERIODS
r={};s1={};c1={};rend={};
lp={'April- RBCC-E','Before RBCC-E',' June-Dec ',' Chamber'};
for i=1:4                              
%%   
if i==1    [gr,h]=group_time(D,datenum(2015,4,[15,65]));
elseif i==2 [gr,h]=group_time(D,datenum(2015,4,[15,45]));  
elseif i==3 [gr,h]=group_time(D,[datenum(2015,6,14),datenum(2016,1,14)]);
else  [gr,h]=group_time(D,[datenum(2016,1,21),datenum(2016,1,26)]);
end
idx=logical(h(:,2));


% remove outliers
[ot,j]=deleteoutliers(F0(idx,1)/F0(1,1),0.1);
f=F0(idx,:);f(j,:)=[];%f0=f(:,[1,3:end]);
t=T(idx);t(j)=[];
d=D(idx);d(j)=[];

%  & linear time fit
p_order=2;
p5=polyfit(d-d(1),mean(f,2),p_order);
p5v=polyval(p5,d-d(1));
n5=p5v/p5v(1);
f=matdiv(f,n5);
%% 

figure
hold on
[md,sd]=grpstats([d,t,matdiv(f,f(1,:))],fix(d*24)/24);
plot(md(:,1),md(:,3:end));
plot(d,matdiv(f,f(1,:)),'.');
legend(mmcellstr(sprintf('Slit #%d|',1:6)),'Orientation','horizontal');

ylabel('normalized log c/s');
xlabel('date');
grid;
set(gca,'Ylim',[0.9985,1.002] );
datetick;
title({'Brewer#185 Ozone Ratio of internal lamp',lp{i}});
%%
figure
plot(t,matdiv(f,f(1,:)),'.');
hold on
[mt,st]=grpstats([d,t,matdiv(f,f(1,:))],t);
plot(mt(:,2),mt(:,3:end))
ylabel('normalized log c/s');
xlabel('temperature');
grid;
set(gca,'Ylim',[0.9985,1.002] );
legend(mmcellstr(sprintf('Slit #%d|',1:6)),'Orientation','horizontal');
title({'Brewer#185 Ozone Ratio of internal lamp',lp{i}});
%%
figure
%subplot(2,1,1)
plot(t,f,'o');
%[r{i},c1{i},s1{i}]=robust_line;

[r{i},c1{i}]=rline;
cf=c1{i}(1,:)

printmatrix(c1{i}([1,4,6],:));

%%
figure
%subplot(2,1,2)
mmplotyy_temp(d,matdiv(f,f(1,:)),'.',[0.9985,1.002],t,'o',[0,50]);

grid;
datetick('x','keeplimits','keepticks')
legend(mmcellstr(sprintf('Slit #%d|',1:6)),'Orientation','horizontal')
set(gca,'Ylim',[0.9985,1.002] );
title({'Brewer#185 Ozone Ratio of internal lamp',lp{i}});

%%

%fc=f+t*c1(2,:);
%load cf
cf=c1{i}(1,:)
fc=f-t*cf;

rc=fc*W;
r0=f*W;

[mt,st]=grpstats([d,t,rc,r0],t);
[md,sd,n]=grpstats([d,t,rc,r0],fix(d*24)/24);

%%
figure
mmplotyy_temp(md(:,1),md(:,8),'.',[300 360],md(:,2),'x:',[10,50]);
hold on
h=errorbar([md(:,1),md(:,1)],md(:,[8,14]),sd(:,[8,14]),'o');
mmplotyy('Temperature');
datetick('keepticks');
grid;
legend(h,'T corr','T uncorr')
ylabel('R6 : Ozone Ratio');
title({'Brewer#185 Ozone Ratio of internal lamp',lp{i}});
set(gca,'YLim',[300 360]);
%%
figure
 hold on;
%errorbar([mt(:,2),mt(:,2)],mt(:,[8,14]),st(:,[8,14]),'o');
errorbar(mt(:,2),mt(:,[8]),st(:,[8]),'ko');
errorbar(mt(:,2),mt(:,[14]),st(:,[14]),'bo');


plot(mt(:,2),mt(:,8),'ko',mt(:,2),mt(:,14),'bo');
legend('corr','uncorr')
[p,rend{i}]=rline;
grid;
set(gca,'YLim',[300 360]);
xlabel('temp');
ylabel('R6 : Ozone Ratio');
title({'Brewer#185 Ozone Ratio of internal lamp',lp{i}});
set(gca,'YLim',[300 360]);
%%

figure
for i=1:6;
    subplot(3,2,i);
    plot(t,rc(:,i),'bx',t,r0(:,i),'ro');
    title(leg_ratios{i});
    xlabel('temperature');
    ylabel(' log cont/sec /C')
    %axis('tight');
    set(gca,'Xlim',[0,50]);
    %set(gca,'YLim',[300 360]);
end

suptitle('Ratios/Double ratios corrected (blue) / uncorrected(red)');





end
% figure
% mt=grpstats([d,t,f],t);
% figure;
% plot(mt(:,2),mt(:,3:end),'o');
% [r,c]=robust_line;

% %%
% figure
% for i=1:6;subplot(3,2,i);
%     plot(T,RC(:,i),'ro',T,R0(:,i),'bx');end
%     suptitle('Ratios corrected (red) / uncorrected(blue)');
% 
% %%
% figure
% for i=1:6;subplot(3,2,i);plot(T,Fcc(:,i),'ro',T,F0(:,i),'bx');end
% suptitle('Counts/second corrected (red) / uncorrected(blue)');
% %%
% figure
% j=T<30;
% for i=1:6;subplot(3,2,i);
%     plot(T(j),RC(j,i),'ro',T(j),Rc(j,i),'bx');
%     rline
% end
% suptitle('Ratios (Temp >38) corrected (red) / uncorrected(blue)');
% %end


%% Chamber 1
% Applty ZD coeff
