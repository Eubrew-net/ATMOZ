[zd1,Fz1,Fz1_0]=readb_zd('bdata185/B02216.185','zd');
[zd2,Fz2,Fz3_0]=readb_zd('bdata185/B02316.185','zd');
[zd3,Fz3,Fz3_0]=readb_zd('bdata185/B02416.185','zd');
[zd4,Fz4,Fz4_0]=readb_zd('bdata185/B02516.185','zd');


%


J=[zd1,zd2,zd3,zd4];
Fz=[Fz1.m;Fz2.m;Fz3.m;Fz4.m];
Fzs=[Fz1.s;Fz2.s;Fz3.s;Fz4.s];

Fz(1:3,:)=[]; %testing
Fz(Fz(:,8)<1E6,3:end)=NaN; %filter #1 to include


%remove bad points
J(:,1:15)=[]; %test
J(3:end,J(end,:)<1E6)=NaN;
figure;
mmplotyy_temp(J(2,:),J(23,:),'o',J(14,:),'x')
datetick('x','dd hh:mm','keeplimits','keepticks')
ylabel('Counts /second');
mmplotyy('temperature');
title(' Ozone obs  (ZD), Direct UV port   ');

%%
figure;
r=matdiv(J([17,19:23],:),J([17,19:23],1));
gscatter(J(13,:),r',fix((J(2,:)-J(2,1))*6),'','o');
grid;
xlabel('temperture')
title(' Ozone obs  (ZD), Direct UV port   ');

%% Monitor correction
% a few temperatures mismatch

load b_monitor.csv

%
figure
plot(b_monitor(:,1),b_monitor(:,6),'.',Fz(:,1),Fz(:,2:4));
legend('Brewer temp monitor','Brw1 meas','Brw2 meas','Brw3 meas')
datetick('x','dd hh:mm','keeplimits','keepticks')
title('Temperatures brewer bfile / brewer monitor file')

Jm=interp1(b_monitor(:,1),b_monitor(:,22),J(2,:),'nearest','extrap');
mz=interp1(b_monitor(:,1),b_monitor(:,22),Fz(:,1),'nearest','extrap');

figure
mmplotyy_temp(J(2,:),Jm,'o',J(14,:),'-x');
datetick('x','dd hh:mm','keeplimits','keepticks')
ylabel('Monitor correction');
mmplotyy('Brewer internal temperature');
title('Monitor correction factor');




%%  monitor correction Counts = Counts / monitor 
Fc=matdiv(Fz(:,[5,7:11]),mz);
Jc=matdiv(J([17,19:23],:),Jm);
r=matdiv(Jc,Jc(:,11));
R=matdiv(Fc,Fc(1,:));
R0=matdiv(Fz(:,[5,7:11]),Fz(1,[5,7:11]))

figure;
plot(Fz(:,2),R); hold on;
gscatter(Fz(:,2),R,fix((Fz(:,1)-Fz(1,1))*6),'','o');
plot(Fz(:,2),R0,'k'); 
ylabel('log counts (normalized to 30? C )')
title(' Ozone (Direct) with/without (black)  monitor correction')
xlabel('Temperature');
grid
grid
%%

%% Ratios calculation
% Weight definition for the seven slits
% slit 0 used for hg calibration slit 1-> dark
O3W=[  0.00      0.00   -1.000    0.500    2.200   -1.700];
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

rleg={'ms4','ms5','ms6','ms7','SO_2','O_3'};

% scale ratios
Fc=round(log10(Fc)*10^4);
Ratios=Fc*W;


%% Temperature correction factor
% Slope of the regresion 
figure;
plot(Fz(:,2),Fc,'x'); hold on;
%gscatter(Fz(:,2),R,fix((Fz(:,1)-Fz(1,1))*6),'','o');
[a,b]=rline;
%% Coeficientes
cf=b(1,:)
cf1=b(4,:);
cf2=b(6,:);
%[a,b,stats]=robust_line;
%cf=b(2,:);
%se=cat(2,stats(:).se);
%load cf
%cf1=cf+se(2,:);
%cf2=cf-se(2,:);

% correction factor
Fcc=Fc;
Fcc=Fc-(Fz(:,2)*cf);
Fc1=Fc-(Fz(:,2)*cf1);
Fc2=Fc-(Fz(:,2)*cf2);

RC=Fcc*W;
RC1=Fc1*W;
RC2=Fc2*W;

% b =
%   Columns 1 through 5
%       -5.0309      -4.8377      -4.5331      -4.4746      -4.7335
%         62252        62366        62510        62971        62899
%       0.70684      0.68803      0.64851      0.64428       0.6631
%   Column 6
%       -4.5481
%         62873
%       0.64089
% %
%%
figure;
for i=1:6;
    subplot(2,3,i);
    plot(Fz(:,2),Fc(:,i),'ro',Fz(:,2),Fcc(:,i),'bx')
    legend('uncorr','temp corr','Location','SouthWest');
    xlabel('temperature');
     ylabel('log counts/sec ')
    axis('tight')
end
suptitle('Counts /sec  for slit (blue temperature corrected)')   
%%
figure
h=plot(Fz(:,2),matdiv(Fc,Fc(1,:)),'-',Fz(:,2),matdiv(Fcc,Fcc(1,:)),':');

    xlabel('temperature');
     ylabel('Normalized log counts/sec ')
    axis('tight')
    legend(mmcellstr(sprintf('Slit #%d|',1:6)),'Orientation','horizontal');
    title({'Brewer#185 Ozone Direct normalized intensity',...
        '(dashed lines after temperature correction)'});
 grid;
 
%% R6 ratio
figure;
%plot(Fz(:,2),Ratios(:,end),'rx');
%robust_line
hold on
%gscatter(Fz(:,2),Ratios(:,end),fix((Fz(:,1)-Fz(1,1))*6),'','o');

%
plot(Fz(:,2),Ratios(:,end),'rx',Fz(:,2),RC(:,end),'bo');
rline
%plot(Fz(:,2),RC1(:,end),':',Fz(:,2),RC2(:,end),':');


legend('O_3 ratio uncorr','O_3 ratio T_{corr}');
%hold on
%gscatter(Fz(:,2),RC(:,end),fix((Fz(:,1)-Fz(1,1))*6),'','o');

grid
title('R6 double ratio uncorrected /corrected');
box on;
box on;
xlabel('Temperature');
ylabel(' O_3 ratio');
%%
figure
for i=1:6;
    subplot(3,2,i);
   
    plot(Fz(:,2),Ratios(:,i),'r-o',Fz(:,2),RC(:,i),'bx:');
    xlabel('temperature');
    ylabel('log counts/sec ');
    title(rleg{i});
end
suptitle({'Direct Ozone Ratios and dobule Ratios',...
    'Temperature corrected (blue x ) / uncorrected(red o)'});

    
