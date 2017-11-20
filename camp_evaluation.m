clear all
close all
l=dir('data_set_1')
date_time= @(x) datetime(x,'ConvertFrom','datenum')
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
plt= {'o','+','*','h','x','s','d','v','>','<','p','+','x','*','x','s'};
colo=num2cell(parula(16),2);
x1=plt(1:15)';
l=dir('data_set_1');
inst={l.name};
inst(1:2)=[]

%%
%load('atmoz_set_1','ref_set_1','osc_set_1','o3_set_1');
%load('atmoz_set_2','ref_set_2','osc_set_2','o3_set_2');

load('atmoz_set_1','tx','t_sync');tx1=tx;
load('atmoz_set_2','tx','t_sync');tx2=tx;
tx=synchronize(tx1,tx2,t_sync,'mean');
o31=tx{:,2:4:64};
o32=tx{:,69:4:end-3};

fecha=datenum(tx.Time);
ref1=o31(:,2:4); ref1=nanmean(ref1,2);
ref2=o32(:,2:4); ref2=nanmean(ref2,2);

% m=(tx{:,4:4:end});
[za,m2,m3]=sza(fecha);
osc1=m2.*ref1;
osc2=m2.*ref2;
%% Dobson
figure;plot(za,100*(o31(:,5:7)-o32(:,5:7))./o31(:,5:7),':.')
grid
legend(inst(5:7))
title('Dobson data set 1 vs data set 2  % diff ratio to operative')
ylabel(' (o3 set 1 - o3 set 2 )/ o3 set 1 %')
figure
boxplot((o31(:,5:7)-o32(:,5:7)),'labels',inst(5:7))
grid
ylabel('DU')
title('Dobson data set 1 - data set 2 ')
xlabel('sza')
%%
figure;plot(za,100*(o31(:,2:4)-o32(:,2:4))./o31(:,2:4),':.');
xlabel('sza')
ylabel(' (o3 set 1 - o3 set 2 )/ o3 set 1 %')
title('Brewer data set 1 - data set 2 ')
grid
legend(inst(1:4))
title('Dobson data set 1 vs data set 2  % diff ratio to operative')

figure;boxplot((o31(:,1:4)-o32(:,1:4)),'labels',inst(1:4))
ylabel('DU')
title('Brewer data set 1 - data set 2 ')
grid
xlabel('sza')
print -clipboard -dbitmap

%%
x=11:14;figure;plot(za,100*(o31(:,x)-o32(:,x))./o31(:,x),':.');legend(inst(x))
%%

%%
load('atmoz_set_2','m_set2');
load('atmoz_set_1','m_set1');

b2=boxplot(m_set2(:,2:end-2,1),'plotstyle','compact','labels',inst,'Color','r','OutlierSize',1,'DataLim',[-10,10]);
hold on
b1=boxplot(m_set1(:,2:end-2,1),'plotstyle','compact','labels',inst,'Color','b','OutlierSize',1,'DataLim',[-10,10]);

legend([b1(1),b2(1)],{'Data Set 1','Data Set 2'})
grid
hline([-1,1])
title('ATMOZ Ratio to Reference ')
ylabel(' o3 - o3_{ref} / o3_{ref} %')

%%
figure;x=6:8;grpstats([m_set1(:,x,1)],fix(za/10)*10,0.05);

hold on;grpstats([m_set2(:,x,1)],fix(za/10)*10,0.05);legend(inst(x-1))
legend('data set 1','data set 2');set(gca,'xlim',[0.5,7]);
title(inst(x-1))
xlabel('sza')