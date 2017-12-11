% Linearity test day 20
%
%addpath(genpath(fullfile('~','CODE','iberonesia','matlab')));
addpath(genpath('~/CODE/rbcce.aemet.es/iberonesia/matlab')) 
path(fullfile(pwd,'matlab'),path)
 set(0,'defaultfigurecolor',[1 1 1])
[jp,F,F_0]=readb_jp('../bdata185/B02016.185');

% 100 ms
%% only for 185
load ../laser_ptb/Monitor_cyclTest.txt

monitor=Monitor_cyclTest;
monitor(1)=NaN; %saturated
F(1:2)=[]; %rubish
F_0(1:2)=[];
brw_cy=reshape(F(1:60),10,[]);  % ten power measurements
figure
plotyy(1:10,monitor/max(monitor),1:10,brw_cy/max(brw_cy))

%% second test
cm=60;
st=[];
mon2=load('../laser_ptb/Mon_1cycle_310_1.txt');
brw2=F(61:60+length(mon2))';
brw20=F_0(61:60+length(mon2))';
r2=brw2./mon2;
r20=brw20./mon2;
cm=cm+length(mon2);
figure
norm_c=r20(30);
semilogy(brw2,r2/norm_c,'o');
st=[st;[brw2,r2/norm_c,brw20,mon2]];
% third test S shape repetition
mon3=load('../laser_ptb/Mon_1cycle_310_2.txt');
brw3=F(cm+1:cm+length(mon3))';
brw30=F_0(cm+1:cm+length(mon3))';
r3=brw3./mon3;
r30=brw30./mon3;
cm=cm+length(mon3);
hold on;
semilogy(brw3,r3/norm_c,'r');
st=[st;[brw3,r3/norm_c,brw30,mon3]];
%% 4  320 global
mon4=load('../laser_ptb/Mon_1cycle_320_(global).txt');
brw4=F(cm+1:cm+length(mon4))';
brw40=F_0(cm+1:cm+length(mon4))';
r4=brw4./mon4;
semilogy(brw4,r4/norm_c,'bo');
cm=cm+length(mon4);
st=[st;[brw4,r4/norm_c,brw40,mon4]];

%% 5  320 direct
mon5=load('../laser_ptb/Mon_1cycle_320_(direct).txt');
brw5=F(cm+1:cm+length(mon5))';
brw50=F_0(cm+1:cm+length(mon5))';
r5=brw5./mon5;
cm=cm+length(mon5);
norm_c=r5(50);
semilogy(brw5,r5/norm_c,'d');
st=[st;[brw5,r5/norm_c,brw50,mon5]];
%% 6  310 direct
mon6=load('../laser_ptb/Mon_1cycle_310_(direct).txt');
brw6=F(cm+1:cm+length(mon6))';
brw60=F_0(cm+1:cm+length(mon6))';
r6=brw6./mon6;
semilogy(brw6,r6/norm_c,'p');
cm=cm+length(mon4);
st=[st;[brw6,r6/norm_c,brw60,mon6]];
%% 7  310 direct whithout Integrating sphere,
[jp1,F1,F1_0]=readb_jp('../bdata185/B02116.185');
brw7=F1';
brw70=F1_0';
mon7=load('../laser_ptb/Monitor_310nm_directPort_noIS-21.1.16.txt');
r7=brw7./mon7;
st=[st;[brw7,r7/norm_c/200,brw70,mon7]];
%%
seto=sortrows(st,1);
semilogy(seto(:,1),seto(:,2),'x')
xlabel('Brewer  Counts/second')
ylabel('Brewer/Monitor (normalized)')
%%  First Region Fit
j=(seto(:,1)>1 & seto(:,1)<11000);
semilogy(seto(j,1),seto(j,2),'rx');
x1=seto(j,1);y1=seto(j,2);ly1=log(y1);

ft1 = fittype( 'rat33' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.StartPoint = [0.712694471678914 0.500471624154843 0.471088374541939 0.0596188675796392 0.681971904149063 0.0424311375007417 0.0714454646006424];

% Fit model to data.
[fit1r, gof] = fit( x1, ly1, ft1, opts );
xt=300:11500;
yt=exp(fit1r(xt));
%  Second Region Fit
j=(seto(:,1)>10000 & seto(:,1)<25000);
semilogy(seto(j,1),seto(j,2),'gx')
x2=seto(j,1);
y2=seto(j,2);
ly2=log(y2);
vt2=0.01:0.0001:1;
p=polyfit(y2,x2,9);
v=polyval(p,vt2);


xt2=10000:20000;
yt2=interp1(v,vt2,xt2,'next');
semilogy(xt2,yt2,'.')
hold on
semilogy(xt,yt,'g.')
% smooth

p2=polyfit(x2,y2,2);
y2s=polyval(p2,xt2);
semilogy(xt2,y2s,'o-');

%% step 3

j=( seto(:,1)>20000);
semilogy(seto(j,1),seto(j,2),'mx')
x3=seto(j,1);y3=seto(j,2);ly3=log(y3);


[p,s,mu]=polyfit(x3,ly3,9);

xt3=20E3:100:6.5E5;
yt3=exp(polyval(p,xt3,[],mu));


semilogy(xt3,yt3,'k-');


%% 
tc_i=[[xt,xt2,xt3]',[yt',y2s,yt3]'];
figure
loglog(st(:,1),st(:,2),'o');
hold on
loglog(tc_i(:,1),tc_i(:,2),'.r');
grid
xlabel(' Brewer counts/second');
ylabel(' Brewer/Monitor normalized');
title(' Brewer nolinearity correction');

%% Smooth curve
J=(tc_i(:,1)<100);
tc_ix=tc_i(~J,1);
tc_iy=tc_i(~J,2);
tc_ly=log(tc_i(~J,2));
figure
plot(tc_ix,tc_ly,'.')
w=1:size(tc_ix);
hold on
a=axis;
plot(seto(:,1),log(seto(:,2)),'ro');
axis(a)

save tc_i tc_i
