close all
clear all

x157 = load('izo16_157.mat');
x183 = load('izo16_183.mat');
x185 = load('izo16_185.mat');

dateFormat=16;

u157_rel = 100*x157.uo3_save./x157.o3_save;
u183_rel = 100*x183.uo3_save./x183.o3_save;
u185_rel = 100*x185.uo3_save./x185.o3_save;

u157_ms9 = 100*x157.ums9_save./x157.ms9_save;
u183_ms9 = 100*x183.ums9_save./x183.ms9_save;
u185_ms9 = 100*x185.ums9_save./x185.ms9_save;

figure(1)
plot(x157.t_save,u157_rel,'+',x183.t_save,u183_rel,'x',x185.t_save,u185_rel,'o')
datetick('x',dateFormat)
set(gca,'Ylim',[0,10])
xlabel('Time')
ylabel('Standard Uncertainty (k=1), %')
legend('157','183','185');
title('Triad Diurnal uncertainty')
grid on

figure(2)
plot(x157.sza_save,u157_rel,'-+',x183.sza_save,u183_rel,'-x',x185.sza_save,u185_rel,'-o')
set(gca,'Ylim',[0,10])
xlabel('Solar Zenith Angle, deg')
ylabel('Standard Uncertainty (k=1), %')
legend('157','183','185');
title('Triad Diurnal uncertainty')
grid on

figure(3)
plot(x157.t_save,x157.uo3_save,'-+',x183.t_save,x183.uo3_save,'-x',x185.t_save,x185.uo3_save,'-o')
datetick('x',dateFormat)
set(gca,'Ylim',[0,30])
xlabel('Time')
ylabel('Standard Uncertainty (k=1), DU')
legend('157','183','185');
title('Triad Diurnal uncertainty')
grid on

figure(5)
plot(x157.sza_save,u157_ms9,'+',x183.sza_save,u183_ms9,'x',x185.sza_save,u185_ms9,'o')
set(gca,'Ylim',[0,0.5])
set(gca,'Xlim',[20,80])
xlabel('Solar Zenith Angle, deg')
ylabel('Standard Uncertainty (k=1), %')
legend('157','183','185');
title('Triad Diurnal MS9 uncertainty')
grid on