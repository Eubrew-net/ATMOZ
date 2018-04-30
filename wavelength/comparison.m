
load(fullfile('o3x','o3_set_45.mat')); % ozone abs (brw, bp, d&m, uip)
o3_set=o3_set_45;
label_o3x={'Brw','B&P','DMB','SGW','SGQ','SG1','BPB'}
label2={'brewer scan','laser quadratic','laser cubic','lamp quadratic','lamps cubic'}
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];
SO2W=[  0  -1.00    0.00    0.00    4.20   -3.20];
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

% ozone comparison Serdyuchenco
SDK='Serdyuchenco';



%%
mic_range=[1020
        2087
        2423
        4174
        6213
        7711];
    sitn=[0,2:6];
%legend(cellstr(num2str([0,2:6]')))    

%%  load laser_scan resutls
load lscan
load brw_scan_cubic
load brw_scan_quad
load lamp_c; 
load lamp_q;
label={'laser brw o3','brw laser q','brw laser c','brw lamp q','brw lamp c'};
label2={'brewer scan','laser quadratic','laser cubic','lamp quadratic','lamps cubic'}
laserscan=lscan;
%%




O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];

ozone_slits=laserscan(1:6,:)
o3abs_ls=-ozone_slits(:,[6:7,end])'*O3W'/log(10)

%% save
for i=1:size(label)
 eval([label{i},' = opo(:,:,',num2str(i),');'])
 eval(['save ',label{i},' ',label{i},' -ascii '])
end



%%
f2=figure
plot(laserscan(:,4),laserscan(:,5),'x')
hold on
plot(brw_scan_cubic(:,3),brw_scan_quad(:,4),'+')
plot(brw_scan_cubic(:,3),brw_scan_cubic(:,4),'p')

plot(lamp_c(:,3),lamp_c(:,4),'s')
plot(lamp_q(:,3),lamp_q(:,4),'o')

legend(label)
grid
xlabel('wavelength [A]')
ylabel('FWHM [A]')

set(f2,'Tag','brewer_scan_cubic_residual')
printfiles_report(f2,'figures','Width',18,'Height',18)
%%
% remove duplicate measurements
[C,IA,IC] = unique(lscan(:,2:3),'rows','sorted');
laserscan=lscan(IA,:)
%% TODO compare repeated measurements
% one measure missing ??
laserscan=insertrows(laserscan,laserscan(1,:)*NaN,17);
% set time, mic, slit, wv, fwhm, sdk
lscan1=laserscan(:,[2:5,7]);
lscan2=laserscan(:,[2:3,8:end]);

opo=cat(3,lscan1,brw_scan_quad(:,[1:4,8]),brw_scan_cubic(:,[1:4,8]),lamp_q(:,[1:4,8]),lamp_c(:,[1:4,8]))
wl=squeeze(opo(:,3,1));
fwl=squeeze(opo(:,4,1));
%%
figure
h=plot(wl,matdiv(squeeze(opo(:,5,:)),squeeze(opo(:,5,1))))
set(h([2,4]),'linewidth',2)
legend(label);
ylabel('ratio')
xlabel('wavelength')
title([SDK,' efective ozone cross section ratio'])

figure
s=[1,3,5];
h=plot(wl,matdiv(squeeze(opo(:,5,s)),squeeze(opo(:,5,1))))
%set(h([3,5]),'linewidth',3)
legend(label(s));
ylabel('ratio')
xlabel('wavelength')
title([SDK,' efective ozone cross section ratio'])
%%
fw=figure
fdwl=matadd(squeeze(opo(:,4,:)),-fwl);fdwl(:,1)=wl; 
fdwl=sortrows(fdwl,1);
h=ploty(fdwl);
set(h(1),'Marker','.')
set(h(2),'Marker','o')
set(h(3),'Marker','x')
set(h(4),'Marker','+')
legend(label2(2:end));
grid
xlabel('wavelength [A]')
ylabel('wavelength [A]')
title('Full Width Half Maximun')
legend(label2(2:end));
hline([-0.1,0.1],'k-')
vline(BWN(3:end)*10)
set(h(1),'LineStyle','-')
set(h(2),'LineStyle',':')
set(h(3),'LineStyle','-.')
set(h(4),'LineStyle','--')
figure(fw)
set(fw,'Tag','fwhm_comparison')
printfiles_report(fw,'figures','Width',12,'Height',12)




%%
cw=figure
dwl=matadd(squeeze(opo(:,3,:)),-wl);dwl(:,1)=wl; dwl=sortrows(dwl,1);
h=ploty(dwl);

set(h(1),'Marker','.')
set(h(2),'Marker','o')
set(h(3),'Marker','x')
set(h(4),'Marker','+')
legend(label2(2:end));
grid
xlabel('wavelength [A] ')
ylabel('wavelength [A]')
title('Central wavelength difference')
legend(label2(2:end),'Location','South') %,'Orientation','Horizontal');
hline([-0.1,0.1],'-k')
vline(BWN(3:end)*10)
set(h(1),'LineStyle','-')
set(h(2),'LineStyle',':')
set(h(3),'LineStyle','-.')
set(h(4),'LineStyle','--')
figure(cw)
set(cw,'Tag','central_comparison')
printfiles_report(cw,'figures','Width',12,'Height',12)


%%
cw=figure
dwl=matadd(squeeze(opo(:,3,:)),-wl);dwl(:,1)=wl; dwl=sortrows(dwl,1);
yyaxis left
h=ploty(dwl);

grid
xlabel('wavelength')
ylabel('wavelength (A)')
title('Central Wavelength')
hline([-0.1,0.1])
vline(BWN(3:end)*10)
set(h(1),'LineStyle','-')
set(h(2),'LineStyle',':')
set(h(3),'LineStyle','-.')
set(h(4),'LineStyle','--')


set(h(1),'Marker','.','Color','b')
set(h(2),'Marker','o','Color','r')
set(h(3),'Marker','x','Color','k')
set(h(4),'Marker','+','Color','g')

yyaxis right
h1=ploty(fdwl);
set(h1(1),'Marker','.')
set(h1(2),'Marker','o')
set(h1(3),'Marker','x')
set(h1(4),'Marker','+')
set(h1(1),'LineStyle','-')
set(h1(2),'LineStyle',':')
set(h1(3),'LineStyle','-.')
set(h1(4),'LineStyle','--')

legend(h,label2(2:end),'Orientation','Horizontal');
hline([-0.1,0,0.1],'b')

set(cw,'Tag','cwfwhm_comparison')
printfiles_report(cw,'figures','Width',12,'Height',12)

%%



%% FWHK

f=figure
h=plot(wl,matadd(squeeze(opo(:,3,:)),-squeeze(opo(:,3,1))))
%set(h([3,5]),'linewidth',3)
legend(label);
ylabel('Angstrom')
xlabel('wavelength')
title([' Central wavelength diff'])
%

figure
s=[1,3,5];
h=plot(wl,matadd(squeeze(opo(:,3,s)),-squeeze(opo(:,3,1))))
%set(h([3,5]),'linewidth',3)
legend(label(s));
ylabel('Angstrom')
xlabel('wavelength')
title([' Central wavelength diff'])

%%
figure
s=[1,3,5];
h=plot(wl,matadd(squeeze(opo(:,4,:)),-squeeze(opo(:,4,1))))
%set(h([3,5]),'linewidth',3)
legend(label);
ylabel('Angstrom')
xlabel('wavelength')
title([' FWHM wavelength diff'])

figure
s=[1,3,5];
h=plot(wl,matadd(squeeze(opo(:,4,s)),-squeeze(opo(:,4,1))))
%set(h([3,5]),'linewidth',3)
legend(label(s));
ylabel('Angstrom')
xlabel('wavelength')
title([' FWHM wavelength diff'])

