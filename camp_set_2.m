clear all
close all

date_time= @(x) datetime(x,'ConvertFrom','datenum')
addpath(genpath(fullfile('D:','CODE','rbcce.aemet.es','iberonesia','matlab')));
plt= {'o','+','*','h','x','s','d','v','>','<','p','+','x','*','x','s'};
colo=num2cell(parula(16),2);
x1=plt(1:15)';
t_set_2={};
%%
l_inst=dir('data_set_2')
inst={l_inst.name};
inst(1:2)=[]
fun_name=strcat('read_',inst,'_2');
figure;hold all

for ii=length(inst):-1:1
    n_inst=ii;
    disp(inst{ii})
    %t_set_2{ii}=[];
    eval(fun_name{ii});
   

    %data from/to 12-30 september
    t_set_2{ii}=t_set_2{ii}(t_set_2{ii}.Date>datenum(2016,9,12) & t_set_2{ii}.Date<datenum(2016,9,30),:);
    plot(t_set_2{ii}{:,2},t_set_2{ii}{:,3});
    
end

%
inst=upper(inst)
legend(inst)
figure
%%
for ii=1:length(inst)
    try
     disp(inst(ii))   
     t_set_2{ii}.Properties.VariableNames{3}='O3_2';
     t_set_2{ii}.Properties.VariableNames{4}='O3_STD';
     t_set_2{ii}.Properties.VariableNames{5}='AIRM';
     t_set_2{ii}(1,:)
    catch
     disp(inst(ii))
    end
end
%%
ts=cellfun(@(x) table2timetable(x),t_set_2,'UniformOutput',false);
s=vartype('numeric');
ts=cellfun(@(x) x(:,s),ts,'UniformOutput',false);
t_sync=[datetime(2016,09,12):minutes(10):datetime(2016,09,30)];
 t_sync=t_sync(hour(t_sync)>5 & hour(t_sync)<21);
 tx=synchronize(ts{:},t_sync,'mean');
 
 fecha=datenum(tx.Time);
 o3=(tx{:,2:4:end});
 ref=o3(:,2:4); 
 ref=nanmean(ref,2);
 m=(tx{:,4:4:end});
 [za,m2,m3]=sza(fecha);
osc=m2.*ref;
ref_set_2=ref;
osc_set_2=osc;
o3_set_2=o3;
 tx.Properties.VariableNames(2:4:end) =...
     strcat(tx.Properties.VariableNames(2:4:end),'_',inst);

 tx.sza=za;
 tx.ref=ref;
 tx.osc=osc;
 writetable(timetable2table(tx),'ATMOZ_SYCN__DATASET_2.xls','Sheet','10_min');

%%
figure
boxplot(100*(o3-ref)./ref,'Labels',inst,'plotstyle','compact')
grid
hline([-1,1]);
set(gca,'Ylim',[-10,7]);
title('ATMOZ Campaing boxplot  10 min simultaneous obs ');
suptitle('ATMOZ  IZO SET 2 ')
ShadePlotForEmpahsis_x([-1,1],'k',0.1)
ylabel('% ratio to reference')

%% OSC table
fecha=datenum(tx.Time);
rat=(o3-ref)./ref;
[za,m2,m3]=sza(fecha);
figure

  [tr,m_set2,set2_stat]=osctable([fecha,100*rat,ref.*m2],[400,600,1000],inst);
  tr=t_table(tr)
  writetable(tr,...
  'ATMOZ_SYCN__DATASET_2.xls','Sheet','osc table (Tsync=10 min) ',...
  'WriteRowNames',true);
grid
hline([-1,1]);
set(gca,'Ylim',[-12,10]);
rotateticklabel(gca,60)
suptitle('ATMOZ Campaing  SET 2 boxplot  10 min simultaneous obs ');
save atmoz_set_2
 %%
figure
ha = tight_subplot(4,4,[.03 .03],[.1 .1],[.1 .1]);
suptitle('ATMOZ  IZO SET 2 (reference in blue)')
for i=1:16
    axes(ha(i));
    plot(fecha,o3(:,i),'.'); hold on;
    plot(fecha,ref,'-');
    title(inst{i});
    set(ha(i),'Xlim',[datenum(2016,9,13),datenum(2016,9,28)]);
    datetick('keeplimits');
    if mod(i,4)~=1
         set(ha(i),'YtickLabel',[]);
    end
    if i<=12 
        set(ha(i),'XtickLabel',[]);
    end
    grid;
end

%%
figure
ha = tight_subplot(4,4,[.1 .1],[.03 .03],[.1 .1]);
for i=1:16
    axes(ha(i));
    histogram(100*(o3(:,i)-ref)./ref,'Normalization','probability');
    title(inst{i});
    %set(ha(i),'Xlim',[-6,3]);
    %datetick('keeplimits');
%     if mod(i,4)~=1
%          set(ha(i),'YtickLabel',[]);
%     end
%     if i<=12 
%         set(ha(i),'XtickLabel',[]);
%     end
end
%% hourly
 th=synchronize(ts{:},'hourly','mean');
 fechah=th{:,1:4:end};
 o3h=th{:,2:4:end};
 refh=o3h(:,2:4); refh=nanmean(refh,2);
 
o3h_2=o3h;

refh_2=refh;
h_set2.fecha=fechah;
h_set2.o3h=o3h;
h_set2.refh=refh;

 figure
ha = tight_subplot(4,4,[.05 .05],[.02 .02],[.1 .1]);
for i=1:16
    axes(ha(i));
    histogram(100*(o3h(:,i)-refh)./refh,'Normalization','probability');
    title(inst{i});
    %set(ha(i),'Xlim',[-6,3]);
    %datetick('keeplimits');
%     if mod(i,4)~=1
%          set(ha(i),'YtickLabel',[]);
%     end
%     if i<=12 
%         set(ha(i),'XtickLabel',[]);
%     end
end
%% 
figure
hr=100*(o3h-refh)./refh;
hr2=hr;
[mh,sh]=grpstats(hr,hour(th.Time));

h_set2.mh=mh;
h_set2.sh=sh;
h_set2.th=th;

save h_set2 h_set2

he=errorbar(repmat(1:24,16,1)',mh,sh);
%set(he(1:2:end),'LineWidth',1)
% pandora
set(he(8:9),'LineStyle',':')
set(he(8:9),'LineWidth',5)
set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-7,7])
 set(gca,'YTick',-7:7)
set(gca,'Ylim',[-7,7]);
set(gca,'Xlim',[7,20]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing hourly ratio ');
xlabel('Hour');
grid 
%%
%%
figure
hr=100*(o3h-refh)./refh;
[md,sd]=grpstats(hr,day(th.Time));
dj=unique(day(th.Time));
he=errorbar(repmat(dj',16,1)',md,sd);
set(he(1:2:end),'LineWidth',2)
set(he(2:4),'LineStyle',':','LineWidth',1)
set(he(5:7),'LineWidth',3)
set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-7,7])
 set(gca,'YTick',-7:7)
set(gca,'Ylim',[-7,7]);
set(gca,'Xlim',[12,30]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing daily ratio ');
xlabel('Hour');
grid 
%%
figure
boxplot(100*(o3h-refh)./refh,'Labels',inst,'plotstyle','compact')
grid
hline([-1,1]);
set(gca,'Ylim',[-10,7]);
title('ATMOZ Campaing boxplot hourly obs ');
%% airmass
 figure
 sza_dep={};
ha = tight_subplot(5,3,[.05 .05],[.05 .05],[.1 .1]);
for i=1:15
    axes(ha(i));
    %sza_dep{i}=mean_smooth_new(za,(o3(:,i)-ref)./ref,0.1,1); 
    plot(za,100*(o3(:,i)-ref)./ref,'x');
    set(ha(i),'Ylim',[-5,5]);
    grid on; hold on;
    title(inst{i});
    set(ha(i),'Xlim',[20,85]);
    sza_dep{i}=mean_smooth_new(za,100*(o3(:,i)-ref)./ref,0.05,0);
    plot(sza_dep{i}(:,1),sza_dep{i}(:,2))
    
    grid on;  
    %datetick('keeplimits');
%     if mod(i,4)~=1
%          set(ha(i),'YtickLabel',[]);
%     end
%     if i<=12 
%         set(ha(i),'XtickLabel',[]);
%     end
end
%%
 x=cell2mat(sza_dep);
 y=reshape(x',7,15,[]);
 figure
za_=squeeze(y(1,:,:))';
 %H=plot(za_,squeeze(y(2,:,:))','Marker','none')
 hold all
 H=plot(za_(1:50:end,:),squeeze(y(2,:,1:50:end))')
 
  set(gca,'YLim',[-7,7])
  set(H,{'Color'},colo(1:15))
  set(H,{'Color'},colo(1:15))
  set(gca,'YTick',-7:7)
  set(H,{'Marker'},plt(1:15)')
 
 hx=legendflex(H,inst,'nrow',3)

 title('ATMOZ Campaing % ratio  (SZA) ');
 xlabel('SZA')
 grid on
 
%%
 %%
 figure
ha = tight_subplot(4,4,[.05 .05],[.02 .02],[.1 .1]);
for i=1:16
    axes(ha(i));
    histogram(100*(o3h(:,i)-refh)./refh,'Normalization','probability');
    title(inst{i});
    %set(ha(i),'Xlim',[-6,3]);
    %datetick('keeplimits');
%     if mod(i,4)~=1
%          set(ha(i),'YtickLabel',[]);
%     end
%     if i<=12 
%         set(ha(i),'XtickLabel',[]);
%     end
end
%% airmass
 figure
 osc_dep={};
ha = tight_subplot(5,3,[.05 .05],[.05 .05],[.1 .1]);
for i=1:15
    axes(ha(i));
    %sza_dep{i}=mean_smooth_new(za,(o3(:,i)-ref)./ref,0.1,1); 
    plot(osc,100*(o3(:,i)-ref)./ref,'x');
    set(ha(i),'Ylim',[-5,7]);
    grid on; hold on;
    title(inst{i});
    set(ha(i),'Xlim',[250,1500]);
    try
    osc_dep{i}=mean_smooth_new(osc,100*(o3(:,i)-ref)./ref,0.1,0);
    plot(osc_dep{i}(:,1),osc_dep{i}(:,2))
    catch
      disp(inst(i))
    end
    grid on;  
    %datetick('keeplimits');
%     if mod(i,4)~=1
%          set(ha(i),'YtickLabel',[]);
%     end
%     if i<=12 
%         set(ha(i),'XtickLabel',[]);
%     end
end
%%
 x=cell2mat(osc_dep);y=reshape(x',7,15,[]);
 figure
 osc=squeeze(y(1,:,:))';
 r=squeeze(y(2,:,:))';
 H=plot(osc,r)
 set(gca,'YLim',[-5,7])
 set(gca,'YTick',-5:7)
 hx=legendflex(H,inst,'nrow',3);
 set(H,{'Color'},colo(1:15))
 title('ATMOZ Campaing ratio against reference ');
 xlabel('Ozone Slant Column')
 grid on
 set(H,{'Marker'},plt(1:15)')
 set(gca,'XLim',[250,1800])
 %% Daily

 td=synchronize(ts{:},'daily','mean');
 fechad=td{:,1:4:end};
 o3d=td{:,2:4:end};
 refd=o3d(:,2:4); refd=nanmean(refd,2);

figure
he=plot(datenum(td.Time),100*(o3d-refd)./refd)
%he=errorbar(repmat(1:24,16,1)',mh,sh);
set(he(1:2:end),'LineWidth',2)
set(he(1:4),'LineStyle',':','LineWidth',1)
set(he(5:7),'LineWidth',3)
set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-5,7])
 set(gca,'YTick',-5:7)
%set(gca,'Xlim',[7,20]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing dayly ratio ');
datetick('keeplimits');
xlabel('Date');
grid 

%% OSC table
fecha=datenum(tx.Time);
rat=(o3-ref)./ref;
[za,m2,m3]=sza(fecha);
  tr=osctable([fecha,100*rat,ref.*m2],[400,600,1000],inst);
  t_table(tr)
  
%% Daily table  
fechah=datenum(th.Time);
rath=(o3h-refh)./refh;
[zah,m2h,m3h]=sza(fechah);
%th_=osctable([fechah,100*rath,refh.*m2h],[400,600,1000],inst);
t_dayly=array2table(100*(o3d-refd)./refd,'VariableNames',inst);
t_dayly.Time=td.Time;
t_dayly=timetable2table(table2timetable(t_dayly))
figure
boxplot(100*(o3d-refd)./refd,'Labels',inst,'plotstyle','compact')
grid
