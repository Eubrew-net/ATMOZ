clear all
close all
l=dir('data_set_1')
date_time= @(x) datetime(x,'ConvertFrom','datenum')
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
plt= {'o','+','*','h','x','s','d','v','>','<','p','+','x','*','x','s'};
colo=num2cell(parula(16),2)
x1=plt(1:15)';
%%
l=dir('data_set_1')
inst={l.name};
inst(1:2)=[]
fun_name=strcat('read_',inst);
figure;hold all
for ii=1:length(inst)
    n_inst=ii;
    disp(inst{ii})
    eval(fun_name{ii});
   

    %  data from/to 12-30 september
    t_set_1{ii}=t_set_1{ii}(t_set_1{ii}.Date>datenum(2016,9,12) & t_set_1{ii}.Date<datenum(2016,9,30),:);
    plot(t_set_1{ii}{:,2},t_set_1{ii}{:,3});
    
end
%
inst=upper(inst)
legend(inst)
figure
%%
for ii=1:length(inst)
    try
     disp(inst(ii))   
     t_set_1{ii}.Properties.VariableNames{3}='O3';
     t_set_1{ii}.Properties.VariableNames{4}='O3_STD';
     t_set_1{ii}.Properties.VariableNames{5}='AIRM';
     t_set_1{ii}(1,:)
    catch
     disp(inst(ii))
    end
end
%%
ts=cellfun(@(x) table2timetable(x),t_set_1,'UniformOutput',false);
s=vartype('numeric');
ts=cellfun(@(x) x(:,s),ts,'UniformOutput',false);
t_sync=[datetime(2016,09,12):minutes(10):datetime(2016,09,30)];
 t_sync=t_sync(hour(t_sync)>5 & hour(t_sync)<21);
 tx=synchronize(ts{:},t_sync,'mean');
 
 fecha=datenum(tx.Time);
 o3=(tx{:,2:4:end});
 ref=o3(:,2:4); ref=nanmean(ref,2);
 m=(tx{:,4:4:end});
 [za,m2,m3]=sza(fecha);
 osc=m2.*ref;
  tx.Properties.VariableNames(2:4:end) =...
     strcat(tx.Properties.VariableNames(2:4:end),'_',inst);

 tx.sza=za;
 tx.ref=ref;
 tx.osc=osc;
 

ref_set_1=ref;
osc_set_1=osc;
o3_set_1=o3;

 writetable(timetable2table(tx),'ATMOZ_SYCN__DATASET_1.xls','Sheet','10 min');

%%
figure
ha = tight_subplot(4,4,[.03 .03],[.1 .1],[.1 .1]);

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
    %grid on;
end
suptitle('ATMOZ  IZO SET 1 (reference in blue)')
%%
figure
ha = tight_subplot(4,4,[.05 .05],[.02 .02],[.1 .1]);
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
suptitle('ATMOZ  IZO SET 1 reference in blue')
%%
%%
figure
boxplot(100*(o3-ref)./ref,'Labels',inst,'plotstyle','compact')
grid
hline([-1,1]);
set(gca,'Ylim',[-10,7]);
title('ATMOZ Campaing boxplot  10 min simultaneous obs ');
suptitle('ATMOZ  IZO SET 1 ')
ShadePlotForEmpahsis_x([-1,1],'k',0.1)
ylabel('% ratio to reference')

%%
%% OSC table
fecha=datenum(tx.Time);
rat=(o3-ref)./ref;
[za,m2,m3]=sza(fecha);
figure

  [tr,m_set1,set1_stat]=osctable([fecha,100*rat,ref.*m2],[400,600,1000],inst);
  tr=t_table(tr)
  writetable(tr,...
  'ATMOZ_SYCN__DATASET_q.xls','Sheet','osc table (Tsync=10 min) ',...
  'WriteRowNames',true)
grid
hline([-1,1]);
set(gca,'Ylim',[-17,10]);
rotateticklabel(gca,60)
suptitle('ATMOZ Campaing  SET 1 boxplot  10 min simultaneous obs ');
save atmoz_set_1
%%
figure
load('atmoz_set_2','m_set2');
b2=boxplot(m_set2(:,2:end-2,1),'plotstyle','compact','labels',inst,'Color','r','OutlierSize',1,'DataLim',[-10,10]);
hold on
b1=boxplot(m_set1(:,2:end-2,1),'plotstyle','compact','labels',inst,'Color','b','OutlierSize',1,'DataLim',[-10,10]);

legend([b1(1),b2(1)],{'Data Set 1','Data Set 2'})
grid
hline([-1,1])
title('ATMOZ Ratio to Reference ')

%% hourly
 th=synchronize(ts{:},'hourly','mean');
 fechah=th{:,1:4:end};
 o3h=th{:,2:4:end};
 refh=o3h(:,2:4); refh=nanmean(refh,2);
   th.Properties.VariableNames(2:4:end) =...
     strcat(th.Properties.VariableNames(2:4:end),'_',inst);


 th.ref=refh;
 writetable(timetable2table(tx),'ATMOZ_SYCN__DATASET_1.xls','Sheet','hour');
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
%%
figure
suptitle('ATMOZ  IZO set 1 comparison')
hr=100*(o3h-refh)./refh;
[mh,sh,nh,hn]=grpstats(hr,hour(th.Time));
writetable(array2table(mh),'ATMOZ_SYCN__DATASET_1.xls','Sheet','hourly means');

he=errorbar(repmat(1:24,16,1)',mh,sh);
%set(he(1:2:end),'LineWidth',1)
set(he(1:4),'LineStyle',':')
set(he(5:7),'LineWidth',3)
set(he(10:11),'LineWidth',3)

set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-10,7])
 set(gca,'YTick',-10:7)
set(gca,'Ylim',[-10,7]);
set(gca,'Xlim',[7,20]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing hourly ratio ');
xlabel('Hour');
grid 
%%
%%
figure
hr=100*(o3h-refh)./refh;
[md,sd,nd,gn]=grpstats(hr,day(th.Time));
dj=unique(day(th.Time));
writetable(array2table([dj,md],...
  'VariableNames',['Day',inst]),... 
  'ATMOZ_SYCN__DATASET_1.xls','Sheet','daily means % ratio');


he=errorbar(repmat(dj',16,1)',md,sd);
%set(he(1:2:end),'LineWidth',2)
set(he(2:4),'LineStyle',':','LineWidth',1)
set(he(5:7),'LineWidth',3)
set(he(16),'LineWidth',3)
set(he(11),'LineWidth',3)
set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-10,7])
 set(gca,'YTick',-10:7)
set(gca,'Ylim',[-10,7]);
set(gca,'Xlim',[12,30]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing daily ratios ');
xlabel('Day');
grid 
%%
figure
boxplot(100*(o3h-refh)./refh,'Labels',inst,'plotstyle','compact')
grid
hline([-1,1]);
set(gca,'Ylim',[-10,7]);
title('ATMOZ Campaing boxplot hourly obs ');
suptitle('ATMOZ  IZO set 1 comparison')
ShadePlotForEmpahsis_x([-1,1],'k',0.1)
ylabel('% ratio to reference')
%% airmass
 figure
 sza_dep={};
ha = tight_subplot(5,3,[.05 .05],[.05 .05],[.1 .1]);
for i=1:15
    axes(ha(i));
    %sza_dep{i}=mean_smooth_new(za,(o3(:,i)-ref)./ref,0.1,1); 
    plot(za,100*(o3(:,i)-ref)./ref,'x');
    set(ha(i),'Ylim',[-10,5]);
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
 suptitle('ATMOZ  IZO set 1 comparison')
za_=squeeze(y(1,:,:))';
 H=plot(za_(1:20:end,:),squeeze(y(2,:,1:20:end))')
  set(gca,'YLim',[-10,7])
  set(H,{'Color'},colo(1:15))
  set(gca,'YTick',-10:7)
  set(H,{'Marker'},plt(1:15)')
  set(H([1,5:7]),'LineWidth',3)
   set(H(11),'LineWidth',3)
   set(H(15),'LineWidth',3)
 hx=legendflex(H,inst,'nrow',3)

 title('ATMOZ Campaing ratio against reference');
 grid on
 xlabel('SZA')
 ylabel('% ratio to reference')
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
 suptitle('ATMOZ  IZO set 1 comparison')
 osc_dep={};
ha = tight_subplot(5,3,[.05 .05],[.05 .05],[.1 .1]);
for i=1:15
    axes(ha(i));
    %sza_dep{i}=mean_smooth_new(za,(o3(:,i)-ref)./ref,0.1,1); 
    plot(osc,100*(o3(:,i)-ref)./ref,'x');
    set(ha(i),'Ylim',[-10,5]);
    grid on; hold on;
    title(inst{i});
    set(ha(i),'Xlim',[250,1500]);
    osc_dep{i}=mean_smooth_new(osc,100*(o3(:,i)-ref)./ref,0.05,0);
    plot(osc_dep{i}(:,1),osc_dep{i}(:,2))
    
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
  suptitle('ATMOZ  IZO SET 1 comparison')
 osc=squeeze(y(1,:,:))';
 r=squeeze(y(2,:,:))';
 H=plot(osc(1:20:end,:),r(1:20:end,:))
 set(gca,'YLim',[-10,7])
 set(gca,'YTick',-10:7)
 hx=legendflex(H,inst,'nrow',3);
 set(H,{'Color'},colo(1:15))
 title('ATMOZ Campaing ratio against reference ');
 xlabel('Ozone Slant Column')
 ylabel('% ratio to reference')
 grid on
 set(H,{'Marker'},plt(1:15)')
 set(gca,'XLim',[250,1800])

 %% Daily

 td=synchronize(ts{:},'daily','mean');
 fechad=td{:,1:4:end};
 o3d=td{:,2:4:end};
 refd=o3d(:,2:4); refd=nanmean(refd,2);
 td.ref=refd;
writetable(timetable2table(td),...
  'ATMOZ_SYCN__DATASET_1.xls','Sheet','daily mean ');


figure
suptitle('ATMOZ  IZO SET 1 comparison')
he=plot(datenum(td.Time),100*(o3d-refd)./refd)
%he=errorbar(repmat(1:24,16,1)',mh,sh);
set(he(1:2:end),'LineWidth',2)
set(he(1:4),'LineStyle',':','LineWidth',1)
set(he(5:7),'LineWidth',3)
set(he,{'Marker'},plt')
set(he,{'Color'},colo)
 set(gca,'YLim',[-10,7])
 set(gca,'YTick',-10:7)
  set(he(11),'LineWidth',3)
   set(he(16),'LineWidth',3)
%set(gca,'Xlim',[7,20]);
legendflex(he,inst,'nrow',4)
title('ATMOZ Campaing daily mean (hourly data) ');
datetick('keeplimits');
xlabel('Date');
grid 

%% OSC table
fecha=datenum(tx.Time);
rat=(o3-ref)./ref;
[za,m2,m3]=sza(fecha);
  tr=osctable([fecha,100*rat,ref.*m2],[400,600,1000],inst);
  tr=t_table(tr)
  writetable(tr,...
  'ATMOZ_SYCN__DATASET_1.xls','Sheet','osc table ',...
  'WriteRowNames',true);

  
%% OSC table  
fechah=datenum(th.Time);
rath=(o3h-refh)./refh;
[zah,m2h,m3h]=sza(fechah);
th_=osctable([fechah,100*rath,refh.*m2h],[400,600,1000],inst);
osc_hour_table=t_table(th_)
writetable(tr,...
  'ATMOZ_SYCN__DATASET_1.xls','Sheet','osc table (hourly data) ',...
  'WriteRowNames',true);

 