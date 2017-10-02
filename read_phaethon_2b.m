fpath='data_set_2/phaethon'
%l=dir(fullfile(fpath,'*.dat'))

t=readtable(fullfile(fpath,'Phaethon_Izana_Bremen(223K).txt'));
t(:,end)=[];
t.Properties.VariableNames=varname(mmcellstr('Doy, Decimal_doy, Decimal_year, time(UT), sza, AMF, Differential_Slant_column_density, Total_ozone_column(TOC), error_TOC',','));
fecha=t.time_UT_/24+datenum(2016,1,0)+t.Doy;
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set2.xls','Sheet','Phaethon_223','WriteRowNames',true);
t_s{1}=t;

t=readtable(fullfile(fpath,'Phaethon_Izana_Bremen(228K).txt'));
t(:,end)=[];
t.Properties.VariableNames=varname(mmcellstr('Doy, Decimal_doy, Decimal_year, time(UT), sza, AMF, Differential_Slant_column_density, Total_ozone_column(TOC), error_TOC',','));
fecha=t.time_UT_/24+datenum(2016,1,0)+t.Doy;
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set2.xls','Sheet','Phaethon_228','WriteRowNames',true);
t_s{2}=t;


t=readtable(fullfile(fpath,'Phaethon_Izana_Bremen(233K).txt'));
t(:,end)=[];
t.Properties.VariableNames=varname(mmcellstr('Doy, Decimal_doy, Decimal_year, time(UT), sza, AMF, Differential_Slant_column_density, Total_ozone_column(TOC), error_TOC',','));
fecha=t.time_UT_/24+datenum(2016,1,0)+t.Doy;
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set2.xls','Sheet','Phaethon_233','WriteRowNames',true);
t_s{3}=t;
%
figure
for i=1:3
t_dep=t_s{i};
t_dep.O3=t_dep.Total_ozone_column_TOC_;
t_dep.O3_STD=t_dep.error_TOC;
t_dep.AIRM=t_dep.AMF;
t_dep.Date=datenum(t_dep.date);
t_dep.Time=t_dep.date;

   t_s{i}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
   plot(t_s{i}{:,2},t_s{i}{:,3});
   hold on;
end
%
inst={'P223','P228','P233','BRW1'}
legend(inst);
x=load('atmoz_set_2.mat','t_set_2');
t_s{4}=x.t_set_2{4};
%%
ts=cellfun(@(x) table2timetable(x),t_s,'UniformOutput',false);
s=vartype('numeric');
ts=cellfun(@(x) x(:,s),ts,'UniformOutput',false);
t_sync=[datetime(2016,09,12):minutes(10):datetime(2016,09,30)];
 t_sync=t_sync(hour(t_sync)>5 & hour(t_sync)<21);
 tx=synchronize(ts{:},t_sync,'mean');
 
 
 fecha=datenum(tx.Time);
 o3=(tx{:,2:4:end}); 
 %ref=nanmean(o3,2);
 ref=o3(:,end);
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

