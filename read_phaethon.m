fpath='data_set_1/phaethon'
%l=dir(fullfile(fpath,'*.dat'))

t=readtable(fullfile(fpath,'Phaethon_Izana_Bass&Paur(228K).txt'));
t(:,end)=[];
t.Properties.VariableNames=varname(mmcellstr('Doy, Decimal_doy, Decimal_year, time(UT), sza, AMF, Differential_Slant_column_density, Total_ozone_column(TOC), error_TOC',','));
fecha=t.time_UT_/24+datenum(2016,1,0)+t.Doy;
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set1.xls','Sheet','Phaeton','WriteRowNames',true);

%%
t_dep=t;
t_dep.O3=t_dep.Total_ozone_column_TOC_;
t_dep.O3_STD=t_dep.error_TOC;
t_dep.AIRM=t_dep.AMF;
t_dep.Date=fecha;
t_dep.Time=t_dep.date;

t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
