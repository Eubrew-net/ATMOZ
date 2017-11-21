fpath='data_set_2/ermis'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'ERMIS_o3_Izana2016_BremenXS.dat'));
fecha=datenum(t{:,1:6});
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set2.xls','Sheet','ERMIS','WriteRowNames',true)

%%
t_dep=t;

t_dep.O3=t_dep.OZONE_DU_;
t_dep.O3_STD=t_dep.Un_OZONE_DU_;
t_dep.AIRM=t_dep.AMF;
t_dep.Time=datetime(datestr(t_dep.date));
t_dep.Date=fecha;
t_dep=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

t_set_2{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
