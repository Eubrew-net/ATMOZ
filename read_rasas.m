fpath='data_set_1/rasas'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'092016_zen.dat'));
fecha=t.DiaFrac+datenum(2016,1,0);

t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set1.xls','Sheet','rasas');
%%
t_dep=t;
t_dep.O3=t_dep.O3vert;
t_dep.O3_STD=NaN*t_dep.O3;
t_dep.AIRM=NaN*t_dep.O3;
t_dep.Date=fecha;
t_dep.Time=t_dep.date;

t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

