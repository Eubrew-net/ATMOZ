fpath='data_set_1/avodor'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'Avodor_all_days.txt'));
fecha=datenum(t{:,1:6});
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set1.xls','Sheet','avodor','WriteRowNames',true);

%%
t_dep=t;

t_dep.O3=t_dep.O3;
t_dep.O3_STD=t_dep.STDOFO3;
t_dep.AIRM=t_dep.MO3;
t_dep.Time=datetime(datestr(t_dep.date));
t_dep.Date=t_dep.MATLABTIME;
t_dep=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});