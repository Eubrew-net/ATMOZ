fpath='data_set_2/ftir'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'O3_izana_ftir_ndaccv4.dat'));

 Adob= 2.687E16 ;% mol/cm^2
 L=2.687E25; %m?3
 o3t=1E18*t.Total_Column_O3total_E18mole_cm_2_/Adob; 
 t.O3_DU=o3t;
 t.Var10=[];
 t.date=datetime(datestr(t.date_matlab));
 writetable(t,'Atmoz_o3_set2.xls','Sheet','FTIR','WriteRowNames',true);
%% 
t_dep=t;
t_dep.O3=t_dep.O3_DU;
t_dep.O3_STD=NaN*t_dep.O3;
t_dep.AIRM=NaN*t_dep.O3;
t_dep.Date=t_dep.date_matlab;
t_dep.Time=datetime(datestr(t_dep.Date));

t_set_2{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

 
