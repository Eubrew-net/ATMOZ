fpath='data_set_1/P101'
l=dir(fullfile(fpath,'P101*.txt'))
t=readtable(fullfile(fpath,l.name));
fecha=datetime(t.Var1,'InputFormat','yyyyMMdd''T''HHmmss''Z''' );
t.Properties.VariableNames={'Date','sza','airm','FW5_flag','FW5_o3','FW6_flag','FW6_o3'};
t.Date_str=t.Date;
t.Date=fecha;
t_dep=t(t.FW5_flag==0,:)
writetable(t_dep,'Atmoz_o3_set1.xls','Sheet','P101');
