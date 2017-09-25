fpath='data_set_1/avodor'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'Avodor_all_days.txt'));
fecha=datenum(t{:,1:6});
t.date=datetime(datestr(fecha));
writetable(t,'Atmoz_o3_set1.xls','Sheet','Avodor','WriteRowNames',true);
