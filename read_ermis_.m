fpath='data_set_1/ermis'
%l=dir(fullfile(fpath,'*.dat'))
t=readtable(fullfile(fpath,'ERMIS_o3_Izana2016.dat'));


t.date=datetime(datestr(datenum(t{:,1:6})));
writetable(t,'Atmoz_o3_set1.xls','Sheet','ermis');
