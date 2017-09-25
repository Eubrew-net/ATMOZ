fpath='data_set_2/B017';
l=dir(fullfile(fpath,'B*.xls'));

file=fullfile(fpath,l(1).name)

B_orig=readtable(file);

writetable(B_orig(:,1:10),'Atmoz_o3_set2.xls','Sheet','B017','WriteRowNames',true);
B_orig.Time=date_time(B_orig.Date);
B_orig.Properties.VariableNames([1,2,11,8,20])
t_set_2{n_inst}=B_orig(:,[1,2,11,8,20]);