fpath='data_set_1/B157';
l=dir(fullfile(fpath,'B*.xls'));

file=fullfile(fpath,l(1).name)

B157_orig=readtable(file);

writetable(B157_orig(:,1:10),'Atmoz_o3_set1.xls','Sheet','B157','WriteRowNames',true);
B157_orig.Time=date_time(B157_orig.Date);
t_set_1{n_inst}=B157_orig(:,[1,2,7,8,4]);