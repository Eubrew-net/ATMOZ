fpath='data_set_1/B183';
l=dir(fullfile(fpath,'B*.xls'));

file=fullfile(fpath,l(1).name)

B183_orig=readtable(file);

writetable(B183_orig(:,1:10),'Atmoz_o3_set1.xls','Sheet','B183','WriteRowNames',true);
B183_orig.Time=date_time(B183_orig.Date);
t_set_1{n_inst}=B183_orig(:,[1,2,7,8,4]);