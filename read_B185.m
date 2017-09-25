fpath='data_set_1/B185';
l=dir(fullfile(fpath,'B*.xls'));

file=fullfile(fpath,l(1).name)

B185_orig=readtable(file);

writetable(B185_orig(:,1:10),'Atmoz_o3_set1.xls','Sheet','B185','WriteRowNames',true);
B185_orig.Time=date_time(B185_orig.Date);
t_set_1{n_inst}=B185_orig(:,[1,2,7,8,4]);