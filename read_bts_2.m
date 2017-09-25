fpath='data_set_2/bts'
l=dir(fullfile(fpath,'*.dat'));
dobson=[];
for i=1:length(l)
   dob=load(fullfile(fpath,l(i).name));
   fecha=datenum(dob(:,1:6));
   dobson=[dobson;[fecha,dob(:,7:end)]];
end
%bts=dobson;   
leg='Date O3 O3_std';
leg=mmcellstr(strrep(leg,' ','|'));

Q_orig=array2table(dobson,'VariableNames' ,leg);
Q_orig.date=datetime(cellstr(datestr(dobson(:,1))));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
writetable(Q_orig,'Atmoz_o3_set2.xls','Sheet','BTS','WriteRowNames',true);

%%
t_dep=Q_orig;
%t_dep.O3=t_dep.Total_ozone_column_TOC_;
t_dep.O3_STD=t_dep.O3_std;
t_dep.AIRM=NaN*t_dep.O3_std;
%t_dep.Date=fecha;
t_dep.Time=t_dep.date;

t_set_2{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

