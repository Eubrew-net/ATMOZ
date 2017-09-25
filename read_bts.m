fpath='data_set_1/bts'
l=dir(fullfile(fpath,'*.dat'))
dobson=[];
for i=1:length(l)
   dob=textread(fullfile(fpath,l(i).name),'','headerlines',2,'emptyvalue',NaN);
   fecha=datenum(dob(:,1:6));
   dobson=[dobson;[fecha,dob(:,7:end)]];
end
bts=dobson;   
leg='Date O3 O3_std';
leg=mmcellstr(strrep(leg,' ','|'));

BTS_orig=array2table(dobson,'VariableNames' ,leg,'RowNames',cellstr(datestr(dobson(:,1))));
BTS_orig.date=datetime(cellstr(datestr(dobson(:,1))));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
writetable(BTS_orig,'Atmoz_o3_set1.xls','Sheet','BTS','WriteRowNames',true);

%%
t_dep=BTS_orig;

t_dep.O3=t_dep.O3;
t_dep.O3_STD=t_dep.O3_std;
t_dep.AIRM=NaN*t_dep.O3;
t_dep.Time=t_dep.date;
%t_dep.Date=Date;
t_dep=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
