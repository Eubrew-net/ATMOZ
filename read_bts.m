fpath='data_set_1/bts'
l=dir(fullfile(fpath,'TOC_BTS*.txt'))
dobson=[];
for i=1:length(l)
   %dob=textread(fullfile(fpath,l(i).name),'%2d.%2d.%4d %2d:%2d:%4d %d','headerlines',3,'emptyvalue',NaN,'delimiter',' ');
   [d,m,a,hh,mm,ss,o3]=textread(fullfile(fpath,l(i).name),'%2d.%2d.%4d %2d:%2d:%4d %d','headerlines',3,'emptyvalue',NaN,'delimiter',' ');
   fecha=datenum(a,m,d,hh,mm,ss);
   dobson=[dobson;[fecha,o3]];
end
bts=dobson;   
leg='Date O3 ';
leg=mmcellstr(strrep(leg,' ','|'));

BTS_orig=array2table(dobson,'VariableNames' ,leg,'RowNames',cellstr(datestr(dobson(:,1))));
BTS_orig.date=datetime(cellstr(datestr(dobson(:,1))));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
writetable(BTS_orig,'Atmoz_o3_set1.xls','Sheet','BTS','WriteRowNames',true);

%%
t_dep=BTS_orig;
t_dep.O3=t_dep.O3;
t_dep.O3_STD=NaN*t_dep.O3;
t_dep.AIRM=NaN*t_dep.O3;
t_dep.Time=t_dep.date;
%t_dep.Date=Date;
t_dep=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

try
    t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
catch
    disp('no set 1');
end
