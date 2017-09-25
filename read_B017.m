fpath='data_set_1/B017';
l=dir(fullfile(fpath,'*.txt'));

file=fullfile(fpath,l(1).name)

[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10]=textread(file,...
    '%02d:%02d:%02d %d %d %04d %f %f %2c %f  ');
    % no temp no o3 std
    %1 2 3   4   5  6   % 7   8    9      10 
    %hora   dia mes year sza airm  type  o3 

%j=find(a1>80);a1(j)=a1(j)+1900;
%j=find(a1<80);a1(j)=a1(j)+2000;


date=datenum(a6,a4,a5,a1,a2,a3);
jds=strmatch('ds',a9);
ty(jds)=1;
jzs=strmatch('zs',a9);
ty(jzs)=2;
ozo_ds=[date(jds),a7(jds),a8(jds),a9(jds),a10(jds)];
%ozo_zs=[date(jzs),a12(jzs),a13(jzs),a7(jzs),a8(jzs),a9(jzs),a11(jzs)];
ozo=[date,a6,a4,a5,a1,a2,a3,a7,a8,a10,NaN*a10,ty'];

leg='Date YY Month AA HH mm ss SZA AIRM O3_ O3_STD DS';
leg=mmcellstr(strrep(leg,' ','|'));
B017_orig=array2table(ozo,'VariableNames' ,varname(leg),'RowNames',cellstr(datestr(date)));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
writetable(B017_orig,'Atmoz_o3_set1.xls','Sheet','B017');
B017_orig.Time=date_time(date);
t_set_1{n_inst}=B017_orig(:,[end,1,10,11,9]);
