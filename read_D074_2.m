fpath='data_set_2/D074';
l=dir(fullfile(fpath,'*.074'));
%/* Date Type Time_C RC NC Ozone_C Time_D RD ND Ozone_D Time_A RA NA Ozone_A Time_AD MU_AD SZA_AD Ozone_AD Time_CD MU_CD SZA_CD Ozone_CD  */
%14.09.2016 DS 8:00:00 165.6 122.7 281.9 8:01:00 104.6  71.9 290.2 8:02:00 278.4 223.3 282.3 8:01:30 3.749  75.3 280.1 8:00:30 3.798  75.5 275.4 
dobson=[];
for i=1:length(l)
   t=readtable(fullfile(fpath,l(i).name),'filetype','text');
   tim=3:4:19;
   tq=t; tq{:,[1,2,tim]}=num2cell(NaN);
   dob=cell2mat(table2cell(tq));
   for i=tim
     tx=([t{:,1},t{:,i}]);
     t(:,i)=cellstr(join(tx));
     dob(:,i)=datenum(t{:,i},'dd.mm.yyyy HH:MM:SS');
   end
   dob(:,1)=datenum(t{:,15},'dd.mm.yyyy HH:MM:SS'); %AD
   dob(:,2)=[];
   dob(:,end)=[];
   dobson=[dobson;dob];
end
   
% error code
%[i,j]=find(dobson==-0.999)
%dobson(i,:)=[];


leg=' Date Time_C RC NC O3_C T_D RD ND O3_D T_A RA NA O3_A Time_AD MU_AD SZA_AD O3_AD T_CD MU_CD SZA_CD O3_CD';
leg=mmcellstr(strrep(leg,' ','|'));

D074_orig=array2table(dobson,'VariableNames' ,leg,'RowNames',cellstr(datestr(dobson(:,1))));
writetable(D074_orig,'Atmoz_o3_set2.xls','Sheet','D074','WriteRowNames',true);

D074_orig.Time=date_time(D074_orig.Date);
AD=D074_orig(:,{'Time','Time_AD','MU_AD','O3_AD'});
CD=D074_orig(:,{'Time','T_CD','MU_CD','O3_CD'});
AD.O3_STD=zeros(size(AD.O3_AD));
x1=AD(AD.MU_AD<3.0,{'Time','Time_AD','O3_AD','O3_STD','MU_AD'});
CD.O3_STD=ones(size(CD.O3_CD));
x2=CD(CD.MU_CD>2.5,{'Time','T_CD','O3_CD','O3_STD','MU_CD'});

x1.O3=x1.O3_AD;
x2.O3=x2.O3_CD;
x1.M=x1.MU_AD;
x2.M=x2.MU_CD;
x1.Time=x1.Time;
x1.Date=x1.Time_AD;
x2.Time=x2.Time;
x2.Date=x2.T_CD;
x2.Properties.RowNames=cellstr(datestr(x2.Date));
t_dob=union(x1(:,{'Time','Date','O3','O3_STD','M'}),x2(:,{'Time','Date','O3','O3_STD','M'}));

t_set_2{n_inst}=sortrows(t_dob,1);


