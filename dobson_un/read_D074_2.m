fpath='data_set_2/D074';
l=dir(fullfile(fpath,'*.074'));
%/* Date Type Time_C RC NC Ozone_C Time_D RD ND Ozone_D Time_A RA NA Ozone_A Time_AD MU_AD SZA_AD Ozone_AD Time_CD MU_CD SZA_CD Ozone_CD  */
%14.09.2016 DS 8:00:00 165.6 122.7 281.9 8:01:00 104.6  71.9 290.2 8:02:00 278.4 223.3 282.3 8:01:30 3.749  75.3 280.1 8:00:30 3.798  75.5 275.4 
dobson=[];
for i=1:length(l)
   t=readtable(fullfile(fpath,l(i).name),'filetype','text','Duration','text');
   tim=3:4:19;
   tq=t(:,1:22); tq{:,[1,2,tim]}=num2cell(NaN);
   %ts=vartype('numeric');
   %dob=t(:,ts);
   dob=cell2mat(table2cell(tq));
   for i=tim
     tx=([t{:,1},t{:,i}]);
     t(:,i)=cellstr(join(tx));
     dob(:,i)=datenum(t{:,i},'dd.mm.yyyy HH:MM:SS');
   end
   dob(:,1)=datenum(t{:,15},'dd.mm.yyyy HH:MM:SS'); %AD
   dob(:,2)=[];
   if size(dob,2)==23 dob(:,end)=[]; end
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

D074_orig.O3_AD=D074_orig.O3_AD/1000;
D074_orig.O3_CD=D074_orig.O3_CD/1000;

Sec_Z=sec(deg2rad(D074_orig.SZA_AD));
D074_orig.M_AD=Sec_Z-0.0018167*(Sec_Z-1)-0.002875*(Sec_Z-1).^2-0.0008083*(Sec_Z-1).^3;

Sec_Z=sec(deg2rad(D074_orig.SZA_CD));
D074_orig.M_CD=Sec_Z-0.0018167*(Sec_Z-1)-0.002875*(Sec_Z-1).^2-0.0008083*(Sec_Z-1).^3;


% error analysis 
er_ad{n_inst} = D074_orig(:,{'Time','Date','M_AD','MU_AD','O3_AD'});
er_cd{n_inst} = D074_orig(:,{'Time','Date','M_CD','MU_CD','O3_CD'});
% t = dobson(:,2);
% m = dobson(:,23);
% mu = dobson(:,22);
% Omega = dobson(:,24);




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


