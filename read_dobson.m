fpath='data_set_1/D074'
l=dir(fullfile(fpath,'*.074'))
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
   dob(:,1)=datenum(t{:,1},'dd.mm.yyyy');
   dob(:,2)=[];
   dob(:,end)=[];
   dobson=[dobson;dob];
end
   
% error code
%[i,j]=find(dobson==-0.999)
%dobson(i,:)=[];

%save D064orig.dat dobson -ascii -double 
%dob=readtable('OD16258Header.064','FileType','text');
%leg='Date HH MM SS Mu_C M_C RC NC O3_C Mu_D M_D RD ND O3_D Mu_A M_A RA NA O3_A Mu_CD M_CD O3_CD Mu_AD M_AD O3_AD';

leg=' Date Time_C RC NC O3_C T_D RD ND O3_D T_A RA NA O3_A Time_AD MU_AD SZA_AD O3_AD T_CD MU_CD SZA_CD O3_CD'
leg=mmcellstr(strrep(leg,' ','|'));

D074_orig=array2table(dobson,'VariableNames' ,leg);
writetable(D074_orig,'D074_orig.xls');
