fpath='data_set_2/D083'
l=dir(fullfile(fpath,'*.txt'));
dobson=[];
for i=1:length(l)
  dob=textread(fullfile(fpath,l(i).name),'','headerlines',8,'emptyvalue',NaN,'whitespace',':\t\b');
  fechas=sscanf(l(i).name,'%04d-%02d-%02d_%*s');
   %fecha=dia_mes(fechas(1),fechas(2));
   fecha=datenum(fechas');
   fecha=fecha(1)+dob(:,1)/24+dob(:,2)/24/60+dob(:,3)/24/60/60;
   
   dobson=[dobson;[fecha,dob]];
end
   
% error code
[i,j]=find(dobson==-0.999);
dobson(i,:)=[];

%save D083orig.dat dobson -ascii -double 
%dob=readtable('OD16258Header.064','FileType','text');
    %Date HH MM SS Mu_C AirM_C R_C N_C O3_C	Mu_D 	AirM_D	R_D	 N_D	 O3_D	 Mu_A	 AirM_A	R_A 	 N_A	  O3_A	 Mu_CD	AirM_CD	O3_CD Mu_AD	 M_AD	 O3_AD
leg='Date HH MM SS Mu_C M_C RC NC O3_C Mu_D M_D RD ND O3_D Mu_A M_A RA NA O3_A Mu_CD M_CD O3_CD Mu_AD M_AD O3_AD';
leg=mmcellstr(strrep(leg,' ','|'));

D083_orig=array2table(dobson,'VariableNames' ,leg,'RowNames',cellstr(datestr(dobson(:,1))));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
writetable(D083_orig,'Atmoz_o3_set2.xls','Sheet','D083');

D083_orig.Time=date_time(D083_orig.Date);
AD=D083_orig(:,{'Time','Date','M_AD','O3_AD'});
CD=D083_orig(:,{'Time','Date','M_CD','O3_CD'});


% error analysis 
er_ad{n_inst} = D083_orig(:,{'Time','Date','M_AD','Mu_AD','O3_AD'});
er_cd{n_inst} = D083_orig(:,{'Time','Date','M_CD','Mu_CD','O3_CD'});
% t = dobson(:,2);
% m = dobson(:,23);
% mu = dobson(:,22);
% Omega = dobson(:,24);



%AD.O3_AD=1000*AD.O3_AD;
AD.O3_STD=zeros(size(AD.O3_AD));
x1=AD(AD.M_AD<3.0,{'Time','Date','O3_AD','O3_STD','M_AD'});
%CD.O3_CD=1000*CD.O3_CD;
CD.O3_STD=ones(size(CD.O3_CD));
x2=CD(CD.M_CD>2.5,{'Time','Date','O3_CD','O3_STD','M_CD'});
x1.O3=x1.O3_AD;
x2.O3=x2.O3_CD;
x1.M=x1.M_AD;
x2.M=x2.M_CD;
x2.Time=x2.Time+datenum(0,0,0,0,30,0);
x2.Properties.RowNames=cellstr(x2.Time);
t_dob=union(x1(:,{'Time','Date','O3','O3_STD','M'}),x2(:,{'Time','Date','O3','O3_STD','M'}));

t_set_2{n_inst}=sortrows(t_dob,1);


