fpath='./dobson_data'
l=dir(fullfile(fpath,'OD16*BPeffcalc.064'));
%l=dir(fullfile(fpath,'OD16*_oriBPnom.064'));
dobson=[];
for i=l'
  dob=textread(fullfile(fpath,i.name),'','headerlines',6,'delimiter',',:', 'emptyvalue',NaN);
  fechas=sscanf(i.name,'%*2c%02d%03d.%03d');
   fecha=dia_mes(fechas(1),fechas(2));
   fecha=fecha(1)+dob(:,1)/24+dob(:,2)/24/60+dob(:,3)/24/60/60;
   
   dobson=[dobson;[fecha,dob]];
end
   
% error code
[i,j]=find(dobson==-0.999);
dobson(i,:)=[];

%save D064orig.dat dobson -ascii -double 
%dob=readtable('OD16258Header.064','FileType','text');
%leg='Date HH MM SS Mu_C M_C RC NC O3_C Mu_D M_D RD ND O3_D Mu_A M_A RA NA O3_A Mu_CD M_CD O3_CD Mu_AD M_AD O3_AD';
%leg=mmcellstr(strrep(leg,' ','|'));
%D064_orig=array2table(dobson,'VariableNames' ,leg,'RowNames',cellstr(datestr(dobson(:,1))));
%t=readtable('Atmoz_o3_set1.xls','Sheet','D083'); t2=table2timetable(t,'RowTimes',datetime(t.Row))
%writetable(D064_orig,'Atmoz_o3_set1.xls','Sheet','D064','WriteRowNames',true);


%%
% dobson equation
%  L (N1-N2) - B1-B2 m p/p0  + S 1/nu = A1-A2  x12 
 % pairs-> mu,m,R,N,O3
 % C pair 4-8
 % D pair 9-13
 % A pair 14-18
 %  Komhyr 1983 xsect
    %A1=[0.800,0.360,1.748];
    %B1=[0.110,0.104,0.116];
PAIR={'C','D','A','CD','AD'};
 
 A1=[0.833,0.374,1.806,0.459,1.432];
 B1=[0.109,0.104,0.114,0.005,0.010];
 i=0:2;
 orx=[8,13,18,21,24];

 %Lamp Corrections:,  -1.3,  -1.3,  -1.9
 %lamp_corr=[  -1.3,  -1.3,  -1.9];
 
 
 % N values
 N=dobson(:,8+i*5);
 N=[N,N(:,1)-N(:,2),N(:,3)-N(:,2)];
 %nu=dobson(:,4+i*5);nu=[nu,mean(nu(:,1:2),2),mean(nu(:,2:3),2)];
 nu=dobson(:,5+i*5);nu=[nu,dobson(:,[20,23])];
 %m=dobson(:,5+i*5);m=[m,mean(m(:,1:2),2),mean(m(:,2:3),2)];
 m=dobson(:,6+i*5);m=[m,dobson(:,[21,24])];
 xf=dobson(:,9+5*i);
 xf=[xf,dobson(:,[22,25])];
 
 % Rayleigh Correction
 RC=matmul(B1,m)*770/1013;
 %  Xsec + airmass
 o3d= matmul(A1,nu);
 % Ozone
 x=(N/100-RC)./o3d;
 
 x1=x;
 %cd 
 % f1 = (((Na / 100) / Amu#) - ((Nd / 100) / Dmu#)) / 1.432 and subsequently:
 %     O3ad = f1 - .007 * Pfact * (Aam# + Dam#) / (Amu# + Dmu#)
 % 0.007 comes from 0.010 / 1.432, Pfact = 770/1013
 
 f1= (N(:,1)./nu(:,1)/100  -N(:,2)./nu(:,2)/100)/A1(4);
 
 f2= (N(:,3)./nu(:,3)/100  -N(:,2)./nu(:,2)/100)/A1(5);
 
 
 x1(:,4)=f1-B1(4)/A1(4)*770/1013* (m(:,1)+m(:,2))./(nu(:,1)+nu(:,2));
 
 x1(:,5)=f2-B1(5)/A1(5)*770/1013* (m(:,3)+m(:,2))./(nu(:,3)+nu(:,2));
 
 figure;
 grpstats(100*(x1-xf)./xf,diaj(dobson(:,1)),0.05); % double pair ?
 title(' % Dobson ozone recalculation - ozone from file / ozone from file');
 legend(PAIR);
 xlabel('day');
 
 
 
 figure;
 grpstats(100*(x1-xf)./xf,fix(dobson(:,5)*4)/4,0.05); % double pair ?
 title(' % Dobson ozone recalculation - ozone from file / ozone from file');
 legend(PAIR)
 xlabel('airm')

 

 figure;
 grpstats(100*(x1-xf)./xf,fix(sza(dobson(:,1))/5)*5,0.05); % double pair ?
 title(' % Dobson ozone recalculation - ozone from file / ozone from file');
 legend(PAIR)
 xlabel('SZA')












% D064_orig.Time=date_time(D064_orig.Date);
% AD=D064_orig(:,{'Time','Date','M_AD','O3_AD'});
% CD=D064_orig(:,{'Time','Date','M_CD','O3_CD'});
% AD.O3_AD=1000*AD.O3_AD;
% AD.O3_STD=zeros(size(AD.O3_AD));
% x1=AD(AD.M_AD<3.0,{'Time','Date','O3_AD','O3_STD','M_AD'});
% CD.O3_CD=1000*CD.O3_CD;
% CD.O3_STD=ones(size(CD.O3_CD));
% x2=CD(CD.M_CD>2.5,{'Time','Date','O3_CD','O3_STD','M_CD'});
% x1.O3=x1.O3_AD;
% x2.O3=x2.O3_CD;
% x1.M=x1.M_AD;
% x2.M=x2.M_CD;
% x2.Time=x2.Time+datenum(0,0,0,0,30,0);
% x2.Properties.RowNames=cellstr(x2.Time);
% t_dob=union(x1(:,{'Time','Date','O3','O3_STD','M'}),x2(:,{'Time','Date','O3','O3_STD','M'}));
% 
% t_set_1{n_inst}=sortrows(t_dob,1);
% 
