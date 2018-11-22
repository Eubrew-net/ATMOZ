%% DOBSON EVALUATION IZO 2012
 clear variables;
 read_dobson
 d=Dobson2012;
 d(:,2)=24*(d(:,1)-fix(d(:,1))); % GMT hour
 td=round(d(:,1)*24*60/5);  % hour to sync
 
  % aod data
 load('aod.mat');
 taod=round(aod(:,1)*24*60/5);
 dobson=scan_join([td,d],[taod,aod(:,[1,8])]);
 dobson(1:find(~isnan(dobson(:,2)),1,'first')-1,:)=[];

%%
%Date;Time Local;Mu_C;M_C;RC;NC;O3_C/1000;Mu_D;M_D;RD;ND;O3_D/1000;Mu_A;M_A;RA;NA;O3_A/1000;Mu_CD;M_CD;O3_CD/1000;Mu_AD;M_AD;O3_AD/1000
%26/09/12; 08:00:00;4,149;4,340;143,1;128,9;0,269;4,090;4,273;88,2;75,1;0,270;4,033;4,208;242,4;231,2;0,267;4,119;4,306;0,268;4,062;4,240;0,267
%26/09/12; 08:03:00;3,978;4,145;138,1;123,7;0,269;3,924;4,084;84,9;71,7;0,269;3,872;4,024;232,3;220,1;0,265;3,951;4,114;0,270;3,898;4,054;0,264


%%
% dobson 
% L (N1-N2) - B1-B2 m p/p0  + S 1/nu = A1-A2  x12 
PAIR={'C','D','A','CD','AD'};
% pairs-> mu,m,R,N,O3
% C pair 4-8      
% D pair 9-13
% A pair 14-18
%  Komhyr 1983
%A1=[0.800,0.360,1.748];
%B1=[0.110,0.104,0.116];
% Komhyr 1993
%  % single pair i=0:2
%  N=dobson(:,7+i*5);
%  nu=dobson(:,4+i*5);
%  m=dobson(:,5+i*5);
%  
 
 A1=[0.833,0.374,1.806,0.459,1.432];
 B1=[0.109,0.104,0.114,0.005,0.010];
 i=0:2;
 orx=[8,13,18,21,24];
 
 N=dobson(:,7+i*5);N=[N,N(:,1)-N(:,2),N(:,3)-N(:,2)];
 %nu=dobson(:,4+i*5);nu=[nu,mean(nu(:,1:2),2),mean(nu(:,2:3),2)];
 nu=dobson(:,4+i*5);nu=[nu,dobson(:,[19,22])];
 %m=dobson(:,5+i*5);m=[m,mean(m(:,1:2),2),mean(m(:,2:3),2)];
 m=dobson(:,5+i*5);m=[m,dobson(:,[20,23])];
 xf=dobson(:,8+5*i);xf=[xf,dobson(:,[21,24])];
 
 %Raleyght
 RC=matmul(B1,m)*770/1013;
 % denominatr
 o3d= matmul(A1,nu);
 x=(N/100-RC)./o3d;
 
 plot(m,100*(x-xf)./xf); % double pair ?
 legend(PAIR)
 
 
 %%
 
 figure
for i=1:5
 subplot(3,2,i) 
 grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
 
 [dd{i},ds{i},dn{i},dg{i}]=grpstats(dobson,grp);
 gscatter(nu(:,i),N(:,i)/100-RC(:,i),grp);
 [h,ld{i}]=rline;
end

figure; hold on; etc_d={};
MK=get(gca,'LineStyleOrder');
CL=jet(8);
for i=1:5
 %figure  
 DETC=ld{i}(1,:); hold all;
 DETC(DETC==0)=NaN;
 DETC(dn{i}(:,7+i)<10)=NaN;
 dt=dd{i}(:,2);
 etc_d{i}=[dt,DETC',dd{i}(:,orx(i)),ds{i}(:,orx(i)),dn{i}(:,orx(i))];
 
 hp(i)=plot(diaj2(dt),DETC,MK(i,:));
 hline(nanmedian(DETC),'-',num2str(round(nanmedian(DETC*1000))/1000));
end
box on;
legend(hp,PAIR);

%%