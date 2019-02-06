 read_dobson
 d=Dobson2012;
 d(:,2)=24*(d(:,1)-fix(d(:,1))); % GMT hour
 td=round(d(:,1)*24*60/5);
 
 %aod data
 %load aod;
 %taod=round(aod(:,1)*24*60/5);
 dobson=scan_join([td,d],[td,d]);
 %dobson(1:find(~isnan(dobson(:,2)),1,'first')-1,:)=[];
 
 %% try to reproduce dobson calculations
 PAIR={'C','D','A'};
 % pairs-> mu,m,R,N,O3
 % C pair 3-7
 % D pair 8-12
 % A pair 13-17
 %  Komhyr 1983
 %A1=[0.800,0.360,1.748];
 %B1=[0.110,0.104,0.116];
 % Komhyr 1993
 
 A1=[0.833,0.374,1.806];
 B1=[0.109,0.104,0.114];
 
 figure; 
 for i=0:2
     i
     subplot(3,1,(i+1));
 
  N=dobson(:,7+i*5);
  nu=dobson(:,4+i*5);
  m=dobson(:,5+i*5);
  RC=B1(i+1)*m*770/1024;
  x=(N+RC)./(A1(i+1)*nu)/100;
  xf=dobson(:,8+5*i);
  plot(dobson(:,2),[x./xf]);
  legend('recalculated','file')
  title(sprintf(' %s alpha= %f beta = %f',PAIR{i+1},A1(i+1),B1(i+1)));
 end
 
 %%  langley
 % pairs-> mu,m,R,N,O3
 % C pair 4-8
 % D pair 9-13
 % A pair 14-18
 % C pair 
 %1983
%  Ac=[1.748,0.491,0.375,0.116,0.0664];
%  Bc=[1.140,0.470,0.357,0.113,0.0991];
%  Cc=[0.800,0.491,0.375,0.116,0.1375];
%  Dc=[0.800,0.491,0.375,0.116,0.1375];
%  

 % single pair A
 % Na - Nd
 L=-dobson(:,12)+dobson(:,17);
 % Nd
 N=dobson(:,17);
 % mu a
 nu=dobson(:,14);
 % m a
 m=dobson(:,15);
 
% por que 0.116??? si es beta(a) - beta(d) seria 0.01. Por que m y no
% m(a)-m(b)??
 RC=0.116*m*770/1013;
 
 % no faltaria dividir por alfa(a)-alfa(d)? por que mu y no mu(a)-mu(d)?
 P=(L-RC)./nu;
 figure
 plot(1./nu,P)
 %%
 figure
 gscatter(nu,L+RC,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 [h,r]=rline;
 [m,s,n,g]=grpstats([dobson(:,[2:3,18]),dobson(:,end)],{diaj(dobson(:,2)),dobson(:,3)>13});
 
 %%  ETC
 figure;
 plot(m(:,1),r(2,:),'ro-');
 hline(median(r(2,:)),'k',num2str(median(r(2,:))));
 datetick('keepticks')
 title(' ETC_{AD}')
 
 
 %% ozone absortion coefficient
 figure;
 plot(m(:,1),r(1,:)'./m(:,3)/100,'ro-');
 hline(1.748);
 datetick('keepticks')
 title(' \alpha AD')
 %% 1/nu
 
 
 figure
 gscatter(1./nu,P,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 xlabel('1/\mu')
 [h,rp]=rline;

 
 
 
 
 
 %%double pairs
 L=-dobson(:,12)+dobson(:,17);
 nu=mean(dobson(:,[10,15]),2);
 m=mean(dobson(:,[11,16]),2);
 RC=0.012*m*770/1013;
 P=(L-RC)./nu;
 %plot(1./dobson(:,4),P)
 figure
 gscatter(nu,L+RC,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 [h,r]=rline;
 [m,s,n,g]=grpstats([dobson(:,[2:3,end-2]),dobson(:,end)],{diaj(dobson(:,2)),dobson(:,3)>13});
 
 %%  ETC
 figure;
 plot(m(:,1),r(2,:),'ro-');
 hline(median(r(2,:)),'k',num2str(median(r(2,:))));
 datetick('keepticks')
 title(' ETC_{AD}')
 
 
 %% ozone absortion coefficient
 figure;
 plot(m(:,1),r(1,:)'./m(:,3)/100,'ro-');
 hline(1.388);
 datetick('keepticks')
 title(' \alpha AD')
 %% 1/nu
 
 
 figure
 gscatter(1./nu,P,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 xlabel('1/\mu')
 
 
 
 
 
 %%double pair
 L=-dobson(:,12)+dobson(:,17);
 nu=mean(dobson(:,[10,15]),2);
 m=mean(dobson(:,[11,16]),2);
 RC=0.012*m*770/1013;
 P=(L-RC)./nu;
 %plot(1./dobson(:,4),P)
 figure
 gscatter(nu,L+RC,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 [h,r]=rline;
 [m,s,n,g]=grpstats([dobson(:,[2:3,end-2]),dobson(:,end)],{diaj(dobson(:,2)),dobson(:,3)>13});
 
 %%  ETC
 figure;
 plot(m(:,1),r(2,:),'ro-');
 hline(median(r(2,:)),'k',num2str(median(r(2,:))));
 datetick('keepticks')
 title(' ETC_{AD}')
 
 
 %% ozone absortion coefficient
 figure;
 plot(m(:,1),r(1,:)'./m(:,3)/100,'ro-');
 hline(1.388);
 datetick('keepticks')
 title(' \alpha AD')
 %% 1/nu
 
 
 figure
 gscatter(1./nu,P,{diaj(dobson(:,2)),dobson(:,3)<13},'','o+');
 xlabel('1/\mu')
 
% diference between RC and NC  -(N are N values, R ?)
gscatter(d(:,1)-fix(d(:,1)),d(:,5),diaj(d(:,1)),'','.'); hold on;
gscatter(d(:,1)-fix(d(:,1)),d(:,10),diaj(d(:,1)),'','+');
gscatter(d(:,1)-fix(d(:,1)),d(:,15),diaj(d(:,1)),'','o');
datetick;

%% Langleys\[\frac{1}{\mu }\]
