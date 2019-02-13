%% Imput data 

clear variables;
 read_dobson  % Ulf obseravtions
 d=Dobson2012;
 d(:,2)=24*(d(:,1)-fix(d(:,1))); % GMT hour
 td=round(d(:,1)*24*60/5);
 
 %aod data
 load('aod.mat');
 taod=round(aod(:,1)*24*60/5);
 dobson=scan_join([td,d],[taod,aod(:,[1,8])]);
 dobson(1:find(~isnan(dobson(:,2)),1,'first')-1,:)=[];

 % brewer data
 load('Izo2012_sumold.mat')
  brd=Izo2012_sumold;
 % ozone comparison
 
 o3d=dobson(:,[2,21,24]);
 o3d(:,2:3)=o3d(:,2:3)*1000; %DU
 o3b=brd.ozone;
 
 %% daily comparison
 c_o3=scan_join(o3d,o3b);
 c_o3=c_o3(~isnan(c_o3(:,2)),:);
 ref=nanmean(c_o3(:,4:end),2);
 [m,s,n,g]=grpstats([diaj(c_o3(:,1)),matdiv(100*matadd(c_o3(:,2:end),...
     -ref),ref)],diaj(c_o3(:,1))');
errorbar(repmat(m(:,1),1,5),m(:,2:end),s(:,2:end))
hline(0)
grid
legend('CD','AD','157','183','185');

title('Daily comparsion whith the mean of the Brewer');

days=unique(fix(m(:,1)))
%% ozone slant path comparison
o3d=dobson(:,[2,19,22,21,24]);
o3d(:,4:end)=o3d(:,4:end)*1000; %DU
o3b=[brd.ozone,brd.airm];
 
c_o3=scan_join(o3d,o3b);
c_o3=c_o3(~isnan(c_o3(:,2)),:);
ref_o3=nanmean(c_o3(:,6:8),2);
 
% [m,s,n,g]=grpstats([diaj(c_o3(:,1)),matdiv(100*matadd(c_o3(:,4:8),...
%     -ref_o3),ref_o3)],diaj(c_o3(:,1))');
% figure;
% errorbar(repmat(m(:,1),1,5),m(:,2:end),s(:,2:end))

ref_osc =nanmean(c_o3(:,10:end),2).*ref_o3;
osc_g=unique(fix(ref_osc/50)*50)';
[m,s,n,g]=grpstats([ref_osc,matdiv(100*matadd(c_o3(:,[4:8]),...
     -ref_o3),ref_o3)],fix(ref_osc/50)*50);
figure;
errorbar(repmat(m(:,1),1,size(m,2)-1),m(:,2:end),s(:,2:end))
hline(0)
title('Brewer Dobson Daily ozone slant path comparison  Brewer');
grid
legend('CD','AD','157','183','185');

brw_str={'157','183','185'}

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
 
 RC=matmul(B1,m)*770/1013;
 o3d= matmul(A1,nu);
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
 grpstats(100*(x1-xf)./xf,diaj(dobson(:,2)),0.05); % double pair ?
 title(' % Dobson ozone recalculation - ozone from file / ozone from file');
 legend(PAIR);
% plot(dobson(:,4),100*(x1-xf)./xf); % double pair ?
 
%  legend(PAIR);
%  title({ '% Ozone recalculate - Ozone file / Ozone file',...
%      ' Double pair ozone calculations  are still wrong (?)'});
%  
%% LANGLEY
%  
% figure
% for i=1:5
%  subplot(3,2,i) 
%  grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
%  
%  [dd{i},ds{i},dn{i},dg{i}]=grpstats(dobson,grp);
%  gscatter(nu(:,i),N(:,i)/100-RC(:,i),grp);
%  [h,ld{i}]=rline;
% end
% 
% figure; hold on; etc_d={};
% MK=get(gca,'LineStyleOrder');
% CL=jet(8);
% for i=1:5
%  %figure  
%  DETC=ld{i}(2,:); hold all;
%  DETC(DETC==0)=NaN;
%  DETC(dn{i}(:,7+i)<10)=NaN;
%  dt=dd{i}(:,2);
%  etc_d{i}=[dt,DETC',dd{i}(:,orx(i)),ds{i}(:,orx(i)),dn{i}(:,orx(i))];
%  
%  hp(i)=plot(diaj2(dt),DETC,MK(i,:));
%  hline(nanmedian(DETC),'-',num2str(round(nanmedian(DETC*1000))/1000));
% end
% box on;
% legend(hp,PAIR);
% title('ETC correction Izaña 2012')



%%
%   DAY LANGLEY PLOT (Brewer ETC ?)
grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
d=unique(grp{1}); 
days=d(~isnan(d));
n=0;B_TABLE=[];B_STABLE=[];

for dj=days'
    %disp(dj)
    n=n+1;
    jd=find(grp{1}==dj);
     grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
   figure;
   for i=1:5
    subplot(2,3,i)
    
    gscatter(nu(jd,i), (N(jd,i)+RC(jd,i)),grp{2}(jd),...
        '','','','off','airm','N-R');
    title(PAIR{i});
    box on;
    % only am
    jam=jd(grp{2}(jd)==0 & nu(jd,i)>=1.2 & nu(jd,i)<=2.4 );
    jpm=jd(grp{2}(jd)==1 & nu(jd,i)>=1.2 & nu(jd,i)<=2.4 );
    plot(nu(jam,i),N(jam,i)/100-RC(jam,i),'rx',...
         nu(jpm,i),N(jpm,i)/100-RC(jpm,i),'c+')
    [B(n,i,1,1:2),BINT1,RAM] = linregress((N(jam,i)/100-RC(jam,i)),nu(jam,i)) ;
    %refline(BINT(:,1))
    %refline(BINT(:,2))
    
    [B(n,i,2,1:2),BINT2,RPM] = linregress((N(jpm,i)/100-RC(jpm,i)),nu(jpm,i)) ;
    %plot(nu(jam,i),RAM,'rx',nu(jpm,i),RPM,'c+')
    title(PAIR(i))
   end
   subplot(2,3,6);
   %plot(nu(jd,i),[xf(jd,:)]);legend(PAIR,-1);
   plot(dobson(jam,2),dobson(jam,end),'r',dobson(jpm,2),dobson(jpm,end),'b');
   
   set(gca,'Ylim',[0,0.1]);title('AOD 340 nm'); axis('tight');datetick('x','keeplimits');
   suptitle(sprintf('Day %d',dj));
   
   
   BETC=[[days(n)+.25,B(n,:,1,2),length(jam),nanmean(dobson(jam,end)),nanstd(dobson(jam,end))]...
       ;[days(n)+.75,B(n,:,2,2)],length(jpm),nanmean(dobson(jpm,end)),nanstd(dobson(jpm,end))];
   BETC(BETC==0)=NaN;
   B_TABLE=[B_TABLE;BETC];
   SETC=[[days(n)+.25,B(n,:,1,1),length(jam),nanmean(dobson(jam,end)),nanstd(dobson(jam,end))]...
       ;[days(n)+.75,B(n,:,2,1)],length(jpm),nanmean(dobson(jpm,end)),nanstd(dobson(jpm,end))];
   SETC(SETC==0)=NaN;
   B_STABLE=[B_STABLE;SETC];
   
   
end



%% ETC

figure;
BTABLE=B_TABLE;
BTABLE(:,2:end-3)=100*B_TABLE(:,2:end-3);
BTABLE=[BTABLE(:,1:6),BTABLE(:,2)-BTABLE(:,3),BTABLE(:,4)-BTABLE(:,3),BTABLE(:,7:end)]
ploty(BTABLE(1:end-1,1:end-3));
hline(nanmean(BTABLE(1:end-1,2:end-3)),'',cellstr(num2str(nanmean(BTABLE(1:end-1,2:end-3))')))
legend([PAIR,'C-D','A-D'],-1)
title(' 100*ETC    (N/100- R) vs nu');
xlabel('day of the year');
  

displaytable(BTABLE,{'Day',PAIR{:},'C-D','A-D','Nobs','AOD 340','AOD340 std'},6,'.2f') 

%% SLOPE
figure;
ploty(B_STABLE(1:end-1,1:end-3));
hline(nanmean(B_STABLE(1:end-1,2:end-3)),'',cellstr(num2str(nanmean(B_STABLE(1:end-1,2:end-3))')))
legend(PAIR,-1)
title(' SLOPE  P(N/100- R) vs nu');
xlabel('day of the year');
%printmatrix(B_STABLE);

%B_STABLE(:,2:end-3)=100*B_STABLE(:,2:end-3);
displaytable(B_STABLE,{'Day',PAIR{:},'Nobs','AOD 340','AOD340 std'},6,'.2f') 
%% Dobson langley 1/nu
figure
G=[];
for i=1:5
 %subplot(3,2,i) 
 figure
 grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
 
 [dd{i},ds{i},dn{i},dg{i}]=grpstats(dobson,grp);
 js=nu(:,i)>=1.2 & nu(:,i)<=2.4;
 %gscatter(1./nu(:,i),(N(:,i)/100+RC(:,i))./nu(:,i),grp,'','o+');
 box on;
 title(PAIR(i));
 hold on;
 plot(1./nu(js,i),(N(js,i)/100-RC(js,i))./nu(js,i),'p')

 [h,G(i,:)]=robust_line;
end


%% DAY PLOT
grp={diaj(dobson(:,2)),24*(dobson(:,2)-fix(dobson(:,2)))>13};
d=unique(grp{1}); 
days=d(~isnan(d));
D=[];n=0; C_TABLE=[];C_STABLE=[];
for dj=days'
    %disp(dj)
    jd=find(grp{1}==dj );
    
    n=n+1;
   figure;
   for i=1:5
    subplot(2,3,i)
    
       % only am
    jam=jd(grp{2}(jd)==0 & nu(jd,i)>=1.2 & nu(jd,i)<=2.4 );
    jpm=jd(grp{2}(jd)==1 & nu(jd,i)>=1.2 & nu(jd,i)<=2.4 );
    
    [D(n,i,1,1:2),DINT,R] = linregress((N(jam,i)/100+RC(jam,i))./nu(jam,i),1./nu(jam,i)) ;
    refline(DINT(:,1))
    refline(DINT(:,2))
    
    [D(n,i,2,1:2),DINT,R] = linregress((N(jpm,i)/100+RC(jpm,i))./nu(jpm,i),1./nu(jpm,i)) ;
    refline(DINT(:,1))
    refline(DINT(:,2))
    if(D(n,i,2,1)==0)D(n,i,2,1:2)=[NaN,NaN];end
    plot(1./nu(jam),(N(jam,i)/100+RC(jam,i))./nu(jam,i),'bo',...
         1./nu(jpm),(N(jpm,i)/100+RC(jpm,i))./nu(jpm,i),'r+');
    title(PAIR{i});
    box on;
    axis('tight');

    
    
    
   end
   subplot(2,3,6);
   %plot(1./nu(jd,i),[x(jd,:)]);legend(PAIR,3);
   %suptitle(sprintf('Day %d',dj));
   plot(dobson(jam,2),dobson(jam,end),'r',dobson(jpm,2),dobson(jpm,end),'b');
  
   set(gca,'Ylim',[0,0.1]);title('AOD 340 nm'); axis('tight');datetick('keeplimits','keeplimits');
   suptitle(sprintf(' Langley P vs 1/nu Day %d',dj));
   
   
   DETC=[[days(n)+.25,D(n,:,1,2),length(jam),nanmean(dobson(jam,end)),nanstd(dobson(jam,end))]...
       ;[days(n)+.75,D(n,:,2,2)],length(jpm),nanmean(dobson(jpm,end)),nanstd(dobson(jpm,end))];
   C_TABLE=[C_TABLE;DETC];
   
   SETC=[[days(n)+.25,D(n,:,1,1),length(jam),nanmean(dobson(jam,end)),nanstd(dobson(jam,end))]...
       ;[days(n)+.75,D(n,:,2,1)],length(jpm),nanmean(dobson(jpm,end)),nanstd(dobson(jpm,end))];
   C_STABLE=[C_STABLE;SETC];
   
end


%%
figure;
ploty(C_TABLE(1:end-1,1:end-3));
hline(nanmean(C_TABLE(1:end-1,2:end-3)),'',cellstr(num2str(nanmean(C_TABLE(1:end-1,2:end-3))')))
legend(PAIR,-1)
title(' ETC P vs 1/nu');
xlabel('day of the year');
  
printmatrix(C_TABLE);

%%
figure;
plot(C_STABLE(1:end-1,1),C_STABLE(1:end-1,2:end-3)*100);
hline(nanmean(C_STABLE(1:end-1,2:end-3))*100,'',cellstr(num2str(100*nanmean(C_STABLE(1:end-1,2:end-3))')))
legend(PAIR,-1)
title(' SLOPE P vs 1/nu');
xlabel('day of the year');



%%
dobson_res; 
figure
PAIR={'C','D','A','CD','AD','C-D','A-D'};
CSTABLE=[C_STABLE(:,1),100*C_STABLE(:,2:6),100*(C_STABLE(:,2)-C_STABLE(:,3)),...
          100*(C_STABLE(:,4)-C_STABLE(:,3)),C_STABLE(:,7:end)]
%Brewer      
displaytable(CSTABLE,{'Day',PAIR{:},'Nobs','AOD 340','AOD340 std'},6,'.2f')  
%Dobson
displaytable(dobson_izo(:,1:end-1),{'Day',PAIR{:}},6,'.2f')  
CSTABLE(end-4,:)=[];     
[di,bj]=ismember(dobson_izo(:,1),CSTABLE(:,1));

aux=[CSTABLE(bj(di),1:end-3),dobson_izo(di,2:end-1)];
aux(:,2:8)=-aux(:,2:8);
displaytable([CSTABLE(bj(di),1:end-3),dobson_izo(di,2:end-1)],{'Day',PAIR{:},PAIR{:}},6,'.2f')  
     
for p=2:8
   
figure;
plot(dobson_izo(1:end-1,1),-dobson_izo(1:end-1,p),'r');
hold on;
plot(CSTABLE(:,1),CSTABLE(:,p),'b');
legend('RDCC','RBCC')

hline(mean(-dobson_izo(2:end,p)),'r',num2str(mean(dobson_izo(:,p))));
hline(nanmean(CSTABLE(:,p)),'b',...
    num2str(nanmean(CSTABLE(:,p))));


title({[ PAIR{p-1}, ' P vs 1/nu Slope*100']  ,...
    sprintf('RDCC  %f    RBCC  %f',...
          mean(-dobson_izo(2:end,p)),... 
          nanmean(CSTABLE(:,p)))});
end

%%

% figure
% %aod340=aod(:,[1,8]);
% %j=ismember(aod340,diaj(aod(:,1)));
% [maod,saod]=grpstats(dobson(:,end-1:end),grp)
% 
% plot(days,D(:,:,2,2),'o',days,D(:,:,1,2),'+');
% hline(nanmean([D(:,:,2,2),D(:,:,2,2)]),,num2str(round(nanmedian(DETC)*100)/100))
%%
% figure
% for i=1:5
%  %figure  
%  DETC=ld{i}(2,:); hold all;
%  DETC(DETC==0)=NaN;
%  DETC(dn{i}(:,7+i)<10)=NaN;
%  dt=dd{i}(:,2);
%  etc_d{i}=[dt,DETC',dd{i}(:,orx(i)),ds{i}(:,orx(i)),dn{i}(:,orx(i))];
%  
%  hp(i)=plot(diaj2(dt),DETC,MK(i,:),'Color',CL(i,:));
%  hline(nanmedian(DETC),'-',num2str(round(nanmedian(DETC)*100)/100));
% end
% box on;
% legend(hp,PAIR);
% 
% f=figure
% for i=1:5
%  %figure  
%  DETC=ld{i}(1,:); hold all;
%  DETC(DETC==0)=NaN;
%  DETC(dn{i}(:,7+i)<10)=NaN;
%  dt=dd{i}(:,2);
%  etc_d{i}=[dt,DETC',dd{i}(:,orx(i)),ds{i}(:,orx(i)),dn{i}(:,orx(i))];
%  figure(f);
%  hp(i)=plot(diaj2(dt),DETC,MK(i,:),'Color',CL(i,:));
%  hline(nanmedian(DETC),'-',num2str(round(nanmedian(DETC)*100)/100));
% end
% box on;
% legend(hp,PAIR);