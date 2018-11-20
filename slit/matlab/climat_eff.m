cd f:\red_brewer\sodankyla\clima


%% ozone climatology
[a,b]=mhdrload('zm_month.n7t');
%elimanmos las medias 
b(37:39:end:39:39:end,:,:)=[];
ozone=b(:,3:end,:);
lat=-87.5:5:87.5;
ozone(ozone==0)=NaN;
%average over years
ozone_climat=nanmean(ozone,3);
figure;
[c,h]=contourf(1:12,lat,ozone_climat,150:5:500);
clabel(c,h,200:50:500);
title('TOMS V8 climatology');
shading flat;

%% Terrain ---->surface pressure-> Not used
% tpdat=zeros(360,720);
% fid_1=fopen('terrain.txt','r+');
% tpdat=fscanf(fid_1,'%f',[360,720]);
% fclose(fid_1);
% long_orig=-179.75:.5:179.75;
% lat_orig=-89.75:.5:89.75;
% pcolor(long_orig,lat_orig,tpdat);shading flat;
% colorbar;colormap bone;
% long=-162:36:162;
% lat=-85:10:85;


%% SBUV-UMKHER a priori profile
%
load climat.txt;
o3=reshape(climat',13,18,12);
% climatological ozone amounts in 13 Umkehr layers
% arranged from layer 12 to layer 0, in 18 10-degree
%  latitude bins from south to north, and
%  12 months, January to December (DU) (input)
%  13 umk layers 
%  18 latitude bands
  umk=12:-1:0;
  p_umk=2.^(-umk); %atm
  h = 7.996;
  z_umk=h*-log(p_umk);
  dz=mean(-diff(z_umk));
  z=z_umk+dz/2; %punto medio del intervalo
  lat=-85:10:85;
  mes=1:12;
% layer 0-> surface to 0.5
% total ozone 
    ozone_sbuv=squeeze(sum(o3,1)); % mean ozone
% nean profile    
    umk_sbuv=squeeze(mean(o3,3));  % mean profile

% Test
    F=figure;
    %ilat=fix(rand*18)+1;
    %imes=fix(rand*12)+1;
    % sodankyla
    ilat=16;
    imes=4;
    oz=o3(:,ilat,imes);
    % Zeff para que sea funcion x= z y= ozono
    zeff=trapz(z,matmul(z,oz'))/trapz(z,oz);
    zeff2=sum(matmul(z,oz'))/sum(oz);
    zeff3=sum(matmul(z_umk,oz'))/sum(oz);
    stairs(oz,z_umk,'o-');
    hold on;
    hline(zeff,'.-b',sprintf(' Trapezoidal rule %.1f Km',zeff));
    hline(zeff2,'.-r',sprintf('Sum aprox %.1f Km',zeff2));
    hline(zeff3,'.:k',sprintf('zo=0 %.1f Km',zeff3));
    plot(oz,z);
    title(sprintf('Ozone SBUV, lat=%.0f Month=%d TO_3 %.1f',lat(ilat),imes,sum(oz)));
    xlabel('ozone (DU)');
    ylabel('height(Km)');
    bold;
 
%effective altitude
 h_sbuv=[];
 %warning to many NaNs 
 for j=1:12
         h_sbuv(j,:)=NaN*ones(1,18);
         h_sbuv(j,:)=trapz(z,matmul(z,o3(:,:,j)')')./trapz(z,o3(:,:,j));
 end        
 
 for j=1:12
         h_sbuv2(j,:)=NaN*ones(1,18);
         h_sbuv2(j,:)=sum(matmul(z,o3(:,:,j)')')./sum(o3(:,:,j));
 end
  %3D
    figure
    [c,h]=contourf(1:12,lat,h_sbuv2',[10:.5:35]);colorbar;
    axis_month;
    clabel(c,h,[20:2.5:35]);shading flat;
    title('Effective heihgt SBUV a priori profile');
  %2D   
    figure
    title('SBUV a priori profile Ozone Efective height');
    [c,h]=contourf(lat,z',umk_sbuv,0:2.5:70);
    hold on;
    clabel(c,h,0:10:70);colorbar;shading flat;
    h_dob = (26 -abs(lat)/10);
    plot(lat,mean(h_sbuv2,1),'py-');
    %plot(lat,mean(h_sbuv,1),':'); 
    hold on;
    plot(lat,h_dob,'g-');
    hline(22,'-w',' H_e_f_f brewer ')
    set(gca,'YLim',[15,45])
    xlabel('Latitude')
    ylabel('Height (Km)')
    title('SBUV a priori profile Ozone Efective height (H_e_f_f)');
    legend('SBUV Ozone (UD)','H_e_f_f SBUV','H_e_f_f  Dobson');
    bold;
  %2D  
%   figure
%    for m=1:12,
%     
%        subplot(3,4,m);
%        title('SBUV a priori profile Ozone Efective height');
%        [c,h]=contourf(lat,z',umk_sbuv,0:5:70);
%        hold on;
%        clabel(c,h,0:10:70);colorbar;shading flat;
%        h_dob = (26 -abs(lat)/10);
%        plot(lat,h_sbuv2(m,:),'py-');
%        hold on;
%        plot(lat,h_dob,'g-');
%        hline(22,'-w',' H_e_f_f brewer ')
%        set(gca,'YLim',[15,45])
%        xlabel('Latitude')
%        ylabel('Height (Km)')
%        bold;
%        title(Month(m,:)); 
%    end
%        suptitle('SBUV a priori profile Ozone Efective height (H_e_f_f) ');
%        legend('SBUV Ozone (UD)','H_e_f_f SBUV','H_e_f_f  Dobson');
%% temperature climatology
 temp=textread('CLIM_TEMP.txt','','whitespace','JanFebMarAprMayJunJulAguSepOctNovDec','emptyvalue',999.00); 
 tmp=reshape(temp',11,18,12);
 % tmp=(umk,lat,time)
 %temp(:,:,12)' ej es la matriz de diciembre del fichero
 %figure;  
 % for m=1:12,
 %  subplot(3,4,m);contourf(lat,z,tmp(:,:,m)-273);colorbar;
 %  title(sprintf('Month= %d ',m));

 %Sbuv layers are inverted  
 o3_sbuv=o3;
 %sum the upper layers
 o3_sbuv(3,:,:,:)=o3_sbuv(1,:,:,:)+o3_sbuv(2,:,:,:)+o3(3,:,:,:);
% inverse the order
 o3_sbuv=o3(13:-1:3,:,:,:); 
 z2=z(end:-1:3);

  figure(F);
     oz_t=o3_sbuv(:,ilat,imes);
%     Zeff para que sea funcion x= z y= ozono
%     zeff=trapz(z2,matmul(z2,oz'))/trapz(z2,oz);
%     zeff2=sum(matmul(z2,oz'))./sum(oz);
%      hp=stairs(oz,z,'_');
%     hold on;
%     hline(zeff,'.:b ',sprintf('\n->%.1f Km',zeff));
%     hline(zeff2,'.-m',sprintf('\n->%.1f Km',zeff2));
      plot(oz_t,z2,'rp','LineWidth',3);
%     title(sprintf('Ozone SBUV, lat=%.0f Month=%d TO_3 %.1f',lat(ilat),imes,sum(oz)));
%     xlabel('ozone (DU)');
%     ylabel('height(Km)');
  
% temperatura efectiva
% int[ t(z)*ozo(z) dz ]/ int[ozo(z) dz]
t_sbuv=[];
for j=1:12  
         t_sbuv(j,:)=NaN*ones(1,18);
         t_sbuv(j,:)=sum(matmul(tmp(:,:,j),o3_sbuv(:,:,j)))./sum(o3_sbuv(:,:,j));
end
 
figure;
[c,h]=contourf(1:12,lat',t_sbuv',180:2.5:250);colorbar;
clabel(c,h,'FontSize',11,'LabelSpacing',172)
shading flat;
axis_month;
title('SBUV Ozone effective temperature')
ylabel('latitude');

%% 2D
figure
gris_line(12);plot(lat,t_sbuv');
hold on;
plot(lat,mean(t_sbuv),'b','LineWidth',3);
xlabel('lat');
ylabel('Temperature Kº');
h=hline(273-46,'r-','Brewer');set(h,'LineWidth',2)
set(gca,'Xlim',[-90,90])
title('SBUV V8 climatology Ozone Efective Temperature  (T_e_f_f)');
legend(strvcat(Meses,'mean'),-1);
bold


 
%% TOMS CLIMATOLOGY  
%
% 0:10 umkher layers ozone 18 lat  10 ozone 12 months
% umk  layers
%
n_umk=0:10;
p_lev=2.^(-n_umk); %atm
%pressure scale height at z = 0 & standard
%temperature (cm) (i.e., p proportional  to exp(-z/h))
h = 7.996; 
z_umk=h*-log(p_lev); 
% half of umkher interval
z=z+dz/2; % interval half point
lat=-85:10:85;

% Read de climatology
%ozone=textread('OZON_CLIM.txt','','whitespace',' Month=Latitude','emptyvalue',999.00);
%ancho de los campos
w=[8,7,7,7,7,7,7,6,6,6,6,6];
ozone=fixed_width_import('OZON_CLIM.txt',1,2376,w);
ozo=ozone;
ozone(1:11:end,:)=[]; % header whith latitude and month
ozone(ozone==999.0)=NaN; % 999.0  to NaN
% fisrt column  total ozone
t_ozone=unique(ozone(~isnan(ozone(:,1)),1));
ozone(:,1)=[];  
o3=reshape(ozone',11,10,18,12);

% o3(:,:,lat,mes)' obtenemos las matrices del fichero.
% o3(ozo,umk,lat,mes)
% INDEX level,t_ozo,lat,time
%% TEST  
%  figure 
%  o3_t=squeeze(nanmean(o3,4));
%  for i=1:12
%      subplot(4,3,i)
%      o3_l=squeeze(nanmean(o3(:,:,:,i),2));
%      contourf(lat,z,o3_l);
%      title(Meses(i,:));
%  end
%%   zo exactly must be topography mean height
%% test 2
 figure;
 oz=o3(:,:,ilat,imes);
 c=jet(12);
 zeff=sum(matmul(oz,z2'))./sum(oz);
 plot(oz,z2,'-');
 hold on;
 for ii=1:length(t_ozone)
     hl=hline(zeff(ii),'-');
     set(hl,'color',c(ii,:));   
 end
 bold
 legend(num2str(t_ozone),-1);
 
 figure
 for ii=1:length(t_ozone)
     subplot(2,5,ii)
     hl=stairs(oz(:,ii),z_umk,'-');
     set(hl,'color',c(ii,:));
     hl=hline(zeff(ii),'-');
     set(hl,'color',c(ii,:));   
 end
 title(sprintf('Ozone SBUV, lat=%.0f Month=%d ',lat(ilat),imes));
    xlabel('ozone (DU)');
    ylabel('height(Km)');
    bold; 
  
%% effective altitude
%  h_eff=[];
%  for i=1:18
%      for j=1:12
%          h_eff(i,j,:)=NaN*ones(1,10);
%          h_eff(i,j,:)=trapz(z,matmul(z,o3(:,:,i,j)')')./trapz(z,o3(:,:,i,j));
%      end
%  end
% % h_eff(lat,time,long)  

%% INTERPOLATION 
% o3
% o3(:,:,lat,mes)' obtenemos las matrices del fichero.
% o3(umk,ozo,lat,mes)
 h3=[];
for i=1:10 %ozone
      for j=1:18 %lat 
          for k=1:12  %month         
          h3(i,j,k)=nansum(z'.*o3(:,i,j,k))./nansum(o3(:,i,j,k));
          %h_eff(i,j,:)=trapz(z,matmul(z,o3(:,:,i,j)')')./trapz(z,o3(:,:,i,j));
          end
      end
  end


% h_eff(ozone,lat,time)
h_eff=h3;    
h_eff_t=squeeze(nanmean(h_eff,3)); % time average
    figure
    contourf(t_ozone,lat,h_eff_t');colorbar;
    title('Effective heihgt TOMS climatology');
%ozone mean-> from SBUV
ozone_mean=ozone_sbuv;
for i=1:18
    for j=1:12
    h_ef(i,j)=interp1(t_ozone,h_eff(:,i,j),ozone_mean(i,j));
    end
end
gris_line(2*18);plot(h_ef');     
legend(num2str(lat'),-1);
title('TOMS V8 climatology  Ozone Efective height');
axis_month;
ylabel('Height (Km)');

figure

 [c,h]=contourf(1:12,lat,h_ef,10:.5:40);
 hold on;
 clabel(c,h,15:2.5:30);colorbar;shading flat;
 title('TOMS a priori profile Ozone Efective height');
 axis_month
 
 
 

figure
title('TOMS V8 climatology  Ozone Efective height');
    [c,h]=contourf(lat,z_umk',umk_sbuv,0:2.5:70);
    hold on;
    clabel(c,h,0:10:70);colorbar;shading flat;

     h = (26 -abs(lat)/10);
     plot(lat,nanmean(h_ef')','yo-');
     plot(lat,mean(h_sbuv2,1),'cp-');

     hold on;
     plot(lat,h,'g-');
     hline(22,'-w','brewer op')

     set(gca,'YLim',[15,45])
     xlabel('Latitude')
     ylabel('Height (Km)')
     title('TOMS V8 climatology Ozone Efective height (H_e_f_f)');
     legend('SBUV Ozone (UD)','H_e_f_f TOMS','H_e_f_f SBUV','H_e_f_f  Dobson');
     bold;
   

%% temperature climatology
 temp=textread('CLIM_TEMP.txt','','whitespace','JanFebMarAprMayJunJulAguSepOctNovDec','emptyvalue',999.00); 
 tmp=reshape(temp',11,18,12);
 % tmp=(umk,lat,time)
 %temp(:,:,12)' ej es la matriz de diciembre del fichero
 %  figure;
 %  contourf(lat,z,tmp(:,:,12)-273);colorbar;
 %  title('temperature dic')

%% temperatura efectiva
 % int[ t(z)*ozo(z) dz ]/ int[ozo(z) dz]
 
 t_eff=[];
 % ozone climatology 
 % temperature climatology 
 
 t_eff=[];
  for i=1:10 %ozone 
      for j=1:18 %lat
          for k=1:12  %month          
          temp=tmp(:,j,k);    
          t_eff(i,j,k)=sum(temp.*o3(:,i,j,k))./sum(o3(:,i,j,k));
          end
      end
  end

%ozone mean-> from SBUV
ozone_mean=ozone_sbuv;
t_ef=[];
 for i=1:18
    for j=1:12
      t_ef(i,j)=interp1(t_ozone,t_eff(:,i,j),ozone_mean(i,j));
    end
end

figure;
[c,h]=contourf((0:11),lat,t_ef,180:2.5:250);colorbar;
clabel(c,h,'FontSize',11,'LabelSpacing',172)
shading flat;
axis_month;
title('Ozone effective temperature')
ylabel('latitude');

%% 2D
figure
gris_line(12);plot(lat,t_ef');
hold on;
plot(lat,mean(t_ef,2),'b','LineWidth',3);
xlabel('lat');
ylabel('Temperature Kº');
h=hline(273-46,'r-','Brewer');set(h,'LineWidth',2)
set(gca,'Xlim',[-90,90])
title('TOMS V8 climatology Ozone Efective Temperature  (H_e_f_f)');
legend(strvcat(Meses,'mean'),-1);
bold
%% Comparison SBUV 