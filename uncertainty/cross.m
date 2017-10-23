function A_unc = cross(brw)


% carga el archivo de la dispersion.
% col1 es el tiempo
% col2 es el brewer
% col3 es el ajuste (2,3)
% col4 a col9 son las longitudes de onda y de la col10  col15 son las fwhm
%brw_str = num2str(brw)
a=strcat('dsp_',brw,'_1516.csv');
data=csvread(a,1,0);
%data=csvread('dsp_' + brw + '_1516.csv',1,0);

%%
% carga las variables de1 1 al 15
dsp_summ=data(:,:);
%bp_r=data(:,[1,2,16:end]);

% sacamos la sigma de las longitudes de onda y de las anchuras
dsp_lambda = std(dsp_summ(:,4:9));
dsp_fwhm = std(dsp_summ(:,10:15));

% vector aleatorio
rnd = randn(1,1000);

% sigmas aleatorias
sigma_lambda1 = rnd*dsp_lambda(1);
sigma_lambda2 = rnd*dsp_lambda(2);
sigma_lambda3 = rnd*dsp_lambda(3);
sigma_lambda4 = rnd*dsp_lambda(4);
sigma_lambda5 = rnd*dsp_lambda(5);
sigma_lambda6 = rnd*dsp_lambda(6);

sigma_lambda = [sigma_lambda1', sigma_lambda2', sigma_lambda3', sigma_lambda4', sigma_lambda5', sigma_lambda6'];

sigma_fwhm1 = rnd*dsp_fwhm(1);
sigma_fwhm2 = rnd*dsp_fwhm(2);
sigma_fwhm3 = rnd*dsp_fwhm(3);
sigma_fwhm4 = rnd*dsp_fwhm(4);
sigma_fwhm5 = rnd*dsp_fwhm(5);
sigma_fwhm6 = rnd*dsp_fwhm(6);

sigma_fwhm = [sigma_fwhm1', sigma_fwhm2', sigma_fwhm3', sigma_fwhm4', sigma_fwhm5', sigma_fwhm6'];
%%

% pesos del o3 y del so2
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];

% longitudes de onda del Brewer
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

%%
% carga las secciones eficaces del so2 (no lo utilizo)
% load so2;
% so2295sm.dat;  % calculated by jim.
%load so2295sm.dat;  % calculated by jim.
%fso2(:,1)=so2295sm(:,6)*10;
%fso2(:,2)=so2295sm(:,7);

% carga las secciones eficaces del o3 y cambia de unidades (cm^-1 y amstrong)
f=liesfile('ozxsec2.dat',1,7);
k=1.38062e-23; % bolzmann
o3x1T=f(:,4)*1.013*1e5/(k*273.1)*1e-6; %in %1/cm
o3x1wl=f(:,1)*10;


%%
% so2 que no me interesa
%aux=find(~isnan(fso2(:,2)));
%o3xf{2}=[fso2(aux,1),fso2(aux,2)];
%aux=find(~isnan(dbm_(:,1)));
%o3xf{2}=[dbm.lamda(aux),dbm_(aux,:)];
%aux=find(~isnan(bpi_(:,1)));

% la celda 1 es la del o3 que es longitud de onda x seccion eficaz
o3xf=[f(:,1)*10,o3x1T];
% o3x->(brewer,temp,slit,o3f)


o3x(size(dsp_summ,1),6,1000)=NaN;
o3a(size(dsp_summ,1))=NaN;
%so2(size(dsp_summ,1),length(T),2)=NaN;
for nb=1:size(dsp_summ,1)
    for ns=1:6
		for nr=1:1000
		%for nr=1:size(sigma_lambda1,2)
			x=o3xf(:,1);
			cs=o3xf(:,2);
			y=trapezoid_brewer(x,dsp_summ(nb,3+ns)+sigma_lambda(nr,ns),(dsp_summ(nb,9+ns)+sigma_fwhm(nr,ns))/2,.87);
			o3x(nb,ns,nr)=trapz(x,cs.*y)/trapz(x,y);
		end
    end
    o3a(nb,1:1000)=-O3W*10*squeeze(o3x(nb,:,1:1000))/log(10);
%        so2(nb,nt,nf)=-SO2W*squeeze(o3x(nb,nt,nf,:))/log(10);
end
disp('t');
%%
dsp_eff=[dsp_summ(:,1),squeeze(o3a)];
aaa = mean(dsp_eff(:,2:1001));
A_unc = std(aaa);
%%
