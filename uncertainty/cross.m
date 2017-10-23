function A_unc = cross(brw)


% load the dispersion file.
% col1 is the time
% col2 is the Brewer number
% col3 is the fit (2,3)
% col4 to col9 are the wavelengths and from col10  col15 are the fwhm
%brw_str = num2str(brw)
a=strcat('dsp_',brw,'_1516.csv');
data=csvread(a,1,0);

%%
% load the variables from 1 to 15
dsp_summ=data(:,:);
%bp_r=data(:,[1,2,16:end]);

% sigma of wavelengths and fwhm
dsp_lambda = std(dsp_summ(:,4:9));
dsp_fwhm = std(dsp_summ(:,10:15));

% random vector
rnd = randn(1,1000);

% random sigmas
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

% weights of o3 and so2
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];

% Brewer's wavelenths
BWN=[302.1,306.3,310.1,313.5,316.8,320.1];

%%
% load the so2 cross sections (not used)
% load so2;
% so2295sm.dat;  % calculated by jim.
%load so2295sm.dat;  % calculated by jim.
%fso2(:,1)=so2295sm(:,6)*10;
%fso2(:,2)=so2295sm(:,7);

% load the o3 cross sections and change the units (cm^-1 and amstrong)
f=liesfile('ozxsec2.dat',1,7);
k=1.38062e-23; % bolzmann
o3x1T=f(:,4)*1.013*1e5/(k*273.1)*1e-6; %in %1/cm
o3x1wl=f(:,1)*10;


%%

% the cell 1 is o3 wavelength x cross section
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
    o3a(nb,1:1000)=-O3W*squeeze(o3x(nb,:,1:1000))/log(10);
%        so2(nb,nt,nf)=-SO2W*squeeze(o3x(nb,nt,nf,:))/log(10);
end
disp('t');
%%
dsp_eff=[dsp_summ(:,1),squeeze(o3a)];
aaa = mean(dsp_eff(:,2:1001));
A_unc = std(aaa);
%%
