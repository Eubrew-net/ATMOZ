fpath='data_set_1/P101'
l=dir(fullfile(fpath,'P*.dat'))
%l=dir(fullfile(fpath,'P101_ATMOZ_Comparison_FWs_*.txt'))
t=readtable(fullfile(fpath,l.name));
fecha=datetime(t.Var1,'InputFormat','yyyyMMdd''T''HHmmss''Z''' );
t.Properties.VariableNames={'Date','mdt','sza','saz','rms','nrms','o3','o3_u','airm','o3_flag','t_int','stray_lev','wv_shift','FW2','fit_flag','o3_eff_t','o3_eff_t_u'};
t.Date_str=t.Date;
t.Date=datenum(fecha);

%% Pandonia Depuration

% #P101 filtering for PanPS format data (P101_IZO_O3_FW5-ATMOZ.txt file), data that agree with these conditions should be filtered:
% # DQP1: Column 13 -> abs(Wavelength shift) > 0.2 [nm]
% # DQP2: Column 3  -> Stray Light, Solar zenith angle > 79.0
% # DQP3: Column 9  -> Stray Light, Air mass > 5.0
% # DQP4: Column 6  -> Scatter, wMRS:normalized root mean squares of the weighted spectral fitting residuals (wRMS) > 0.015
% # DQP5: Column 17 -> Cloud sreening: Uncertainty in the vertical column amount > 1 [DU]
% # DQP6: Column 15 -> Fitting result index: 1,2 = ok,  2=error

t_dep=t(abs(t.wv_shift)<0.2 &  t.sza<79.0 & t.airm<5 & t.nrms<0.015 &  t.o3_u<1 & t.fit_flag<=2 ,:);


%writetable(t_dep,'Atmoz_o3_set1.xls','Sheet','P101');

%% flag==0
t_dep.O3=t_dep.o3;
t_dep.O3_STD=t_dep.o3_u;
t_dep.AIRM=t_dep.airm;
t_dep.Time=datetime(datestr(t_dep.Date));
t_dep=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});

t_set_1{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});