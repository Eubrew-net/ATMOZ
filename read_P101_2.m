fpath='data_set_2/P101'
l=dir(fullfile(fpath,'P101*.txt'))
t=readtable(fullfile(fpath,l.name));
fecha=datetime(t.Var1,'InputFormat','yyyyMMdd''T''HHmmss''Z''' );
% #------------------------------------------------------------------------------------------
% #Column 1: Datetimes
% #Column 2: SZA
% #Column 3: AMF
% #Column 4: P121_CF121_20160415_v2_Daumont4TGOME_225K_FW5_filter_code
% #Column 5: P121_CF121_20160415_v2_Daumont4TGOME_225K_FW5_gas_vc
% #Column 6: P121_CF121_20160415_v2_Daumont4TGOME_225K_FW6_filter_code
% #Column 7: P121_CF121_20160415_v2_Daumont4TGOME_225K_FW6_gas_vc
% #Column 8: P121_CF121_20160415_v2_Harmonics2013_227K_FW5_filter_code
% #Column 9: P121_CF121_20160415_v2_Harmonics2013_227K_FW5_gas_vc
% #Column 10: P121_CF121_20160415_v2_Harmonics2013_227K_FW6_filter_code
% #Column 11: P121_CF121_20160415_v2_Harmonics2013_227K_FW6_gas_vc


t.Properties.VariableNames={'Date','sza','airm','FW5_flag','FW5_o3','FW6_flag','FW6_o3',...
                                                'FW5b_flag','FW5b_o3','FW6b_flag','FW6b_o3'};
t.Date_str=t.Date;
t.Date=datenum(fecha);
t_dep=t(t.FW5b_flag==0,:);
writetable(t_dep,'Atmoz_o3_set2.xls','Sheet','P101');


t_dep.O3=t_dep.FW6b_o3;
t_dep.O3_STD=t_dep.FW6b_flag;
t_dep.AIRM=t_dep.airm;
t_dep.Time=datetime(datestr(t_dep.Date));

t_set_2{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});




