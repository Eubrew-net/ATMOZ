fpath='data_set_2/P121'
l=dir(fullfile(fpath,'Pandora121s1*.txt'))
t=readtable(fullfile(fpath,l.name),'HeaderLines',63);
fecha=datetime(t.Var1,'InputFormat','yyyyMMdd''T''HHmmss''Z''' );
t.Date=datenum(fecha);
t=t(month(fecha)==9,:);



%t.Properties.VariableNames={'Date','sza','airm','FW5_flag','FW5_o3','FW6_flag','FW6_o3',...
%                                                'FW5b_flag','FW5b_o3','FW6b_flag','FW6b_o3'};
% ---------------------------------------------------------------------------------------
% Column 1: UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)
% Column 2: Fractional days since 1-Jan-2000 UT midnight for center of measurement
% Column 3: Effective duration of measurement in seconds
% Column 4: Solar zenith angle for center of measurement in degree
% Column 5: Solar azimuth for center of measurement in degree, 0=north, increases clockwise
% Column 6: Lunar zenith angle for center of measurement in degree
% Column 7: Lunar azimuth for center of measurement in degree, 0=north, increases clockwise
% Column 8: Ozone total vertical column amount [Dobson Units], -9e99=retrieval not successful
% Column 9: Uncertainty of ozone total vertical column amount [Dobson Units] based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved
% Column 10: Ozone effective temperature [K], -8=retrieval not successful, 0=cross section is zero in this wavelength range or differential optical depth is too small to retrieve the temperature
% Column 11: Uncertainty of ozone effective temperature [K] based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -2=differential optical depth is too small to retrieve the temperature, -3=spectral fitting was done, but no uncertainty could be retrieved
% Column 12: Direct ozone air mass factor
% Column 13: L2 data quality flag for ozone: 0=high quality, 1=medium quality, 2=low quality
% Column 14: Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large
% Column 15: Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ2 limit (same parameters as for DQ1)
% 

t.Date_str=t.Var1;
t_dep=t(t.Var13==0,:);
t_dep=t_dep(t_dep.Var14==0,:);






t_dep.O3=t_dep.Var8;
t_dep.O3_STD=t_dep.Var9;
t_dep.AIRM=t_dep.Var12;
t_dep.Time=datetime(datestr(t_dep.Date));

t_set_2{n_inst}=t_dep(:,{'Time','Date','O3','O3_STD','AIRM'});
writetable(t_dep,'Atmoz_o3_set2.xls','Sheet','P121');




