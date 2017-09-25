% options_pub.outputDir=fullfile(pwd,'latex','017','html'); options_pub.showCode=true;
% close all; publish(fullfile(pwd,'cal_report_017c2.m'),options_pub);
%% Brewer Evaluation
%%
clear all;
file_setup='atmoz2016_setup';

eval(file_setup);     % configuracion por defecto
Cal.n_inst=find(Cal.brw==017);
Cal.file_latex=fullfile('.','latex',Cal.brw_str{Cal.n_inst});
Cal.dir_figs=fullfile('latex',filesep(),Cal.brw_str{Cal.n_inst},...
                              filesep(),[Cal.brw_str{Cal.n_inst},'_figures'],filesep());
mkdir(Cal.dir_figs);
%Cal.file_save='atmoz.mat'
try
 save(Cal.file_save,'-Append','Cal'); %sobreescribimos la configuracion guardada.
 load(Cal.file_save);
catch
    disp('clean');
    save(Cal.file_save);
end
%% Brewer setup
%%
%Cal.Date.CALC_DAYS=[130:157,165:170]
%Cal.Date.CALC_DAYS=[240:280]
Cal.n_ref=[1:3];%find(Cal.brw==185);

% READ Brewer Summaries
 for i=[Cal.n_inst,Cal.n_ref]
    dsum{i}={};       ozone_raw{i}={};   hg{i}={};
    ozone_sum{i}={};  ozone_raw0{i}={};  bhg{i}={};
    config{i}={};     sl{i}={};          log{i}={};
    ozone_ds{i}={};   sl_cr{i}={};       missing{i}=[];

    [ozone,log_,missing_]=read_bdata(i,Cal);

    dsum{i}=ozone.dsum;
    ozone_sum{i}=ozone.ozone_sum;
    config{i}=ozone.config;
    ozone_ds{i}=ozone.ozone_ds;
    ozone_raw{i}=ozone.raw;
    ozone_raw0{i}=ozone.raw0;
    sl{i}=ozone.sl; %first calibration/ bfiles
    sl_cr{i}=ozone.sl_cr; %recalculated with 2? configuration
    hg{i}=ozone.hg;
    bhg{i}=ozone.bhg;
    log{i}=cat(1,log_{:});
    missing{i}=missing_';
    disp(log{i});
 end

save(Cal.file_save,'-APPEND','ozone_raw','dsum','ozone_sum','ozone_ds','config','sl','sl_cr','log','missing');
matrix2latex(log{Cal.n_inst},fullfile(Cal.file_latex,['tabla_fileprocess_',Cal.brw_str{Cal.n_inst},'.tex']));
%% configuration files
% Changes detected on B files configuration:
%%
changes=write_excel_config(config,Cal,Cal.n_inst);
%% 
% Configuration :

[config_ref,TCref,DTref,ETCref,A1ref,ATref,leg]  =read_icf(Cal.brw_config_files{Cal.n_ref(1),2});
[config_def,TCdef,DTdef,ETCdef,A1def,ATdef]      =read_icf(Cal.brw_config_files{Cal.n_inst,2});
[config_orig,TCorig,DTorig,ETCorig,A1orig,ATorig]=read_icf(Cal.brw_config_files{Cal.n_inst,1});
% latex output
latexcmd(fullfile(Cal.file_latex,['cal_config_',Cal.brw_str{Cal.n_inst}]),...
                     DTref,'\ETCref',ETCref(1),'\Aref',A1ref(1),...
                     DTdef,'\ETCdef',ETCdef(1),'\Adef',A1def(1),...
                     DTorig,'\ETCorig',ETCorig(1),'\Aorig',A1orig(1));
str_leg=cellstr(leg);

[a b c]=fileparts(Cal.brw_config_files{Cal.n_inst,1}); confini=strcat(b,c);
[a b c]=fileparts(Cal.brw_config_files{Cal.n_inst,2}); confend=strcat(b,c);
columnlabels={sprintf('%s %c%s%c','Initial','(',confini,')'),...
              sprintf('%s %c%s%c','Final','(',confend,')')};
matrix2latex_config([config_orig(2:52),config_def(2:52)],...
                     fullfile(Cal.file_latex,['table_config_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                              'rowlabels',str_leg(2:end),'columnlabels',columnlabels,...
                                              'size','footnotesize');

%makeHtmlTable([config_orig,config_def],[],cellstr(leg),[Cal.brw_config_files(Cal.n_inst,1),Cal.brw_config_files(Cal.n_inst,2)])  
makeHtml_Table([config_orig,config_def],[],cellstr(leg),[Cal.brw_config_files(Cal.n_inst,1),Cal.brw_config_files(Cal.n_inst,2)])  

%% 
% Configuration Table
%%
[p,fref]=fileparts(Cal.brw_config_files{Cal.n_ref(1),2});  
[p,f1]=fileparts(Cal.brw_config_files{Cal.n_inst(1),1});  
%[p,f2]=fileparts(Cal.brw_config_files{Cal.n_inst(1),2});  
t_config=dispTable([config_orig(1:52),config_def(1:52),config_ref],[f1,columnlabels(2),fref],cellstr(leg));
t_config(1:21,:)
%% DT analysis
%%
 DT_analysis(Cal,ozone_raw0,config,'plot_flag',0);

 figure(findobj('tag','DT_comp')); set(get(gca,'title'),'FontSize',8); 
 printfiles_report(gcf,Cal.dir_figs,'LockAxes',0,'Height',6);
 
%figure(findobj('tag','raw_counts_187')); title('Day 18711'); set(findobj(gcf,'Tag','legend'),'FontSize',8); printfiles_report(gcf,Cal.dir_figs,'Width',8,'Height',6,'LockAxes',0);
%figure(findobj('tag','raw_counts_191')); title('Day 19111'); ylabel(''); set(findobj(gcf,'Tag','legend'),'FontSize',8); printfiles_report(gcf,Cal.dir_figs,'Width',8,'Height',6,'LockAxes',0);
%close all
%% 
% This shows the ratio between two configurations analyzed.
%% SL Report
% Standad Lamp analysis.
% 
% New TC introduced, old TC are not bad shows  4 ETC/10º but can be improved.
%%
close all;
for ii=[Cal.n_ref,Cal.n_inst]
    sl_mov_o{ii}={}; sl_median_o{ii}={}; sl_out_o{ii}={}; R6_o{ii}={};
    sl_mov_n{ii}={}; sl_median_n{ii}={}; sl_out_n{ii}={}; R6_n{ii}={};
    try
    if ii==Cal.n_inst
% old instrumental constants
      [sl_mov_o{ii},sl_median_o{ii},sl_out_o{ii},R6_o{ii}]=sl_report_jday(ii,sl,Cal.brw_str,...
                               'date_range',datenum(Cal.Date.cal_year,1,Cal.calibration_days{Cal.n_inst,1}([1 end])),...
                               'outlier_flag',1,'hgflag',1,'fplot',1);
      fprintf('%s Old constants: %5.0f +/-%5.1f\n',Cal.brw_name{ii},...
                  nanmedian(R6_o{ii}(diajul(R6_o{ii}(:,1))>=Cal.calibration_days{ii,3}(1),2)),...
                  nanstd(   R6_o{ii}(diajul(R6_o{ii}(:,1))>=Cal.calibration_days{ii,3}(1),2)));
% new instrumental constants
      [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
                               'date_range',datenum(Cal.Date.cal_year,1,Cal.calibration_days{Cal.n_inst,1}([1 end])),...
                               'outlier_flag',1,'hgflag',1,'fplot',1);
      fprintf('%s New constants: %5.0f +/-%5.1f\n',Cal.brw_name{ii},...
                  nanmedian(R6_n{ii}(diajul(R6_n{ii}(:,1))>=Cal.calibration_days{ii,3}(1),2)),...
                  nanstd(   R6_n{ii}(diajul(R6_n{ii}(:,1))>=Cal.calibration_days{ii,3}(1),2)));
    else
      [sl_mov_n{ii},sl_median_n{ii},sl_out_n{ii},R6_n{ii}]=sl_report_jday(ii,sl_cr,Cal.brw_str,...
                               'date_range',datenum(Cal.Date.cal_year,1,Cal.calibration_days{Cal.n_inst,1}([1 end])),...
                               'outlier_flag',0,'fplot',0);
    end
    catch exception
          fprintf('%s, brewer: %s\n',exception.message,Cal.brw_str{ii});
    end
end
% anadimos el sl
save(Cal.file_save,'-APPEND','sl_mov_o','sl_median_o','sl_out_o','R6_o',...
                             'sl_mov_n','sl_median_n','sl_out_n','R6_n');

sl_report{Cal.n_inst}.sl_mov_n=sl_mov_n{Cal.n_inst};         sl_report{Cal.n_inst}.sl_mov_o=sl_mov_o{Cal.n_inst};
sl_report{Cal.n_inst}.sl_median_n=sl_median_n{Cal.n_inst};   sl_report{Cal.n_inst}.sl_median_o=sl_median_o{Cal.n_inst};
sl_report{Cal.n_inst}.sl_out_n=sl_out_n{Cal.n_inst};         sl_report{Cal.n_inst}.sl_out_o=sl_out_o{Cal.n_inst};
sl_report{Cal.n_inst}.R6_n{Cal.n_inst}=R6_n{Cal.n_inst};     sl_report{Cal.n_inst}.R6_o{Cal.n_inst}=R6_o{Cal.n_inst};
save(Cal.file_save,'-APPEND','sl_report');
%% SL report plot
%
ix=sort(findobj('-regexp','Tag','SL_R\w+'));
%if length(ix)>2 for i=ix, set( findobj(i,'Tag','legend'),'FontSize',8);end 
% set(ix(1),'tag','SL_R6_report_old');
% set(ix(2),'tag','SL_R5_report_old');
% Width=8; Height=6; 
%else
%    Width=14; Height=6.5; 
%end
%printfiles_report(ix',Cal.dir_figs,'Width',Width,'Height',Height);
%figure(double(findobj('tag','SL_I5_report'))); printfiles_report(gcf,Cal.dir_figs,'LockAxes',0,'no_export');
%close all
%% Hg Report
% % Hg failure

 hg_fail=read_custom(fullfile(Cal.path_root,sprintf('bdata%s',Cal.brw_str{Cal.n_inst})),sprintf('B*.%s',Cal.brw_str{Cal.n_inst}),...
     'hg: lamp not detected','printfile',0,...
     'date_range',datenum(Cal.Date.cal_year,1,1)); 
 if ~isempty(hg_fail) [mean_hg_fail N_hg_fail] = grpstats(fix(hg_fail),fix(hg_fail),{'nanmean','numel'});
  displaytable( N_hg_fail,{'N(Hg failure)'},13,'d',cellstr(datestr(mean_hg_fail(:,1),1)))
 else 
     fprintf('Ho Hg failures recorded for Brewer %s\r\n',Cal.brw_name{Cal.n_inst});
 end % Hg during campaign hg_report(Cal.n_inst,hg,Cal,'outlier_flag',1,'date_range',[]);
%% READ Configuration

close all
[A,ETC,SL_B,SL_R,F_corr,cfg]=read_cal_config_new(config,Cal,{sl_median_o,sl_median_n});
% F_corr is read from setup file (or matrix)
try
r=[ETC.old,ETC.new,A.old,A.new,fix(SL_B.old),fix(SL_B.new),SL_R.old,SL_R.new];
 cfgt=reshape(r,[],Cal.n_brw+1,8);
 tabla_sl=[squeeze(cfgt(:,1,1)),squeeze(cfgt(:,Cal.n_inst+1,:)),F_corr{Cal.n_inst}.old(:,5:end),F_corr{Cal.n_inst}.new(:,5:end)];
 Var=varname({'Date','ETC_1','ETC_2','A1_1','A1_2','SL_B1','SL_B2','SL_R1','SL_R2','Fc1#3','Fc1#4','Fc1#5','Fc2#3','Fc2#4','Fc2#5'});
 tabla_sl=array2table(tabla_sl,'VariableNames',Var);
 tabla_sl.Date=datetime(datestr(squeeze(cfgt(:,1,1))))
 %cell2csv(fullfile(Cal.file_latex,'tabla_config.csv'),tabla_sl',';');
 writetable(tabla_sl,'cfg_table.xls','Filetype','spreadsheet','Sheet',Cal.brw_name{Cal.n_inst})
catch
  l=lasterror;
  disp(l.message)
  disp('tabla de configuracion posiblemente incompleta')
end
%% Data recalculation for summaries and individual observations
%
for i=[Cal.n_inst Cal.n_ref]
    cal{i}={}; summary{i}={}; summary_old{i}={};
    if i==Cal.n_inst
    [cal{i},summary{i},summary_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,'flag_sl',1);
    else
    [cal{i},summary{i},summary_old{i}]=test_recalculation(Cal,i,ozone_ds,A,SL_R,SL_B,'flag_sl',1);
    end
end
% outiler
%id_out=find(summary{Cal.n_inst}(:,6)>400 | summary{Cal.n_inst}(:,6)<200);
%disp(unique(diaj(summary{Cal.n_inst}(id_out,1))));
%summary{Cal.n_inst}(id_out,:)=[];

summary_orig=summary; summary_orig_old=summary_old;
save(Cal.file_save,'-APPEND','summary_old','summary_orig_old','summary','summary_orig');
%% Filter distribution
%%
close all
ozone_filter_analysis_mi(summary,Cal);
%%
printfiles_report(findobj('Tag','FILTER_DISTRIBUTION'),Cal.dir_figs,'Width',13,'Height',7.5);
printfiles_report(findobj('Tag','Ozone_diff_filter_rel'),Cal.dir_figs,'Width',13,'Height',7.5);

close all
%% filter correction
% Si queremos eliminar algun filtro CORREGIR a NaN
%
for ii=[Cal.n_ref Cal.n_inst]
   [summary_old_corr summary_corr]=filter_corr(summary_orig,summary_orig_old,ii,A,F_corr{ii});
   summary_old{ii}=summary_old_corr; summary{ii}=summary_corr;
end
 figure; plot(summary{Cal.n_inst}(:,1),summary{Cal.n_inst}(:,6),'r.',...
              summary{Cal.n_ref(1)}(:,1),summary{Cal.n_ref(1)}(:,6),'k.');%,...
%             summary{Cal.n_ref(2)}(:,1),summary{Cal.n_ref(2)}(:,6),'m.');
% legend(gca,'inst','IZO#183','IZO#185','Location','NorthEast'); grid;
% datetick('x',26,'KeepLimits','KeepTicks');
%ozone_filter_analysis_mi(summary,Cal);

%% ATMOZ output
t_atmoz=write_summary_cfg(Cal.n_inst,2016,summary_old,summary,SL_B,SL_R,A,ETC);
t_017=table2timetable(t_atmoz{Cal.n_inst,2},'RowTimes',date_time(t_atmoz{Cal.n_inst,2}.Date));
writetable(timetable2table(t_017),'B017_atmoz.xls')

% test
%o3x=(t_017.ms9_corr-t_017.ETC)./(t_017.A1.*t_017.m2*10);
%plot(o3x-t_017.O3_1)
meteo=table2timetable(readtable('atmoz_mteo_eff.xls'));
ts=synchronize(t_017,meteo,'first','nearest','EndValues',NaN);
m_heff= (6370./(6370+ts.Heff_ecc)).* (ts.m2/0.99656);
m3=ts.m2*1.0027;
t_017.m_heff=m_heff;
BE=[4870,4620,4410,4220,4040];
O3W=[ 0 1.00    -0.50    -2.20   +1.70];
Bates=BE*O3W'
% brewer software Fi=Fio + Bi 
% Bi= B(i)p/po*m3 
% the ratios are calculated with -O3W (indirectly)
% but the effective paramters  use O3W-> B=O3W*Bi
% MS9=MS90 + Beta 
% -O3W*Fi = -O3W*Fio + -O3W*Bi p/po*m3
%  MS9 =MS90 - B (p/po*m3)
%  MS90= MS9 + B (p/po*m3)
sum_wiBi=Bates; % orii
MS9c= sum_wiBi*ts.pressure_hPa_.*ts.m2*1.0027/1013.24;
MS90=ts.ms9_corr+MS9c;

%MS9= MS90 - sum(wiBi) 770/1013.25 * m2*1.0027
sum_wiB=10; 
MS9=MS90-sum_wiB*ts.pressure_hPa_.*ts.m2*1.0027/1013.24;
t_017.MS9_2=MS9;
t_017.p=ts.pressure_hPa_;
t_017.T_eff=ts.Teff_ecc;
%o3x=(ts.ms9_corr-ts.ETC)./(ts.A1.*ts.m2*10);
o3c=(MS9-ts.ETC)./(ts.A1.*ts.m2*10);
figure
 histogram((o3c-ts.O3_1));
 xlabel('DU');
 title('Ozone difference O3 (Nicolet)- O3 Bates(Fixed)')
Ci=[-4.9188E-8,2.8781E-5 3.4591E-1, 3.4452E-1];
%adjust to brewer value
%A1=A.new(:,Cal.n_inst+1);
%A1_new=unique(A1(~isnan(A1)))
A1_new=  0.3450;  
To=polyval(Ci(1:end-1),-45)-A1_new;
Ci(3)=Ci(3)-To;
Aop=polyval(Ci(1:end-1),-45);
At=polyval(Ci(1:end-1),ts.Teff_ecc-273.15);
o3temp=(ts.ms9_corr-ts.ETC)./(At.*ts.m2*10);
o3to=(ts.ms9_corr-ts.ETC)./(Aop.*ts.m2*10);
t_017.AT=At;
t_017.ATo=Aop*ones(size(At));
figure
histogram((o3temp-ts.O3_1))
title('Ozone difference O3 (B&P)- O3 (Bremen))')
disp('DU')


%all 
o3t=(MS9-ts.ETC)./(At.*m_heff*10);
figure
histogram(100*(o3t-ts.O3_1)./o3t)
title('% Ozone difference O3 BREWER - O3 (ATMOZ))')
xlabel('%')
figure
plot(ts.Time,100*(o3t-ts.O3_1)./ts.O3_1)
t_017.O3_2=o3t;

writetable(timetable2table(t_017),'B017_atmoz.xls')






%% Langley Calculation
%%
summ_old{Cal.n_inst}=summary_orig{Cal.n_inst};
summ{Cal.n_inst}=summary{Cal.n_inst};
airm_rang=[1.15 3.75];
Fcorr={[],[0,0,0,0,0,0],[],[0,0,0,0,0,0]}; % F_corr Usaremos F's de la configuracion
N_data=12;
O3_std=2.5;

AOD_file=fullfile('Triad','aod','150101_161231_Izana.lev15');
CLOUD_file=fullfile('..','cloudScreening.txt');

%% ---- langley from Indiv. Measurements ----
%
[ozone_lgl{Cal.n_inst},cfg_indv,leg,ozone_lgl_sum{Cal.n_inst}] = langley_data_cell(ozone_raw{Cal.n_inst},ozone_ds{Cal.n_inst},config{Cal.n_inst});
save(Cal.file_save,'-APPEND','ozone_lgl','cfg_indv','leg','ozone_lgl_sum');

  % Filters regression (aplicamos filtros de AOD y O3 stricto)
   cfgs=1;      % Operativa
   lgl_filt{Cal.n_inst}=langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
                                'airmass',airm_rang,'O3_hday',2,'N_hday',N_data);%,...
                                %'AOD',AOD_file);%,'date_range',datenum(2014,1,[66 80]));

   brw_indv_{Cal.n_inst} = langley_analys_filter(lgl_filt,Cal.n_inst,...
                                         'res_filt',1,'plot_flag',0);
% Filter
%
   nd0_=cat(1,brw_indv_{Cal.n_inst}(:,[1 2],cfgs),brw_indv_{Cal.n_inst}(:,[1 6],cfgs)); nd0=sortrows(nd0_,1);
   nd3_=cat(1,brw_indv_{Cal.n_inst}(:,[1 3],cfgs),brw_indv_{Cal.n_inst}(:,[1 7],cfgs)); nd3=sortrows(nd3_,1);
   nd4_=cat(1,brw_indv_{Cal.n_inst}(:,[1 4],cfgs),brw_indv_{Cal.n_inst}(:,[1 8],cfgs)); nd4=sortrows(nd4_,1);

   % Hist
   figure; set(gcf,'Tag',sprintf('Lag%s_residuals',Cal.brw_str{Cal.n_inst}));
   aux={nd0(:,2),nd3(:,2),nd4(:,2)};
   nhist(aux,'box','smooth','samebins','ylabel','Pdf (Langley fit residuals)'); grid; box on;
   set(findobj(gca,'Type','Line'),'Marker','None','LineWidth',2)
   
   legendflex({'ND#0,1,2','ND#3','ND#4'},'anchor', [2 6], 'buffer',[0 -20], ...
                'nrow',1,'fontsize',10,'box','on','xscale',.7,...
                'title',sprintf('Brewer %s: ND#3 Corr.=%d, ND#4 Corr.=%d',Cal.brw_name{Cal.n_inst},...
                        round(nanmedian(nd3(:,2))-nanmedian(nd0(:,2))),...
                        round(nanmedian(nd4(:,2))-nanmedian(nd0(:,2)))));
%% Final Langley  With Filter correction applied
[ozone_lgl_dep{Cal.n_inst} s_day]=langley_filter_lvl1(ozone_lgl{Cal.n_inst},'plots',0,...
                       'F_corr',Fcorr{Cal.n_inst},'airmass',airm_rang,'O3_hday',O3_std);%,...
                       % 'AOD',AOD_file,'lgl_days',1,'plots',1);%,...
%                   
[brw_indv{Cal.n_inst} dbs_indv{Cal.n_inst} st_brw{Cal.n_inst} st_dbs{Cal.n_inst}] = langley_analys(ozone_lgl_dep,Cal.n_inst,Cal,...
                       'res_filt',1,'plot_flag',0);

write_langley(Cal.brw(Cal.n_inst),Cal.Date.cal_year,brw_indv(Cal.n_inst));                   
                   
%

cfgs=1; [a b]=fileparts(Cal.brw_config_files{Cal.n_inst,cfgs});

 figure; ha=tight_subplot(2,1,.08,.1,.075); hold all;
 axes(ha(1)); set(gca,'XTicklabel',[],'box','on','YTickLabelMode','auto'); grid; hold on;
 axes(ha(2)); set(gca,'box','on','YTickLabelMode','auto'); grid; hold on;

 plot(ha(1),brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,2,cfgs),'-bd',brw_indv{Cal.n_inst}(:,1,cfgs),brw_indv{Cal.n_inst}(:,3,cfgs),':rd');
 plot(ha(2),dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,2,cfgs),'-bd',dbs_indv{Cal.n_inst}(:,1,cfgs),dbs_indv{Cal.n_inst}(:,3,cfgs),':rd');
 set(findobj(ha,'Type','Line'),'MarkerSize',4);
 legend({'AM','PM'},'Location','North','Orientation','Horizontal','FontSize',7);
%
figure
boxplot(dbs_indv{Cal.n_inst}(:,2:3,1))
%% 
% 
%%
med1=round(trimmean(brw_indv{Cal.n_inst}(:,2,1),20))
med2=round(trimmean(brw_indv{Cal.n_inst}(:,2,2),20))
me1=round(trimmean(reshape(dbs_indv{Cal.n_inst}(:,2:3,1),1,[]),20))
me2=round(trimmean(reshape(dbs_indv{Cal.n_inst}(:,2:3,2),1,[]),20))
%%
figure
boxplot(reshape(dbs_indv{Cal.n_inst}(:,2:3,1:2),[],2),'Labels',{'Op_config','Alt_config'})
grid

hline([me1,me2]);
title({sprintf('ETC correction =%d / %d',[me1,me2]),sprintf('ETC=%d / %d',[med1,med2])})
%%
med1=round(trimmean(brw_indv{Cal.n_inst}(:,2,1),20))
med2=round(trimmean(brw_indv{Cal.n_inst}(:,2,2),20))
figure
scatterhist(brw_indv{Cal.n_inst}(:,1,2),brw_indv{Cal.n_inst}(:,2,1));
hold on
plot(brw_indv{Cal.n_inst}(:,1,2),brw_indv{Cal.n_inst}(:,2,2),'x');
legend('initial config','new config')
datetick;

hline([med1,med2]);
title(sprintf('ETC=%d / %d',[med1,med2]))
%% Ozone Calibration
% Reference Brewer #185
%
 close all
 Cal.n_ref=1:3;
 n_ref=Cal.n_ref; % Cuidado con cual referencia estamos asignando
 n_inst=Cal.n_inst
 brw=Cal.brw; brw_str=Cal.brw_str;
 % Artificial reference
 summary{3}(summary{3}(:,6)<200,:)=[];
  ref=[]; ref_o=[];
  for n_ref=Cal.n_ref     
    jday_ref=findm((diaj(summary{n_ref}(:,1))),Cal.calibration_days{Cal.n_inst,1},0.5);
    ref_=summary{n_ref}(jday_ref,:);
    ref_old=summary_old{n_ref}(jday_ref,:);
    ref=[ref;ref_]; 
    ref_o=[ref_o;ref_old]; 
  end
%ref=summary{3};
 summary{5}=sortrows(ref,1);
 summary_old{5}=sortrows(ref_o,1);
% summary{4}=summary_old{4}; % still not calibrated
Cal.brw_name=[Cal.brw_name, 'REF'];
%  old config
 summ=summary_old;
 summary_old{1}=summary{1};
 reference_brw=5; 
 analyzed_brewer=[1 2 3 4 5];
 osc_interval=[400,700,1000,1200];
 Cal.analyzed_brewer=analyzed_brewer; 
 Cal.sl_c=[0,0,0,0,0];
% 
 Cal.brw=[157 183 185 017 100 ];
 Cal.brw_str=[Cal.brw_str; 'REF'];
 brw_str=Cal.brw_str;
 Cal.n_ref=5;
%% osc table
%%
 [reft,ratio_ref]=join_summary(Cal,summary_old,reference_brw,analyzed_brewer,5);
 [ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
  set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
  [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
  set(gca,'Ylim',[-5,2])
%%
   [reft,ratio_ref]=join_summary(Cal,summary,reference_brw,analyzed_brewer,5);
   [ratio_osc_table,osc_matrix,osc_stats]=osc_table(Cal,ratio_ref,osc_interval);
  set(gcf,'Tag','osc_box_plot'); set(findobj(gcf,'Tag','legend'),'Location','SouthWest'); grid
  [f_hist,f_ev,f_sc,f_smooth]=ratio_ref_plots(Cal,ratio_ref,'plot_smooth',1);
  set(gca,'Ylim',[-5,2])
%% 
% 
%% Original config
%%
ETC_SUG(1).NEW=NaN; o3c_SUG=NaN;  osc_smooth_sug=NaN*ones(1,7);
blinddays=Cal.calibration_days{Cal.n_inst,2};
%
n_ref=5;Cal.n_ref;
n_inst=Cal.n_inst;
ref=summary_old{n_ref};
blinddays=Cal.Date.CALC_DAYS;
jday=findm(diaj(summary_old{Cal.n_inst}(:,1)),blinddays,0.5);
inst1_b=summary_old{Cal.n_inst}(jday,:);

A1=A.old(ismember(Cal.Date.CALC_DAYS,blinddays),Cal.n_inst+1); 
A1_old=unique(A1(~isnan(A1)))
osc_range=0.6;
[ETC_SUG,o3c_SUG,m_etc_SUG]=ETC_calibration_C(Cal,summary_old,A1_old,...
                   Cal.n_inst,n_ref,5,osc_range,0.02,blinddays);

data_tabl=[nanmean(ETC.old(:,n_inst+1)),round([ETC_SUG(1).NEW,ETC_SUG(1).TP(1)]),...
           nanmean(A.old(:,n_inst+1)),round(ETC_SUG(1).TP(2)/10000,4),A1_old
%         solo el rango seleccionado
           nanmean(ETC.old(:,n_inst+1)),round([ETC_SUG(2).NEW,ETC_SUG(2).TP(1)]),...
           nanmean(A.old(:,n_inst+1)),round(ETC_SUG(2).TP(2)/10000,4),A1_old];
%         todo el rango
display_table(data_tabl,{'ETCorig','ETCnew 1p','ETCnew 2p','O3Abs (ICF)','O3Abs 2p','O3Abs sug.'},...
                     11,{'d','d','d','.4f','.4f','.4f'},{sprintf('OSC < %.2f DU',osc_range),'All OSC range'})
                 
%% 
%% STRAY OLD
%%
close all
A1_old=0.3416
osc_range=2;
[ETC_S,o3c_S,m_etc_S]=ETC_calibration_C(Cal,summary_old,A1_old,...
                   Cal.n_inst,n_ref,5,osc_range,0.02,blinddays);

close all;                                                    
straylight_{n_inst}=straylight_model(o3c_S,A1_old,ETC_S(1).NEW,Cal);
save(Cal.file_save,'-APPEND','straylight_');

%

a=straylight_{n_inst}.coeff.R(1);
b=straylight_{n_inst}.coeff.R(2);
etc=straylight_{n_inst}.coeff.R(3);
%etc=ETC_NEW(1).NEW
inst2=inst1_b;
o3=inst2(:,6);
m =inst2(:,3);
o3r= o3 -  a/A1_old*(m.*o3/1000).^b./m/10;
%o3r1= o3 -  a/A1_new*(m.*o3r).^b./m/10;

inst2(:,10)=o3r;
[x,r,rp,ra,dat,ox,osc_stray]=ratio_min_ozone(...
        inst2(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
        5,brw_str{n_inst},brw_str{n_inst},'plot_flag',0);


o3r1= o3 -  a/A1_old*(m.*o3r/1000).^b./m/10;

inst2(:,10)=o3r1;
[x,r,rp,ra,dat,ox,osc_stray2]=ratio_min_ozone(...
        inst2(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
        5,brw_str{n_inst},brw_str{n_inst},'plot_flag',0);
    
% inst2(:,10)=o3r1;
 [x,r,rp,ra,dat,ox,osc_orig]=ratio_min_ozone(...
         inst1_b(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
         5,brw_str{n_inst},brw_str{n_inst},'plot_flag',0);    
    
%   f=figure;
% set(f,'tag','RATIO_ERRORBAR_stray');
% h=plot_smooth_(osc_smooth_inisl,osc_smooth_fin,osc_stray2);
% legend(h,'Config. initial (SL corrected)','Final Configuration','Stray Light corrected ');
 set(gca,'YLim',[-4 2]);
% title([ brw_str{n_inst},' - ',brw_str{n_ref},' / ',brw_str{n_ref}]);
%% 
% inst1b original constant inst2 stray light corrected
%%
%inst2=summary{Cal.n_inst}(jday & jlim ,:);
inst1=inst1_b;
o3r= (inst1(:,8)-ETC_SUG(1).NEW)./(A1_old*inst2(:,3)*10);
inst1(:,10)=o3r;


    [x,r,rp,ra,dat,ox,osc_smooth_ini]=ratio_min_ozone(...
       inst1_b(:,[1,6,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
       5,brw_str{n_inst},Cal.brw_str{n_ref},'plot_flag',0);% original config
   
   
    [x,r,rp,ra,dat,ox,osc_smooth_fin]=ratio_min_ozone(...
       inst1(:,[1,12,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
       5,Cal.brw_str{n_ref},Cal.brw_str{n_ref},'plot_flag',0);
    
       
    
 f=figure;
 set(f,'tag','RATIO_ERRORBAR_stray');
 h=plot_smooth_(osc_smooth_ini,osc_smooth_fin,osc_stray2);
 legend(h,'Config. initial (SL corrected)','Final Configuration','Stray Light corrected ');
 set(gca,'YLim',[-4 2]);
 title([ brw_str{n_inst},' - ',brw_str{n_ref},' / ',brw_str{n_ref}]);
%% Final Period
%%
close all
finaldays=Cal.Date.CALC_DAYS;%=130:170



A1=A.new(ismember(Cal.Date.CALC_DAYS,finaldays),Cal.n_inst+1);
A1_new=unique(A1(~isnan(A1)))
osc_range=[0.9];
[ETC_NEW,o3c_NEW,m_etc_NEW]=ETC_calibration_C(Cal,summary,A1_new,n_inst,n_ref,...
                                                                5,osc_range,0.01,finaldays);
hidx=ismember(Cal.Date.CALC_DAYS,finaldays);
cc=nanmean([ETC.old(hidx,Cal.n_inst+1),10000*A.new(hidx,Cal.n_inst+1),10000*A.old(hidx,Cal.n_inst+1)]);

tableform({'ETCorig','ETCnew 1p','ETCnew 2p','O3Abs (ICF)','O3Abs 2p','O3Abs dsp'},...
          [round([cc(1),ETC_NEW(1).NEW,ETC_NEW(1).TP(1), cc(3),ETC_NEW(1).TP(2),cc(2)])
%         solo el rango seleccionado
           round([cc(1),ETC_NEW(2).NEW,ETC_NEW(2).TP(1), cc(2),ETC_NEW(2).TP(2),cc(3)])]);
%         todo el rango
%latexcmd(fullfile(['>',Cal.file_latex],['cal_etc_',brw_str{n_inst}]),'\ETCfin',num2str(round(ETC_NEW(1).NEW)),...
%    '\ETCSOfin',config_def(12),'\SOABSfin',config_def(10),'\OtresSOABSfin',config_def(9),...
%    '\NobsCalfin',size(o3c_NEW,1));
%etc{Cal.n_inst}.new=ETC_NEW;
save(Cal.file_save,'-APPEND','etc');
%%
jday=ismember(diaj(summary{Cal.n_inst}(:,1)),fix(finaldays));
jlim=(diaj2(summary{Cal.n_inst}(:,1))>finaldays(1) & ...    % 2st set the limits
      diaj2(summary{Cal.n_inst}(:,1))<finaldays(end)); 
inst2=summary{Cal.n_inst}(jday & jlim ,:);
inst1=summary_orig{Cal.n_inst}(jday & jlim ,:);

  o3r= (inst2(:,8)-ETC_NEW(1).NEW)./(A1_new*inst2(:,3)*10);
% o3r= (inst2(:,8)-ETC_NEW(1).NEW)./(A1_new*inst2(:,3)*10);
inst2(:,10)=o3r;

      [x,r,rp,ra,dat,ox,osc_smooth_ini]=ratio_min_ozone(...
         inst1(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
         5,brw_str{n_inst},Cal.brw_str{n_ref},'plot_flag',0);

    [x,r,rp,ra,dat,ox,osc_smooth_fin]=ratio_min_ozone(...
       inst2(:,[1,6,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
       5,brw_str{n_inst},Cal.brw_str{n_ref},'plot_flag',0);

figure(maxf(findobj('tag','RATIO_FILTER_INST')));
set(findobj(gcf,'Marker','o'),'MarkerSize',4); ax=findobj(gcf,'Tag','legend');
set(ax,'FontSize',6,'Orientation','Vertical'); 
printfiles_report(gcf,Cal.dir_figs,'aux_pattern',{'corr'});
   
figure; set(gcf,'tag','RATIO_ERRORBAR_all');
h=plot_smooth(osc_smooth_ini,osc_smooth_fin);
legend(h,'Config. initial','Final configuration'); set(gca,'YLim',[-4 2]);
title([ brw_str{n_inst},' - ',brw_str{n_ref},' / ',brw_str{n_ref}]);
%% 2P calibration
%%
   o3r= (inst2(:,8)-ETC_NEW(1).TP(1))./(ETC_NEW(1).TP(2)/10000*inst2(:,3)*10);
  inst2(:,10)=o3r;
       [x,r,rp,ra,dat,ox,osc_smooth_2P]=ratio_min_ozone(...
         inst2(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
          5,brw_str{n_inst},Cal.brw_str{n_ref},'plot_flag',0);
 
osc_smooth{Cal.n_inst}.fin=osc_smooth_fin;
osc_smooth{Cal.n_inst}.twoP=osc_smooth_2P;
 save(Cal.file_save,'-APPEND','osc_smooth');
%%
figure(maxf(findobj('tag','CAL_2P_SCHIST')));
ax=findobj(gca,'type','text');
set(ax(2),'FontSize',8,'Backgroundcolor','w'); set(ax(3),'FontSize',8,'Backgroundcolor','w');
printfiles_report(gcf,Cal.dir_figs,'aux_pattern',{'fin'});

figure(maxf(findobj('tag','RATIO_ERRORBAR')));
printfiles_report(gcf,Cal.dir_figs,'aux_pattern',{'fin'});

figure(findobj('tag','RATIO_ERRORBAR_all'));
set(gca,'YLim',[-2 7]); printfiles_report(gcf,Cal.dir_figs);

close all
%% 
% 
%% STRAY
%%
close all

osc_range=2;
[ETC_NEWo,o3c_NEWo,m_etc_NEWo]=ETC_calibration_C(Cal,summary_old,A1_new,n_inst,n_ref,...
                                                                10,osc_range,0.05,finaldays);

close all;                                                    
straylight_{n_inst}=straylight_model(o3c_NEWo,A1_new,ETC_NEW(1).NEW,Cal);
save(Cal.file_save,'-APPEND','straylight_');
%
a=straylight_{n_inst}.coeff.R(1);
b=straylight_{n_inst}.coeff.R(2);
etc=straylight_{n_inst}.coeff.R(3);
%etc=ETC_NEW(1).NEW
o3=inst2(:,6);
m =inst2(:,3);
%o3r= o3 -  a/0.3341*(m.*o3/1000).^b./m/10;
o3r1= o3 -  a/A1_new*(m.*o3r).^b./m/10;

inst2(:,10)=o3r;
[x,r,rp,ra,dat,ox,osc_stray]=ratio_min_ozone(...
        inst2(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
        5,brw_str{n_inst},brw_str{n_inst},'plot_flag',0);


o3r1= o3 -  a/A1_new*(m.*o3r/1000).^b./m/10;

inst2(:,10)=o3r1;
[x,r,rp,ra,dat,ox,osc_stray2]=ratio_min_ozone(...
        inst2(:,[1,10,3,2,8,9,4,5]),ref(:,[1,6,3,2,8,9,4,5]),...
        5,brw_str{n_inst},brw_str{n_inst},'plot_flag',0);
    
    
  f=figure;
set(f,'tag','RATIO_ERRORBAR_stray');
h=plot_smooth_(osc_smooth_fin,osc_smooth_fin,osc_stray2);
legend(h,'Config. initial (SL corrected)','Final Configuration','Stray Light corrected ');
set(gca,'YLim',[-4 2]);
title([ brw_str{n_inst},' - ',brw_str{n_ref},' / ',brw_str{n_ref}]);
%% Final days table
% m_etc_(lo que sea) son los datos filtrados por Tsync, OSC y sza !!! -> salida 
% de ETC_calibration En ste caso se usa summary, as? que la configuraci?n original 
% ser? 10, y la final 6
%%
m_etc=[m_etc_NEW(:,[1:4 7:8]),NaN*m_etc_NEW(:,1),m_etc_NEW(:,5:6)];
m_etc(:,7)=100*(m_etc(:,5)-m_etc(:,2))./m_etc(:,2);
m_etc(:,10)=100*(m_etc(:,8)-m_etc(:,2))./m_etc(:,2);
% La diferencia relativa anterior se refiere al ozono recalculado con la ETC que resulta de los c?lculos
% Esto podr?a ser una incoherencia: por ejemplo para el 072 sale 3179, pero mantenemos 3185.

%m_etc(:,[end-1,end])=[];
M=cell(size(m_etc,1)+1,size(m_etc,2)+1);
label={'Date','Day',['O3#',brw_str{n_ref}],'O3std','N',...
         ['O3#',brw_str{n_inst}],'O3 std',...
         ['%(',brw_str{n_inst},'-',brw_str{n_ref},')/',brw_str{n_ref}]...
         ['O3(*)#',brw_str{n_inst}],'O3std',...
         ['(*)%(',brw_str{n_inst},'-',brw_str{n_ref},')/',brw_str{n_ref}]};
M(1,:)=label;
M(2:end,1)=cellstr(datestr(fix(m_etc(:,1))));
m_etc(:,1)=diaj(m_etc(:,1));
M(2:end,2:end)=num2cell(round(m_etc*10)/10);

disp(M);
matrix2latex_ctable(M(2:end,2:end),fullfile(Cal.file_latex,['table_ETCdatafin_',Cal.brw_str{Cal.n_inst},'.tex']),...
                                   'rowlabels',M(2:end,1),...
                                   'columnlabels',label(2:end),'alignment','c','resize',0.9)

% So2 calibration from icf file
% Pendiente
% latexcmd(['>',file_latex,'_etc_',brw_str{n_inst}],'\ETCSOdos',config_def(12),'\SOdosABS',config_def(10),'\OtresSOdosABS',config_def(9));
%% Summary Final Comparison all data detailed for osc ranges
%%
ETC_SUG=ETC_NEW
blinddays=Cal.calibration_days{Cal.n_inst,2}
A1=A.old(ismember(Cal.Date.CALC_DAYS,blinddays),Cal.n_inst+1); 
A1_old=unique(A1(~isnan(A1)))
 


TIME_SYNC=5;
caption=strcat('Ozone Summary Report. Mean daily ozone, grouped by ozone slant path ranges,',...
               ' with original and final configuration (with an asterisk)');
tags_={'osc$>$1500' '1500$>$osc$>$1000' '1000$>$osc$>$700' '700$>$osc$>$400' 'osc$<$400'};
label_={'Day','osc range',['O3\#',brw_str{n_ref}],'O3std','N',...
                          ['O3\#',brw_str{n_inst}],'O3 std','\%diff',...
                          ['(*)O3\#',brw_str{n_inst}],'O3 std','(*)\%diff'};

ozone_osc_sum=o3_daily_osc(Cal,TIME_SYNC,n_ref,ETC_SUG,A1_old,summary_orig_old,summary_old,summary);
dat=cat(2,num2cell(ozone_osc_sum(:,1)),tags_(ozone_osc_sum(:,end))',num2cell(ozone_osc_sum(:,2:end-1)));

displaytable(dat,label_,12);
matrix2latex_ctable(dat,fullfile(Cal.file_latex,['table_summarydetailed_',brw_str{n_inst},'.tex']),...
                                  'size','tiny','columnlabels',label_,'alignment','c',...
                                  'format',{'%.0f','%s','%.0f','%.1f','%.0f','%.0f','%.1f','%.1f','%.0f','%.1f','%.1f'});
%% 
% [m,s,n,grpn]=grpstats(ozone_osc_sum,{ozone_osc_sum(:,1)},{'mean','std','numel','gname'}); 
% ozone_day_sum=round([m(:,1),m(:,2),s(:,2),m(:,4),m(:,8),s(:,8),100*(m(:,8)-m(:,2))./m(:,2)]*10)/10;
% 
% |makeHtmlTable(ozone_day_sum,[],cellstr(datestr(ozone_day_sum(:,1)+datenum(Cal.Date.cal_year,1,0))),...|
% 
% |       {'Day',['O3 #',brw_str{n_ref}],'O3 std','N obs',['O3 #',brw_str{n_inst}],...|
% 
% |        'O3 std',[' % ',brw_str{n_ref},'-',brw_str{n_inst},'/',brw_str{n_ref}]})|
%% Plot daily summary
%%
close all
figure; set(gcf,'Tag','_GlobalPlot_');
plot(summary{n_ref}(:,1),summary{n_ref}(:,6),'g*','MarkerSize',5);
hold on; plot(summary{Cal.n_inst}(:,1),summary{Cal.n_inst}(:,6),'b.','MarkerSize',12);
legend(gca,Cal.brw_name{n_ref},Cal.brw_name{Cal.n_inst},'Location','Best',...
                                                               'Orientation','Horizontal');
ylabel('Total Ozone (DU)'); xlabel('Date'); grid;
title(Cal.campaign);
datetick('x',6,'keepLimits','KeepTicks');

% blind days: suggested config.
f0=figure;
if ~Cal.no_maint(Cal.n_inst)     % if the instrument has changed due to maintenance
   inst1_b      =summary_old{Cal.n_inst}(findm(diaj(summary_old{Cal.n_inst}(:,1)),blinddays,0.5),:);
   inst1_b(:,10)=(inst1_b(:,8)-ETC_SUG(1).NEW)./(A1_old*inst1_b(:,3)*10);     
   for dd=1:length(blinddays)
      j=find(diajul(floor(inst1_b(:,1)))==blinddays(dd));
      j_=find(diajul(floor(ref(:,1)))==blinddays(dd));
      if (isempty(j) || length(j)<4), continue; end
         f=figure; 
         set(f,'Tag',sprintf('%s%s','DayPlot_',num2str(blinddays(dd))));
         plot(ref(j_,1),ref(j_,6),'g-s','MarkerSize',6,'MarkerFaceColor','g');
         hold on; plot(inst1_b(j,1),inst1_b(j,10),'b-d','MarkerSize',7,'MarkerFaceColor','b');
         plot(inst1_b(j,1),inst1_b(j,6),'r:.','MarkerSize',9);
         title(sprintf('Suggested configuration \n  Blind day %d%d',blinddays(dd), Cal.Date.cal_year-2000))
         legend(gca,Cal.brw_name{n_ref},Cal.brw_name{Cal.n_inst},[Cal.brw_name{Cal.n_inst} ' old config.'],'Location','SouthEast',...
                                                               'Orientation','Horizontal');
         ylabel('Total Ozone (DU)'); grid;
         datetick('x',15,'keepLimits','KeepTicks');
         set(gca,'YLim',[min(inst1_b(j,10))-8 max(inst1_b(j,10))+8])
   end
end
%% final days: final config.
%%
jday=findm(diaj(summary{Cal.n_inst}(:,1)),finaldays,0.5);% quiero mostrar la primera config. con sl
inst2=summary{Cal.n_inst}(jday,:);
f0=double(figure);
for dd=1:length(finaldays)
j=find(diajul(floor(inst2(:,1)))==finaldays(dd));
j_=find(diajul(floor(ref(:,1)))==finaldays(dd));
if (isempty(j) || length(j)<4), continue; end
f1=double(figure); 
set(f1,'Tag',sprintf('%s%s','DayPlot_',num2str(finaldays(dd))));
plot(ref(j_,1),ref(j_,6),'g-s','MarkerSize',6,'MarkerFaceColor','g');
hold on;
plot(inst2(j,1),inst2(j,6),'b--d','MarkerSize',7,'MarkerFaceColor','b');
plot(inst2(j,1),inst2(j,10),'r:.','MarkerSize',9);
title(sprintf('Final configuration \nFinal day %d%02d',finaldays(dd), Cal.Date.cal_year-2000))
legend(gca,Cal.brw_name{n_ref},Cal.brw_name{Cal.n_inst},[Cal.brw_name{Cal.n_inst} ' old config.'],...
                                                      'Location','SouthEast','Orientation','Horizontal');
ylabel('Total Ozone (DU)'); grid;
datetick('x',15,'keepLimits','KeepTicks');
set(gca,'YLim',[min(inst2(j,10))-20 max(inst2(j,10))+20])
end

figure(maxf(findobj('tag','_GlobalPlot_')));
printfiles_report(gcf,Cal.dir_figs,'Height',7,'Width',13);

printfiles_report(f0:f1,...
                              Cal.dir_figs,'Width',14.5,'Height',7.5);
graph2latex(Cal.file_latex,{'summaryplot','DayPlot_'},brw_str{n_inst},'scale',0.8);