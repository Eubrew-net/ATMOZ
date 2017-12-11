% Este script lee las variables ###_yy_ddd.mat que se generan automaticamente cada vez que se procesa
% un test de dispersion y produce un fichero de texto con la siguiente informacion (labels):
% 
%     'date','Brw','idx','wl_0','wl_2','wl_3','wl_4','wl_5','wl_6',... 
%     'fwhm_0','fwhm_2','fwhm_3','fwhm_4','fwhm_5','fwhm_6','cal_ozonepos','ozonepos',... 
%     'o3_0','o3_2','o3_3','o3_4','o3_5','o3_6'
%  
% Se ejecuta desde el mismo directorio donde se encuentran los test de dispersion clasificados segun 
% el convenio conocido (###_yy_ddd/orig). 
% No es necesario pasarle el fichero setup (el script genera de foma automatica la variable Cal.brw).
% No es necesario tener en el path de MATLAB el directorio de funciones /matlab (para eso se ponen en el 
% mismo directorio los scripts brewer_data.m y diames.m. 
% 
% Para procesar los test de dispersion ver dsp_summary_all.m (notar que la version de dspreport.m 
% es ligeramente diferente de la operativa).

%% DSP data
clear all;

l_all=dir(fullfile('.')); 
l_all=l_all(cellfun(@(x) ~isempty(x),cellfun(@(x) regexp(x,'^\d\d\d.'),extractfield(l_all,'name'),'UniformOutput' ,0)));
Cal.brw=unique(cellfun(@(x) sscanf(x,'%d%d%d%*s'),extractfield(l_all,'name')));
Cal.n_brw=length(Cal.brw);

Cal.brw_str=cellfun(@(x) num2str(x,'%03d'),num2cell(Cal.brw),'UniformOutput',0);
Cal.brw_name=cellfun(@(x) strcat('#',x),Cal.brw_str,'UniformOutput' ,0);

labels={'date','Brw','idx','wl_0','wl_2','wl_3','wl_4','wl_5','wl_6',... 
        'fwhm_0','fwhm_2','fwhm_3','fwhm_4','fwhm_5','fwhm_6',... 
        'cal_ozonepos','ozonepos','o3_0','o3_2','o3_3','o3_4','o3_5','o3_6',...
        'r_0','r_2','r_3','r_4','r_5','r_6'};

%% dsp data
close all

wv_matrix=NaN*ones(length(l_all),length(labels));
idx=1;
for brwi=1:Cal.n_brw
    %%
    l=dir(fullfile('.',[Cal.brw_str{brwi},'*'])); 
    ldsp=cellstr(cat(1,l.name));
    for indx=1:length(ldsp)  
        info=sscanf(ldsp{indx},'%d_%d_%d'); info_=brewer_date(str2double(sprintf('%03d%02d',info(3),info(2))));
        try
           load(fullfile(pwd,ldsp{indx},strcat(ldsp{indx},'.mat')));
           % date Brw idx wl_0 wl_2 wl_3 wl_4 wl_5 wl_6 fwhm_0 fwhm_2 fwhm_3 fwhm_4 fwhm_5 fwhm_6 
           % cal_ozonepos ozonepos o3_0 o3_2 o3_3 o3_4 o3_5 o3_6
           % rayleight   r_0 r_2 r_3 r_4 r_5 r_6
           wv_matrix(idx,:)=cat(2,info_(1),dsp_sum.brewnb,indx,dsp_sum.salida.QUAD{end-1}.thiswl,dsp_sum.salida.QUAD{end-1}.fwhmwl/2,...
                                  dsp_sum.salida.QUAD{end-1}.cal_ozonepos,dsp_sum.salida.QUAD{end-1}.ozone_pos,dsp_sum.salida.QUAD{end-1}.o3coeff,dsp_sum.salida.QUAD{end-1}.raycoeff);        
        catch exception
           wv_matrix(idx,:)=cat(2,info_(1),info(1),indx,NaN*ones(1,26)); 
        end
        idx=idx+1;
    end
end

%% Ordenamos en tuplas brewer - fecha 
dsps=NaN*ones(size(wv_matrix,1),size(wv_matrix,2)); 

dsps(:,[1:2 4:end])=sortrows(wv_matrix(:,[1:2 4:end]),[2 1]);
dsps(:,3)=wv_matrix(:,3); 

fid=fopen(fullfile('.','dsp_summ_are2015.txt'),'w');
fprintf(fid,'%%date Brw idx wl_0 wl_2 wl_3 wl_4 wl_5 wl_6 fwhm_0 fwhm_2 fwhm_3 fwhm_4 fwhm_5 fwhm_6 cal_ozonepos ozonepos o3_0 o3_2 o3_3 o3_4 o3_5 o3_6 r_0 r_2 r_3 r_4 r_5 r_6\n\r');
for l=1:size(dsps,1)
    fprintf(fid,'%f %d %d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %d %d %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n\r',dsps(l,:));
end
fclose(fid);
