l=dir('data_set_2')
date_time= @(x) datetime(x,'ConvertFrom','datenum')
addpath(genpath(fullfile('~','CODE','rbcce.aemet.es','iberonesia','matlab')));
plt= {'o','+','*','h','x','s','d','v','>','<','p','+','x','*','x','s'};
colo=num2cell(colorcube(16),2)
x1=plt(1:15)';
%%
l=dir('data_set_2')
inst={l.name};
inst(1:2)=[]
fun_name=strcat('read_',inst)

figure;hold all

for ii=1:length(inst)
    n_inst=ii;
    disp(inst{ii})
    eval(fun_name{ii});
   
     %t = dobson(:,2);
     %m = dobson(:,23);
     %mu = dobson(:,22);
     %Omega = dobson(:,24);
    %  data from/to 12-30 september
    %t_set_1{ii}=t_set_1{ii}(t_set_1{ii}.Date>datenum(2016,9,12) & t_set_1{ii}.Date<datenum(2016,9,30),:);
    %plot(t_set_1{ii}{:,2},t_set_1{ii}{:,3});
    
    % error analysis 
    %er_ad{n_inst} = D064_orig(:,{'Time','Date','M_AD','Mu_AD','O3_AD'});
    %er_cd{n_inst} = D064_orig(:,{'Time','Date','M_CD','Mu_CD','O3_CD'});

    t_set_2{ii}=  [table2array(er_ad{n_inst}(:,2:end));  table2array(er_cd{n_inst}(:,2:end))];
    
    ploty(t_set_2{ii}(:,[3,4]),plt{ii})
    
end
xlabel('air mass')
ylabel('ozone cm')
total_set_2=cell2mat(t_set_2');


