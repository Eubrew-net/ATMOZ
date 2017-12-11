function [jp,F,F0]=readb_j3(bfile,measure)
%  Input : fichero B, configuracion (ICF en formato ASCII)
%  Output:
%
%
% Alberto Redondas 2014
%if nargin==1 
    DT=2.9E-8
%end

if nargin==1
measure='j1';
end

dt=[]; rs=[]; jp=[];
fmtdt='dto3  %c%c%c %f/ %f %d:%d:%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f  ';
fmtrs='rso3  %c%c%c %f/ %f %d:%d:%d ';
fmt_js=['js %*s %d %f %d %d %d %d %d %d %d %d %d %d rat %d %d %d %d'];
fmt_ju=['ju %*s %d %f %d %d %d %d %d %d %d %d %d %d rat %d %d %d %d'];
%fmt_xx=['xx %*s %d %f %d %d %d %d %d %d %d %d %d %d %d rat %f %f %f %f'];
fmt_xx=['xx %*s %d %f %d %d %d %d %d %d %d %d %d %d  rat %f %f %f %f'];
fmt_jp=['jp %c %d %f %d %d %d %d %d rat %f %f %f %f'];
fmt_ds=['ds %*s %d %f %d %d %d %d %d %d %d %d %d %d rat %d %d %d %d'];
fmtsum=['summary %d:%d:%d %c%c%c %f/ %f %f %f %f %c%c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f']; % summary format

fmt=strrep(fmt_xx,'xx',measure);

%leemos el fichero en memoria
f=fopen(bfile);
if f < 0
    disp(['Error',bfile]);
    return
end
    s=fread(f);
    fclose(f);
    s=char(s)';
    [PATHSTR,NAME,EXT] = fileparts(bfile);
    fileinfo=sscanf([NAME,EXT],'%c%03d%02d.%03d');
    datefich=datejul(fileinfo(3),fileinfo(2));
%     datestr(datefich(1));

    l=mmstrtok(s,char(10));
    jjp=strmatch(measure,l);
   % jjps=strfind(l(jjp),'summary');
   % j_pss=find(~cellfun(@isempty,jjps));
    
%     jds=strfind(l,'ds');
%     j_ds=find(~cellfun(@isempty,jds));
%     jdss=strfind(l(j_ds),'summary');
%     j_dss=find(~cellfun(@isempty,jdss));
%     j_dss=j_ds(j_dss);
%% js is preceded by a ds measurement 
%  Temperature and airmass calculation can be take from previous ds
%  we need the DSP to have the wavelength
%  then we have 20 measurements taken at slit 1 wavelength 
%  
    
    jp=[];
    if ~isempty(jjp)
        jp_str=l(jjp);
        for i=1:length(jjp)
           dt_=sscanf(l{jjp(i)},fmt);
           try
             jp=[jp,dt_];
           catch
            %solo las ultimas son buenas   
            jp=[];
            jp=[jp,dt_];
            disp(l{jjp(i)})
           end
        end         
    else
        jp=[];
    end
    

x=jp;
%j_c=find(diff(x(2,:))<-1);  % date change
%jp(2,:)=datefich(1)+jp(2,:)/60/24; % we have to correct for LAT
%jp(2,j_c:end)=jp(2,j_c:end)+1;
IT=0.1147;
if ~isempty(x)
 F.time=x(2,:);
 F.dark=x(7,:);
F.mic=x(end,:);
F.temp=x(13:14,:);
F.wv=x(15,:);
F.CY=x(5,:);

CY=x(5,:);
F.raw=x(6:12,:);
F_0=matadd(F.raw,-F.dark);
%F_dark=x(7,1);
 Fcs= 2*matdiv(matadd(F.raw,-F.dark),F.CY)/IT;
  Fcs(Fcs<0)=2;
  Fcs(Fcs>1E07)=1E07;
  
  % dead time correction
  F0=Fcs;
  for j=1:9
    Fcs=F0.*exp(Fcs*DT);   
  end
jp=[jp;Fcs];
F.cs=Fcs;

%[i,j,k]=unique(F.mic);
% meds=[0,cumsum(diff(k)==-3)'];
%[F.m,F.s]=grpstats([(datefich(1)+F.time/60/24)',mean(F.temp)',F.mic',F.cs'],{meds,F.mic})
end
%%  