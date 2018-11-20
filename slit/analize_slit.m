% analize slits

%l=dir('slits037/*.dat')
fpath='lazer_singles_2015';
l=dir(fullfile(fpath,'*.dat'))
slit=[];s=[];
wv=280:0.05:400;
for i=1:length(l)
  
     
   slit(i).raw=load(fullfile(fpath,l(i).name));
   sinfo=sscanf(l(i).name,"%03d_slit%d.dat")
   slit(i).brw=sinfo(1)
   slit(i).nslit=sinfo(2) 
   f=figure;
   %title(fullfile(fpath,l(i).name))
   semilogyy(slit(i).raw)
    set(gca,'Xlim',[3000,3500]);
    set(gca,'Ylim',[1E-5,1.1]);
   [slit(i).wl,slit(i).fwhm]=slitfit(slit(i).raw ,0,1,double(f));
   hold all
   slit(i).slr=min(slit(i).raw(:,2));
   hline(min(slit(i).raw(:,2)),'-',num2str(slit(i).slr));
   vline(slit(i).wl,'-k',num2str(slit(i).fwhm));
   % in nm
   slit(i).fit=get_slit_fit2(slit(i).fwhm/10,slit(i).wl/10,slit(i).slr,wv);
   semilogy(10*wv,slit(i).fit,'x')
   title({fullfile(fpath,l(i).name),sprintf(' %f %f %f ',slit(i).wl,slit(i).fwhm,slit(i).slr)})
   set(gca,'Xlim',[3000,3500])
   
   s=scan_join(s,slit(i).raw);
   
   
   %% Scan tipo
   figure
   load slit_are2015
   semilogyy(slit_are2015);
   hold on
   semilogy(wv*10,get_slit_fit2(5.3/10,325.03,0.3E-4,wv),'.r')
   semilogy(wv*10,get_slit_fit2(5.3/10,325.03,0.5E-4,wv),'.r')
   semilogy(wv*10,get_slit_fit2(5.3/10,325.03,1.3E-4,wv),'.-b')
   semilogy(wv*10,get_slit_fit2(5.3/10,325.03,0.15E-4,wv),'.k')
   %set(gca,'Xlim',[3100,3400]);
   set(gca,'Ylim',[1E-5,1]);
   
   %% Scan tipo
   figure
   ploty(slit_are2015);
   hold on
   plot(wv*10,get_slit_fit2(5.3/10,325.03,0.3E-4,wv),'.-r','linewidth',4)
   plot(wv*10,get_slit_fit2(5.3/10,325.03,0.5E-4,wv),'.r')
   plot(wv*10,get_slit_fit2(5.3/10,325.03,1.3E-4,wv),'.-b','linewidth',3)
   plot(wv*10,get_slit_fit2(5.5/10,325.03,0.15E-4,wv),'-k','linewidth',4)
   set(gca,'Xlim',[3230,3270]);
   set(gca,'Ylim',[1E-5,1]);
   
   
end

% load xsections


  