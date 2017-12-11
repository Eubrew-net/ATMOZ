load('/Users/aredondas/Google Drive/ATMOZ/Laboratory/TemperatureChamber/ZD.185');

figure
ZD(ZD(:,17)<1E6,17)=NaN;
for i=2:5
   
  si=(ZD(:,10)==i);
  norm=ZD(si,17); 
  if ~isempty(norm)
   plot(ZD(si,5),ZD(si,17)/norm(100),'o')
  end
   hold on
end



%datetick('x','ddd hh:mm');