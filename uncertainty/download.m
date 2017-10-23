function x = download(url)

%url='"http://rbcce.aemet.es/eubrewnet/data/get/DS?brewerid=185&date=2016-09-12&enddate=2016-09-30&format=text"';
[a,b]=system(['curl -s --connect-timeout 120 ',url]);
if a==0
    x=textscan(b,'','headerlines',1,'delimiter',',TZa','commentstyle','matlab','TreatAsEmpty','None');
else
    disp('error');
    disp(b)
end
%x=cell2mat(data_sc);