l=dir('*.xls')
dobson=[];
for i=1:length(l)
   t=readtable(l(i).name); 
   dob=table2array(t);
   dobson=cat(3,dob);
   head{i}=l(i).name;
end
dobson64=dobson;
save('dobson64.mat','dobson64','head');