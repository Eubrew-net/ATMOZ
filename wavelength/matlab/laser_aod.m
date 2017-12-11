%laser ozone mode99999u
J={};
for i=1:10
   J{i}=readb_j3('bdata185/B02716.185',sprintf('y%d',i));
   if i==10
   J{i}=readb_j3('bdata185/B02716.185','ya');
   end
end
for i=11:24    
    ch=char(87+i);
    J{i}=readb_j3('bdata185/B02816.185',sprintf('y%c',ch));
end

%% measurement check
s=cellfun(@(x) size(x,2)',J)
%%


% j2 repetido
aux1=J{2}(:,1:77);
aux2=J{2}(:,78:end);
J{2}=aux1;
J{end+1}=aux2;
% elimianmos las pruebas

m=cat(2,J{:});
%% Raw counts
figure
semilogy(m(15,:),m(6:12,:))
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer Raw counts')
title('Brewer ozone mode')

%& Brewer counts second
figure
norf=max(max(m(17:end,:)));
semilogy(m(15,:),m(17:end,:)/norf,'.')
legend(cellstr(num2str([0:6]')))
grid
xlabel('wavelength (set)')
ylabel('Brewer  counts/sec')
title('Brewer ozone mode  (counts/sec)')
hold on
%% no linearity correction
nl=m(17:end,:);
for i=1:7
  nl(i,:)=monitor_brewer(nl(i,:));
end
lc=m(17:end,:)./nl;
norf=max(max(lc));

semilogy(m(15,:),lc/norf,'o-')
% linearity correction
% for counts/sec >12000