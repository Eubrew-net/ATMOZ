function [data,head]=read_ptb_ab(file)
s=fileread(file);
t1=textscan(s,'','Delimiter','\t','HeaderLines',1);
 s=strtok(s, char(10)); 
 head=strsplit(s,'\t');
%files are not ended by crlf 
% t1{end}=[t1{end};NaN];
 data=cell2mat(t1);