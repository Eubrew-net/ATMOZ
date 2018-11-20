function [ o3xsec ] = o3xsec_temp( bp,temp )
%function [ o3xsec ] = o3xsec_temp( bp,temp )
%   Detailed explanation goes here
  j=find(fix(bp.temp)==temp); 
  if ~isempty(j)
      o3xsec=[bp.lamda,bp.o3x(:,j)];
  else
      o3xsec=[bp.lamda,NaN*bp.lamda];
      %o3xsec=[bp.lamda,interp1(bp.temp,bp.o3x',temp')'];
  end
  

