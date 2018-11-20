function y=naninterp2(X)
% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.

x0=X(:,1);
x=X(:,2:end);
y=x;

for n=1:size(x,2)
    x_=x(:,n);
    y(isnan(x_),n) = interp1(x0(find(~isnan(x_))), x_(~isnan(x_)), x0(isnan(x_)),'linear',0);
end

y=[x0,y];

