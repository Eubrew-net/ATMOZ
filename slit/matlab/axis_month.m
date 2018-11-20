function axis_month(h)

Meses = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';
        'Aug';'Sep';'Oct';'Nov';'Dic'];

if nargin==0
    h=gca;
end
set(h,'Xtick',1:12);
set(h,'XtickLabel',Meses);
set(h,'XLim',[1,12]);
