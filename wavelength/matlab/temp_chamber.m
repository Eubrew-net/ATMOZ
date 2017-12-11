leg={'Date'
    'AAAAMMDDhhmmss'
    'BFILE'
    'min (GMT)'
    'B temp1'
    'B temp2'
    'B temp2'
    'Measure'
    'Port'
    'Slit'
    'Fiilter#1'
    'Filter#2'
    'Micro'
    'Wavelengt'
    'Dark'
    'CY'
    'counts'
    'counts/sec'
    ' T_Almemo'
    'H_Almemo'
    'T_PT100'
    'Monitor'
    'Monitor_u'
    'Hgtest'
    'Hgtest'};

load CC_B185.txt
ploty(CC_B185(:,[1,5,6,7,19,21]));
legend(leg([5,6,7,19,21]));
datetick
xlabel('Chamber')
ylabel('Brewer')


figure
plot(CC_B185(:,19),CC_B185(:,5:7),'.')
rline
xlabel('Chamber')
ylabel('Brewer')
legend(leg([5,6,7]))

figure
plot(CC_B185(:,21),CC_B185(:,5:7),'.')
xlabel('Chamber PT100')
ylabel('Brewer')
legend(leg([5,6,7]))
rline

%% Monitor vs Brewer temperture
ft=fix(CC_B185(:,1)*6-CC_B185(1,1)*6);% every 4 hours
h=gscatter(CC_B185(:,5),CC_B185(:,22),ft);
set(h,'Marker','o');
title('Monitor vs Brewer Temperature (color indicates time ');
ylabel('Ration %');
xlabel('Brewer Temperature');
grid

%% Interpolation

m=grpstats(CC_B185,CC_B185(:,1));
