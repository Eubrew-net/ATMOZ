function [x_up,x_dwn,fwhm,resup]=dsp_brewer(a,s,cutmin,cutmax,pl);
% function [x_up,x_down,fwhm]=dsp(a,s,cutmin,cutmax,pl);
% 2 12 96 julian
% calculates disperison constants for brewer 119
% day is julian day.
% a and s obtained from alldsp.

% 4 12 96 corrects for deadtime, and converts to cnts per sec.
% 2 5 97 julian add symmetric conditioning
% 7 10 97 julian discover error in backward scan.
% 17 10 97 julian add %changes in point adding...
% 11 12 97 julian Add test for not complete line.
% 19 12 97 julian change how dspchi is called up.
% 5 3 98 julian new dspchi3 withour fmins is much faster now.

resup=0;
if nargin<2,s='';end
if nargin<5,pl=[];end
if isempty(pl),pl=0;end  % default no plot
if nargin<3,cutmin=[];end
if isempty(cutmin),cutmin=0.2;end
if nargin<4,cutmax=[];end
if isempty(cutmax),cutmax=0.8;end



indmax=find(a(:,1)==max(a(:,1)));indmax=indmax(1);
up=a(1:indmax,2);
%down=a(indmax+1:size(a,1),2);down=down(size(down,1):-1:1);
if indmax<length(a)
    down=a(size(a,1):-1:(indmax+1),2);
else
    down=up;
end
steps=a(1:indmax,1);  % get steps
try
cnts=[up down];
catch
    warning('DSP: Not complete measuremnt ');
    x_up=0;%steps(indmax(1));
    x_dwn=0;%steps(indmax(2));
    fwhm=nan;
    return
end
[maxcnts,indmax]=max(cnts);  % get max value and position.
ind1=find((cnts(:,1)<(cutmax*maxcnts(1))) & (cnts(:,1)>(cutmin*maxcnts(1))));
ind2=find((cnts(:,2)<(cutmax*maxcnts(2))) & (cnts(:,2)>(cutmin*maxcnts(2))));

% remove background here: 9 3 98 julian
%indback1=find(cnts(:,1)<(cutmin*maxcnts(1)));
%indback2=find(cnts(:,2)<(cutmin*maxcnts(1)));
%if ~isempty(indback1),cnts(:,1)=cnts(:,1)-mean(cnts(indback1,1));end
%if ~isempty(indback2),cnts(:,2)=cnts(:,2)-mean(cnts(indback2,2)); end


bufa=ind1(ind1<indmax(1));   %pup1=polyfit(steps(bufa),cnts(bufa,1),1)
bufb=ind1(ind1>indmax(1));   %pup2=polyfit(steps(bufb),cnts(bufb,1),1)
if length(bufa)<2 | length(bufb)<2 %| length(bufa)>6 | length(bufb)>6,
    warning('DSP: Not enough or too many points on either side of the maximum of the line');
    x_up=0;%steps(indmax(1));
    x_dwn=0;%steps(indmax(2));
    fwhm=nan;
    return
end

if length(bufa)~=length(bufb), % try to find one more point closest to cut
    if length(bufa)>length(bufb),
        iminus=min([max(bufb)+1 length(steps)]);  % get two new points around cut. changed 12 11 97
        iplus=max([min(bufb)-1 1]);
        if (abs(cnts(iminus,1)-cutmin*maxcnts(1))/(cutmin*maxcnts(1)))<(abs(cnts(iplus,1)-cutmax*maxcnts(1))/(cutmax*maxcnts(1))),
            bufb=[bufb;iminus];
        else
            bufb=[iplus;bufb];
        end
    else
        iminus=max([min(bufa)-1 1]);
        iplus=min([max(bufa)+1 length(steps)]);
        if (abs(cnts(iminus,1)-cutmin*maxcnts(1))/(cutmin*maxcnts(1)))<(abs(cnts(iplus,1)-cutmax*maxcnts(1))/(cutmax*maxcnts(1))),
            bufa=[iminus; bufa];
        else
            bufa=[bufa; iplus];
        end
    end
    ind1=[bufa;bufb];
end
%pup=fmins('dspchi2',[1e3 1e3 1],[],[],steps(bufa),cnts(bufa,1),steps(bufb),cnts(bufb,1))
pup=dspchi3(steps(bufa),cnts(bufa,1),steps(bufb),cnts(bufb,1)); % 5 3 98 julian

resupa=(polyval([pup(3) pup(1)],steps(bufa))-cnts(bufa,1))/maxcnts(1);
resupb=(polyval([-pup(3) pup(2)],steps(bufb))-cnts(bufb,1))/maxcnts(1);


%(jul(2)-jul(1))/(2*jul(3))
bufa=ind2(ind2<indmax(2));    %pdwn1=polyfit(steps(buf),cnts(buf,2),1);
bufb=ind2(ind2>indmax(2));    %pdwn2=polyfit(steps(buf),cnts(buf,2),1);
if length(bufa)<2 | length(bufb)<2,
    warning('DSP: Not enough points on either side of the maximum of the line');
    x_up=0;%steps(indmax(1));
    x_dwn=0;%steps(indmax(2));
    fwhm=nan;
    return
end

if length(bufa)~=length(bufb), % try to find one more point closest to cut
    if length(bufa)>length(bufb),
        iminus=min([max(bufb)+1 length(steps)]);  % get two new points around cut.
        iplus=max([min(bufb)-1 1]);
        if (abs(cnts(iminus,2)-cutmin*maxcnts(2))/(cutmin*maxcnts(2)))<(abs(cnts(iplus,2)-cutmax*maxcnts(2))/(cutmax*maxcnts(2))),
            bufb=[bufb;iminus];
        else
            bufb=[iplus;bufb];
        end
    else
        iminus=max([min(bufa)-1 1]);
        iplus=min([max(bufa)+1 length(steps)]);
        if (abs(cnts(iminus,2)-cutmin*maxcnts(2))/(cutmin*maxcnts(2)))<(abs(cnts(iplus,2)-cutmax*maxcnts(2))/(cutmax*maxcnts(2))),
            bufa=[iminus; bufa];
        else
            bufa=[bufa; iplus];
        end
    end
    ind2=[bufb;bufa];
end
%pdwn=fmins('dspchi2',[1e3 1e3 1],[],[],steps(bufa),cnts(bufa,2),steps(bufb),cnts(bufb,2));
pdwn=dspchi3(steps(bufa),cnts(bufa,2),steps(bufb),cnts(bufb,2)); % 5 3 98 julian

resdwna=(polyval([pdwn(3) pdwn(1)],steps(bufa))-cnts(bufa,2))/maxcnts(2);
resdwnb=(polyval([-pdwn(3) pdwn(2)],steps(bufb))-cnts(bufb,2))/maxcnts(2);

x_up=(pup(2)-pup(1))/(2*pup(3));
x_dwn=(pdwn(2)-pdwn(1))/(2*pdwn(3));

fwhm=(pup(2)+pup(1)-maxcnts(1))/pup(3);
%fwhm=2*(x_up-(maxcnts(1)/2-pup(1))/pup(3));% 6 6 97 julian gives same results!!!

fwhm_up=(pup(2)+pup(1)-(pup(1)+pup(3)*x_up))/pup(3); % 19 8 97 julian fwhm from ideal triangle
fwhm_dw=(pdwn(2)+pdwn(1)-(pdwn(1)+pdwn(3)*x_dwn))/pdwn(3); % 19 8 97 julian fwhm from ideal triangle
fwhm=[fwhm_up,fwhm_dw];
%(pup(1)+pup(3)*x_up)


% OR:
% calculate intersection of fitted lines.
%x_up=(pup2(2)-pup1(2))/(pup1(1)-pup2(1))
%x_dwn=(pdwn2(2)-pdwn1(2))/(pdwn1(1)-pdwn2(1))


% get data for plotting
%x_up=steps(ind1);x_dwn=steps(ind2);
%y_up1=polyval(pup1,x_up);y_dwn1=polyval(pdwn1,x_dwn);


% 21 8 2000: Julian plot differences to linefit
%subplot(2,1,1);
if pl
    %figure;
    plot(steps,cnts(:,1),'r',steps,cnts(:,2),'b',steps(ind1),cnts(ind1,1),'rx',steps(ind2),cnts(ind2,2),'bo',[x_up;x_up],[0;maxcnts(1)*1.2],'r',[x_dwn;x_dwn],[0;maxcnts(2)*1.2],'b');
    %ax=axis;
    %subplot(2,1,2);
    %plot(steps(ind1),[resupa;resupb],'rx',steps(ind2),[resdwna;resdwnb],'bx');
    %axis;axis([ax(1:2) ans(3:4)]);
    %subplot(2,1,1);

    g=text(steps(indmax(1))+5,maxcnts(1),sprintf('%5.2f',x_up));set(g,'color','red','fontsize',9);
    g=text(steps(indmax(2))+5,maxcnts(2)*0.9,sprintf('%5.2f',x_dwn));set(g,'color','blue','fontsize',9);
    g=text(steps(indmax(2))-10,maxcnts(2)*1.1,sprintf('FWHM:%3.2f',fwhm(1)));set(g,'color','red','fontsize',9);
    g=text(steps(indmax(2))+10,maxcnts(2)*1.1,sprintf('FWHM:%3.2f',fwhm(2)));set(g,'color','blue','fontsize',9);
    %buf=axis;axis([steps(indmax(1))-200 steps(indmax(1))+200 buf(3:4)]);

    %title(filename);
    xlabel('steps');
    %ylabel('photons per sec');
    title(s);
    %ylabel(sprintf('%f',sum(abs([resupa;resupb]))));
end

resup=sum(abs([resupa;resupb]));



function f=dspchi3(x1,y1,x2,y2,sig1,sig2)
% function f=dspchi3(x1,y1,x2,y2,sig1,sig2)
% f=[x0 x1 x2], x0=b1, x1=b2, x2=a   for b+ax
% calculates slopes for isosceles triangle for brewer dsp.
% sig1 sig2 are uncertainities.
% 5 3 98 julian gives same result as fmins(dspchi2)

if nargin<6,sig2=[];end
if nargin<5,sig1=[];end

if isempty(sig1),sig1=ones(size(x1));end
if isempty(sig2),sig2=ones(size(x2));end


%X=[x1(:) x2(:)];
%Y=[y1(:) y2(:)];
%SIG=[sig1(:) sig2(:)];

S(1)=sum(1./sig1.^2);S(2)=sum(1./sig2.^2);

Sx(1)=sum(x1./sig1.^2);Sx(2)=sum(x2./sig2.^2);

Sy(1)=sum(y1./sig1.^2);Sy(2)=sum(y2./sig2.^2);

Sxx(1)=sum(x1.^2./sig1.^2);Sxx(2)=sum(x2.^2./sig2.^2);

Sxy(1)=sum(x1.*y1./sig1.^2);Sxy(2)=sum(x2.*y2./sig2.^2);


M=[sum(Sxx) Sx(1) -Sx(2);...
    -Sx(1)   -S(1)   0   ;...
    Sx(2)     0   -S(2) ];

C=[Sxy(1)-Sxy(2);-Sy(1);-Sy(2)];
A=M\C; % A(1) is slope, A(2)=b1, A(2)=b2;

f=[A(2) A(3) A(1)];


