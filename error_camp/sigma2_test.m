function [unc_bp,unc_sg] = sigma2_test(brw)


bp = load('BP.dat');
sg = load('SG.dat');

T = -45.0;

lambda_bp = bp(:,1); 

c_0_bp = bp(:,2);
c_1_bp = bp(:,3);
c_2_bp = bp(:,4);

sec_bp = (1e-20)*(c_0_bp + c_1_bp*T + c_2_bp*T^2);

bp_0 = bp(:,5);
bp_1 = bp(:,6);
bp_2 = bp(:,7);
bp_3 = bp(:,8);
bp_4 = bp(:,9);

lambda_sg = sg(:,1);

c_0_sg = sg(:,2);
c_1_sg = sg(:,3);
c_2_sg = sg(:,4);

sec_sg = (1e-20)*(c_0_sg + c_1_sg*T + c_2_sg*T^2);

sg_0 = sg(:,5);
sg_1 = sg(:,6);
sg_2 = sg(:,7); 
sg_3 = sg(:,8);
sg_4 = sg(:,9);

sigma_bp = (bp_0 + bp_1*T + bp_2*T^2 + bp_3*T^3 + bp_4*T^4)*1e-20;

sigma_sg = (sg_0 + sg_1*T + sg_2*T^2 + sg_3*T^3 + sg_4*T^4)*1e-20;

a=strcat('dsp_',brw,'_17.csv');
data=csvread(a,1,0);

dsp_summ=data(:,:);

% pesos del o3
O3W=[   0   0.00   -1.00    0.50    2.20   -1.70];

s_bp(size(dsp_summ,1),6)=NaN;
sig_bp(size(dsp_summ,1),1)=NaN;

s_sg(size(dsp_summ,1),6)=NaN;
sig_sg(size(dsp_summ,1),1)=NaN;

for nb=1:size(dsp_summ,1)
    for ns=1:6
        y_bp=trapezoid_brewer(lambda_bp,dsp_summ(nb,3+ns)/10,dsp_summ(nb,9+ns)/20,.87);
        s_bp(nb,ns)=trapz(lambda_bp,sigma_bp.*y_bp);
	    y_sg=trapezoid_brewer(lambda_sg,dsp_summ(nb,3+ns)/10,dsp_summ(nb,9+ns)/20,.87);
        s_sg(nb,ns)=trapz(lambda_sg,sigma_sg.*y_sg);
    end
    sig_bp(nb,1)=-O3W*(2.687e+19)*squeeze(s_bp(nb,:)')/log(10);
	sig_sg(nb,1)=-O3W*(2.687e+19)*squeeze(s_sg(nb,:)')/log(10);
end

bp_bp=[dsp_summ(:,1),squeeze(sig_bp)];
sg_sg=[dsp_summ(:,1),squeeze(sig_sg)];

unc_bp = mean(bp_bp(:,2));
unc_sg = mean(sg_sg(:,2));