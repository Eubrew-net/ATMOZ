[jl1,Fl1,Fl1_0]=readb_jl('bdata185/B02216.185');
[jl2,Fl2,Fl3_0]=readb_jl('bdata185/B02316.185');
[jl3,Fl3,Fl3_0]=readb_jl('bdata185/B02416.185');
[jl4,Fl4,Fl4_0]=readb_jl('bdata185/B02516.185');

% summaries

F=[Fl1.m;Fl2.m;Fl3.m;Fl4.m];
Fs=[Fl1.s;Fl2.s;Fl3.s;Fl4.s];
J=[jl1,jl2,jl3,jl4];

save jl_sum.txt F -ascii
save jl_ssum.txt Fs -ascii
save jl_obs.text J -ascii


% remove bad temperature
F(F(:,2)>100,:)=[];
plot(F(:,1),F(:,2),'o');
