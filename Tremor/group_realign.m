clear all
load ('C:\Users\creis\OneDrive - Nexus365\Periph_tremor_data\F_group.mat')

LS1=squeeze(mean(LS,2));
for i=1:size(S,1);
phase_peak(1,i)=find(S(i,:)==max(S(i,:)));
if phase_peak(i)==1;
    s_re_alg(i,:)=S(i,:);
%     ns_re_alg(i,:)=NS(i,:);
    ls_re_alg(i,:)=LS1(i,:);
else
s_re_alg(i,:)=[S(i,phase_peak(i):end) S(i,1:phase_peak(i)-1)];
% ns_re_alg(i,:)=[NS(i,phase_peak(i):end) NS(i,1:phase_peak(i)-1)];
ls_re_alg(i,:)= [LS1(i,phase_peak(i):end) LS1(i,1:phase_peak(i)-1)];
end

% [p,h1]=ttest(s_re_alg,ns_re_alg);
% s_ns=find(h1<(0.05/size(s_re_alg,1)))
[p,h2]=ttest(s_re_alg,ls_re_alg);
s_rs=find(h2<(0.05/size(s_re_alg,1)))
end
% 
% 
% figure()
% subplot(1,2,1)
% bar(mean(s_re_alg))
% ylim([-0.1 0.3])
% subplot(1,2,2)
% bar(mean(ns_re_alg))
% ylim([-0.1 0.3])