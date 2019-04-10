% clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\Juxta SUA_act_mat\mat')
load('output_count_bz.mat')
%
%  m=output_count{1,1}(1,:);
time=1:401;
%  plot(time(m==1),m(m==1),'r.')

for f=1:size(output_count,1)
    for i=1:size(output_count{f,1},1)
        m=[];
        m=output_count{f,1}(i,:);
        idx=find(m==1);
        m(idx)=steps(i);
        figure(1)
        plot(time(idx),m(idx),'r.')
        hold on
    end
    xlim ([0 400])
    close all
end
