function [xx]=pha_sta(d,s,iii)

mm=find(d.ds_dall(2,:)>4.9);

% t=1:length((d.ds_dall(2,:)));
% plot(t,d.ds_dall(2,:))
% hold on
% plot(t(mm),d.ds_dall(2,mm),'.')

nn=s.phase{iii,1}(match(cond,iii),mm);  

xx=circ_r(nn');
end