clear; close
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','out')

for iii=1:10
d=out.fchange{iii,2}{1,1}; 
data= d(:) ; data(isnan(data)) = []; data=abs(data);

thresh=prctile(abs(out.fchange{iii,1}(1,:)),99.7917);

above_th(iii,1)=(numel(find(data>thresh))/numel(data))*100;

end
bar(above_th)
ylim([0 100])
ylabel({'percentage of stim freq abs change','>95prctl non-stim freq abs change'})
box('off')
xlabel('patients')

% for iii=1:10
% d=out.fchange{iii,2}{1,1}; 
% data= d(:) ; data(isnan(data)) = []; data=(data);
% 
% thresh1=prctile((out.fchange{iii,1}(1,:)),95);
% thresh2=prctile((out.fchange{iii,1}(1,:)),5);
% 
% above_th2(iii,1)=((numel(find(data>thresh1))+numel(find(data<thresh2)))/numel(data))*100
% 
% end



for iii=1:10
d=out.change_c{iii,2}{1,1}; 
data= d(:) ; data(isnan(data)) = []; data=abs(data);

thresh=prctile(abs(out.change_c{iii,1}(1,:)),95);

above_th(iii,1)=(numel(find(data>thresh))/numel(data))*100;

end
bar(above_th)
ylim([0 100])
ylabel({'percentage of stim amp abs change','>95prctl non-stim amp abs change'})
box('off')
xlabel('patients')