clear; close
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/aux_out.mat','x_all')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/cluster_out_mc.mat','clust')
load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/mean_amplitude_5sec.mat','norm_5env_t')


for ii=1:10
    stim_coef{ii,1}=x_all{ii,1}(50001:end,:);
end
clear x_all

 sub=find(~isnan(clust.win));

for ii=1:length(sub)
    cl=clust.C{sub(ii),2};
      data=stim_coef{sub(ii),1};
    % % data=(m_env5_time{sub(ii),1})'./(10*9.81/0.5);
%      data=(norm_5env_t{sub(ii),1})';
    f1=figure(1)
    subplot(length(sub),1,ii)
    for m=1:size(data,1)
        x=1:size(data,1);
        colo=[[1 0 0];[0 0 1];[0 0 0]];
        for a =1:3
            if cl(m)==1
                plot(x(m),data(m,a),'.','MarkerSize',4,'Color',colo(a,:))
            else
                plot(x(m),data(m,a),'o','MarkerSize',4,'Color',colo(a,:))
            end
            hold on
        end
                   ylim([-1 1])
        box('off')
        title(['patient',num2str(sub(ii))])
    end
    clear x data cl a 
    f1= set(f1,'color','w');
end



%     colo=[[1 0 0];[0 0 0]];
%      x=1:size(env,1);
%   for m=1:size(env,1)
%       for a=1:3
%           subplot(3,1,a)
%       plot(x(m),env(m,a)/(10*9.81/0.5),'.','Color',colo(cl(m),:))
%       hold on
%       ylim([0 0.1])
%       box('off')
%       end
%   end


% % % 
% % % clear
% % % load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/mean_amplitude.mat','m_env5_time','mean_amp_ns','std_amp_ns')
% % % for sub=1:10
% % %             mean_stim(sub,1)=nanmean(m_env5_time{sub,1}(1,:));
% % %             std_stim(sub,1)=std(m_env5_time{sub,1}(1,:));
% % % end
% % % 
% % % m=mean([[mean_stim]';[mean_amp_ns]],1);
% % % std=mean([[std_stim]';[std_amp_ns]],1);