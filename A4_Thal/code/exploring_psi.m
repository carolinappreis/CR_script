% clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
%
% %cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% region=char({'BZ';'BZ.'});
% load(region(2,:))
%
% for ik=1
% %     :size(BZ.filt_thal,1);
%     for ct=1:size(BZ.filt_thal{ik,1},1)
%         non_norm=BZ.phase_ctx(ik,1)-BZ.phase_thal{ik,1}(ct,:);
%         for x =1:size(non_norm,2)
%             if non_norm(1,x)>pi
%                 non_norm(1,x)=(non_norm(1,x))-(2.*pi);
%             elseif non_norm(1,x)<-pi
%                 non_norm(1,x)=(non_norm(1,x))+(2.*pi);
%             else
%                 non_norm(1,x)= non_norm(1,x);
%             end
%         end
%
%         ref=BZ.onset_phase_al{1,1}{2,1}(1);
%         steps=[50:5:550];
%         for i=1:length(steps)
%             output(ct,i)=abs(mean(exp(sqrt(-1).*(non_norm(ref:ref+steps(i)).*hanning(length(ref:ref+steps(i)))'))));
%         end
%     end
% end
%
% plot(steps,output')
% figure()
% boxplot(output')
% %
% test=[BZ.phase_thal{1,1}(1,:);BZ.phase_ctx(1,:)];
% testa=abs(mean(exp(sqrt(-1).*(test(1,:)-test(2,:)))));
% test1=[BZ.filt_thal{1,1}(1,:);BZ.filt_ctx(1,:)];
% n=computePPC(test);
% n1=computePPC(test1);



clear all
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
load('BZ_opt.mat');
%  load('BZ.mat');
short=1;long=2;


for ik=1:size(BZ.filt_thal,1);
    ref=[];
    for ct=1:size(BZ.filt_thal{ik,1},1)
        non_nomr=[];epochs_idx=[];epochs_t=[];
        ref=BZ.onset_raw{1,ik}{long,1};
        non_norm=BZ.phase_ctx(ik,:)-BZ.phase_thal{ik,1}(ct,:);
        for x =1:size(non_norm,2)
            if non_norm(1,x)>pi
                non_norm(1,x)=(non_norm(1,x))-(2.*pi);
            elseif non_norm(1,x)<-pi
                non_norm(1,x)=(non_norm(1,x))+(2.*pi);
            else
                non_norm(1,x)= non_norm(1,x);
            end
        end
        
        el=200;
        for ii=1:length(ref)
            if ref(ii)>el
                epochs_idx(ii,:)=ref(ii)-el:ref(ii)+el;
                epochs_t(ii,:)=non_norm(ref(ii)-el:ref(ii)+el);
            end
        end
        
        epochs_idx = epochs_idx(any(epochs_idx,2),:);
        epochs_t = epochs_t(any(epochs_t,2),:);
        
        for n=1:(length(non_norm)/1000)
            idx_sur=randi([el+1,(length(non_norm)-el)],1,1);
            epochs_idx_sur(n,:)= idx_sur-el:idx_sur+el;
            epochs_t_sur(n,:)= non_norm(idx_sur-el:idx_sur+el);
        end
        
        ol=50;
        for z= 1:size(epochs_idx,1)
            for w=1:size(epochs_idx,2)
                ep_b(ct,z,w)=abs(mean(exp(sqrt(-1).*(non_norm(epochs_idx(z,w):epochs_idx(z,w)+ol)))));
                ep_t(ct,w,:)=abs(mean(exp(sqrt(-1).*(epochs_t(:,w)))));
            end
        end
        
        for z= 1:size(epochs_idx_sur,1)
            for w=1:size(epochs_idx_sur,2)
                if  epochs_idx_sur(z,w)+ol<length(non_norm)
                    ep_b_s(ct,z,w)=abs(mean(exp(sqrt(-1).*(non_norm(epochs_idx_sur(z,w):epochs_idx_sur(z,w)+ol)))));
                    ep_t_s(ct,w,:)=abs(mean(exp(sqrt(-1).*(epochs_t_sur(:,w)))));
                end
            end
        end
        
    end
    if size(ep_b,2)==1
        ep_b=ep_b';
        ep_t=ep_t';
    end
    if size(ep_b,2)==1
        ep_b_s=ep_b_s';
        ep_t_s=ep_t_s';
    end
    if size(ep_b,1)==1
        BZ.across_b_long{ik,1}=squeeze(mean(ep_b,2))';
%         BZ.across_b_surr{ik,1}=squeeze(mean(ep_b_s,2))';
    else
        
        BZ.across_b_long{ik,1}=squeeze(mean(ep_b,2));
%         BZ.across_b_surr{ik,1}=squeeze(mean(ep_b_s,2));
        % BZ.across_t_short{ik,:}=ep_t;
        % BZ.across_t_surr{ik,:}=ep_t_s;
    end
        clearvars ep_b ep_t ep_b_s ep_t_s
    
end

clearvars -except BZ
% cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% save 'BZ_opt'
