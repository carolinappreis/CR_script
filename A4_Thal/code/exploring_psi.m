% clear all
% % cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% region=char({'SNR';'BZ.'});
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
%         ref=BZ.ctxb_phase_al{1,5}{2,1}(5);
%         steps=[50:25:550];
%         for i=1:length(steps)
%             output(ct,i)=abs(mean(exp(sqrt(-1).*(non_norm(ref:ref+steps(i))))));
%         end
%     end
% end

% plot(steps,output')
% figure()
% boxplot(output')
%
% test=[SNR.phase_thal{1,1}(1,:);SNR.phase_ctx(1,:)];
% testa=abs(mean(exp(sqrt(-1).*(test(1,:)-test(2,:)))));
% test1=[SNR.filt_thal{1,1}(1,:);SNR.filt_ctx(1,:)];
% n=computePPC(test);
% n1=computePPC(test1);



clear all
% cd('C:\Users\creis\Documents\GitHub\CRcode\codes_thal')
cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal')

% load('BZ.mat');
load('SNR.mat');


for ik=1:size(SNR.filt_thal,1);
    ref=[];
    for ct=1:size(SNR.filt_thal{ik,1},1)
        non_nomr=[];epochs_idx=[];epochs_t=[];
        ref=SNR.onset_raw{1,ik}{2,1};
        non_norm=SNR.phase_ctx(ik,:)-SNR.phase_thal{ik,1}(ct,:);
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
        epochs_idx = epochs_idx(any(epochs_idx,2),:);
        
        ol=50;
        for z= 1:size(epochs_idx,1)
            for w=1:size(epochs_idx,2)
                ep_b(ct,z,w)=abs(mean(exp(sqrt(-1).*(non_norm(epochs_idx(z,w):epochs_idx(z,w)+ol)))));
                ep_t(ct,w,:)=abs(mean(exp(sqrt(-1).*(epochs_t(:,w)))));
            end
        end
    end
    SNR.across_b_long{ik,1}=squeeze(mean(ep_b,2));
    SNR.across_t_long{ik,:}=ep_t;
end

clearvars -except BZ
cd('/Users/Carolina/Documents/GitHub/CRcode/codes_thal')
save 'BZ'
