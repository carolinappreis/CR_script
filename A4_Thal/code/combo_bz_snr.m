% 
% 
% 
% clear all
% % cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal')
% load('BZ.mat');
% 
% % idx=[size(SNR.bua_coh,1)-1 size(SNR.bua_coh,1)];
%     for i=size(BZ.bua_coh,1)-1
% %         : size(BZ.bua_coh,1);
% %       plot(log(SNR.power_ctx(i,1:80)))
%  
%       for m=1:size(BZ.power_thal{i,1},1)
%       plot(log(BZ.power_thal{i,1}(m,1:80)))
%       
%       end
%     end

    
load('SNR.mat');
load('BZ.mat');
data1=BZ.phase_thal{11,1};
data2=SNR.phase_thal{13,1}([3 4],:);
env=SNR.env_ctx(13,:);
ref=SNR.onset_raw_all{13,1};
    m=1;
for a=1:size(data1,1)
    for b=1:size(data2,1)
        non_norm=data1(a,:)-data2(b,:);
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
                ep_b(m,z,w)=abs(mean(exp(sqrt(-1).*(non_norm(epochs_idx(z,w):epochs_idx(z,w)+ol)))));
                ep_t(m,w,:)=abs(mean(exp(sqrt(-1).*(epochs_t(:,w)))));
            end
        end
        
        for z= 1:size(epochs_idx_sur,1)
            for w=1:size(epochs_idx_sur,2)
                if  epochs_idx_sur(z,w)+ol<length(non_norm)
                    ep_b_s(m,z,w)=abs(mean(exp(sqrt(-1).*(non_norm(epochs_idx_sur(z,w):epochs_idx_sur(z,w)+ol)))));
                    ep_t_s(m,w,:)=abs(mean(exp(sqrt(-1).*(epochs_t_sur(:,w)))));
                end
            end
        end
        
      m=m+1;  
      non_norm=[];
    end
   
end

time=1:401;
plot(time,mean(ep_t,1))


ep_t1=reshape(ep_t,size(ep_t,1)*size(ep_t,2),size(ep_t,3));
ep_ts1=reshape(ep_t_s,size(ep_t_s,1)*size(ep_t_s,2),size(ep_t_s,3));

    b_l_rat=mean(ep_t,1);
    
    b_s_rat=mean(SNR.across_b_short{q,1},1);
    t_l_rat=mean(SNR.across_t_long{q,1},1);
    t_s_rat(q,:)=mean(SNR.across_t_short{q,1},1);
    b_surr_rat(q,:)=mean(SNR.across_b_surr{q,1},1);
    t_surr_rat(q,:)=mean(SNR.across_t_surr{q,1},1);



cd('C:\Users\creis\Documents\GitHub\CR_script\A4_Thal\code')
st=NaN(1,401);
clear A; A=t_l_rat; %b1{f,1};
clear B; B=t_surr_rat; %s1{f,1}(1:size(A,1),:);
hayriye_c; st(1,:)=stats.prob; st2(1,:)=stats.posclusterslabelmat;
beg=find(st(1,:)<0.01 & st2(1,:)~=0);
if ~isempty(beg)
    sig_rise_all=[beg(1) beg(end)]; 
end

y2=zscore(mean(t_l_rat)); 
y1=zscore(mean(t_l_rat)+std(t_l_rat)./sqrt(size(t_l_rat,1))); 
y3=zscore(mean(t_l_rat)-std(t_l_rat)./sqrt(size(t_l_rat,1)));
p2=plot(time, y2,'LineStyle','-', 'LineWidth',1.5,'Color',color_b)
patch([time fliplr(time)], [y1 fliplr(y2)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
patch([time fliplr(time)], [y2 fliplr(y3)],[color_b],'FaceAlpha',[0.2],'EdgeColor','none')
hold on
patch([sig_rise_all(1) sig_rise_all(2) sig_rise_all(2) sig_rise_all(1)],[2.5 2.5 2.6 2.6],color_b,'EdgeColor','none')

legend([p1 p2],{'short burst','long burst'})
ylabel ('PSI bursts (zscore)')
xlabel ('Time (msec)')
