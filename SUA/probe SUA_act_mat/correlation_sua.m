clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\SUA\probe SUA_act_mat')
load('SNR_sua_skrate.mat') ;
load('BZ_sua_skrate.mat') ;
subj= SNR.beta_rats(ismember(SNR.beta_rats,BZ.beta_rats));

for i =1:size(subj,1)
    for ii=1:size(SNR.beta_idx{subj(i),1},2)
        hp=SNR.beta_idx{subj(i),1}(1,ii);
        s{i,1}(ii,:)=SNR.sua{subj(i),1}(hp,:);
        clear hp
    end
end

for i =1:size(subj,1)
    for ii=1:size(BZ.beta_idx{subj(i),1},2)
        hp=BZ.beta_idx{subj(i),1}(1,ii);
        b{i,1}(ii,:)=BZ.sua{subj(i),1}(hp,:);
        clear hp
    end
end

clearvars -except b s


% b_2=(find(b{1,1}(2,:)==1))';
% s_2=(find(s{1,1}(1,:)==1))';
m=1;
for rat=1:size(b,1);
    for i=1:size(b{rat,1},1)
        for ii=1:size(s{rat,1},1)
            
            b_2=(find(b{rat,1}(i,:)==1))';
            s_2=(find(s{rat,1}(ii,:)==1))';
            
            seg_pwr=8;
            lag_tot=150;
            sec_tot=100;
            samp_rate=1000;
            freq=100;
            lag_neg=50;
            ch_max=0.1;
            opt_str='';
            
%             [f,t,cl]=sp2_m1(0,s_2,b_2,sec_tot,samp_rate,seg_pwr,opt_str);
%             psp_ch1(f,cl,freq,ch_max)
            hold on
            % m=m+1;
            % cor(m,:)=cl;
        end
        clear b_2 s_2
    end
    [f,t,cl]=sp2_m1(0,s_2,b_2,sec_tot,samp_rate,seg_pwr,opt_str);
    psp2(f,t,cl,freq,lag_tot,lag_neg,ch_max)
    
    % psp_ph1(f4,cl4,freq)
end