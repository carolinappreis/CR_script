clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
B=char('BZ_ctx_probe.mat ','CZ_ctx_probe.mat ','SNR_ctx_probe.mat');
load 'animal_lesion_nolesion.mat'
lesion_numb=unique(str2num(A(lesion,4:6)));

for xxx=1:size(lesion_numb,1)
    for xx=1
        %1:size(B,1)
        name=B(xx,1:(find(B(xx,:)=='.')-1));
        load(name)
        
        idx=[];
        for i=1:size(WaveData_DCall,1)
            if size(WaveData_DCall{i,1},1)~=1
                idx=[idx i];
            end
        end
            animals{xx,1}=str2num(A(lesion(idx),4:6));
%         animal_rep{xxx,xx}=data(str2num(A(lesion(idx),4:6))==lesion_numb(xxx));
    end
end
for xx=1:size(B,1)
    name=B(xx,1:(find(B(xx,:)=='.')-1));
    load(name)
    
    idx=[];
    for i=1:size(WaveData_DCall,1)
        if size(WaveData_DCall{i,1},1)~=1
            idx=[idx i];
        end
    end
    animals{xx,1}=idx;
    
end




% visually pic rep animals
bcs=animal_rep(1,:);
bc=animal_rep (2:3,:);
bs=animal_rep(7,:);
% bs(cellfun('isempty',bs)) = []

clearvars -except bcs bc bs 
% load ('data_all.mat','freq')
% % save 'rep_animal_bua'

clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('rep_animal_bua.mat')
data1=bc;

for i=1:size(data1,1)
data{:,i}=cat(1, bc{:,i})
end

dif_angs=zeros(size(freq,1),size(filtBZ,2),size(filtSNR,2),size(filtSNR,3));
for t=1:size(freq,1)
    Fs=1000;
    for i=1:size(data,2)
        for r=1:size(data{i,:},1)
            [b,a]=butter(2,[(freq(t)-5)/(0.5*Fs) (freq(t)+5)/(0.5*Fs)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*Fs) 100/(0.5*Fs)],'bandpass');
            end
            filtall{i,:}(r,:)=filtfilt(b,a,data{i,:}(r,:));
            non_norm=squeeze(angle(hilbert(filtall{i,1}(1,:)))-angle(hilbert(filtall{i,1}(r,:))))';
            for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=non_norm(1,x)-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=non_norm(1,x)+(2.*pi);
                end
            end
            dif_angs{i,1}(r,:)=non_norm;
            
        end
        dif_angs{i,1}=dif_angs{i,1}(2:end,:);
        for r=1:size(dif_angs{i,:},1)
            euler{i,1}(r,:)=(sum(exp(sqrt(-1)*(dif_angs{i,1}(r,:))))./(length(non_norm)));
        end
        psi_rec{t,xx,1}(1,i)=abs(mean(euler{i,1}));
        ang_rec{t,xx,1}(1,i)=angle(mean(euler{i,1}));
    end
    
    euler1=cell2mat(euler);
    psi_all(t,xx,1)=abs(mean(euler1));
    ang_all{t,xx,:}=angle(euler1);
    
end


