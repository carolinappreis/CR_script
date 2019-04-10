clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
load ('data_all.mat','freq')
A=char('AV_ctx_probe.mat','BZ_ctx_probe.mat ','CZ_ctx_probe.mat ','SNR_ctx_probe.mat','Rt_ctx_probe.mat ');
for xx=1:size(A,1)
    name=A(xx,1:(find(A(xx,:)=='.')-1));
    load(name)
   for t=1:size(freq,1)
       
    Fs=1000;
    for i=1:size(data,1)
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
end

    
bar(psi_all)
xticklabels({'5-15','16-26','27-37','38-48','49-100'})
legend ('VA-CTX','BZ-CTX','CZ-CTX','SNR-CTX','RT-CTX',"FontSize", 12)
xlabel ('Frequencies(Hz)',"FontSize", 14)
ylabel ('PSI',"FontSize", 14)
box off
% for i=1:size(psi_rec,1)
% subplot(size(psi_rec,1),1,i)
% bar([psi_rec{i,1}])
% end

% figure()
% for ii=1:3
% for i=1:4
%     figure(ii)
%     subplot(1,size(ang_all,1),i)
%     polarhistogram(ang_all{i,ii},12)
% end
% end


   