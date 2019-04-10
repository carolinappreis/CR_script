clear all
%cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'SNR_ctx_probe'

Fs=1000;
for t=1:size(freq,1)
    for ii=1:size(data,1)
        for r=1:size(data{ii,:},1)
            [b,a]=butter(2,[(freq(t)-5)/(0.5*Fs) (freq(t)+5)/(0.5*Fs)],'bandpass');
            if t==length(freq)
                [b,a]=butter(2,[49/(0.5*Fs) 60/(0.5*Fs)],'bandpass');
            end
            non_norm=angle(hilbert(filtfilt(b,a,data{ii,:}(1,:))))-(angle(hilbert(filtfilt(b,a,data{ii,:}(r,:)))));
            for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=(non_norm(1,x))-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=(non_norm(1,x))+(2.*pi);
                else
                    non_norm(1,x)= non_norm(1,x);
                end
            end
            dif_angs1(r,:)=non_norm;
            clearvars non_norm;
        end
        dif_angs{ii,1}=dif_angs1(2:end,:);
        
        
        %---- bursts
        
        env=abs(hilbert(filtfilt(b,a,data{ii,:}(1,:))));
        threshold=prctile(env,75);
        tt(size(env,1):size(env,2))=threshold;
        indexexceed=find(env>threshold);
        diffindex=diff(indexexceed);
        pnts=find(diffindex>1);
        begin=indexexceed(pnts+1);
        ending=indexexceed(pnts);
        begin2=[indexexceed(1) begin];
        ending2=[ending indexexceed(end)];
        
        ind_b=[];
        for i=1:(length(begin2))
            if (ending2(i)-begin2(i))>=100
                ind_b=[ind_b i];
            end
        end
        
        if ~isempty (ind_b)
            begin3=begin2(ind_b);
            ending3=ending2(ind_b);
            
            space_betb=200;
            ind_b1=[];
            for i=1:(length(begin3)-2)
                if (begin3(i+1)-ending3(i))>=space_betb && (begin3(i+2)-ending3(i+1))>=space_betb
                    ind_b1=[ind_b1 i+1];
                end
            end
            if (begin3(2)-ending(1))>=space_betb
                ind_b1= [1 ind_b1];
            end
            if (begin3(length(begin3))-ending3(length(begin3)-1))>=space_betb
                ind_b1=[ind_b1 length(begin3)];
            end
            onset1=begin3(ind_b1);
            offset1=ending3(ind_b1);
            %-------
            epoch=500;
            
            for s=1:size(dif_angs{ii,:},1)
                for j=1:length(onset1)
                    if onset1(j)>epoch && onset1(j)+epoch<length(env)
                        dumo=dif_angs{ii,1}(s,onset1(j)-epoch:onset1(j));
                        dumi=dif_angs{ii,1}(s,onset1(j):onset1(j)+(round((offset1(j)-onset1(j))/2)));
                        pre_bi{t,ii}(s,j,:)= (sum(exp(sqrt(-1)*(dumi))))./(size(dumi,2));
                        pre_bo{t,ii}(s,j,:)= (sum(exp(sqrt(-1)*(dumo))))./(size(dumo,2));
                    end
                end
                
                for n=1:(length(env)/100)
                    idx_sur=randi([epoch+1,(length(env)-2*epoch)],1,1);
                    pre_sur{t,ii,1}(s,n,:)= dif_angs{ii,1}(s,idx_sur+2*epoch);
                end
            end
        end
    end
end

time= [1:2*epoch+1];
for f=1:size(pre_bi,1)
    bi=[];
    bo=[];
    s=[];
    for i=1:size(pre_bi,2)
        bi=[bi ; reshape(pre_bi{f,i},[size(pre_bi{f,i},1)*size(pre_bi{f,i},2) size(pre_bi{f,i},3)])];
        bo=[bo ; reshape(pre_bo{f,i},[size(pre_bo{f,i},1)*size(pre_bo{f,i},2) size(pre_bo{f,i},3)])];
        s=[s ;reshape(pre_sur{f,i},[size(pre_sur{f,i},1)*size(pre_sur{f,i},2) size(pre_sur{f,i},3)])];
        stats(f,1)=ranksum(angle(bi),angle(bo));
        b1{f,1}=bi;
        b2{f,1}=bo;
    end
end

for i=1:2
    for f=1:size(pre_bi,1)
        if i==1
            figure(i)
            subplot(size(b1,1),1,f)
            polarhistogram(angle(b1{f,1}))
        else
            figure(i)
            subplot(size(b2,1),1,f)
            polarhistogram(angle(b2{f,1}))
        end
    end
end




