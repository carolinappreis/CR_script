
clear all

% cd('/Users/Carolina/OneDrive - Nexus365/BNDU_computer/Documents/Carolina_code/codes_thal/A3_Thal/mat')
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A3_Thal\mat')
load ('data_all' , 'freq')
load 'BZ_ctx_probe'

samprate=1000;

m=1;
for r=1:size(data,1);
    for rr=2:size(data{r,1},1)
        [power,f]=pwelch(data{r,1}(rr,:),1000,[],1000,1000);
        [Pxx_ind,F_ind]=mscohere(data{r,1}(1,:),data{r,1}(rr,:),samprate,[],samprate,samprate);
        filtrange=freq(3);
        [b,a]=butter(2,[(filtrange-5)/(0.5*samprate) (filtrange+5)/(0.5*samprate)],'bandpass');
        filtc=filtfilt(b,a,data{r,:}(1,:));
        filthal=filtfilt(b,a,data{r,:}(rr,:));

        if sum(Pxx_ind(filtrange-5:filtrange+5))./sum(Pxx_ind(1:end))>0.1

        non_norm=squeeze(angle(hilbert(filtc))-angle(hilbert(filthal)))';
            for x =1:size(non_norm,2)
                if non_norm(1,x)>pi
                    non_norm(1,x)=non_norm(1,x)-(2.*pi);
                elseif non_norm(1,x)<-pi
                    non_norm(1,x)=non_norm(1,x)+(2.*pi);
                end
            end
            dif_angs(m,:)=non_norm;
            ctx_filt(m,:)=filtc;
            p(m,:)=sum(Pxx_ind(filtrange-5:filtrange+5))./sum(Pxx_ind(1:end));
            m=m+1;

        end
    end
end

for ii=1:size(ctx_filt,1)
    
    env=abs(hilbert(ctx_filt(ii,:)));
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
        epoch=200;
    end
    
    for j=1:length(onset1)
        if onset1(j)>epoch && onset1(j)+epoch<length(env)
            pre_b{ii,1}(j,:)=dif_angs(ii,onset1(j)-epoch:onset1(j)+epoch);
        end
    end
    
    for n=1:(length(env)/100)
        idx_sur=randi([epoch+1,(length(env)-epoch)],1,1);
        pre_sur{ii,1}(n,:)= dif_angs(ii,idx_sur-epoch:idx_sur+epoch);
    end
end

b1=vertcat(pre_b{:});
s1=vertcat(pre_sur{:});

psi_b=abs(sum(exp(sqrt(-1)*(b1)),1)./size(b1,1));
psi_s=abs(sum(exp(sqrt(-1)*(s1)),1)./size(s1,1));
time= [1:2*epoch+1]; 
  
plot(time,psi_b)
hold on
plot (time,psi_s)
