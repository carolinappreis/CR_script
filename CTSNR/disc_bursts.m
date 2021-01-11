function [data,nr_b,dur_b,dur_nb]=disc_bursts(dat,data,co1,iii,nr_b,dur_b,dur_nb,filtrange)
samprate=1000;
pre_b_in=zeros(size(dat,1),size(dat,2));
    for ct=1:size(dat,1)
%         [Pxx_ind,F_i]=pwelch(dat(ct,:),samprate,[],samprate,samprate);
%         frange=find(F_i==15):find(F_i==35);
%         Pxx_ind_beta=Pxx_ind(frange);
%         filtrange=14+find(Pxx_ind_beta==max(Pxx_ind_beta));
        [b,a]=butter(2,[(filtrange(iii,1)-5)/(0.5*samprate) (filtrange(iii,1)+5)/(0.5*samprate)],'bandpass');
        flt=filtfilt(b,a,dat(ct,:));
%          flt=dat(ct,:);
        env=abs(hilbert(flt));
        [onset1,offset1]=bursts(env);
        onset1=horzcat(onset1{:});
        offset1=horzcat(offset1{:});
        for jj=1:length(onset1)
            pre_b_in(ct,onset1(jj):offset1(jj))=1;
            n_b(1,ct)=numel(onset1);
            d_b(1,jj)=numel(onset1(jj):offset1(jj));
        end
        dur_b{co1,iii}(1,ct)=sum(d_b);
        dur_nb{co1,iii}(1,ct)=length(env)-sum(d_b);
        clear d_b flt onset1 offset1
    end
    data{co1,iii}=pre_b_in;
    nr_b{co1,iii}=n_b;
end
