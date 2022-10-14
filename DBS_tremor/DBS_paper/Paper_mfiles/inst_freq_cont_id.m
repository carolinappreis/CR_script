
function[dev_freq, sev]=inst_freq_cont_id(s,start,ending,tp_s,tp_e,match_ax,pt)

%%%% you need 1 second (1000 ponits segments)

peaki=[NaN 4 NaN 3.5];

dum2=s.freq{1,3};
dum3=s.env_acc{1,3};
freqi=dum2(match_ax(2,pt,1),:);
ch_fr=abs(freqi-peaki(pt));
dev_freq =[mean(ch_fr) std(ch_fr)];

[ref]=pt_post(dum2, start, ending,tp_s,tp_e); %%% uncomment plot to see breaks in data

r=[];

%%% stick data again across trials and pauses
for t=1:size(ref,1)
    pr=[];
    for j = 2:2:length(ref{t,1})
        if ~isempty(length(ref{t,1}))
            period=ref{t,1}(j-1):ref{t,1}(j);
            r=[r period];
            pr=[pr period];
        end
        clear period tempo
    end
        trial_menv(t,:)=[mean(dum3(match_ax(2,pt,1),pr)) std(dum3(match_ax(2,pt,1),pr))];
end

sev=mean(trial_menv);


n=ch_fr(r);
if max(n)<4
    p=[];
    for ii=1:length(n)
        if n(ii)<0.5
            p(1,ii)=1;
        elseif n(ii)>=0.5 && n(ii)<1
            p(1,ii)=2;
        elseif n(ii)>=1 && n(ii)<2
            p(1,ii)=3;
        else
            p(1,ii)=4;
        end
    end
else
    error('outlier')
end


for ii=1:size(p,1)
    total=[numel(find(p==1))./length(p) numel(find(p==2))./length(p) numel(find(p==3))./length(p) numel(find(p==4))./length(p)];
end

f1=figure(2)
pie(total)
box('off')
set(f1,'color','w')
legend('{<0.5HZ}','{0.5-1Hz}','{1-2Hz}','{>2Hz}')
legend('boxoff')
title(sprintf('patient %d',(pt)))


end






