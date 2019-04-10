
clear all
cd('C:\Users\creis\OneDrive - Nexus365\BNDU_computer\Documents\Carolina_code\codes_thal\A4_Thal\mat')
load('BZ_beta_clustb.mat')
d=clust_b;
for i=1:size(d,1)
    [m,e(i)]=max(d(i,180:220));
end
for i=1:size(d,1)
    dd(i,:)=(d(i,:)-mean(d(i,:)))./std(d(i,:));
    dd_max(i)=dd(i,e(i)+179);
end
%P=sum(Pxx(:,16:36),2)./sum(Pxx(:,1:501),2);
sub=find(dd_max'>1);
for k=1:length(sub)
    dum_var2=diff(smooth(d(sub(k),:),101));
    dum_var4=dum_var2(1800:e(sub(k))+1799);
    x0=find(dum_var4>0);
    if ~isempty(find(diff(x0)>1))
        XXXstr(k,1)=x0(max(find(diff(x0)>1))+1)-199;
    elseif isempty(find(diff(x0)>1)) && ~isempty(x0)
        XXXstr(k,1)=x0(1)-199;
        onset=180;
        while x0(1)==1
            dum_var4=dum_var2(onset-5:e(sub(k))+179);
            x0=find(dum_var4>0);
            onset=onset-5;
            if ~isempty(find(diff(x0)>1))
                XXXstr(k,1)=x0(max(find(diff(x0)>1))+1)-(200-onset+1);
            elseif isempty(find(diff(x0)>1)) && ~isempty(x0)
                XXXstr(k,1)=x0(1)-(200-onset+1);
            end
            
        end
    elseif isempty(x0)
        clear dum_var4
        dum_var4=dum_var2(180:220);
        x00=find(dum_var4>0);
        if ~isempty(x00)
            XXXstr(k,1)=x00(1)-19;
        elseif isempty(x00)
            XXXstr(k,1)= NaN;
        end
    end
    clear dumvar2 dumvar4 x0
end
% mean(XXXstr)