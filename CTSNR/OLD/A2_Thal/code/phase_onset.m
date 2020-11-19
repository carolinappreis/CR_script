    
env=abs(hilbert(muafilt(1,:)));
pha=(angle(hilbert(muafilt)));
apha=abs(angle(hilbert(muafilt)));
threshold=prctile(env,75);
tt(size(env,1):size(env,2))=threshold;

indexexceed=find(env>threshold);
diffindex=diff(indexexceed);
pnts=find(diffindex>1);
begin=indexexceed(pnts+1);
ending=indexexceed(pnts);
begin2=[indexexceed(1) begin];
ending2=[ending indexexceed(end)];

m=1
e_burst=[];
b_burst=[];
for i=2:(length(begin2)-1)
    if  (ending2(i)-begin2(i))>= 50 &&(begin2(i+1)-ending2(i))>=150 && (begin2(i)-ending2(i-1))>=150
        b_burst(1,m)=begin2(i)
        e_burst(1,m)=ending2(i)
        m=m+1;
    end
end


[maxvalM,maxidxM] = findpeaks(muafilt);

pre_onset=[];
for b = 1:length(b_burst)
    for p=1:length(maxidxM)
        if min(abs(b_burst(b)-maxidxM(p)))<=30;
          pre_onset=[pre_onset p];
        end
    end
end

onset=maxidxM(pre_onset);

plot(time,muafilt)
hold on
plot(time(b_burst),env(e_burst),'b.')
plot(time(onset),muafilt(onset),'.')
plot(time,tt)
plot(time,env)

clearvars z output output out_evoked out_ctx
output=[];
surrogate=[];
for z=1:size(thal_raw,1)
    for j=1:length(onset)
        if onset(j)-2000>0 && onset(j)+2000<length(thal_raw)
             idx_sur=randi([2001,(length(env)-2000)],1,1);
            surrogate (z,j,1:4001)= (thal_raw(z,idx_sur-2000:idx_sur+2000)-median(thal_raw(z,idx_sur-2000:idx_sur+2000)))./median(thal_raw(z,idx_sur-2000:idx_sur+2000));
            output(z,j,1:4001)= (thal_raw(z,onset(j)-2000:onset(j)+2000));
        end
    end
end
output=surrogate;
out_evoked(1:size(output,1),1:size(output,3))=(sum(output,2))./size(output,2); % taking the mean of segments per channel

if isempty (out_evoked)==1
    output_ctx=[];
else  
    for j=1:length(onset)
        if onset(j)-2000>0 && onset(j)+2000<length(env)
            output_ctx(j,1:4001)=(env(onset(j)-2000:onset(j)+2000)-median(env(onset(j)-2000:onset(j))))./median(env(onset(j)-2000:onset(j)));
        else
            output_ctx=[];
        end
    end
end




