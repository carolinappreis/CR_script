%From this code we geto to:
%a)generate beta bursts from  filtered EEG
%b)exclude bursts less than 50msec
%c)slect the ones isolates for 150msec before and after
%d)phase align them to the beat peak with allowing for a lag of max 30msec
%e)get trials 4sec long centred around phase-aligned burst onset


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

ind_b=[];
for i=1:(length(begin2))
    if (ending2(i)-begin2(i))>50
        ind_b=[ind_b i];
    end
end

begin3=begin2(ind_b);
ending3=ending2(ind_b);

duration=ending3-begin3;
median_b=median(duration);
SD_b=std(duration);

ind_b1=[];
for i=2:(length(begin3)-1)
    if (begin3(i+1)-ending3(i))>=150 && (begin3(i)-ending3(i-1))>=150
        ind_b1=[ind_b1 i];
    end
end

onset1=begin3(ind_b1);
offset1=ending3(ind_b1);

% plot(time,env)
% hold on
% plot(time,tt)
% plot(time(begin3),env(begin3),'r.')
% plot(time(onset1),env(onset1),'ro')

[maxvalM,maxidxM] = findpeaks(muafilt);
Inv = 1.01*max(muafilt) - muafilt;
[minval,minidx] = findpeaks(Inv);

pre_onset=cell(1,1);
for b = 1:length(onset1)
    for p=1:length(minidx)
        if min(abs(onset1(b)-minidx(p)))<=30;
          pre_onset{1,b}=p;
        end
    end
end
pre_onset=cell2mat(pre_onset);
onset=minidx(pre_onset);

plot(time,muafilt)
hold on
plot(time(onset1),env(offset1),'bo')
plot(time(onset),env(onset),'r.')
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
%output=surrogate;
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




