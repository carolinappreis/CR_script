%a)generate beta bursts from  filtered EEG
%b)exclude bursts less than 50msec
%c)select the ones isolates for 150msec before and after
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

thr=150;
ind_b1=[];
for i=1:(length(begin3)-2)
    if (begin3(i+1)-ending3(i))>=thr && (begin3(i+2)-ending3(i+1))>=thr
        ind_b1=[ind_b1 i+1];
    end
end

if (begin3(2)-ending(1))>=thr
    ind_b1= [1 ind_b1];
end

if (begin3(length(begin3))-ending3(length(begin3)-1))>=thr
    ind_b1=[ind_b1 length(begin3)];
end
    
onset1=begin3(ind_b1);
offset1=ending3(ind_b1);

% plot(time,env)
% hold on
% plot(time,tt)
% plot(time(begin2),env(begin2),'r.')
% plot(time(begin3),env(begin3),'y*')
% plot(time(onset1),env(onset1),'bo')

[maxvalM,maxidxM] = findpeaks(muafilt);

pre_onset=cell(1,1);
for b = 1:length(onset1)
    for p=1:length(maxidxM)
        if min(abs(onset1(b)-maxidxM(p)))<=30;
          pre_onset{1,b}=p;
        end
    end
end
pre_onset=cell2mat(pre_onset);
onset=maxidxM(pre_onset);

% plot(time,muafilt)
% hold on
% plot(time(onset1),env(onset1),'bo','MarkerSize', 5)
% plot(time(onset),env(onset),'r.','MarkerSize', 10)
% plot(time,tt)
% plot(time,env)

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
out_evoked(1:size(output,1),1:size(output,3))=(sum(output,2))./size(output,2);

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




