clearvars -except time muafilt A nfile_i state_mat nn x med_thal_cell med_ctx chanofinterest B WaveData WaveData_DC ind_p

env=abs(hilbert(muafilt(1,:)));
threshold=prctile(env,75);
tt(size(env,1):size(env,2))=threshold;
foward=500;
backward=500;

% V_transition_pnts script
indexexceed=find(env>threshold);
diffindex=diff(indexexceed);
pnts=find(diffindex>1);
begin=indexexceed(pnts+1);
ending=indexexceed(pnts);

begin2=[indexexceed(1) begin];
ending2=[ending indexexceed(end)];

duration=ending2-begin2;
mean_duration=mean(duration);
% histfit(duration, length(duration))
% mean(duration)
% plot(time,env,time,threshold*ones(length(env),2))
% hold on
% plot(time(onset),env(onset),'r.')
% plot(time(offset),env(offset),'b.')

m=1
offset=[];
onset=[];
for i=2:(length(begin2)-1)
    if  (ending2(i)-begin2(i))>200 && (begin2(i+1)-ending2(i))>=150 && (begin2(i)-ending2(i-1))>=150   % ending2(i)-begin2(i))>=200 ending2(i)-begin2(i)>150 && ending2(i)-begin2(i)<=250 % ending2(i)-begin2(i)>=250
        %if (ending2(i)-begin2(i))>=50 && (ending2(i)-begin2(i))<=200 && (begin2(i+1)-ending2(i))>=150 && (begin2(i)-ending2(i-1))>=150  
        onset(1,m)=begin2(i)
        offset(1,m)=ending2(i)
        m=m+1;
    end
end

% plot(time,env)
% hold on
% plot(time,tt)
% plot(time(onset),env(onset),'b*')
% plot(time(offset),env(offset),'r*')

%check2, check3

thal=(size(muafilt,1));

for z=2:thal
    thal_env(z-1,:)=abs(hilbert(muafilt(z,:)));
end

clearvars z output surrogate output1 med_thal_sub ouput_ctx med_ctx_sub
start1=[];
end1=[];

for i=1:length(onset)
    if (onset(i)-500)>0
    start1=[start1 i];
    end
    if (onset(i)+500)<length(env)
        end1=[end1 i];
    end
end
start=start1(1);
stop=end1(end);

onset=onset(start:stop);
    
    for j=1:length(onset)
            output_ctx(j,1:1001)=(env(onset(j)-500:onset(j)+500)-median(env(onset(j)-500:onset(j))))./median(env(onset(j)-500:onset(j)));
    end

output_ctx( ~any(output_ctx,2), : ) = [];

output=[];
surrogate_c=[];
surrogate_t=[];
for z=1:(thal-1) 
    for j=1:length(onset)
        if onset(j)-500>0 && onset(j)+500<length(thal_env)
            idx_sur=randi([501,(length(env)-500)],1,1);
            surrogate_t (z,j,1:1001)= (thal_env(z,idx_sur-500:idx_sur+500)-median(thal_env(z,idx_sur-500:idx_sur+500)))./median(thal_env(z,idx_sur-500:idx_sur+500));
            surrogate_c (j,1:1001)= (env(idx_sur-500:idx_sur+500)-median(env(idx_sur-500:idx_sur+500)))./median(env(idx_sur-500:idx_sur+500));
            output(z,j,1:1001)= ((thal_env(z,onset(j)-500:onset(j)+500)-median(thal_env(z,onset(j)-500:onset(j))))./median(thal_env(z,onset(j)-500:onset(j))));
        end
    end
end

output=surrogate_t;       
%organized by:probe1, segments68, probe2 , segments68 (...)
m=1;
output1=zeros(size(output,1).*size(output,2),size(output,3));
for i=1:size(output,1)
    for h=1:size(output,2)
        output1(m,:)=output(i,h,:);
        m=m+1;
    end
end
output1( ~any(output1,2), : ) = [];
med_thal_sub=median(output1);
%check4 - (size(muafilt,1)-1).*length(onset)==size(output1,1)


%median of probes'signal for each bursts (organized by segments - segment 1
%of cortex and median o all probes during that time)
m=1;
median_probes=zeros(size(onset,2),1001);
for i=1:size(output,2)
    index_b=i:size(output,2):size(output1,1)
    median_probes(m,:)=median(output1(index_b,:));
    m=m+1;
end
%check4
median_probes( ~any(median_probes,2), : ) = [];

if isempty (output_ctx) | isempty (output1)
    med_ctx_sub=zeros(1,1001);
elseif ~isempty (output_ctx) && size(output_ctx,1)==1
    med_ctx_sub=(output_ctx);
else ~isempty (output_ctx) && size(output_ctx,1)>=2
    med_ctx_sub=median(output_ctx);
end


%check5


