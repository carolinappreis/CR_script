clearvars -except surr time subjects muafilt A nfile_i state_mat nn x med_thal_cell med_ctx chanofinterest B WaveData WaveData_DC ind_p nbursts

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

ind_b1=cell(2,1);
for i=1:(length(begin3)-1)
    if (ending3(i)-begin3(i))<=median_b && (begin3(i+1)-ending3(i))>=200
        ind_b1{1}(1,i)=i;
    end
    if (ending3(i)-begin3(i))>=median_b && (begin3(i+1)-ending3(i))>=200
        ind_b1{2}(1,i)=i;
    end
end


for i=1:size(ind_b1,1)
    onset1{i,1}=begin3(ind_b1{i}(ind_b1{i}(1,:)~=0));
    offset1{i,1}=ending3(ind_b1{i}(ind_b1{i}(1,:)~=0));
end

duration1=cell(2,1);
for i=1:2;
    duration1{i,1}=offset1{i}-onset1{i};
end

% plot(time,env,time,threshold*ones(length(env),2))
% hold on
% plot(time(onset1{1}),env(onset1{1}),'r.')
% plot(time(offset1{1}),env(offset1{1}),'b.')
% plot(time(onset1{2}),env(onset1{2}),'y.')
% plot(time(offset1{2}),env(offset1{2}),'k.')
% plot(time,env)
% hold on
% plot(time,tt)
% plot(time(onset1),env(onset1),'b*')
% plot(time(offset1),env(offset1),'r*')

%check2, check3

if size(muafilt,1)~=1
    
    thal=(size(muafilt,1));
    
    for z=2:thal
        thal_env(z-1,:)=abs(hilbert(muafilt(z,:)));
    end
    
    clearvars z start1 end1 output surrogate med_thal ouput_ctx med_ctx_array surrogate1
    start1=cell(2,1);
    end1=cell(2,1);
    for j=1:2;
        for i=1:size(onset1{j},2)
            if (onset1{j}(1,i)-500)>0
                start1{j}(1,i)=i;
            end
            if (onset1{j}(1,i)+500)<length(env)
                end1{j}(1,i)=i;
            end
        end
    end
    
    start=zeros(2,1);
    stop=zeros(2,1);
    for i =1:2;
        start(i,1)=min((start1{i}(start1{i}(1,:)~=0)));
        stop(i,1)=max((end1{i}(end1{i}(1,:)~=0)));
    end
    
    for j =1:2;
        onset{j,1}=onset1{j}(start(j,1):stop(j,1));
    end
    
    for k=1:2;
        for i =1:length(onset{k});
            figure(1)
            subplot(2,1,k)
            plot(env(onset{k}(i)-500:onset{k}(i)+500))
            hold on
        end
    end
    
    
    for k=1:2;
        for j=1:length(onset{k})
            output_ctx{k,1}(j,1:1001)=(env(onset{k}(j)-500:onset{k}(j)+500)-median(env(onset{k}(j)-500:onset{k}(j))))./median(env(onset{k}(j)-500:onset{k}(j)));
        end
    end
    
    
    output=[];
    surrogate=[];
    for k=1:2;
        for z=1:(thal-1)
            for j=1:length(onset{k})
                if onset{k}(j)-500>0 && onset{k}(j)+500<length(thal_env)
                    idx_sur=randi([501,(length(env)-500)],1,1);
                    surrogate{k,1}(z,j,1:1001)= (thal_env(z,idx_sur-500:idx_sur+500)-median(thal_env(z,idx_sur-500:idx_sur+500)))./median(thal_env(z,idx_sur-500:idx_sur+500));
                    output{k,1}(z,j,1:1001)= ((thal_env(z,onset{k}(j)-500:onset{k}(j)+500)-median(thal_env(z,onset{k}(j)-500:onset{k}(j))))./median(thal_env(z,onset{k}(j)-500:onset{k}(j))));
                end
            end
        end
    end
    
    if size(surrogate{1},2) >= size(surrogate{2},2)
        surrogate=surrogate{1};
    else
        surrogate=surrogate{2};
    end
    
    
    %output{1}(1,1,:) 1st probe, 1st segment , condition 1 (short)
    
    
    %     med_thal=vertcat(median(reshape(output{1},size(output{1},1).*size(output{1},2),size(output{1},3))),...
    %         median(reshape(output{2},size(output{2},1).*size(output{2},2),size(output{2},3))));
    
    surrogate1=[];
    for k=1:2
        for i=1:size(output{k},1)
            med_thal{k,1}(i,1:1001)=median(output{k}(i,:,1:1001))
            surrogate1(i,1:1001)=median(surrogate(i,:,1:1001))
        end
    end
    
    
    figure(2)
    subplot(2,1,1)
    plot(med_thal{1,:}')
    subplot(2,1,2)
    plot(med_thal{2,:}')
    hold on
    
    if isempty (output_ctx) | isempty (med_thal)
        med_ctx_array=zeros(1,1001);
    elseif ~isempty (output_ctx) && size(output_ctx,1)==2
        med_ctx_array=zeros(2,1001);
        for i=1:2;
            med_ctx_array(i,:)=median(output_ctx{i})
        end
    end
    
else
    med_thal=[];
    med_ctx_array=[];
end

close all