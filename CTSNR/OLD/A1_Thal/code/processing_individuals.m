clearvars -except time muafilt A nfile_i state_mat nn x med_thal_cell med_ctx chanofinterest B WaveData WaveData_DC ind_p n median_indv beta_coh beta_coh_all

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
    if  (ending2(i)-begin2(i))>=150 && (begin2(i+1)-ending2(i))>=150 && (begin2(i)-ending2(i-1))>=150   % ending2(i)-begin2(i))>=200 ending2(i)-begin2(i)>150 && ending2(i)-begin2(i)<=250 % ending2(i)-begin2(i)>=250
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
output=[];
surrogate=[];
for z=1:(thal-1)
    for j=1:length(onset)
        if onset(j)-500>0 && onset(j)+500<length(thal_env)
            idx_sur=randi([501,(length(env)-500)],1,1);
            surrogate (z,j,1:1001)= (thal_env(z,idx_sur-500:idx_sur+500)-median(thal_env(z,idx_sur-500:idx_sur+500)))./median(thal_env(z,idx_sur-500:idx_sur+500));
            output(z,j,1:1001)= ((thal_env(z,onset(j)-500:onset(j)+500)-median(thal_env(z,onset(j)-500:onset(j))))./median(thal_env(z,onset(j)-500:onset(j))));
        end
    end
end

median_chan=[];
for i=1:size(output,1)
median_chan(i,1:size(output,3))=median(reshape(output(i,1:size(output,2),1:size(output,3)),[size(output,2) , size(output,3)]));
end