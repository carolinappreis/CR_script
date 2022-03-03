function [f1]=spirals_cond(tremor_ds,samplerate,start,ending)

[b,a]=butter(2,[0.1/(2*samplerate) 1/(2*samplerate)],'bandpass');
% [b,a]=butter(2,[4/(2*samplerate) 8/(2*samplerate)],'bandpass');

for ix=1:size(tremor_ds,1)
    to_pc(ix,:)=filtfilt(b,a,tremor_ds(ix,:));
end

dx = [to_pc(1,:); to_pc(2,:);to_pc(3,:)];
[pc, score, latent, tsquare, explained] = pca(dx');
pc_sig= to_pc(1,:).*pc(1,1)+ to_pc(2,:).*pc(2,1)+ to_pc(3,:).*pc(3,1);



for pp=1:size(start,2)
    segment=start(pp):ending(pp);
%     plot(pc_sig(segment))
end

handup = [];
for ix = 1:length(start)
    handup = [handup start(ix):ending(ix)]; %#ok<*AGROW>
end
handup = sort(handup,'ascend');

time=1:length(pc_sig);
f1=figure(5);
subplot(2,1,1)
plot(time,pc_sig)
hold on
plot(time(start),pc_sig(start),'k.')
plot(time(ending),pc_sig(ending),'ro')
subplot(2,1,2)
plot(time,tremor_ds(1,:))
hold on
plot(time(start),tremor_ds(1,start),'k.')
plot(time(ending),tremor_ds(1,ending),'ro')
close
        
end
