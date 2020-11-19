function [env_var]=burst_var(ctx)

for i=1:size(ctx,1)
    env=abs(hilbert(ctx(i,:)));
    [onset1,offset1]=bursts(env);
    onset=horzcat(onset1{:}); offset=horzcat(offset1{:});
    figure(1)
    subplot(size(ctx,1),1,i)
    time=1:length(env);
    for ii=1:length(onset)
        if onset(ii)>150
            change(ii,1)=(median(env(onset(ii)-150:onset(ii)))-median(env(onset(ii):offset(ii))))./(median(env(onset(ii):offset(ii))));
            plot(time,env,'Color', [0.5 0.5 0.5])
            hold on
            plot(time(onset(ii):offset(ii)),env(onset(ii):offset(ii)),'r.')
            box('off')
        end
    end
    xlim([0 10000])
    
    figure(2)
    subplot(1,size(ctx,1),i)
    histogram(change)
    box('off')
    env_var(i,1)=var(env);
    clear env
end
figure(3)
bar(env_var)
box('off')
ylabel('envelope variance')
xlabel('patients')
end
