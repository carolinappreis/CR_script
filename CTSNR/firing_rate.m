function [srate]=firing_rate(data_region)
data_a=cell2mat(data_region);
srn=1000;
    for ii=1:size(data_a,1)
        data=data_a(ii,:);

        spkrate_1=[];
        for i =1:srn:(length(data)-srn);
            spkrate_1=[spkrate_1 numel(find(data(i:i+srn)==1))];
        end
        srate(ii,1)=mean(spkrate_1);
        clear data
    end

end