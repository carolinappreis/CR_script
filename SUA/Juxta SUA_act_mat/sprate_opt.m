spikerate=[];
m=1;
for i =1:srn:(length(data)-srn);
    spikerate(1,m)=numel(find(data(i:i+srn)==1));
    m=m+1;
end



m=1;
for i =1:length(onset1)-1;
    spkr_inb(1,m)=numel(find(data(onset1(i):offset1(i))==1)).*srn./(length(time(onset1(i):offset1(i))));
    spkr_outb(1,m)=numel(find(data(offset1(i):onset1(i+1))==1)).*srn./(length(time(offset1(i):onset1(i+1))));
    m=m+1;
end

% find(spkr_inb==0)
% find(spkr_outb==0)
% 
% mean(nonzeros(spkr_inb))
% mean(nonzeros(spkr_outb))

ISI=diff(time(data==1));

for j=1:length(onset1)
    ISI_inb1{:,j}=diff(time(find(data(onset1(j):offset1(j))==1)));
    output_in{:,j}= (onset1(j)+find(data(onset1(j):offset1(j))==1))-1;
    
    if offset1(j)==offset1(end)
       ISI_outb1{:,j}=diff(time(find(data(offset1(j):end)==1)));
       output_out{:,j}= (offset1(j)+find(data(offset1(j):end)==1))-1;
    else
        ISI_outb1{:,j}=diff(time(find(data(offset1(j):onset1(j+1))==1)));
        output_out{:,j}= (offset1(j)+find(data(offset1(j):onset1(j+1))==1))-1;
    end
end

 for n=1:(length(data)/4)
            idx_sur=randi([thr+1,(length(data)-thr)],1,1);
            output_sur{:,n}=idx_sur+find(data(idx_sur:idx_sur+thr)==1);
end

ISI_inb= cell2mat(ISI_inb1);
ISI_outb=[diff(time(find(data(1:onset1(1))==1))) cell2mat(ISI_outb1)];



% Beta ECoG phases of thalamic spiking 

idx_firingin= cell2mat(output_in);
idx_firingsur= cell2mat(output_sur);
idx_firingout= [find(data(1:onset1(1))==1) cell2mat(output_out)];
[~, pos] = ismember(idx_firingin,idx_firingout);
idx_firingout(nonzeros(pos))=[];



if length(data_ones)==(length(idx_firingin)+length(idx_firingout))
    spkang_inb=ang(idx_firingin); % ang= phases of filtered beta
    spkang_outb=ang(idx_firingout);
    spkang_sur=ang(idx_firingsur);
end


