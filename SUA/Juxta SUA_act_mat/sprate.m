spikerate=[];
m=1;
for i =1:srn:(length(data)-srn);
    spikerate(1,m)=numel(find(data(i:i+srn)==1));
    m=m+1;
end

epoch=srn/5;
nspike_bin=[];
m=1;
for i =1:epoch:(length(data)-epoch);
    nspike_bin(1,m)=numel(find(data(i:i+epoch)==1));
    m=m+1;
end

m=1;
for i =1:length(onset1);
    spikes_inb(1,m)=numel(find(data(onset1(i):onset1(i)+epoch)==1));
    if onset1(i)>epoch
        spikes_outb(1,m)=numel(find(data(onset1(i)-epoch:onset1(i))==1));
        m=m+1;
    end
end


ISI=diff(time(data==1));

ISI_inb1=cell(1,1);
for i=1:length(onset1)
    if length(time(find(data(onset1(i):onset1(i)+epoch)==1)))>1
        ISI_inb1{:,i}=diff(time(find(data(onset1(i):onset1(i)+epoch)==1)));
    end
end

ISI_outb1=cell(1,1);
for i=1:length(onset1)
    if onset1(i)>epoch && length(time(find(data(onset1(i)-epoch:onset1(i))==1)))>1
        ISI_outb1{:,i}=diff(time(find(data(onset1(i)-epoch:onset1(i))==1)));
    end
end

ISI_inb=cell2mat(ISI_inb1);
ISI_outb=cell2mat(ISI_outb1);


% phases thalamic phase locking to Beta ECoG

output_in=cell(1,1);
output_out=cell(1,1);
for j=1:length(onset1)
    if onset1(j)>epoch
    output_in{:,j}= onset1(j)+find(data(onset1(j):onset1(j)+epoch)==1);
    output_out{:,j}= onset1(j)+find(data(onset1(j)-epoch:onset1(j))==1);
    end
end


idx_firingin= cell2mat(output_in);
idx_firingout= cell2mat(output_out);

spkang_inb=ang(idx_firingin);
spkang_outb=ang(idx_firingout);

