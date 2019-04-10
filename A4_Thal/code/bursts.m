
function [onset1]=bursts(env)

threshold=prctile(env,75);
indexexceed=find(env>threshold);
diffindex=diff(indexexceed);
pnts=find(diffindex>1);
begin=indexexceed(pnts+1);
ending=indexexceed(pnts);
begin2=[indexexceed(1) begin];
ending2=[ending indexexceed(end)];
space_b=200;

ind_b=[];
for i=1:(length(begin2))-1
    if (ending2(i)-begin2(i))>=50 && begin2(i+1)-ending2(i)>space_b% min duration of bursts
        ind_b=[ind_b i];
    end
end

if ~isempty (ind_b)
    begin3=begin2(ind_b);
    ending3=ending2(ind_b);
    
    duration=ending3-begin3;
    median_b=median(duration);
    SD_b=std(duration);
    
   
    bs=[];
    bl=[];
    for i=2:(length(begin3)-1)
        if (ending3(i)-begin3(i))<=median_b     % && (begin3(i)-ending3(i-1))>=space_b && (begin3(i+1)-ending3(i))>=space_b  min duration of bursts
            bs=[bs i];
        elseif (ending3(i+1)-begin3(i))>median_b             %&& (begin3(i)-ending3(i-1))>=space_b && (begin3(i+1)-ending3(i))>=space_b
            bl=[bl i];
        end
    end
    
    
    if (ending3(1)-begin3(1))<=median_b    %&& (begin3(2)-ending3(1))>=space_b
        bs= [1 bs];
    elseif (ending3(1)-begin3(1))>median_b    %&& (begin3(2)-ending3(1))>=space_b
        bl= [1 bl];
    end
    
    if (ending3(end)-begin3(end))<=median_b   %&& (begin3(end-1)-ending3(end))>=space_b
        bs= [bs length(begin3)];
    elseif (ending3(end)-begin3(end))>median_b    %&& (begin3(end-1)-ending3(end))>=space_b
        bl= [bl length(begin3)];
    end
    
    b_short=[begin3(bs);ending3(bs)];
    b_long=[begin3(bl);ending3(bl)];
    
    onset1{1,1}=b_short(1,:);
    onset1{2,1}=b_long(1,:);
    offset1{1,1}=b_short(2,:);
    offset1{2,1}=b_long(2,:);
end

