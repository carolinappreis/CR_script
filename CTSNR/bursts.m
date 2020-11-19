
function [onset1,offset1]=bursts(env)

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
    if (ending2(i)-begin2(i))>=50
        ind_b=[ind_b i];
    end
end


begin3=begin2(ind_b);
ending3=ending2(ind_b);

for r=1:length(ending3)-1
    if begin3(r+1)-ending3(r)<space_b
        begin3(1,r+i)=NaN;
        ending3(1,r)=NaN;
    end
end


begin3=begin2(ind_b);
ending3=ending2(ind_b);
for r=1:2:length(ending3)-2
    
    if begin3(r+1)-ending3(r)<space_b
        begin3(1,r+1)=NaN;
        ending3(1,r+1)=NaN;
    end
    
    if begin3(r+2)-ending3(r)<space_b
        begin3(1,r+2)=NaN;
        ending3(1,r+2)=NaN;
    end
    
    if ~isnan(begin3(1,r+1)) && ~isnan(begin3(1,r+2))
        if begin3(r+2)-ending3(r+1)<space_b
            begin3(1,r+2)=NaN;
            ending3(1,r+2)=NaN;
        end
    end
    
end


st=begin3(~isnan(begin3)); clear begin3
ed=ending3(~isnan(ending3)); clear ending3


% time=1:length(env);
% plot(time,env)
% hold on
% plot(time(begin2(ind_b)),env(begin2(ind_b)),'r.')
% plot(time(ending2(ind_b)),env(ending2(ind_b)),'b.')
% plot(time(st),env(st),'go')
% plot(time(ed),env(ed),'bo')
% find((ed(2:end)-st(1:end-1))<200)


if ~isempty (st)
    begin3=st;
    ending3=ed;
    
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

