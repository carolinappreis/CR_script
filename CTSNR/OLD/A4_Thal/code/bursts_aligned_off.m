
function [offset]=bursts_aligned(env,Ecogfiltered,phase)
samprate=1000;
threshold=prctile(env,75);
tt(size(env,1):size(env,2))=threshold;
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
    if (ending2(i)-begin2(i))>=50 && begin2(i+1)-ending2(i)>200% min duration of bursts
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
    
    %-----
    
    [maxvalM,maxidxM] = findpeaks(Ecogfiltered);
    
    for hh=1:size(offset1,1)
        for b = 1:length(offset1{hh,1})
            for p=1:length(maxidxM)
                if min(abs(offset1{hh,1}(b)-maxidxM(p)))<=30;
                    pre_offset{hh,1}(b,:)=p;
                end
            end
        end
        offset{hh,1}=maxidxM(nonzeros(pre_offset{hh,1}));
    end
    
    for hh=1:size(onset1,1)
        for b = 1:length(onset1{hh,1})
            for p=1:length(maxidxM)
                if min(abs(onset1{hh,1}(b)-maxidxM(p)))<=30;
                    pre_onset{hh,1}(b,:)=p;
                end
            end
        end
        onset{hh,1}=maxidxM(nonzeros(pre_onset{hh,1}));
    end
    % -----------------------------------
    %     onset_all=horzcat(onset{:});
    
    %     for i=1:length(onset_all)
    %     bursts_singular(i,:)= phase(1,onset_all(i)-200:onset_all(i)+200);
    %     end
    %
%     samprate=1000;
%     time=0:1/samprate:(size(Ecogfiltered,2)-1)/samprate;
%     plot(time,env,'LineWidth',1)
%     hold on
%     plot(time,Ecogfiltered,'LineWidth',1,'Color',[0.5 0.5 0.5])
%     plot(time,tt,'LineWidth',2,'Color',[0.6 0 0])
%     plot(time(onset1{2,1}),env(onset1{2,1}),'b.','MarkerSize', 20)
%     plot(time(onset{2,1}),Ecogfiltered(onset{2,1}),'ko')
%     
    
    
%         plot(time(n1(5):n2(5)),env(n1(5):n2(5)),'LineWidth',2,'Color',
%     n1=double(onset1{2,1});
%     n2=double(offset1{2,1});
%     
%     n3=double(onset1{1,1});
%     n4=double(offset1{1,1});
%     
%     m1=double(onset{2,1});
%     m2=double(onset{1,1});
%         plot(time,env,'LineWidth',1,'Color',[0.5 0.5 0.5])
%         hold on
%         plot(time,Ecogfiltered,'LineWidth',1,'Color',[0.5 0.5 0.5])
%         plot(time(n1(5):n2(5)),env(n1(5):n2(5)),'LineWidth',2,'Color',[0 0 0])
%         plot(time(n3(5):n4(5)),env(n3(5):n4(5)),'LineWidth',2,'Color',[0 0 0])
%         plot(time(n1(5):n2(5)),Ecogfiltered(n1(5):n2(5)),'LineWidth',2,'Color',[0 0 0])
%         plot(time(n3(5):n4(5)),Ecogfiltered(n3(5):n4(5)),'LineWidth',2,'Color',[0 0 0])
%     
%         plot(time,tt,'LineWidth',2,'Color',[0.6 0 0])
%         plot(time([n1 n3]),env([n1 n3]),'b.','MarkerSize', 20)
%         plot(time([m1 m2]),env([m1 m2]),'bo','MarkerSize', 10)
%         xlim ([time(n1(5))-0.5 time(n1(5))+0.5])
%         xticks([time(n1(5))-0.5 :0.1:time(n1(5))+0.5])
%         xticklabels ({'-500','-400','-300','-200','-100','0','100','200','300','400','500'})
%         box('off')
%     %
%     plot(time,Ecogfiltered,'LineWidth',1.5,'Color','k')
%     xlim ([time(n1(5))-0.5 time(n1(5))+0.5].*1000)
%     xticks([time(n1(5))-0.5 :0.1:time(n1(5))+0.5].*1000))
%     xticklabels ({'-500','-400','-300','-200','-100','0','100','200','300','400','500'})
%     box('off')
    %
    
end

