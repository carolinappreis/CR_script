% cd('/Users/Carolina/Documents/GitHub/CR_script/DBS_tremor/DBS_paper')

clear all
close all
% toplot_spirals1
load(strcat('/Users/Carolina/OneDrive - Nexus365/Phasic_DBS/patient data/DBS_DATA/clean_SH_spirals/spiral_4bins.mat'));
%%% from 'organise_bins.m'

ad.metrics{1,1,1}=[ad.metrics{1,1,1}; ad.metrics{1,1,2}];
ad.metrics{2,1,1}=[ad.metrics{2,1,1}; ad.metrics{2,1,2}];
met=cell(3,3);

for i=1:3
    for ii=1:3
        met{i,ii}=ad.metrics{i,ii,1}; %%%% met{pt,cond}(n_spirals,n_quadrants,metrics)
    end
end



% % %ANOVA
me=1;

for ss=2
%     1:3
    for qq=3
%         1:size(met{1,1},2)
        for cc=1:3
            a1{cc,1}= met{ss,cc}(:,qq,me);
        end

        sze=cellfun('size',a1,1);
        dum=find(sze<max(sze));
        for f=1:size(dum,1)
            adding(f,:)=max(sze)-size(a1{dum(f),1},1);
            for i = 1:adding(f)
                a1{dum(f),1}(end+1,1)=NaN;
            end
        end
        
        
%         [P,ANOVATAB,STATS] = anova1(cell2mat(a1')); 
%         close
%         [c,m,h,nms] = multcompare(STATS);
%         close
        join{ss,1}(:,qq)=m(:,1);
        if ~isempty(find(c(:,end)>0.05))
            results{ss,qq}=find(c(:,end)>0.05);
%             results{ss,qq}=c(:,end);

        else
            results{ss,qq}=NaN;
        end
        clear sze dum adding a1 c m h nms P ANOVATAB STATS

    end
    
    figure(ss)
    bar(join{ss,1}')
end

%%TTEST

% me=3;
% % 
% % for ss=1:3
% %     for cc=1:3
% %         for qq=1:2
% %             a1{cc,qq}= [met{ss,cc}(:,qq,me) ; met{ss,cc}(:,qq+2,me)] ;
% %         end
% %     end
% %     sze=cellfun('size',a1,1);
% %     sze=sze(:,1);
% %     dum=find(sze<max(sze));
% %     for f=1:size(dum)
% %         adding(f,:)=max(sze)-size(a1{dum(f),1},1);
% %         for i = 1:adding(f)
% %             a1{dum(f),1}(end+1,1)=NaN;
% %             a1{dum(f),2}(end+1,1)=NaN;
% %         end
% %     end
% %     
% %     for qua=1:2
% %         for i=1:3
% %             mat_test(:,i)=a1{i,qua};
% %         end
% %         
% %         [P,ANOVATAB,STATS] = anova1(mat_test);
% %         close
% %         [c,m,h,nms] = multcompare(STATS);
% %         close
% %         
% %         if ~isempty(find(c(:,end)>0.05))
% %             results{ss,qua}=find(c(:,end)>0.05);
% %             %         results{ss,qua}=c(:,end);
% %             
% %         else
% %             results{ss,qua}=NaN;
% %         end
% %         clear sze dum adding c m h nms P ANOVATAB STATS mat_test
% %     end
% %     
% % end
