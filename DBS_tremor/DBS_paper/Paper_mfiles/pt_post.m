function [ref]=pt_post(dum2, start, ending, tp_s,tp_e)

for i=1:3
    if length(tp_s{i,1})~=2
        ref{i,1}= [start{1,3}(i) tp_s{i,1}(1) tp_e{i,1}(1) ending{1,3}(i)];
    else
        ref{i,1}= [start{1,3}(i) tp_s{i,1}(1) tp_e{i,1}(1) tp_s{i,1}(2) tp_e{i,1}(2) ending{1,3}(i)];
    end
end

% plot(dum2(1,:))
% hold on
% for m=1:3
%     for i=1:length(ref{m,1})
%         xline(ref{m,1}(i))
%     end
% end

for i=1:3
    breaks(i,:)=tp_e{i,1}-tp_s{i,1};
end
end


