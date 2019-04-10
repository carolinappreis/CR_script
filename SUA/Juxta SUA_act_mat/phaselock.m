

output=cell(1,1);
for j=1:length(onset1)
    output{j,:}= onset1(j)+find(data(onset1(j):onset1(j)+200)==1);
end

% output_c=[];
% m=1;
% for j=1:size(output,1)
%     output_c(1,m)=size(output{j,:},2)
%      m=m+1;
% end
% sum(output_c)


angles=cell(1,1)
for j=1:size(output,1)
    angles{:,j}=ang(output{j}(1,:));
end

spike_angles=cell2mat(angles);

% output_a=[];
% m=1;
% for j=1:size(angles,2)
%     output_a(1,m)=size(angles{1,j},2)
%      m=m+1;
% end
% sum(output_a)
