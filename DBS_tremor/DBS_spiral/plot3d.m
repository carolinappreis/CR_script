function [f]=plot3d(lf_data,st,et)



timeor=1:size(lf_data,2);

for f=1:length(st)
    figure(1)
    subplot(length(st),1,f)
    epoch=st(f):et(f);
    plot3(lf_data(1,epoch),lf_data(2,epoch), lf_data(3,epoch),'color',rand(1,3));
    hold on
    n(f,:)=et(f)-st(f);
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    for ii=1:3
    figure(2)
    subplot(3,1,ii)
    plot(timeor(1,epoch),lf_data(ii,epoch));
    hold on
    end
   
    figure(3)
    plot3(timeor(1,epoch),lf_data(2,epoch), lf_data(3,epoch),'color',rand(1,3));
    xlabel('time')
    ylabel('y axis')
    zlabel('z axis')
    clear epoch
end


%%%%%%%%%%%%%%% animation drawing spiral in time
% figure()
% xlabel('x axis')
% ylabel('y axis')
% zlabel('z axis')
% for f=1:length(st)
%     epoch=st(f):et(f);
%     h = animatedline;
%     numpoints = length(epoch);
%     % x = tremorx(1,epoch);
%     % y = tremory(1,epoch);
%     % z= tremorz(1,epoch);
%     x = lf_data(1,epoch);
%     y = lf_data(2,epoch);
%     z= lf_data(3,epoch);
%     a = tic; % start timer
%     for k = 1:numpoints
%         addpoints(h,x(k),y(k),z(k));
%         b = toc(a); % check timer
%         if b > (1/10)
%             drawnow % update screen every 1/30 seconds
%             a = tic; % reset timer after updating
%             k;
%         end
%     end
%     h = animatedline(x,y,z,'color',rand(1,3));
%     hold on
% end

end
