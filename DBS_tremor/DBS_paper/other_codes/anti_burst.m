function [ab_2,norm_bd,begin2,nr]=anti_burst(env,thr)

idx=find(env <thr);
% %                         figure(iii+1)
%                         plot(data_e(1,:))
%                         hold on
%                         xline(start(1))
%                         xline(ending(1))
%                         yline(thr)
d_idx=diff(idx);
pnts=find(d_idx>1);
begin=idx(pnts+1);
ending=idx(pnts);

begin2=[idx(1) begin];
ending2=[ending idx(end)];

                                
if length(begin2)==1
    ab_2=0;
    norm_bd=ending2-1;
    nr=1/length(env);
else
    r=0;
    for ii=1:length(begin2)
        if  ii+1 <length (begin2)
            ab_2(ii,:)=begin2(ii+1)-ending2(ii);
            %                 plot(time(ending2(ii):begin2(ii+1)),env(ending2(ii):begin2(ii+1)),'k.')
            r=r+1;
        else
            ab_2(ii,:)=NaN;
        end
    end

    norm_bd=(ending2-begin2);
    nr=r./length(env);

end
end


% %                    time=1:length(env);
% %                    plot(time,env)
% %                    hold on
% %                    plot(time(begin2),env(begin2),'ro')
% %                    plot(time(ending2),env(ending2),'bo')
% %                    yline(thr)
% %                    
% %                    [m,n]=butter(1,[1.5/(0.5*d.samplerateold)],'low');
% %                    for i=1:3
% %                    lf_data(i,:)=filtfilt(m,n,d.data_ds(i,:));
% %                    end
% % 
% % 
% %                    lf_data=s.filt{3,1};
% %                    figure(1)
% %                    time=1:length(env);
% %                    plot(time,env)
% %                    
% %                    plot3(time,lf_data(1,handup), lf_data(3,handup));
% % 
% %                    x=26400
% %                    figure(2)
% %                    plot3(time(1:x),lf_data(1,handup(1:x)), lf_data(2,handup(1:x)));
% %                    hold on
% %                    for i=1:1
% %                    epoch=begin2(i):ending2(i);
% % %                    figure(1)
% % %                    plot(time(epoch),env(epoch),'r.')
% %                    figure(2)
% %                    plot3(time(1,epoch),lf_data(1,epoch), lf_data(2,epoch),'r.');
% %                    hold on
% %                    end
% %                   
% %                    

    