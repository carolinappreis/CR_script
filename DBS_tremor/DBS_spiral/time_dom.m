function [g1]=time_dom(out,clust,iii)

p=[[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250]];
tt={'NS';'RS'};

for co=1:2
    figure(10+co)
    
    data=out.filt{iii,co};
    
    time=1:length(data(1,:));
    plot3(time(1,out.h_up{iii,co}),data(2,out.h_up{iii,co}),data(3,out.h_up{iii,co}),'Color',[0.5 0.5 0.5])
    %     plot3(data(1,out.h_up{iii,co}),data(2,out.h_up{iii,co}),data(3,out.h_up{iii,co}),'Color',[0.5 0.5 0.5])
    xlabel('time')
    ylabel('y axis')
    zlabel('z axis')
    hold on
    for cl=1:3
        id=[clust.idx{iii,co}{cl,1}];
        if ~isempty(id)
            for ii=1:length(id)
                if co==1
                    dummy=out.h_up{iii,co}(out.start{iii,co}(id(ii),1):out.ending{iii,co}(id(ii),1));
                else
                    dummy=(out.start{iii,co}(1,id(ii)):out.ending{iii,co}(1,id(ii)));
                end
                plot3(time(dummy),data(2,(dummy)),data(3,(dummy)),'.','Color',p(cl,:))
                %                 plot3(data(1,(dummy)),data(2,(dummy)),data(3,(dummy)),'.','Color',p(cl,:))                
                hold on
                clear dummy
            end
            clear id
        end
    end
    view(-30,15)
    title(sprintf(tt{co,1}))
    clear time   data  
end

end

