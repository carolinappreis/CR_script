gcf=figure(1)

NS2=[];

kk=1;

for i=1:size(NS1,1)

    if NS1(i,1)~=-1 && NS1(i,2)~=-1 && NS1(i,3)~=-1

        NS2(kk,:)=NS1(i,:);

        kk=kk+1;

    end

end

 

set(gcf,'color','white')

xlim([1 1000])

ylim([-300 1000])

set(gca,'xtick',[]);

set(gca,'xcolor',[1 1 1])

set(gca,'ytick',[]);

set(gca,'ycolor',[1 1 1])

time=NS2(:,3);

amplitude1=NS2(:,1);

amplitude2=NS2(:,2);

v = VideoWriter('pres4','MPEG-4');

open(v);

for k = 1:10:2000

    

    plot(amplitude1(1:k),amplitude2(1:k),'r');

    set(gca,'xtick',[]);

    set(gca,'xcolor',[1 1 1])

    set(gca,'ytick',[]);

    set(gca,'ycolor',[1 1 1])

    xlim([1 1000])

    ylim([-300 1000])

    hold on

    frame = getframe(gcf);

    writeVideo(v,frame);

end

 

close(v);