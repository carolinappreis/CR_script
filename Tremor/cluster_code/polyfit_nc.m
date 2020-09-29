function[out]=polyfit_nc(out)
win_m_ns=zeros(10,100,2);
for iii=1:10
    m_ax=1;
    y1 = squeeze(out.ns_arc{iii,m_ax});
%     parpool(10)

    parfor uu=1:100
       
        y=y1(uu,:); x = 1:12;
        y=y';x=x';
        
        [win]=aid_fx(iii,uu,x,y)
        win_m_ns(iii,uu,:)=win;
    end
 figure(1)
    subplot(2,5,iii)
    histogram(win_m_ns(iii,:,1))
end

out.win_ns=win_m_ns;
end

% for iii=1:10
%    rat(iii,:)=[numel(find( win_m_ns(iii,:,1)~=1))];
%  r=[(win_m_ns(iii,:,1)); (win_m_ns(iii,:,2))];
% avg(iii,1)=mean(r(2,find(r(1,:)~=1))); clear r
% end


