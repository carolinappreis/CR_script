x_all{iii,1}=[pc_trials_ns; pc_trials];
pc1_exp{iii,1}(1:3,:)=[(explained_ns(:,1:3))' (explained_rs(:,1:3))'];


st=out.start_c{iii,co};
et=out.ending_c{iii,co};
yyt=out.yy{iii,co};

for ax=1:3
    tremor_or2=NaN(1,length(st));
    tremor_or3=NaN(1,length(st));
    tremor_or4=NaN(1,length(st));
    tremor_f2=NaN(20,5001);
    tremor_f22=NaN(20,5001);
    z_seg1=NaN(length(st),5000);
    e_env1=NaN(length(st),5000);
    env_5=NaN(1,length(st));
    
    
    for i=1:length(st)
        if (~isnan(st(i)))
            tremor_or2(1,i)=(mean(envelope(ax,et(i)-1000:et(i)))-mean(envelope(ax,st(i)-1000:st(i))))/mean(envelope(ax,st(i)-1000:st(i)));
            tremor_or3(1,i)=mean(envelope(ax,st(i)-1000:st(i)));
            tremor_or4(1,i)=mean(env_acc(ax,et(i)-1000:et(i)));
            
            env_5(1,i)=mean(((envelope(ax,st(i):st(i)+5000-1))./(max(envelope(ax,st(i):st(i)+5000-1)))));
            
            tremor_f2(i,1:(et(i)-st(i)+1))=unwrap(phase(ax,st(i):et(i)));
            tremor_f22(i,1:(et(i)-st(i)+1))=(phase(ax,st(i))+(0:1:(et(i)-st(i)))*2*pi/(1000./mean(freqi(ax,st(i)-1000:st(i)))));
            tremor_k(i,1)= (tremor_f2(i,(et(i)-st(i)+1))-tremor_f22(i,(et(i)-st(i)+1)))/(2*pi*0.001*(et(i)-st(i))); %mean(frequency(et(i)-1000:et(i)));%
            
            z_seg1(i,1:5000)=z_sig(m_ax,st(i):st(i)+5000-1);
            e_env1(i,1:5000)=envelope(m_ax,st(i):st(i)+5000-1);
            e_env1(i,1:5000)=zenv(m_ax,st(i):st(i)+5000-1);
            
            
            
        else
            env_5(1,i)=NaN;
            tremor_or2(1,i)=NaN;
            tremor_or3(1,i)=NaN;
            tremor_or4(1,i)=NaN;
            
            tremor_f22(i,1:5001)=NaN;
            tremor_f2(i,1:5001)=NaN;
            tremor_k(i,1)=NaN;
            
            
            z_seg1(i,1:5000)=NaN;
            e_env1(i,1:5000)=NaN;
            yyt(i)= NaN;
        end
    end
    
    out.all_amp{iii,co}=tremor_or3;
    out.env5_time{iii,1}(ax,:)=env_5;
    
    amp_1=NaN(2,round(size(tremor_or3,2)./2));
    ch_a1=NaN(2,round(size(tremor_or3,2)./2));
    pha_idx=NaN(2,round(size(tremor_or3,2)./2));
    m=1;
    n=1;
    for i=1:length(st)
        if tremor_or3(i)<=nanmedian(tremor_or3)
            amp_1(1,n)= tremor_or3(i);
            ch_a1(1,n)= tremor_or2(i);
            pha_idx(1,n)=yyt(i);
            n=n+1;
            
        else
            amp_1(2,m)= tremor_or3(i);
            ch_a1(2,m)= tremor_or2(i);
            pha_idx(2,m)=yyt(i);
            m=m+1;
            
        end
    end
    
    env5_p=NaN(20,12);
    
    tt=NaN(20,12);
    ee=NaN(20,12);
    
    tt1=NaN(20,12);
    amp1=NaN(10,12);
    tt2=NaN(20,12);
    amp2=NaN(10,12);
    ttf=NaN(20,12);
    
    for i=1:12
        
        env5_p(1:sum(yyt==i),i)=env_5(1,find(yyt==i));
        tt(1:sum(yyt==i),i)=tremor_or2(1,find(yyt==i));
        ee(1:sum(yyt==i),i)=tremor_or4(1,find(yyt==i));
        tt(tt==0)=NaN;
        ttf(1:sum(yyt==i),i)=tremor_k(find(yyt==i));
        
        out.z_seg{iii,i}= z_seg1(find(yyt==i),:);
        out.env_seg{iii,i}= e_env1(find(yyt==i),:);
        
        
        tt1(1:sum(pha_idx(1,:)==i),i)=ch_a1(1,find(pha_idx(1,:)==i));
        amp1(1:sum(pha_idx(1,:)==i),i)=amp_1(1,find(pha_idx(1,:)==i));
        tt2(1:sum(pha_idx(2,:)==i),i)=ch_a1(2,find(pha_idx(2,:)==i));
        amp2(1:sum(pha_idx(2,:)==i),i)=amp_1(2,find(pha_idx(2,:)==i));
    end
    
    out.env5_p{iii,co}{ax,1}=env5_p;
    
    out.arc1{iii,co}{ax,1}=tt1;
    out.arc2{iii,co}{ax,1}=tt2;
    out.change_c{iii,co}{ax,1}=tt;
    out.end_env{iii,co}{ax,1}=ee;
    out.fchange{iii,co}{ax,1}=ttf;
    
    clear tremor_or2 tremor_or3 tt tt1 tt2 ttf tremor_f2 tremor_f22 env5_p env_5
end

else
    
    sig=[2 4 5 7 10];
    if ismember(iii,sig)
        load('/Users/Carolina/OneDrive - Nexus365/Periph_tremor_data/hu_hd_pls.mat','st')
        st=st{find(sig==iii),1};
        et=out.ending_c{iii,co};
        z_seg1=NaN(length(st),3);
        
        for i=1:length(st)
            if (~isnan(st(i))) && (st(i)-5000)>0
                ghost=[st(i)+5000:et(i)-5000];
                duma1=st(i)+5000+round(length(ghost)/2);
                z_seg1(i,:)=[mean(zenv(m_ax,st(i):st(i)+5000)) mean(zenv(m_ax,duma1-5000:duma1)) mean(zenv(m_ax,et(i)-5000:et(i)))];
                clear ghost dum dum1
            end
        end
        out.pls3(iii,:)=median(z_seg1);
        out.pls_hu{iii,1}=z_sig(m_ax,h_up{iii,co});
    else
        out.pls3(iii,:)=NaN(1,3);
        out.pls_hu{iii,1}=[];
    end
end