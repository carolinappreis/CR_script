function  [pc_t,out]=stitch_xall(out,iii,co,pc_t)

f_data=out.filt{iii,co};
lf_data=out.l_filt{iii,co};
st=out.b_trials{iii,co};
et=out.e_trials{iii,co};

% [f]=plot3d(lf_data,st,et);

x_signal=f_data;

if co==1
    
    L = floor((5000./out.samplerateold)*out.samplerate);
    for o=1:3
    [v_new,out]= split_vector_length(x_signal(o,out.h_up{iii,co}), L,out,iii,co,o);
    out.ns_sseg{iii,o}=v_new; clear v_new
    end
    
    for j = 1:size(out.ns_sseg{iii,o},1)
        x = [out.ns_sseg{iii,1}(j,:); out.ns_sseg{iii,2}(j,:); out.ns_sseg{iii,3}(j,:)];
        [pc, score, latent, tsquare, explained] = pca(x');
        pc_trials_ns(j, 1:3) = pc(1:3, 1);
        explained_ns(j, 1:3) = explained;
    end
    
    figure(4)
    subplot(2,1,co)
    plot(pc_trials_ns(:,1), 'r.')
    hold on
    plot(pc_trials_ns(:,2), 'b.')
    plot(pc_trials_ns(:,3), 'k.')
    xlim([0, size(pc_trials_ns, 1)])
    title('NS')
    
    
    pc_t{iii,co}=pc_trials_ns;
    
    
else
    
    bg=out.start{iii,co};
    en=out.ending{iii,co};
    time=1:length(x_signal);
    figure;
%     plot(time,x_signal(1,:))
    hold on
    test=[]
    for j = 1:length(bg)
        if (~isnan(bg(j)))
            epoch=bg(j):en(j);
            test=[test epoch];
            x = [x_signal(1,epoch); x_signal(2,epoch); x_signal(3,epoch)];
%             plot(time(1,epoch),x_signal(1,epoch))
%             hold on
            clear epoch 
            [pc, score, latent, tsquare, explained] = pca(x');
            pc_trials(j, 1:3) = pc(1:3, 1);
            explained1(j, 1:3) = explained;
            clear x
        end
    end
    
    figure(4)
    subplot(2,1,co)
    plot(pc_trials(:,1), 'r.')
    hold on
    plot(pc_trials(:,2), 'b.')
    plot(pc_trials(:,3), 'k.')
    xlim([0, size(pc_trials, 1)])
    title('RS')
    
    pc_t{iii,co}=pc_trials;
    
end
end


