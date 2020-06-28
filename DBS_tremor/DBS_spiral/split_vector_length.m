function [v_new ,out,e_new] = split_vector_length(v, e, L, out, iii, co)
n = floor(length(v)./ L); % Number of segments
t=1:length(v);
idx_st = NaN(n, 1);
idx_st(1, 1) = 1;

idx_end = NaN(n, 1);
idx_end(1, 1) = L;

v_new = NaN(n, L);
for i = 2:n
    idx_st(i) = idx_st(i-1) + L;
    idx_end(i) = idx_end(i-1) + L;
end

for i = 1:n
    v_new(i, :) = v(1, idx_st(i, 1) : idx_end(i, 1));
    t_new(i, :) = t(1, idx_st(i, 1) : idx_end(i, 1));
    e_new(i,:)= e(1,idx_st(i, 1) : idx_end(i, 1));
    %     figure(1)
    %     plot(t_new(i,:),v_new(i,:))
    %     hold on
    
end


out.start{iii,co}=idx_st;
out.ending{iii,co}=idx_end;
end