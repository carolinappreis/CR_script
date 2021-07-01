function[match_ax]=link_ax(spiral)
    if spiral==0
        match_ax=NaN(2,4,3);
        match_ax(1,1:4,1:3)=[repmat([3 2 1],3,1); [1 2 3]];
        match_ax(2,1:4,1:3)=([[3 2 1];[3 2 1];[3 2 1];[1 2 3]]);
    else
        match_ax=NaN(2,4,3);
        match_ax(1,1:4,1:3)=[repmat([3 2 1],3,1); [2 1 3]];
        match_ax(2,1:4,1:3)=([[3 2 1];[3 2 1];[3 2 1];[2 1 3]]);
    end

end
