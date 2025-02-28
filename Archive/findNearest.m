function min_idx = findNearest(t_ref, t_compare)
    min_diff = 1;
    min_idx = 0;
    for idx=1:length(t_compare)
        curr_t = t_compare(idx);
        if abs(curr_t-t_ref) < min_diff
            min_diff = abs(curr_t-t_ref);
            min_idx = idx;
        end
    end
end