function [PPG_max_index, PPG_min_index] = FindSegMaxMin(total_seg_idx, PPG_filtered, Fs)
    segment = PPG_filtered(total_seg_idx, :);
    %find maxima and minima of current segment
    PPG_max = islocalmax(segment ,"MinProminence",0.5,"FlatSelection","all",...
        "MinSeparation", ceil(Fs/5));
    %tries to find a peak close after minima
    PPG_min = islocalmin(segment ,"MinProminence",0.5, "ProminenceWindow",[0 ceil(Fs/5)], "FlatSelection","all",...
        "MinSeparation", ceil(Fs/5));
    %exception if it is the last minima
    PPG_last_min = islocalmin(flip(segment(end - ceil(Fs/5):end)), "MaxNumExtrema", 1);

    PPG_max_index = find(PPG_max);
    PPG_min_index = find(PPG_min);
    PPG_last_min_index = find(flip(PPG_last_min)) + (length(segment) - ceil(Fs/5));

    if ~isempty(PPG_last_min_index) && ~isempty(PPG_min_index) && PPG_last_min_index(1) ~= PPG_min_index(end)
        PPG_min_index(end + 1) = PPG_last_min_index(1);
    end
end