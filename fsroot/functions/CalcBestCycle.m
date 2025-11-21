function [index, quality] = CalcBestCycle(start_index, end_index, PPG)
    index = 1;
    quality = skewness(PPG(start_index(1):end_index(1)));
    
    for i = 2:length(start_index)-1
        curr_quality = skewness(PPG(start_index(i):end_index(i)));
        if curr_quality > quality
            index = i;
            quality = curr_quality;
        end
    end
    return;
end