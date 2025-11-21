function [index, quality] = CalcBestCycle(start_index, end_index, RFs, PPG)
    index = 0;
    quality = 0;
    seg_len = length(PPG);
    
    for i = 1:length(start_index)
        
        p1 = start_index(i) - RFs;
        p2 = end_index(i) + RFs;

        if p1 < 1
            p1 = 1;
        end
        if p2 > seg_len
            p2 = seg_len;
        end
        %Get skewness for cycle plus, upto, 1 seconds before and after
        curr_quality = skewness(PPG(p1:p2));
        if curr_quality > quality
            index = i;
            quality = curr_quality;
        end
    end
    return;
end