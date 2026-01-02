function [index, corr_quality, skew_quality, seg_quality] = CalcBestCycle(start_index, end_index, peak_index, RFs, PPG)
    % How many cycles at the start/end to cut off
    start_buffer = 0;
    end_buffer = 0;

    index = 0;
    corr_sum = 0;
    quality = -100;
    max_corr_quality = 0;
    max_skew_quality = 0;
    seg_len = length(PPG);

    num_cycle = length(start_index);

    mean_cycle = 0;

    if num_cycle > 2
        % Calculate the median beat-to-beat interval
        if (~mod(num_cycle,2))
            median_interval = ceil(peak_index(floor(num_cycle / 2) + 1) - peak_index(floor(num_cycle / 2)));
        else
            median = ceil(num_cycle / 2);
            median_interval = ceil(((peak_index(median) - peak_index(median - 1)) + ...
                                    (peak_index(median + 1) - peak_index(median))) / 2);
        end
        centered_cycles = zeros(num_cycle, median_interval);

        for i = 1:length(start_index)
            % Center cycle to peak with length of mean interval
            i1 = peak_index(i) - ceil(median_interval / 2);
            i2 = peak_index(i) + floor(median_interval / 2) - 1;

            % See of centered cycle fits in segment
            if i1 < 1
                start_buffer = 1;
                % Remove first column
                centered_cycles(1,:) = [];
                num_cycle = num_cycle - 1;
                continue;
            end
            if i2 > seg_len
                end_buffer = 1;
                % Remove last column
                centered_cycles(num_cycle - start_buffer,:) = [];
                num_cycle = num_cycle - 1;
                continue;
            end
            
            centered_cycles(i,:) = PPG(i1:i2);
        end

        % Take the mean of all centered cycles in segment
        mean_cycle = mean(centered_cycles);
    end

    
    
    for i = 1 + start_buffer:num_cycle - end_buffer
        
        % Add a second before and after cycle
        p1 = start_index(i) - RFs;
        p2 = end_index(i) + RFs;

        if p1 < 1
            p1 = 1;
        end
        if p2 > seg_len
            p2 = seg_len;
        end

        % Get correlation ccoefficient between cycle and median cycle
        if num_cycle > 2
            corr_coef = corrcoef(centered_cycles(i, :), mean_cycle);
            corr_quality = corr_coef(1, 2);
            corr_sum = corr_sum + corr_quality;
        else
            corr_quality = 1;
        end

        %Get skewness for cycle plus, upto, 1 seconds before and after
        skew_quality = skewness(PPG(p1:p2));

        if (corr_quality * (skew_quality / 2)) > quality
            index = i;
            quality = (corr_quality * (skew_quality / 2));
            max_corr_quality = corr_quality;
            max_skew_quality = skew_quality;
        end
    end

    corr_quality = max_corr_quality;
    skew_quality = max_skew_quality;
    
    seg_quality = corr_sum / (num_cycle - start_buffer - end_buffer);
    return;
end