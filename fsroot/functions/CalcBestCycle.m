function [index, corr_quality, skew_quality, seg_corr_quality, seg_skew_quality, quality] = CalcBestCycle(start_index, end_index, peak_index, RFs, PPG)
    % How many cycles at the start/end to cut off
    start_buffer = 0;
    end_buffer = 0;

    index = NaN;
    corr_sum = 0;
    skew_sum = 0;
    quality = NaN;
    max_corr_quality = NaN;
    max_skew_quality = NaN;
    corr_quality = NaN;
    skew_quality = NaN;
    seg_corr_quality = NaN;
    seg_skew_quality = NaN;

    seg_len = length(PPG);

    num_cycle = length(start_index);

    mean_cycle = 0;

    %ecg_corr = NaN;

    if num_cycle > 2
        %if length(ecg_max) > 1
        %    ecg_corr = length(peak_index) / length(ecg_max);
        %end

        % Calculate the median beat-to-beat interval
        intervals = zeros(length(peak_index) - 1, 1);
        for i = 1:length(peak_index) - 1
            intervals(i) = ceil(peak_index(i + 1) - peak_index(i));
        end
        median_interval = ceil(median(intervals));

        centered_cycles = zeros(num_cycle, median_interval);

        for i = 1:length(start_index)
            % Center cycle to peak with length of mean interval
            i1 = peak_index(i) - ceil(median_interval / 2);
            i2 = peak_index(i) + floor(median_interval / 2) - 1;

            % See of centered cycle fits in segment
            if i1 < 1
                start_buffer = 1;

                if length(centered_cycles) < 1
                    return
                end

                % Remove first column
                centered_cycles(1,:) = [];
                num_cycle = num_cycle - 1;
                continue;
            end
            if i2 > seg_len
                end_buffer = 1;

                if (num_cycle - start_buffer) < 1
                    return
                end

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

        %p1 = start_index(i);
        %p2 = end_index(i);


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
        skew_sum = skew_sum + skew_quality;

        if skew_quality > 1
            skew_quality = 1;
        end

        %if (corr_quality + (0.5* skew_quality)) > quality || isnan(quality)
        if isnan(quality)
            index = i;
            quality = (corr_quality * (skew_quality / 2));
            max_corr_quality = corr_quality;
            max_skew_quality = skew_quality;
        end
    end

    corr_quality = max_corr_quality;
    skew_quality = max_skew_quality;
    
    seg_corr_quality = corr_sum / (num_cycle - start_buffer - end_buffer);
    seg_skew_quality = skew_sum / (num_cycle - start_buffer - end_buffer);

    return;
end