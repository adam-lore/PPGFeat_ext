function [index, quality] = CalcBestCycle(start_index, end_index, peak_index, RFs, PPG)
    % How many cycles at the start/end to cut off
    start_buffer = 0;
    end_buffer = 0;

    index = 0;
    quality = 0;
    seg_len = length(PPG);

    num_cycle = length(start_index);

    mean_cycle = 0;

    if num_cycle > 1
        % Calculate the mean beat-to-beat interval
        mean_interval = ceil((peak_index(num_cycle) - peak_index(1)) / (num_cycle - 1));

        centered_cycles = zeros(num_cycle, mean_interval);

        for i = 1:length(start_index)
            % Center cycle to peak with length of mean interval
            i1 = peak_index(i) - ceil(mean_interval / 2);
            i2 = peak_index(i) + floor(mean_interval / 2) - 1;

            % See of centered cycle fits in segment
            if i1 < 1
                start_buffer = start_buffer + 1;
                % Remove first column
                centered_cycles(1,:) = [];
                continue;
            end
            if i2 > seg_len
                end_buffer = end_buffer + 1;
                % Remove flast column
                centered_cycles(length(start_index),:) = [];
                continue;
            end
            
            centered_cycles(i,:) = PPG(i1:i2);
        end

        % Take the mean of all centered cycles in segment
        mean_cycle = mean(centered_cycles);
    end

    
    
    for i = 1:length(start_index)
        
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
        if num_cycle > 1
            corr_coef = corrcoef(centered_cycles(i, :), mean_cycle);
            curr_quality = corr_coef(1, 2);
        end

        %{
        %Get skewness for cycle plus, upto, 1 seconds before and after
        curr_quality = skewness(PPG(p1:p2));
        %}

        if curr_quality > quality
            index = i;
            quality = curr_quality;
        end
    end
    return;
end