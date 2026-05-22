function [start_index, end_index, peak_index] = FindCycles(total_seg_idx, PPG_filtered, Fs)
    [index_PPG_max, index_PPG_min] = FindSegMaxMin(total_seg_idx, PPG_filtered, Fs);

    % Set start, end and peak index arrays to zero of size index_PPG_min
    [start_index, end_index, peak_index] = deal(zeros(size(index_PPG_min)));

    threshold = 40 * Fs/100;

    cycle_index = 0;

    % Find two minima that are sufficiently far apart and has a peak in between
    % Loop over all minima indices
    for i = 1:length(index_PPG_min)-1
        
        % Check if the next minimum index is greater than the current one by
        % the threshold value
        if index_PPG_min(i+1) > index_PPG_min(i) + threshold
            % If yes, check if there is any maximum index between the two minima
            maxima_between_minima = index_PPG_max(index_PPG_max > index_PPG_min(i) & index_PPG_max < index_PPG_min(i+1));
            next_peak = index_PPG_max(find(index_PPG_max > index_PPG_min(i+1), 1));
            if ~isempty(maxima_between_minima) & ~isempty(next_peak)
                cycle_index = cycle_index + 1;

                start_index(cycle_index) = index_PPG_min(i);
                end_index(cycle_index) = index_PPG_min(i+1);

                peak_index(cycle_index) = maxima_between_minima(1);
                peak_value = PPG_filtered(maxima_between_minima(1));

                % Find the highest peak in the cycle
                if length(maxima_between_minima) > 1
                    for j = 2:length(maxima_between_minima)
                        if PPG_filtered(maxima_between_minima(j)) > peak_value
                            peak_index(cycle_index) = maxima_between_minima(j);
                            peak_value = PPG_filtered(maxima_between_minima(j));
                        end
                    end
                end
            end
        end
    end

    start_index = resize(start_index, cycle_index);
    end_index = resize(end_index, cycle_index); 
    peak_index = resize(peak_index, cycle_index);
end