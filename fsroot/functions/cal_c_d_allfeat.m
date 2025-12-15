function [N, D, c, d,e, f ] = ...
          cal_c_d_allfeat(minapg, maxapg, minjpg, maxjpg, z_apg, z_jpg, Jpg, T2_5)
    %disp("cal_c_d_allfeat is running now");
    % Initialize outputs
    [N, D, c, d, e, f] = deal(0);
    

    if length(maxapg) < 2 || length(minjpg) < 2
        warning('Not enough APG extrema for reliable c/d detection.');
        return
    end

    if length(z_jpg) < 4
        warning('Not enough JPG zero crossings for reliable c/d detection.');
        return
    end
    
 
    % --- First Condition: Jpg(minjpg(2)) < 0
    if Jpg(minjpg(2)) < 0
        if maxjpg(1) > minjpg(1)
            c = maxjpg(1);
            d = z_apg(2);
        elseif length(maxjpg) >= 2 && maxjpg(2) > minjpg(1)
            c = maxjpg(2);
            d = z_apg(2);
        elseif length(maxjpg) >= 3 && maxjpg(3) > minjpg(1)
            c = maxjpg(3);
            d = z_apg(2);
        else
            warning('No suitable maxjpg found for case Jpg(minjpg(2)) < 0');
        end

    % --- Second Condition: Jpg(minjpg(2)) > 0
    elseif Jpg(minjpg(2)) > 0
        if maxapg(2) < 0
            c = maxapg(2);
            d = minapg(2);
        elseif maxapg(2) > 0
            c = minjpg(2) - T2_5;
            d = minjpg(2) + T2_5;
        else
            warning('Unexpected APG_maxima(2) value.');
        end
    end

    % Assign shared outputs
    e = z_jpg(3);
    f = z_jpg(4);
    N = e;
    D = f;
end
