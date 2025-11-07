function [c_d_presence, c, d, e, f, N, D] = ...
    APG_c_d_test(APG, APG_maxima, APG_minima, ...
    Jpg, maxjpg, minjpg, z_apg, z_jpg, T2_5, Fs)

% Initialize all outputs
c_d_presence = 1;
[c, d, e, f, N, D] = deal(0);


% Ensure sufficient indices
if length(APG_maxima) < 3 || length(APG_minima) < 3
    warning('Not enough APG extrema for reliable c/d detection.');
    return;
end

% Case 1: c and d are NOT present (positive c)
if APG(APG_maxima(2)) > 0
    c_d_presence = 0;
    disp("c and d is NOT present in APG");

    % Compute fallback using cal_c_d_allfeat
    [N, D, c, d,e, f ] = ...
        cal_c_d_allfeat(APG_minima, APG_maxima, ...
        minjpg, maxjpg, z_apg, z_jpg, Jpg, T2_5);

    return;
end

% Case 2: c and d ARE present and close enough
if APG(APG_maxima(2)) < 0 && ...
        (APG_maxima(2) - APG_minima(2)) < ceil(Fs / 10)

    c_d_presence = 1;
    disp("c and d is present in APG");

    c = APG_maxima(2);
    d = APG_minima(2);
    e = APG_maxima(3);
    f = APG_minima(3);
    N = e;
    D = f;

    return;
else
    % Case 3: Fallback — use cal_c_d_allfeat
    [N, D, c, d,e, f ] = ...
        cal_c_d_allfeat(APG_minima, APG_maxima, ...
        minjpg, maxjpg, z_apg, z_jpg, Jpg, T2_5);
end



end
