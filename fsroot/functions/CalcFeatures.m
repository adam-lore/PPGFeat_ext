function [feature] = CalcFeatures(OnSpDnDpOff, uxvw, next_u, abcdef, OnSpDnDpOff_time, uxvw_time, abcdef_time, ppg, vpg, apg, sample_freq)%store filtered data
    [val_on, val_sp, val_dn, val_dp, val_off] = deal(OnSpDnDpOff(1), OnSpDnDpOff(2), OnSpDnDpOff(3), OnSpDnDpOff(4), OnSpDnDpOff(5));
    [val_u, val_x, val_v, val_w] = deal(uxvw(1), uxvw(2), uxvw(3), uxvw(4));
    [val_a, val_b, val_c, val_d, val_e, val_f] = deal(abcdef(1), abcdef(2), abcdef(3), abcdef(4), abcdef(5), abcdef(6));

    [time_on, time_sp, time_dn, time_dp, time_off] ...
     = deal(OnSpDnDpOff_time(1), OnSpDnDpOff_time(2), OnSpDnDpOff_time(3), OnSpDnDpOff_time(4), OnSpDnDpOff_time(5));
    [time_u, time_x, time_v, time_w] = deal(uxvw_time(1), uxvw_time(2), uxvw_time(3), uxvw_time(4));
    [time_a, time_b, time_c, time_d, time_e, time_f] ...
     = deal(abcdef_time(1), abcdef_time(2), abcdef_time(3), abcdef_time(4), abcdef_time(5), abcdef_time(6));

    sample_time = 1000/sample_freq;

    feature.fiducial_value(1) = val_on;
    feature.fiducial_value(2) = val_sp;
    feature.fiducial_value(3) = val_dn;
    feature.fiducial_value(4) = val_dp;
    feature.fiducial_value(5) = val_off;
    feature.fiducial_value(6) = val_u;
    feature.fiducial_value(7) = val_x;
    feature.fiducial_value(8) = val_v;
    feature.fiducial_value(9) = val_w;
    feature.fiducial_value(10) = val_a;
    feature.fiducial_value(11) = val_b;
    feature.fiducial_value(12) = val_c;
    feature.fiducial_value(13) = val_d;
    feature.fiducial_value(14) = val_e;
    feature.fiducial_value(15) = val_f;

    feature.fiducial_time(1) = time_on * sample_time;
    feature.fiducial_time(2) = time_sp * sample_time;
    feature.fiducial_time(3) = time_dn * sample_time;
    feature.fiducial_time(4) = time_dp * sample_time;
    feature.fiducial_time(5) = time_off * sample_time;
    feature.fiducial_time(6) = time_u * sample_time;
    feature.fiducial_time(7) = time_x * sample_time;
    feature.fiducial_time(8) = time_v * sample_time;
    feature.fiducial_time(9) = time_w * sample_time;
    feature.fiducial_time(10) = time_a * sample_time;
    feature.fiducial_time(11) = time_b * sample_time;
    feature.fiducial_time(12) = time_c * sample_time;
    feature.fiducial_time(13) = time_d * sample_time;
    feature.fiducial_time(14) = time_e * sample_time;
    feature.fiducial_time(15) = time_f * sample_time;
    
    % ts means time span
    ts_on_sp = (time_sp - time_on) * sample_time;
    ts_on_dn = (time_dn - time_on) * sample_time;
    ts_on_dp = (time_dp - time_on) * sample_time;
    ts_on_u = (time_u - time_on) * sample_time;
    ts_on_v = (time_v - time_on) * sample_time;
    ts_on_a = (time_a - time_on) * sample_time;
    ts_on_b = (time_b - time_on) * sample_time;
    ts_on_c = (time_c - time_on) * sample_time;
    
    % Need to ask about this, if I understood correctly, needs next cycle to calculate
    ts_u_next_u = (next_u - time_u) * sample_time;

    ts_sp_c = (time_c - time_sp) * sample_time;
    ts_sp_d = (time_d - time_sp) * sample_time;
    ts_sp_e = (time_e - time_sp) * sample_time;
    ts_sp_dp = (time_dp - time_sp) * sample_time;

    ts_dn_dp = (time_dp - time_dn) * sample_time;

    % Need to ask about this, if I understood correctly, needs next cycle to calculate
    %Tm_bb2 = 0;

    ts_b_c = (time_c - time_b) * sample_time;
    ts_b_d = (time_d - time_b) * sample_time;

    ts_u_sp = (time_sp - time_u) * sample_time;
    ts_u_w = (time_w - time_u) * sample_time;
    ts_u_b = (time_b - time_u) * sample_time;
    ts_u_c = (time_c - time_u) * sample_time;
    ts_u_d = (time_d - time_u) * sample_time;
    
    ts_a_c = (time_c - time_a) * sample_time;

    feature.timespan(1) = ts_on_sp;
    feature.timespan(2) = ts_on_dn;
    feature.timespan(3) = ts_on_dp;
    feature.timespan(4) = ts_on_u;
    feature.timespan(5) = ts_on_v;
    feature.timespan(6) = ts_on_a;
    feature.timespan(7) = ts_on_b;
    feature.timespan(8) = ts_on_c;
    feature.timespan(9) = ts_u_next_u;
    feature.timespan(10) = ts_sp_c;
    feature.timespan(11) = ts_sp_d;
    feature.timespan(12) = ts_sp_e;
    feature.timespan(13) = ts_sp_dp;
    feature.timespan(14) = ts_dn_dp;
    %feature.timespan(15) = Tm_bb2;
    feature.timespan(15) = ts_b_c;
    feature.timespan(16) = ts_b_d;
    feature.timespan(17) = ts_u_sp;
    feature.timespan(18) = ts_u_w;
    feature.timespan(19) = ts_u_b;
    feature.timespan(20) = ts_u_c;
    feature.timespan(21) = ts_u_d;
    feature.timespan(22) = ts_a_c;

    % am means amplitude
    am_on_sp = val_sp - val_on;
    am_on_dn = val_dn - val_on;
    am_on_dp = val_dp - val_on;

    am_on_u = nanCalculate(time_u, @(arr, in) (arr(in(1)) - in(2)), ppg, [time_u val_on]);
    %am_on_u = ppg(time_u) - val_on;
    am_on_v = nanCalculate(time_v, @(arr, in) (arr(in(1)) - in(2)), ppg, [time_v val_on]);
    %am_on_v = ppg(time_v) - val_on;

    am_on_a = nanCalculate(time_a, @(arr, in) (arr(in(1)) - in(2)), ppg, [time_a val_on]);
    %am_on_a = ppg(time_a) - val_on;
    am_on_b = nanCalculate(time_b, @(arr, in) (arr(in(1)) - in(2)), ppg, [time_b val_on]);
    %am_on_b = ppg(time_b) - val_on;
    am_on_c = nanCalculate(time_c, @(arr, in) (arr(in(1)) - in(2)), ppg, [time_c val_on]);
    %am_on_c = ppg(time_c) - val_on;

    am_on_off = time_off - time_on;

    am_dn_sp = val_sp - val_dn;

    % ar mean amplitude ratio
    ar_on_dn__on_sp = am_on_dn / am_on_sp;
    ar_on_dp__on_sp = am_on_dp / am_on_sp;
    ar_dn_sp__on_sp = am_dn_sp / am_on_sp;
    ar_dp_sp__on_sp = (val_sp - val_dp) / am_on_sp;

    feature.amplitude(1) = am_on_sp;
    feature.amplitude(2) = am_on_dn;
    feature.amplitude(3) = am_on_dp;
    feature.amplitude(4) = am_on_u;
    feature.amplitude(5) = am_on_v;
    feature.amplitude(6) = am_on_a;
    feature.amplitude(7) = am_on_b;
    feature.amplitude(8) = am_on_c;
    feature.amplitude(9) = am_on_off;
    feature.amplitude(10) = am_dn_sp;
    feature.amplitude(11) = ar_on_dn__on_sp;
    feature.amplitude(12) = ar_on_dp__on_sp;
    feature.amplitude(13) = ar_dn_sp__on_sp;
    feature.amplitude(14) = ar_dp_sp__on_sp;


    val_vpg_c = nanCalculate(time_c, @(arr, in) (arr(in(1))), vpg, time_c);
    %val_vpg_c = vpg(time_c);
    val_vpg_d = nanCalculate(time_d, @(arr, in) (arr(in(1))), vpg, time_d);
    %val_vpg_d = vpg(time_d);

    % r means ratio
    r_w_u = val_w / val_u;
    r_v_u = val_v / val_u;
    r_val_vpg_c__u = val_vpg_c / val_u;
    r_val_vpg_d__u = val_vpg_d / val_u;

    r_b_a = val_b / val_a;
    r_c_a = val_c / val_a;
    r_d_a = val_d / val_a;
    r_e_a = val_e / val_a;
    r_bcde_a = (val_b - val_c - val_d - val_e) / val_a;
    r_bcd_a = (val_b - val_c - val_d) / val_a;

    feature.vpg_apg(1) = val_vpg_c;
    feature.vpg_apg(2) = val_vpg_d;
    feature.vpg_apg(3) = r_w_u;
    feature.vpg_apg(4) = r_v_u;
    feature.vpg_apg(5) = r_val_vpg_c__u;
    feature.vpg_apg(6) = r_val_vpg_d__u;
    feature.vpg_apg(7) = r_b_a;
    feature.vpg_apg(8) = r_c_a;
    feature.vpg_apg(9) = r_d_a;
    feature.vpg_apg(10) = r_e_a;
    feature.vpg_apg(11) = r_bcde_a;
    feature.vpg_apg(12) = r_bcd_a;

    % wa means waveform area
    wa_on_off = nanCalculate([time_on time_off], @(arr, in) (trapz(arr(in(1) : in(2)))), ppg, [time_on time_off]);
    %wa_on_off = trapz(ppg(time_on : time_off));
    wa_on_sp = nanCalculate([time_on time_sp], @(arr, in) (trapz(arr(in(1) : in(2)))), ppg, [time_on time_sp]);
    %wa_on_sp = trapz(ppg(time_on : time_sp));
    wa_on_c = nanCalculate([time_on time_c], @(arr, in) (trapz(arr(in(1) : in(2)))), ppg, [time_on time_c]);
    %wa_on_c = trapz(ppg(time_on : time_c));
    wa_on_dn = nanCalculate([time_on time_dn], @(arr, in) (trapz(arr(in(1) : in(2)))), ppg, [time_on time_dn]);
    %wa_on_dn = trapz(ppg(time_on : time_dn));

    feature.waveform_area(1) = wa_on_off;
    feature.waveform_area(2) = wa_on_sp;
    feature.waveform_area(3) = wa_on_c;
    feature.waveform_area(4) = wa_on_dn;

    % pa means power area
    pa_on_sp_ppg = nanCalculate([time_on time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), ppg, [time_on time_sp]);
    %pa_on_sp_ppg = sumsqr(ppg(time_on : time_sp));
    pa_u_sp_ppg = nanCalculate([time_u time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), ppg, [time_u time_sp]);
    %pa_u_sp_ppg = sumsqr(ppg(time_u : time_sp));
    pa_sp_c_ppg = nanCalculate([time_sp time_c], @(arr, in) (sumsqr(arr(in(1) : in(2)))), ppg, [time_sp time_c]);
    %pa_sp_c_ppg = sumsqr(ppg(time_sp : time_c));
    pa_sp_d_ppg = nanCalculate([time_sp time_d], @(arr, in) (sumsqr(arr(in(1) : in(2)))), ppg, [time_sp time_d]);
    %pa_sp_d_ppg = sumsqr(ppg(time_sp : time_d));

    pa_on_sp_vpg = nanCalculate([time_on time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), vpg, [time_on time_sp]);
    %pa_on_sp_vpg = sumsqr(vpg(time_on : time_sp));
    pa_u_sp_vpg = nanCalculate([time_u time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), vpg, [time_u time_sp]);
    %pa_u_sp_vpg = sumsqr(vpg(time_u : time_sp));
    pa_sp_c_vpg = nanCalculate([time_sp time_c], @(arr, in) (sumsqr(arr(in(1) : in(2)))), vpg, [time_sp time_c]);
    %pa_sp_c_vpg = sumsqr(vpg(time_sp : time_c));
    pa_sp_d_vpg = nanCalculate([time_sp time_d], @(arr, in) (sumsqr(arr(in(1) : in(2)))), vpg, [time_sp time_d]);
    %pa_sp_d_vpg = sumsqr(vpg(time_sp : time_d));

    pa_on_sp_apg = nanCalculate([time_on time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), apg, [time_on time_sp]);
    %pa_on_sp_apg = sumsqr(apg(time_on : time_sp));
    pa_u_sp_apg = nanCalculate([time_u time_sp], @(arr, in) (sumsqr(arr(in(1) : in(2)))), apg, [time_u time_sp]);
    %pa_u_sp_apg = sumsqr(apg(time_u : time_sp));
    pa_sp_c_apg = nanCalculate([time_sp time_c], @(arr, in) (sumsqr(arr(in(1) : in(2)))), apg, [time_sp time_c]);
    %pa_sp_c_apg = sumsqr(apg(time_sp : time_c));
    pa_sp_d_apg = nanCalculate([time_sp time_d], @(arr, in) (sumsqr(arr(in(1) : in(2)))), apg, [time_sp time_d]);
    %pa_sp_d_apg = sumsqr(apg(time_sp : time_d));

    pa_on_off_ppg = nanCalculate([time_on time_off], @(arr, in) (sumsqr(arr(in(1) : in(2)))), ppg, [time_on time_off]);
    %pa_on_off_ppg = sumsqr(ppg(time_on : time_off));
    pa_on_off_vpg = nanCalculate([time_on time_off], @(arr, in) (sumsqr(arr(in(1) : in(2)))), vpg, [time_on time_off]);
    %pa_on_off_vpg = sumsqr(vpg(time_on : time_off));
    pa_on_off_apg = nanCalculate([time_on time_off], @(arr, in) (sumsqr(arr(in(1) : in(2)))), apg, [time_on time_off]);
    %pa_on_off_apg = sumsqr(apg(time_on : time_off));

    feature.power_area(1) = pa_on_sp_ppg;
    feature.power_area(2) = pa_u_sp_ppg;
    feature.power_area(3) = pa_sp_c_ppg;
    feature.power_area(4) = pa_sp_d_ppg;
    feature.power_area(5) = pa_on_sp_vpg;
    feature.power_area(6) = pa_u_sp_vpg;
    feature.power_area(7) = pa_sp_c_vpg;
    feature.power_area(8) = pa_sp_d_vpg;
    feature.power_area(9) = pa_on_sp_apg;
    feature.power_area(10) = pa_u_sp_apg;
    feature.power_area(11) = pa_sp_c_apg;
    feature.power_area(12) = pa_sp_d_apg;
    feature.power_area(13) = pa_on_off_ppg;
    feature.power_area(14) = pa_on_off_vpg;
    feature.power_area(15) = pa_on_off_apg;

    % r means ratio (again)
    r_ts_on_a__ts_u_next_u = ts_on_a / ts_u_next_u;
    r_ts_on_u__ts_u_next_u = ts_on_u / ts_u_next_u;
    r_ts_on_b__ts_u_next_u = ts_on_b / ts_u_next_u;
    r_ts_on_sp__ts_u_next_u = ts_on_sp / ts_u_next_u;
    r_ts_on_c__ts_u_next_u = ts_on_c / ts_u_next_u;
    r_ts_on_v__ts_u_next_u = ts_on_v / ts_u_next_u;
    r_ts_on_dn__ts_u_next_u = ts_on_dn / ts_u_next_u;
    r_ts_u_w__ts_u_next_u = ts_u_w / ts_u_next_u;
    r_ts_sp_dp__ts_u_next_u = ts_sp_dp / ts_u_next_u;
    %r_Tm_bb2_Tss = 0; % Needs future b (Tm_bb2)

    r_am_on_a__am_on_sp = am_on_a / am_on_sp;
    r_am_on_u__am_on_sp = am_on_u / am_on_sp;
    r_am_on_b__am_on_sp = am_on_b / am_on_sp;
    r_am_on_c__am_on_sp = am_on_c / am_on_sp;
    r_am_on_v__am_on_sp = am_on_v / am_on_sp;
    r_am_on_off__am_on_sp = am_on_off / am_on_sp; 

    wa_dn_off = nanCalculate([time_dn time_off], @(arr, in) (trapz(arr(in(1) : in(2)))), ppg, [time_dn time_off]);
    %wa_dn_off = trapz(ppg(time_dn : time_off));
    r_wa_dn_off__wa_on_dn = wa_dn_off / wa_on_dn;

    r_sp_on = val_sp / val_on;

    r_wa_on_sp__wa_on_off = wa_on_sp / wa_on_off;
    r_wa_on_c__wa_on_off = wa_on_c / wa_on_off;
    r_wa_on_dn__wa_on_off = wa_on_dn / wa_on_off;

    r_pa_on_sp_ppg__pa_on_off_ppg = pa_on_sp_ppg / pa_on_off_ppg;
    r_pa_u_sp_ppg__pa_on_off_ppg = pa_u_sp_ppg / pa_on_off_ppg;
    r_pa_sp_c_ppg__pa_on_off_ppg = pa_sp_c_ppg / pa_on_off_ppg;
    r_pa_sp_d_ppg__pa_on_off_ppg = pa_sp_d_ppg / pa_on_off_ppg;

    r_pa_on_sp_vpg__pa_on_off_vpg = pa_on_sp_vpg / pa_on_off_vpg;
    r_pa_u_sp_vpg__pa_on_off_vpg = pa_u_sp_vpg / pa_on_off_vpg;
    r_pa_sp_c_vpg__pa_on_off_vpg = pa_sp_c_vpg / pa_on_off_vpg;
    r_pa_sp_d_vpg__pa_on_off_vpg = pa_sp_d_vpg / pa_on_off_vpg;

    r_pa_on_sp_apg__pa_on_off_apg = pa_on_sp_apg / pa_on_off_apg;
    r_pa_u_sp_apg__pa_on_off_apg = pa_u_sp_apg / pa_on_off_apg;
    r_pa_sp_c_apg__pa_on_off_apg = pa_sp_c_apg / pa_on_off_apg;
    r_pa_sp_d_apg__pa_on_off_apg = pa_sp_d_apg / pa_on_off_apg;

    feature.ratio(1) = r_ts_on_a__ts_u_next_u;
    feature.ratio(2) = r_ts_on_u__ts_u_next_u;
    feature.ratio(3) = r_ts_on_b__ts_u_next_u;
    feature.ratio(4) = r_ts_on_sp__ts_u_next_u;
    feature.ratio(5) = r_ts_on_c__ts_u_next_u;
    feature.ratio(6) = r_ts_on_v__ts_u_next_u;
    feature.ratio(7) = r_ts_on_dn__ts_u_next_u;
    feature.ratio(8) = r_ts_u_w__ts_u_next_u;
    feature.ratio(9) = r_ts_sp_dp__ts_u_next_u;
    %feature.ratio(10) = r_Tm_bb2_Tss;
    feature.ratio(10) = r_am_on_a__am_on_sp;
    feature.ratio(11) = r_am_on_u__am_on_sp;
    feature.ratio(12) = r_am_on_b__am_on_sp;
    feature.ratio(13) = r_am_on_c__am_on_sp;
    feature.ratio(14) = r_am_on_v__am_on_sp;
    feature.ratio(15) = r_am_on_off__am_on_sp;
    feature.ratio(16) = r_wa_dn_off__wa_on_dn;
    feature.ratio(17) = r_sp_on;
    feature.ratio(18) = r_wa_on_sp__wa_on_off;
    feature.ratio(19) = r_wa_on_c__wa_on_off;
    feature.ratio(20) = r_wa_on_dn__wa_on_off;
    feature.ratio(21) = r_pa_on_sp_ppg__pa_on_off_ppg;
    feature.ratio(22) = r_pa_u_sp_ppg__pa_on_off_ppg;
    feature.ratio(23) = r_pa_sp_c_ppg__pa_on_off_ppg;
    feature.ratio(24) = r_pa_sp_d_ppg__pa_on_off_ppg;
    feature.ratio(25) = r_pa_on_sp_vpg__pa_on_off_vpg;
    feature.ratio(26) = r_pa_u_sp_vpg__pa_on_off_vpg;
    feature.ratio(27) = r_pa_sp_c_vpg__pa_on_off_vpg;
    feature.ratio(28) = r_pa_sp_d_vpg__pa_on_off_vpg;
    feature.ratio(29) = r_pa_on_sp_apg__pa_on_off_apg;
    feature.ratio(30) = r_pa_u_sp_apg__pa_on_off_apg;
    feature.ratio(31) = r_pa_sp_c_apg__pa_on_off_apg;
    feature.ratio(32) = r_pa_sp_d_apg__pa_on_off_apg;

    % s means slope

    s_sp_c_ppg = (nanCalculate(time_c, @(arr,in) (arr(in(1))), ppg, time_c) - val_sp) / ((time_c - time_sp) * sample_time);
    %s_sp_c_ppg = (ppg(time_c) - val_sp) / ((time_c - time_sp) * sample_time);
    s_sp_d_ppg = (nanCalculate(time_d, @(arr,in) (arr(in(1))), ppg, time_d) - val_sp) / ((time_d - time_sp) * sample_time);
    %s_sp_d_ppg = (ppg(time_d) - val_sp) / ((time_d - time_sp) * sample_time);
    s_b_sp_ppg = (val_sp - nanCalculate(time_b, @(arr,in) (arr(in(1))), ppg, time_b)) / ((time_sp - time_b) * sample_time);
    %s_b_sp_ppg = (val_sp - ppg(time_b)) / ((time_sp - time_b) * sample_time);
    s_b_c_ppg = nanCalculate([time_c time_b], @(arr,in) (arr(in(1)) - arr(in(2))), ppg, [time_c time_b]) / ((time_c - time_b) * sample_time);
    %s_b_c_ppg = (ppg(time_c) - ppg(time_b)) / ((time_c - time_b) * sample_time);
    s_b_d_ppg = nanCalculate([time_d time_b], @(arr,in) (arr(in(1)) - arr(in(2))), ppg, [time_d time_b]) / ((time_d - time_b) * sample_time);
    %s_b_d_ppg = (ppg(time_d) - ppg(time_b)) / ((time_d - time_b) * sample_time);
    s_u_sp_ppg = (val_sp - nanCalculate(time_u, @(arr,in) (arr(in(1))), ppg, time_u)) / ((time_sp - time_u) * sample_time);
    %s_u_sp_ppg = (val_sp - ppg(time_u)) / ((time_sp - time_u) * sample_time);
    s_on_sp_ppg = (val_sp - val_on) / ((time_sp - time_on) * sample_time);
    s_a_b_ppg = nanCalculate([time_b time_a], @(arr,in) (arr(in(1)) - arr(in(2))), ppg, [time_b time_a]) / ((time_b - time_a) * sample_time);
    %s_a_b_ppg = (ppg(time_b) - ppg(time_a)) / ((time_b - time_a) * sample_time);
    
    s_a_b_apg = (val_b - val_a) / ((time_b - time_a) * sample_time);
    s_b_sp_apg = (nanCalculate(time_sp, @(arr,in) (arr(in(1))), apg, time_sp) - val_b) / ((time_sp - time_b) * sample_time);
    %s_b_sp_apg = (apg(time_sp) - val_b) / ((time_sp - time_b) * sample_time);
    s_b_c_apg = (val_c - val_b) / ((time_c - time_b) * sample_time);
    s_b_d_apg = (val_d - val_b) / ((time_d - time_b) * sample_time);
    s_b_e_apg = (val_e - val_b) / ((time_e - time_b) * sample_time);
    s_sp_c_apg = (val_c - nanCalculate(time_sp, @(arr,in) (arr(in(1))), apg, time_sp)) / ((time_c - time_sp) * sample_time);
    %s_sp_c_apg = (val_c - apg(time_sp)) / ((time_c - time_sp) * sample_time);
    s_u_sp_apg = nanCalculate([time_sp time_u], @(arr,in) (arr(in(1)) - arr(in(2))), apg, [time_sp time_u]) / ((time_sp - time_u) * sample_time);
    %s_u_sp_apg = (apg(time_sp) - apg(time_u)) / ((time_sp - time_u) * sample_time);
    s_on_sp_apg = nanCalculate([time_sp time_on], @(arr,in) (arr(in(1)) - arr(in(2))), apg, [time_sp time_on]) / ((time_sp - time_on) * sample_time);
    %s_on_sp_apg = (apg(time_sp) - apg(time_on)) / ((time_sp - time_on) * sample_time);

    feature.slope(1) = s_sp_c_ppg;
    feature.slope(2) = s_sp_d_ppg;
    feature.slope(3) = s_b_sp_ppg;
    feature.slope(4) = s_b_c_ppg;
    feature.slope(5) = s_b_d_ppg;
    feature.slope(6) = s_u_sp_ppg;
    feature.slope(7) = s_on_sp_ppg;
    feature.slope(8) = s_a_b_ppg;
    feature.slope(9) = s_a_b_apg;
    feature.slope(10) = s_b_sp_apg;
    feature.slope(11) = s_b_c_apg;
    feature.slope(12) = s_b_d_apg;
    feature.slope(13) = s_b_e_apg;
    feature.slope(14) = s_sp_c_apg;
    feature.slope(15) = s_u_sp_apg;
    feature.slope(16) = s_on_sp_apg;

    feature.total = [feature.fiducial_value feature.fiducial_time feature.timespan ...
                     feature.amplitude feature.vpg_apg feature.waveform_area ...
                     feature.power_area feature.ratio feature.slope];
end

% calculates function if all variables in check are not NaN, otherwise returns NaN
function out = nanCalculate(check, func, arr, inputs)
    for i = 1:length(check)
        if isnan(check(i))
            out = NaN;
            return
        end
    end
    out = func(arr, inputs);
end