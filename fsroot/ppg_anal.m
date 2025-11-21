classdef ppg_anal < handle

    properties (Access = public)
        loaded_data % Description  
        seg %to store ppg segment
        PPG_SEG % for PPG segment
        data_idx % data index
        entry_idx % index of current entry
        seg_idx % index of current segment in entry
        total_seg_idx % index of current segment
        Ssqi % Description
        Sub_ID % Description
        PPG_filtered % Description
        size_data % Description
        SEG_min_max % Description
        VPG % VPG segment
        APG % APG segment
        APG_SEG %to store APG segment
        APG_maxima %store all maxima of APG
        APG_minima %store all minima of APG
        c_d_APG %to store result of c and d presence
        OFs %original sampling freq
        RFs %resampled sampling freq
        FL %filter low freq
        FH %filter high freq
        is_seg %if entries should be divided into smaller segments
        seg_len %max length of segments (ms)
        num_seg %number of segments each entry is divided into
        T2_5 %2.5 % of total length of selected segment

        % For PPG Segments
        on % onset
        sp % systolic peak
        dn % dicrotic notch
        dp % diastolic peak

        % For VPG Segments
        u % global maxima in systolic phase
        x % local maxima in systolic phase
        v % global minima in systolic phase
        w % first local maxima in diastolic phase

        % For APG Segments
        a % early systolic positive peak
        b % early systolic negative peak
        c % late systolic re-increasing wave
        d % eate systolic re-decreasing wave
        e % early diastolic positive wave
        f % diastolic negative wave

        c_d_pres % presence of c and d

        %Vectors for segment

        OnSpDnDp %for PPG
        uxvw %for VPG
        abcde %for APG
        feature %store feature table

        OnSpDnDp_time % dt values of PPG
        uxvw_time % dt values of VPG
        abcde_time % dt values of APG

        entry_and_seg_id % id of the entry and the segment of the entry
    end

    properties (Access = private)
        is_dir % wheter data is loaded from a directory or a file
        dir_files % to store files in directory
        num_files % number of files in directory
        dir_ext % extension of files in directory to load
        data_col % what column the data is in
        filter_sos % filter 
        g % storing filter
    end
  
methods
    function obj = ppg_anal()
        addpath('fsroot/functions');
    end

    function res = LoadPPG(obj, OFs, RFs, FL, FH, is_seg, seg_len)
        %data is loaded from file
        obj.is_dir = false;

        %start the counter for reading all PPG data
        obj.entry_idx = 1;
        obj.seg_idx = 1;
        obj.data_idx = 1;
        obj.total_seg_idx = 1;

        %set initial parameters for filter
        obj.OFs = OFs;
        obj.RFs = RFs;
        obj.FL = FL;
        obj.FH = FH;
        obj.is_seg = is_seg;
        obj.seg_len = seg_len;

        %REMEMBER TO CHANGE
        [file, path]= uigetfile({'*.txt';'*.csv;*.xlsx';'*.mat'}, 'Load PPG File', 'D:\Research\Examensarbete\Datasets\');  %read CVS file
        if (path == 0)
            res = false;
            return
        end

        [~, ~, file_extension] = fileparts(file);
        %load all RAW data depending on file format
        switch file_extension
            case {'.txt' , '.csv', '.xlsx'}
                obj.loaded_data = importdata([path file]); 
            case '.mat'
                mat_var_arr = who('-file', [path file]);

                if (size(mat_var_arr) == 1)
                    mat_var = mat_var_arr{1};
                else
                    [var_idx, tf] = listdlg('SelectionMode', 'single', 'ListString', mat_var_arr);
                    if (tf == false)
                        res = false;
                        return
                    end
                    mat_var = mat_var_arr{var_idx};
                end
                loaded_struct = load([path file], mat_var);
                obj.loaded_data = loaded_struct.(mat_var);
            otherwise
                res = false;
                return
        end

        obj.size_data = size(obj.loaded_data); %read size of loaded data

        if (is_seg)
            obj.num_seg = ceil(obj.size_data(2) / ((seg_len * OFs) / 1000)); %max number of segments each entry is divided into
        else
            obj.num_seg = 1;
        end

        if (obj.num_seg == 1)
            obj.PPG_filtered = zeros(obj.size_data(1) * obj.num_seg, ceil(obj.size_data(2) * (RFs / OFs))); %Create matrix to store resampled PPG
        else
            obj.PPG_filtered = zeros(obj.size_data(1) * obj.num_seg, ceil(seg_len * (RFs / OFs))); %Create matrix to store resampled PPG
        end

        obj.SEG_min_max = zeros(obj.size_data(1) * obj.num_seg, 3); %variable to store location of min and max with data id
        filtered = zeros(obj.size_data); %Create matrix to store filtered PPG 
        
        
        obj.PPG_SEG = zeros(obj.size_data(1) * obj.num_seg, ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (RFs / OFs) / obj.num_seg)); %create matrix to store segment
        obj.APG_SEG = obj.PPG_SEG; %APG database

        %Create Filter 
        % Read all PPG RAW data and apply filters and store
        [A,B,C,F] = cheby2(4,20,[obj.FL obj.FH]/(obj.OFs/2));
        [obj.filter_sos , obj.g] = ss2sos(A,B,C,F);

        for i = 1:obj.size_data(1)
            if (anynan(obj.loaded_data(i)))
                obj.loaded_data(i, isnan(obj.loaded_data)) = 0; %replace NaN with 0
            end

            %apply filter to all data
            filtered_data = filtfilt(obj.filter_sos, obj.g, obj.loaded_data(i,:));

            %normalizing the data
            m = mean(filtered_data);
            st = std(filtered_data);
            filtered(i, :) = (filtered_data - m)/st; %store filtered data
        end

        if (RFs ~= OFs)
            %resample the data from OFs to RFs
            if (obj.num_seg == 1)
                for i = 1:obj.size_data(1)
                    obj.PPG_filtered(i, :) = resample(filtered(i, :), RFs, OFs);
                end
            else    
                num_per_seg = ceil(seg_len * RFs / 1000); %how many points per segment
                for i = 1:obj.size_data(1)
                    resampled = resample(filtered(i, :), RFs, OFs);
                    for j = 0:obj.num_seg - 1
                        start_seg = 1 + j * (num_per_seg);
                        end_seg = start_seg + num_per_seg - 1;
                        if (end_seg > ceil(obj.size_data(2) * (RFs / OFs)))
                            end_seg = ceil(obj.size_data(2) * (RFs / OFs));
                        end
                        obj.PPG_filtered(1 + ((i - 1) * obj.num_seg) + j, :) = resize(resampled(start_seg:end_seg), num_per_seg);
                        obj.entry_and_seg_id(1 + ((i - 1) * obj.num_seg) + j, :) = [i j + 1];
                    end
                end
            end
        else
            if (obj.num_seg == 1)
                obj.PPG_filtered = filtered;
            else
                num_per_seg = ceil((seg_len * OFs) / 1000); %how many points per segment

                for i = 1:obj.size_data(1)
                    %divide entry into segments of length num_per_seg
                    for j = 0:obj.num_seg - 1
                        start_seg = 1 + j * (num_per_seg);
                        end_seg = start_seg + num_per_seg - 1;
                        if (end_seg > obj.size_data(2))
                            end_seg = obj.size_data(2);
                        end
                        obj.PPG_filtered(1 + ((i - 1) * obj.num_seg) + j, :) = resize(filtered(i, start_seg:end_seg), num_per_seg);
                        obj.entry_and_seg_id(1 + ((i - 1) * obj.num_seg) + j, :) = [i j + 1];
                    end
                end
            end
        end

        %initialize the variables
        obj.OnSpDnDp = zeros(1,4);
        obj.uxvw = zeros(1,4);
        obj.abcde = zeros(1,5);
        obj.feature = zeros(obj.size_data(1) * obj.num_seg, 30);
        obj.c_d_APG = zeros(obj.size_data(1) * obj.num_seg, 1);
        obj.entry_and_seg_id = zeros(obj.size_data(1) * obj.num_seg, 2);

        % Create empty Ssqi values
        A = ones(obj.size_data(1),1);
        B = 1:1:obj.size_data(1);
        B = B';
        obj.Ssqi = [A B];

        obj.Sub_ID = obj.Ssqi(obj.entry_idx,2);

        res = true;
    end

    function res = LoadDirectory(obj, OFs, RFs, FL, FH, is_seg, seg_len)
        %data is loaded from directory
        obj.is_dir = true;

        %start the counter for reading all PPG data
        obj.entry_idx = 1;
        obj.seg_idx = 1;
        obj.data_idx = 1;
        obj.total_seg_idx = 1;

        %set initial parameters for filter
        obj.OFs = OFs;
        obj.RFs = RFs;
        obj.FL = FL;
        obj.FH = FH;
        obj.is_seg = is_seg;
        obj.seg_len = seg_len;

        %REMEMBER TO CHANGE
        path = uigetdir('D:\Research\Examensarbete\Datasets');  %read CVS file
        if (path == 0)
            res = false;
            return
        end
        allowed_extensions = {'.mat', '.csv', '.xlsx', '.txt'};
        [idx, tf] = listdlg('PromptString', {'Select file format'}, ...
                                'SelectionMode', 'single', 'ListString', allowed_extensions);
        if (tf == false)
            res = false;
            return
        end

        obj.dir_ext = allowed_extensions{idx};
        obj.dir_files = dir([path, '\*' obj.dir_ext]);

        obj.size_data = [size(obj.dir_files), 1]; %read number of entries
        obj.num_files = obj.size_data(1);

        if (obj.num_files == 0)
            res = false;
            return
        end

        first_file = importdata([obj.dir_files(1).folder '\' obj.dir_files(1).name]);

        if (isa(first_file, 'struct'))
            if (strcmp(obj.dir_ext, '.mat'))
                col_arr = fieldnames(first_file);
            else
                col_arr = first_file.colheaders;
            end
            if (size(col_arr) == 1)
                obj.data_col = 1;
            else
                [obj.data_col, tf] = listdlg('PromptString', {'Select column'}, 'SelectionMode','single', 'ListString', col_arr);
            end
        end
        

        if (tf == false)
            res = false;
            return
        end

        obj.SEG_min_max = zeros(obj.size_data(1), 3); %variable to store location of min and max with data id

        if (is_seg)
            num_per_seg = ceil(seg_len * RFs / 1000); %how many points per segment
            obj.size_data(2) = num_per_seg;
        else
           num_per_seg = 1; %placeholder 
        end

        obj.PPG_filtered = zeros(obj.size_data(1), num_per_seg); %Create matrix to store resampled PPG
        obj.PPG_SEG = obj.PPG_filtered;
        obj.APG_SEG = obj.PPG_filtered;

        %initialize the variables
        obj.OnSpDnDp = zeros(1,4);
        obj.uxvw = zeros(1,4);
        obj.abcde = zeros(1,5);
        obj.feature = zeros(obj.size_data(1), 30);
        obj.c_d_APG = zeros(obj.size_data(1), 1);
        obj.entry_and_seg_id = zeros(obj.size_data(1), 2);

        % Create empty Ssqi values
        A = ones(1,1);
        B = 1:1:1;
        B = B';
        obj.Ssqi = [A B];

        obj.Sub_ID = obj.Ssqi(obj.entry_idx,2);

        obj.LoadFile();

        res = true;
    end

    function LoadSsqi(obj)
        %REMEMBER TO CHANGE
        [file1,path1]= uigetfile('D:\Research\Examensarbete\Datasets\*.csv');
        obj.Ssqi= importdata([path1 file1]);
    end

    function res = Next(obj)
        if (obj.seg_idx < obj.num_seg)
            obj.seg_idx = obj.seg_idx + 1;
            obj.total_seg_idx = obj.total_seg_idx + 1;
            res = true;
            return
        elseif (obj.entry_idx < obj.size_data(1))

            obj.entry_idx = obj.entry_idx + 1;
            obj.total_seg_idx = obj.total_seg_idx + 1;
            obj.seg_idx = 1;

            if (obj.is_dir)
                obj.LoadFile();
            else
                obj.data_idx = obj.data_idx + 1;
            end

            obj.Sub_ID = obj.Ssqi(obj.entry_idx,2);

            res = true;
        else
            res = false;
        end
    end

    function [PPG_max, PPG_min] = FindSegMaxMin(obj)

        segment = obj.PPG_filtered(obj.total_seg_idx ,:);
        %find maxima and minima of current segment
        PPG_max = islocalmax(segment ,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/5));
        PPG_min = islocalmin(segment ,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/10));
    end

    function [start_idx, end_idx] = FindBestCycle(obj)
        [start_index, end_index] = obj.FindCycles();
        [cycle_idx, quality] = CalcBestCycle(start_index, end_index, obj.RFs, obj.PPG_filtered(obj.total_seg_idx ,:));
        disp(quality);

        if cycle_idx == 0
            start_idx = 0;
            end_idx = 0;
            return
        end

        start_idx = start_index(cycle_idx);
        end_idx = end_index(cycle_idx);
    end

    function res = CalculateFiducial(obj, min1, min2)
        res = true;

        %select the segment from the filtered ppg
        obj.seg = obj.PPG_filtered(obj.total_seg_idx, (min1 - ceil(obj.RFs/67)):(min2 + ceil(obj.RFs/67)));
        %plot(obj.seg);
        %calculate 2.5% of T
        obj.T2_5 = floor((((min2 - min1))/100)*2.5);

        obj.VPG = diff(obj.seg)*1000;
        obj.VPG = smoothdata(obj.VPG, "movmean", ceil(obj.RFs/20)); %data smoothing using 50 ms window at 1000Hz
        %plot(obj.VPG);

        obj.APG = diff(obj.VPG)*1000;
        obj.APG = smoothdata(obj.APG, "movmean", ceil(obj.RFs/15)); %data smoothing using 65 ms window at 1000Hz
        %obj.APG = smoothdata(obj.APG, "movmean", 65); %data smoothing using 65 ms window at 1000Hz
        %plot(obj.APG);

        %create matrix for segment information
        obj.SEG_min_max(obj.total_seg_idx,1) = obj.Sub_ID;
        obj.SEG_min_max(obj.total_seg_idx,2) = min1 - ceil(obj.RFs/67);
        obj.SEG_min_max(obj.total_seg_idx,3) = min2 + ceil(obj.RFs/67);

        %find maxima and minima of current ppg segment
        PPG_SEG_max = islocalmax(obj.seg,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/100) ,"MaxNumExtrema",2);
        PPG_SEG_min = islocalmin(obj.seg,"MinProminence",0.01,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/100) ,"MaxNumExtrema",3);

        index_max = find(PPG_SEG_max); %index of maxima points
        value_max = obj.seg(PPG_SEG_max);  %value at maxima point

        index_min = find(PPG_SEG_min) ; %index of maxima points
        value_min = obj.seg(PPG_SEG_min) ;  %value at maxima point
        
        %find maxima and minima of current vpg segment
        PPG_vpg_max = islocalmax(obj.VPG,"MinProminence",0.2,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",3);
        PPG_vpg_min = islocalmin(obj.VPG,"MinProminence",0.2,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",2);

        index_max_vpg = find(PPG_vpg_max); %index of maxima points
        value_max_vpg = obj.VPG(PPG_vpg_max);  %value at maxima point

        index_min_vpg = find(PPG_vpg_min) ; %index of maxima points
        value_min_vpg = obj.VPG(PPG_vpg_min) ;  %value at maxima point

        %find maxima and minima of current apg segment
        PPG_apg_max = islocalmax(obj.APG,"MinProminence",1.5,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/100) ,"MaxNumExtrema",5);
        PPG_apg_min = islocalmin(obj.APG,"MinProminence",1.5,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/100) ,"MaxNumExtrema",5);

        index_max_apg = find(PPG_apg_max); %index of maxima points
        value_max_apg = obj.APG(PPG_apg_max);  %value at maxima point

        index_min_apg = find(PPG_apg_min) ; %index of maxima points
        value_min_apg = obj.APG(PPG_apg_min) ;  %value at maxima point

        % new code on 15th MARCH 2023 to check two points close to each other
        if abs(index_min_apg(1) - index_max_apg(1)) <= ceil(obj.RFs/20)
            disp('The first minimum and maximum are close to each other');
            % Shift the index_min_apg vector by 1
            index_min_apg = index_min_apg(2:end);
            value_min_apg = value_min_apg(2:end);
        end

        %corrected datatips and indexes of APG
        obj.APG_maxima = index_max_apg;
        obj.APG_minima = index_min_apg;

        %read all the feature points from PPG VPG and APG
        obj.on = index_min(1);
        obj.sp = index_max(1);

        %check the presence of c and d
        %obj.APG_c_d_test();

        z_apg = zerocrossing(obj.APG);
        Jpg = diff(obj.APG)*1000;
        Jpg = smoothdata(Jpg,"movmean", ceil(obj.RFs/12)); %data smoothing using 85 ms window at 1000Hz

        max_index_jpg = islocalmax(Jpg,"MinProminence",40,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",5);
        min_index_jpg = islocalmin(Jpg,"MinProminence",50,"FlatSelection","all",...
            "MinSeparation",ceil(obj.RFs/20) ,"MaxNumExtrema",5);

        max_value_jpg = find(max_index_jpg); %first point of jpg is not detected
        min_value_jpg  = find(min_index_jpg);
        z_jpg = zerocrossing(Jpg);

        [obj.c_d_pres, obj.c, obj.d, obj.e, obj.f, obj.dn, obj.dp] = APG_c_d_test(obj.APG, obj.APG_maxima, obj.APG_minima, ...
                                                    Jpg, max_value_jpg, min_value_jpg, z_apg, z_jpg, obj.T2_5, obj.RFs);

        seg_size = size(obj.PPG_SEG);
        %store segments by zero padding remaining values
        obj.PPG_SEG(obj.total_seg_idx,:) = [obj.seg zeros(1,(seg_size(2))-length(obj.seg))];
        obj.APG_SEG(obj.total_seg_idx,:) = [obj.APG zeros(1,(seg_size(2))-length(obj.APG))];

        if (obj.c == 0 | obj.d == 0 | obj.e == 0 | obj.f == 0 | obj.dn == 0 | obj.dp == 0 | ...
            size(index_max_apg) < 1 | size(index_min_apg) < 1)
            warning('Cannot calculate PPG and APG fiducial points');
            res = false;

            obj.OnSpDnDp = [obj.seg(obj.on) obj.seg(obj.sp) 0 0];

            obj.a = 0;
            obj.b = 0;
            obj.abcde = [0 0 0 0 0];
            F_Value = 0;
        else
            obj.OnSpDnDp = [obj.seg(obj.on) obj.seg(obj.sp) obj.seg(obj.dn) obj.seg(obj.dp)];
            obj.a = index_max_apg(1);
            obj.b = index_min_apg(1);
            obj.abcde = [obj.APG(obj.a) obj.APG(obj.b) obj.APG(obj.c) obj.APG(obj.d) obj.APG(obj.e)];
            F_Value = obj.APG(obj.f);
        end

        if (size(index_max_vpg) < 2 | size(index_max) < 1 | size(index_min_vpg) < 1)
            warning('Cannot calculate VPG fiducial points');
            res = false;
            obj.u = 0;
            obj.x = 0;
            obj.v = 0;
            obj.w = 0;
            obj.uxvw = [0 0 0 0];
        else
            obj.u = index_max_vpg(1);
            obj.x = index_max(1);
            obj.v = index_min_vpg(1);
            obj.w = index_max_vpg(2);
            obj.uxvw = [obj.VPG(obj.u) obj.VPG(obj.x) obj.VPG(obj.v) obj.VPG(obj.w)];
        end

        %time variables of all waveform
        obj.OnSpDnDp_time  = [obj.on obj.sp obj.dn obj.dp];
        obj.uxvw_time  = [obj.u obj.x obj.v obj.w];
        obj.abcde_time = [obj.a obj.b obj.c obj.d obj.e];

        F_t = obj.f;
        %for On+1 values of PPG
        O_next_t = (min2 - min1);
        O_next = obj.seg(O_next_t);

        %feature table
        obj.feature(obj.total_seg_idx,:) = [obj.OnSpDnDp O_next obj.uxvw obj.abcde F_Value obj.OnSpDnDp_time O_next_t obj.uxvw_time obj.abcde_time F_t];

        %c and d presence
        obj.c_d_APG(obj.total_seg_idx) = obj.c_d_pres;
    end

    function GenerateOutput(obj)

        %save segmented data of PPG for zero padded values
        arr_size = size(obj.PPG_SEG);
        obj.PPG_SEG = resize(obj.PPG_SEG, [obj.total_seg_idx arr_size(2)]);
        writematrix(obj.PPG_SEG, 'PPG_Segments.xlsx')
        PPG = obj.PPG_SEG;

        %save segmented data of APG for zero padded values
        arr_size = size(obj.APG_SEG);
        obj.APG_SEG = resize(obj.APG_SEG, [obj.total_seg_idx arr_size(2)]);
        writematrix(obj.APG_SEG, 'APG_Segments.xlsx')
        APG_s = obj.APG_SEG;

        %save location of min and max with segment
        obj.SEG_min_max = resize(obj.SEG_min_max, [obj.total_seg_idx 3]);
        writematrix(obj.SEG_min_max, 'ID_min1_min2.xlsx')
        SMM = obj.SEG_min_max;

        %save filtered data of the whole input PPG 219 x 2100
        writematrix(obj.PPG_filtered, 'PPG_Filtered_HighSQI.xlsx')
        PPG_fil = obj.PPG_filtered;

        %save feature table
        obj.feature = resize(obj.feature, [obj.total_seg_idx 30]);
        writematrix(obj.feature, 'PPG_features.xlsx')
        P_feat = obj.feature;

        %save c and d presence table
        obj.c_d_APG = resize(obj.c_d_APG, [obj.total_seg_idx 1]);
        writematrix(obj.c_d_APG, 'c_d_presence.xlsx')
        C_D = obj.c_d_APG;

        %save entry and segment index table
        obj.entry_and_seg_id = resize(obj.entry_and_seg_id, [obj.total_seg_idx 2]);
        writematrix(obj.entry_and_seg_id, 'entry_and_seg_id.xlsx')
        index = obj.entry_and_seg_id;

        save 'Results30Jan' 'PPG' 'APG_s' 'SMM' 'P_feat' PPG_fil 'C_D' 'index','-mat';
    end

end

methods (Access = private)
    function LoadFile(obj)
        filepath = [obj.dir_files(obj.entry_idx).folder '\' obj.dir_files(obj.entry_idx).name];

        
        obj.loaded_data = importdata(filepath);

        if (isa(obj.loaded_data, 'struct'))

            if (strcmp(obj.dir_ext, '.mat'))
                mat_var_arr = who('-file', filepath);
                mat_var = mat_var_arr{obj.data_col};
                obj.loaded_data = load(filepath, mat_var);
                obj.loaded_data = obj.loaded_data.(mat_var);

                if (isa(obj.loaded_data, 'timeseries'))
                    obj.loaded_data = rot90(obj.loaded_data.Data);
                else
                    obj.loaded_data = rot90(obj.loaded_data);
                end

            else
                obj.loaded_data = rot90(obj.loaded_data.data(:, obj.data_col));
            end
        else
            if (~strcmp(obj.dir_ext, '.txt'))
                obj.loaded_data = rot90(obj.loaded_data);
            end
        end

        if (anynan(obj.loaded_data))
            obj.loaded_data(isnan(obj.loaded_data)) = 0; %replace NaN with 0
        end

        size_file = size(obj.loaded_data);

        if (obj.is_seg)
            obj.num_seg = ceil(size_file(2) / ((obj.seg_len * obj.OFs) / 1000)); %max number of segments each entry is divided into
        else
            obj.num_seg = 1;
        end

        %if this entry is the currently longest entry, risize arrays to fit it
        if (obj.is_seg == false && size_file(2) > obj.size_data(2))
            obj.size_data(2) = size_file(2);
            obj.PPG_filtered = resize(obj.PPG_filtered, [obj.size_data(1) ceil(obj.size_data(2) * (obj.RFs / obj.OFs))]);
            obj.PPG_SEG = resize(obj.PPG_SEG, [obj.size_data(1) ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (obj.RFs / obj.OFs))]);
            obj.APG_SEG = resize(obj.APG_SEG, [obj.size_data(1) ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (obj.RFs / obj.OFs))]);
        
            %if the assumed total number of segemnts after division is larger than previously assumed, resize array to fit assumtion
        elseif (obj.is_seg == true && obj.total_seg_idx + (obj.num_files - obj.entry_idx + 1) * obj.num_seg > obj.size_data(1))
            obj.size_data(1) = obj.total_seg_idx + (obj.num_files - obj.entry_idx + 1) * obj.num_seg;
            obj.SEG_min_max = resize(obj.SEG_min_max, [obj.size_data(1) 3]);
            obj.PPG_filtered = resize(obj.PPG_filtered, [obj.size_data(1) obj.size_data(2)]);
            obj.PPG_SEG = resize(obj.PPG_SEG, [obj.size_data(1) obj.size_data(2)]);
            obj.APG_SEG = resize(obj.APG_SEG, [obj.size_data(1) obj.size_data(2)]);
            obj.feature = resize(obj.feature, [obj.size_data(1) 30]);
            obj.c_d_APG = resize(obj.c_d_APG, [obj.size_data(1) 1]);
            obj.entry_and_seg_id = resize(obj.entry_and_seg_id, [obj.size_data(1) 2]);
        end

        %Create Filter 
        % Read all PPG RAW data and apply filters and store
        [A,B,C,F] = cheby2(4,20,[obj.FL obj.FH]/(obj.OFs/2));
        [obj.filter_sos , obj.g] = ss2sos(A,B,C,F);

        %apply filter to all data
        filtered_data = filtfilt(obj.filter_sos, obj.g, obj.loaded_data(1,:));

        %normalizing the data
        m = mean(filtered_data);
        st = std(filtered_data);
        filtered(obj.entry_idx, :) = (filtered_data - m)/st; %store filtered data

        if (obj.RFs ~= obj.OFs)
            %resample the data from OFs to RFs   
            if (obj.num_seg == 1) 
                obj.PPG_filtered(obj.entry_idx, :) = resample(filtered(obj.entry_idx, :), obj.RFs, obj.OFs);
            else
                resampled = resample(filtered(obj.entry_idx, :), obj.RFs, obj.OFs);
                %divide entry into segments of length obj.size_data(2)
                for j = 0:obj.num_seg - 1
                    start_seg = 1 + j * (obj.size_data(2));
                    end_seg = start_seg +  obj.size_data(2) - 1;
                    if (end_seg > ceil(size_file(2) * (obj.RFs / obj.OFs)))
                        end_seg =  ceil(size_file(2) * (obj.RFs / obj.OFs));
                    end
                    obj.PPG_filtered(obj.total_seg_idx + j, :) = resize(resampled(start_seg:end_seg), obj.size_data(2));
                    obj.entry_and_seg_id(obj.total_seg_idx + j, :) = [obj.entry_idx j + 1];
                end
            end
        else
            if (obj.num_seg == 1)
                obj.PPG_filtered = filtered;
            else
                %divide entry into segments of length obj.size_data(2)
                for j = 0:obj.num_seg - 1
                    start_seg = 1 + j * (obj.size_data(2));
                    end_seg = start_seg + obj.size_data(2) - 1;
                    if (end_seg > size_file(2))
                        end_seg = size_file(2);
                    end
                    obj.PPG_filtered(obj.total_seg_idx + j, :) = resize(filtered(obj.entry_idx, start_seg:end_seg), obj.size_data(2));
                    obj.entry_and_seg_id(obj.total_seg_idx + j, :) = [obj.entry_idx j + 1];
                end
            end
        end

        A = ones(obj.size_data(1),1);
        B = 1:1:obj.size_data(1);
        B = B';
        obj.Ssqi = [A B];
    end

    function [start_index, end_index] = FindCycles(obj)
        [PPG_max, PPG_min] = FindSegMaxMin(obj);

        index_PPG_max = find(PPG_max);
        index_PPG_min = find(PPG_min);

        [start_index, end_index] = deal(zeros(size(index_PPG_min)));

        threshold = 40 * obj.RFs/100;

        cycle_index = 0;

        % Find two minima that are sufficiently far apart and has a peak in between
        % Loop over all minima indices
        for i = 1:length(index_PPG_min)-1
            
            % Check if the next minimum index is greater than the current one by
            % the threshold value
            if index_PPG_min(i+1) > index_PPG_min(i) + threshold
                % If yes, check if there is any maximum index between the two minima
                maxima_between_minima = index_PPG_max(index_PPG_max > index_PPG_min(i) & index_PPG_max < index_PPG_min(i+1));
                if ~isempty(maxima_between_minima)
                    cycle_index = cycle_index + 1;

                    start_index(cycle_index) = index_PPG_min(i);
                    end_index(cycle_index) = index_PPG_min(i+1);
                end
            end
        end

        start_index = resize(start_index, cycle_index);
        end_index = resize(end_index, cycle_index); 
    end
end


end