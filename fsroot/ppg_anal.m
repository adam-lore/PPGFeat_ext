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
        JPG % JPG segment
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

        c_d_pres % presence of c and d

        %Vectors for segment

        OnSpDnDpOff %for PPG
        uxvw %for VPG
        abcdef %for APG

        fiducial %store fiducial table
        feature %store feature table

        OnSpDnDpOff_time % dt values of PPG
        uxvw_time % dt values of VPG
        abcdef_time % dt values of APG

        next_peak % systolic peak of next cycle

        entry_and_seg_id % id of the entry and the segment of the entry

        quality_arr % different quality metrics for segment/cycle 


        dir_files % to store files in directory
        dir_folder % to store folder name
    end

    properties (Access = private)
        is_dir % wheter data is loaded from a directory or a file
        %dir_files % to store files in directory
        %dir_folder % to store folder name
        num_entries % number of files in directory
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
        [file, path]= uigetfile({'*.*';'*.txt';'*.csv;*.xlsx';'*.mat'}, 'Load PPG File', 'D:\Research\Examensarbete\Datasets\');  %read CVS file
        if (path == 0)
            res = false;
            return
        end

        [~, ~, file_extension] = fileparts(file);
        %load all RAW data depending on file format
        switch file_extension
            case '.txt' 
                obj.loaded_data = importdata([path file]); 
            case {'.csv', '.xlsx'}
                loaded_file = importdata([path file]);
                if ~isa(loaded_file, 'struct')
                    obj.loaded_data = loaded_file;
                else
                    obj.loaded_data = loaded_file.data;

                    obj.size_data = size(obj.loaded_data); %read size of loaded data

                    % The data may be orientated differently when there is only one entry
                    % (And anyone who organized each entry by column is crazy)
                    if obj.size_data(1) > obj.size_data(2)
                        col_arr = loaded_file.colheaders;
                        if (size(col_arr) == 1)
                            ppg_col = 1;
                        else
                            [ppg_col, tf] = listdlg('PromptString', {'Select column'}, 'SelectionMode','single', 'ListString', col_arr);
                            if (tf == false)
                                res = false;
                                return
                            end
                        end
                        obj.loaded_data = rot90(obj.loaded_data(:, ppg_col));
                    end
                end
            case '.mat'
                loaded_struct = importdata([path file]);

                % Get field containing the data
                [field, res] = SelectField(loaded_struct);
                if res == false
                    return
                elseif isempty(field)
                    obj.loaded_data = loaded_struct;
                else
                    obj.loaded_data = getfield(loaded_struct, field{:});
                end

                % Find field containing the data if each entry is a struct
                [entry_field, res] = SelectField(obj.loaded_data(1));
                if res == false
                    return
                elseif ~isempty(entry_field)
                    %old_ecg_data = obj.loaded_data;
                    obj.loaded_data = ExtractDotPath(obj.loaded_data, entry_field);
                    %obj.loaded_ecg = ExtractDotPath(old_ecg_data, {'ekg', 'v'});
                end
                obj.loaded_data = rot90(obj.loaded_data);
                %obj.loaded_ecg = rot90(obj.loaded_ecg, 3);
            otherwise
                res = false;
                return
        end

        obj.size_data = size(obj.loaded_data); %read size of loaded data

        % The data may be orientated differently when there is only one entry
        % (And anyone who organized each entry by column is crazy)
        if obj.size_data(1) > obj.size_data(2)
            obj.loaded_data = rot90(obj.loaded_data);

            [row_num, tf] = listdlg('PromptString', {'Select the column with PPG data'}, ...
                                        'SelectionMode', 'single', 'ListString', string(1:obj.size_data(2)));
            if (tf == false)
                res = false;
                return
            end
            obj.loaded_data = obj.loaded_data(row_num, :);
            obj.size_data = size(obj.loaded_data);
        end

        obj.num_entries = obj.size_data(1);

        if (is_seg)
            obj.num_seg = ceil(obj.size_data(2) / ((seg_len * OFs) / 1000)); %max number of segments each entry is divided into
        else
            obj.num_seg = 1;
        end

        if (obj.num_seg == 1)
            obj.PPG_filtered = zeros(obj.size_data(1) * obj.num_seg, ceil(obj.size_data(2) * (RFs / OFs))); %Create matrix to store resampled PPG
        else
            obj.PPG_filtered = zeros(obj.size_data(1) * obj.num_seg, ((seg_len * RFs) / 1000)); %Create matrix to store resampled PPG
            %obj.ECG_Filtered = zeros(obj.size_data(1) * obj.num_seg, ((seg_len * RFs) / 1000));
        end

        total_seg = obj.size_data(1) * obj.num_seg;

        obj.entry_and_seg_id = zeros(total_seg, 2);
        obj.quality_arr = zeros(total_seg, 6);

        obj.SEG_min_max = zeros(obj.size_data(1) * obj.num_seg, 2); %variable to store location of min and max with data id
        filtered = zeros(obj.size_data); %Create matrix to store filtered PPG 
        %ecg_filtered = zeros(obj.size_data); %Create matrix to store filtered ECG
        
        
        obj.PPG_SEG = zeros(obj.size_data(1) * obj.num_seg, ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (RFs / OFs) / obj.num_seg)); %create matrix to store segment
        obj.APG_SEG = obj.PPG_SEG; %APG database

        %Create Filter 
        % Read all PPG RAW data and apply filters and store
        [A,B,C,F] = cheby2(4,20,[obj.FL obj.FH]/(obj.OFs/2));
        [obj.filter_sos , obj.g] = ss2sos(A,B,C,F);

        for i = 1:obj.size_data(1)
            %if (anynan(obj.loaded_data(i)) || anynan(obj.loaded_ecg(i)))
            if (anynan(obj.loaded_data(i)))
                obj.loaded_data(i, isnan(obj.loaded_data)) = 0; %replace NaN with 0
                %obj.loaded_ecg(i, isnan(obj.loaded_ecg)) = 0; %replace NaN with 0
            end

            %apply filter to all data
            filtered_data = filtfilt(obj.filter_sos, obj.g, obj.loaded_data(i,:));
            %filtered_ecg = filtfilt(obj.filter_sos, obj.g, obj.loaded_ecg(i,:));

            %normalizing the data
            m = mean(filtered_data);
            st = std(filtered_data);
            filtered(i, :) = (filtered_data - m)/st; %store filtered data

            %normalizing the data
            %m = mean(filtered_ecg);
            %st = std(filtered_ecg);
            %ecg_filtered(i, :) = (filtered_ecg - m)/st; %store filtered data
        end

        if (RFs ~= OFs)
            %resample the data from OFs to RFs
            if (obj.num_seg == 1)
                for i = 1:obj.size_data(1)
                    obj.PPG_filtered(i, :) = resample(filtered(i, :), RFs, OFs);
                    obj.entry_and_seg_id(i, :) = [i 1];
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
                obj.entry_and_seg_id = [(1:obj.size_data(1))', ones(obj.size_data(1),1)];
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
                        %obj.ECG_Filtered(1 + ((i - 1) * obj.num_seg) + j, :) = resize(obj.loaded_ecg(i, start_seg:end_seg), num_per_seg);
                        obj.entry_and_seg_id(1 + ((i - 1) * obj.num_seg) + j, :) = [i j + 1];
                    end
                end
            end
        end

        %initialize the variables
        obj.OnSpDnDpOff = NaN(1,4);
        obj.uxvw = NaN(1,4);
        obj.abcdef = NaN(1,5);

        obj.fiducial = struct('OnSpDnDpOff_value', NaN(obj.size_data(1),5), 'uxvw_value', NaN(obj.size_data(1),4), ...
                              'abcdef_value', NaN(obj.size_data(1),6), 'OnSpDnDpOff_index', NaN(obj.size_data(1),5), ...
                              'uxvw_index', NaN(obj.size_data(1),4), 'abcdef_index', NaN(obj.size_data(1),6));
        obj.feature = struct('total', NaN(total_seg,145), 'fiducial_value', NaN(total_seg,15), 'fiducial_time', NaN(total_seg,15), ...
                                'timespan', NaN(total_seg,22), 'amplitude', NaN(total_seg,14), ...
                                'vpg_apg', NaN(total_seg,12), 'waveform_area', NaN(total_seg,4), ...
                                'power_area', NaN(total_seg,15), 'ratio', NaN(total_seg,32), 'slope', NaN(total_seg,16));

        obj.c_d_APG = zeros(total_seg, 1);

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
        directory = dir([path, '\*' obj.dir_ext]);

        obj.dir_files = {directory.name};
        obj.dir_files = natsortfiles(obj.dir_files);
        first_file = directory(1);
        obj.dir_folder = first_file.folder;

        obj.size_data = [length(obj.dir_files), 1]; %read number of entries
        obj.num_entries = obj.size_data(1);

        if (obj.num_entries == 0)
            res = false;
            return
        end

        first_file = importdata(strjoin([obj.dir_folder '\' obj.dir_files(1)], ''));

        if (isa(first_file, 'struct'))
            if (strcmp(obj.dir_ext, '.mat'))
                col_arr = fieldnames(first_file);
            else
                if isfield(first_file, 'colheaders')
                    col_arr = first_file.colheaders;
                elseif isfield(first_file, 'textdata')
                    col_arr = first_file.textdata(1, :);
                else
                    res = false;
                    return
                end
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

        obj.SEG_min_max = zeros(obj.size_data(1), 2); %variable to store location of min and max with data id

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
        obj.OnSpDnDpOff = NaN(1,4);
        obj.uxvw = NaN(1,4);
        obj.abcdef = NaN(1,5);

        obj.fiducial = struct('OnSpDnDpOff_value', NaN(obj.size_data(1),5), 'uxvw_value', NaN(obj.size_data(1),4), ...
                              'abcdef_value', NaN(obj.size_data(1),6), 'OnSpDnDpOff_index', NaN(obj.size_data(1),5), ...
                              'uxvw_index', NaN(obj.size_data(1),4), 'abcdef_index', NaN(obj.size_data(1),6));
        obj.feature = struct('total', NaN(obj.size_data(1),145), 'fiducial_value', NaN(obj.size_data(1),15), 'fiducial_time', NaN(obj.size_data(1),15), ...
                                'timespan', NaN(obj.size_data(1),22), 'amplitude', NaN(obj.size_data(1),14), ...
                                'vpg_apg', NaN(obj.size_data(1),12), 'waveform_area', NaN(obj.size_data(1),4), ...
                                'power_area', NaN(obj.size_data(1),15), 'ratio', NaN(obj.size_data(1),32), 'slope', NaN(obj.size_data(1),16));
        
        obj.c_d_APG = zeros(obj.size_data(1), 1);
        obj.entry_and_seg_id = zeros(obj.size_data(1), 2);
        obj.quality_arr = zeros(obj.size_data(1), 6);

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
        elseif (obj.entry_idx < obj.num_entries)

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

    function [PPG_max_index, PPG_min_index] = FindSegMaxMin(obj)

        %if (if_ecg)
        %    segment = obj.ECG_Filtered(obj.total_seg_idx ,:);
        %else
        %    segment = obj.PPG_filtered(obj.total_seg_idx ,:);
        %end

        segment = obj.PPG_filtered(obj.total_seg_idx ,:);
        %find maxima and minima of current segment
        PPG_max = islocalmax(segment ,"MinProminence",0.5,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/5));
        %tries to find a peak close after minima
        PPG_min = islocalmin(segment ,"MinProminence",0.5, "ProminenceWindow",[0 ceil(obj.RFs/5)], "FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/5));
        %exception if it is the last minima
        PPG_last_min = islocalmin(flip(segment(end - ceil(obj.RFs/5):end)), "MaxNumExtrema", 1);

        PPG_max_index = find(PPG_max);
        PPG_min_index = find(PPG_min);
        PPG_last_min_index = find(flip(PPG_last_min)) + (length(segment) - ceil(obj.RFs/5));

        if ~isempty(PPG_last_min_index) && ~isempty(PPG_min_index) && PPG_last_min_index(1) ~= PPG_min_index(end)
            PPG_min_index(end + 1) = PPG_last_min_index(1);
        end
    end

    function [start_idx, end_idx, corr_qulity, skew_quality, seg_corr_quality, seg_skew_quality, quality, seg_quality] = FindBestCycle(obj)
        [start_index, end_index, peak_index] = obj.FindCycles();
        %[ecg_max, ecg_min] = obj.FindSegMaxMin();

        %index_ECG_max = find(ecg_max);

        [cycle_idx, corr_qulity, skew_quality, seg_corr_quality, seg_skew_quality, quality, seg_quality] = CalcBestCycle(start_index, end_index, peak_index, obj.RFs, obj.PPG_filtered(obj.total_seg_idx ,:));


        if isnan(cycle_idx)
            start_idx = NaN;
            end_idx = NaN;
            return
        end

        obj.quality_arr(obj.total_seg_idx, :) = [corr_qulity, skew_quality, seg_corr_quality, seg_skew_quality, quality, seg_quality];

        start_idx = start_index(cycle_idx);
        end_idx = end_index(cycle_idx);
    end

    function res = UpdateFeatures(obj)
        min2 = obj.SEG_min_max(obj.total_seg_idx,2);

        if isempty(obj.next_peak) || length(obj.PPG_filtered) < obj.next_peak || min2 >= obj.next_peak
            next_u = NaN;
        else

            %calculate u of the next cycle which is needed to calculate extra features
            seg_next_max = obj.PPG_filtered(obj.total_seg_idx, min2:obj.next_peak);
            seg_next_VPG = diff(seg_next_max)*1000;
            seg_next_VPG = smoothdata(seg_next_VPG, "movmean", ceil(obj.RFs/20)); %data smoothing using 50 ms window at 1000Hz
            
            PPG_vpg_max = islocalmax(seg_next_VPG,"MinProminence",0.2,"FlatSelection","all",...
                "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",3);

            next_u = find(PPG_vpg_max, 1); %index of maxima points

            if isempty(next_u)
                next_u = NaN;
            end
        end
        obj.fiducial.OnSpDnDpOff_value(obj.total_seg_idx, :) = obj.OnSpDnDpOff;
        obj.fiducial.uxvw_value(obj.total_seg_idx, :) = obj.uxvw;
        obj.fiducial.abcdef_value(obj.total_seg_idx, :) = obj.abcdef;
        obj.fiducial.OnSpDnDpOff_index(obj.total_seg_idx, :) = obj.OnSpDnDpOff_time;
        obj.fiducial.uxvw_index(obj.total_seg_idx, :) = obj.uxvw_time;
        obj.fiducial.abcdef_index(obj.total_seg_idx, :) = obj.abcdef_time;

        feature_struct  = CalcFeatures(obj.OnSpDnDpOff, obj.uxvw, next_u, obj.abcdef, obj.OnSpDnDpOff_time, obj.uxvw_time, obj.abcdef_time, obj.seg, obj.VPG, obj.APG, obj.RFs);

        obj.feature.total(obj.total_seg_idx,:) = feature_struct.total;
        obj.feature.fiducial_value(obj.total_seg_idx,:) = feature_struct.fiducial_value;
        obj.feature.fiducial_time(obj.total_seg_idx,:) = feature_struct.fiducial_time;
        obj.feature.timespan(obj.total_seg_idx,:) = feature_struct.timespan;
        obj.feature.amplitude(obj.total_seg_idx,:) = feature_struct.amplitude;
        obj.feature.vpg_apg(obj.total_seg_idx,:) = feature_struct.vpg_apg;
        obj.feature.waveform_area(obj.total_seg_idx,:) = feature_struct.waveform_area;
        obj.feature.ratio(obj.total_seg_idx,:) = feature_struct.ratio;

        if anynan(obj.feature.total)
            res = false;
        else
            res = true;
        end
    end

    function res = CalculateFiducial(obj, min1, min2)
        res = true;

        if isnan(min1) || isnan(min2)
            res = false;
            return;
        end

        %create matrix for segment information
        if min1 - ceil(obj.RFs/50) < 1
            obj.SEG_min_max(obj.total_seg_idx,1) = 1;
        else
            obj.SEG_min_max(obj.total_seg_idx,1) = min1 - ceil(obj.RFs/50);
        end
        if min2 + ceil(obj.RFs/50) > length(obj.PPG_filtered)
            obj.SEG_min_max(obj.total_seg_idx,2) = length(obj.PPG_filtered);
        else
            obj.SEG_min_max(obj.total_seg_idx,2) = min2 + ceil(obj.RFs/50);
        end

        %select the segment from the filtered ppg
        obj.seg = obj.PPG_filtered(obj.total_seg_idx, obj.SEG_min_max(obj.total_seg_idx,1):obj.SEG_min_max(obj.total_seg_idx,2));
        %plot(obj.seg);
        %calculate 2.5% of T
        obj.T2_5 = floor((((min2 - min1))/100)*2.5);

        onset = min1 - obj.SEG_min_max(obj.total_seg_idx,1);
        offset = min2 - obj.SEG_min_max(obj.total_seg_idx,1);

        obj.VPG = diff(obj.seg)*1000;
        obj.VPG = smoothdata(obj.VPG, "movmean", ceil(obj.RFs/20)); %data smoothing using 50 ms window at 1000Hz
        %plot(obj.VPG);

        obj.APG = diff(obj.VPG)*1000;
        obj.APG = smoothdata(obj.APG, "movmean", ceil(obj.RFs/15)); %data smoothing using 65 ms window at 1000Hz
        %obj.APG = smoothdata(obj.APG, "movmean", 65); %data smoothing using 65 ms window at 1000Hz
        %plot(obj.APG);

        obj.JPG = diff(obj.APG)*1000;
        obj.JPG = smoothdata(obj.JPG,"movmean", ceil(obj.RFs/12)); %data smoothing using 85 ms window at 1000Hz

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

        if length(index_min) < 1 || length(index_max) < 1
            warning('Cannot calculate PPG');
            sp = NaN;
        else
            %read all the feature points from PPG VPG and APG
            sp = index_max(1);
        end

        z_apg = zerocrossing(obj.APG);

        max_index_jpg = islocalmax(obj.JPG,"MinProminence",40,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",5);
        min_index_jpg = islocalmin(obj.JPG,"MinProminence",50,"FlatSelection","all",...
            "MinSeparation",ceil(obj.RFs/20) ,"MaxNumExtrema",5);

        index_max_jpg = find(max_index_jpg); %first point of jpg is not detected
        index_min_jpg  = find(min_index_jpg);
        z_jpg = zerocrossing(obj.JPG);

        [obj.c_d_pres, c, d, e, f, dn, dp] = APG_c_d_test(obj.APG, obj.APG_maxima, obj.APG_minima, ...
                                                    obj.JPG, index_max_jpg, index_min_jpg, z_apg, z_jpg, obj.T2_5, obj.RFs);

        seg_size = size(obj.PPG_SEG);
        %store segments by zero padding remaining values
        %obj.PPG_SEG(obj.total_seg_idx,:) = [obj.seg zeros(1,(seg_size(2))-length(obj.seg))];
        %obj.APG_SEG(obj.total_seg_idx,:) = [obj.APG zeros(1,(seg_size(2))-length(obj.APG))];

        % nomralize vpg
        m_vpg = mean(obj.VPG);
        st_vpg = std(obj.VPG);
        obj.VPG = (obj.VPG - m_vpg)/st_vpg;
        % nomralize apg
        m_apg = mean(obj.APG);
        st_apg = std(obj.APG);
        obj.APG = (obj.APG - m_apg)/st_apg;

        if (any(isnan([c, d, e, f, dn, dp, onset, sp, offset])) | length(index_max_apg) < 1 | length(index_min_apg) < 1)
            warning('Cannot calculate PPG and APG fiducial points');
            res = false;

            if any(isnan([onset, sp, offset]))
                obj.OnSpDnDpOff = [NaN NaN NaN NaN NaN];
            else
                obj.OnSpDnDpOff = [obj.seg(onset) obj.seg(sp) NaN NaN obj.seg(offset)];
            end

            a = NaN;
            b = NaN;
            obj.abcdef = [NaN NaN NaN NaN NaN NaN];
        else
            obj.OnSpDnDpOff = [obj.seg(onset) obj.seg(sp) obj.seg(dn) obj.seg(dp) obj.seg(offset)];
            a = index_max_apg(1);
            b = index_min_apg(1);
            obj.abcdef = [obj.APG(a) obj.APG(b) obj.APG(c) obj.APG(d) obj.APG(e) obj.APG(f)];
        end

        if (size(index_max_vpg) < 2 | size(index_max) < 1 | size(index_min_vpg) < 1)
            warning('Cannot calculate VPG fiducial points');
            res = false;
            u = NaN;
            x = NaN;
            v = NaN;
            w = NaN;
            obj.uxvw = [NaN NaN NaN NaN];
        else
            u = index_max_vpg(1);
            x = index_max(1);
            v = index_min_vpg(1);
            w = index_max_vpg(2);
            obj.uxvw = [obj.VPG(u) obj.VPG(x) obj.VPG(v) obj.VPG(w)];
        end

        %time variables of all waveform
        obj.OnSpDnDpOff_time  = [onset sp dn dp offset];
        obj.uxvw_time  = [u x v w];
        obj.abcdef_time = [a b c d e f];
        
        %c and d presence
        obj.c_d_APG(obj.total_seg_idx) = obj.c_d_pres;

        if length(obj.PPG_filtered(obj.total_seg_idx, :)) < min2 + ceil(obj.RFs/2)
            next_cycle_seg = obj.PPG_filtered(obj.total_seg_idx , min2:end);
        else
            next_cycle_seg = obj.PPG_filtered(obj.total_seg_idx , min2:(min2 + ceil(obj.RFs/2)));
        end
        %find systolic peak of next cycle
        %find maxima and minima of current segment
        next_max = islocalmax(next_cycle_seg ,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/5));
        obj.next_peak = find(next_max, 1) + min2;
        UpdateFeatures(obj);
    end

    function GenerateOutput(obj, if_mat)

        disp("generating output");

        date_string = string(datetime, 'MM-dd_HH-mm');

        %save segmented data of PPG for zero padded values
        %arr_size = size(obj.PPG_SEG);
        %obj.PPG_SEG = resize(obj.PPG_SEG, [obj.total_seg_idx arr_size(2)]);
        %writematrix(obj.PPG_SEG, 'PPG_Segments.xlsx')
        %PPG = obj.PPG_SEG;

        %save segmented data of APG for zero padded values
        %arr_size = size(obj.APG_SEG);
        %obj.APG_SEG = resize(obj.APG_SEG, [obj.total_seg_idx arr_size(2)]);
        %writematrix(obj.APG_SEG, 'APG_Segments.xlsx')
        %APG_s = obj.APG_SEG;

        %save filtered data of the whole input PPG 219 x 2100
        %writematrix(obj.PPG_filtered, 'PPG_Filtered_HighSQI.xlsx')
        %PPG_fil = obj.PPG_filtered;

        %save fiducial table
        obj.fiducial.OnSpDnDpOff_value = resize(obj.fiducial.OnSpDnDpOff_value, [obj.total_seg_idx 5]);
        obj.fiducial.uxvw_value = resize(obj.fiducial.uxvw_value, [obj.total_seg_idx 4]);
        obj.fiducial.abcdef_value = resize(obj.fiducial.abcdef_value, [obj.total_seg_idx 6]);
        obj.fiducial.OnSpDnDpOff_index = resize(obj.fiducial.OnSpDnDpOff_index, [obj.total_seg_idx 5]);
        obj.fiducial.uxvw_index = resize(obj.fiducial.uxvw_index, [obj.total_seg_idx 4]);
        obj.fiducial.abcdef_index = resize(obj.fiducial.abcdef_index, [obj.total_seg_idx 6]);

        fiducial_table = array2table([obj.fiducial.OnSpDnDpOff_value obj.fiducial.uxvw_value obj.fiducial.abcdef_value ...
                                      obj.fiducial.OnSpDnDpOff_index obj.fiducial.uxvw_index obj.fiducial.abcdef_index], ...
                                      'VariableNames', {'val_on', 'val_sp', 'val_dn', 'val_dp', 'val_off', 'val_u', 'val_x', ...
                                      'val_v', 'val_w', 'val_a', 'val_b', 'val_c', 'val_d', 'val_e', 'val_f', 'index_on', ...
                                      'index_sp', 'index_dn', 'index_dp', 'index_off', 'index_u', 'index_x', 'index_v', ...
                                      'index_w', 'index_a', 'index_b', 'index_c', 'index_d', 'index_e', 'index_f'});

        %save total feature table
        obj.feature.total = resize(obj.feature.total, [obj.total_seg_idx 145]);
        %feature_table = array2table(obj.feature.total, 'VariableNames', {'val_on', 'val_sp', 'val_dn', 'val_dp', 'val_off', 'val_u', 'val_x', 'val_v', 'val_w', 'val_a', 'val_b', 'val_c', 'val_d', 'val_e', 'val_f', 'time_on', 'time_sp', 'time_dn', 'time_dp', 'time_off', 'time_u', 'time_x', 'time_v', 'time_w', 'time_a', 'time_b', 'time_c', 'time_d', 'time_e', 'time_f', 'ts_on_sp', 'ts_on_dn', 'ts_on_dp', 'ts_on_u', 'ts_on_v', 'ts_on_a', 'ts_on_b', 'ts_on_c', 'ts_u_next_u', 'ts_sp_c', 'ts_sp_d', 'ts_sp_e', 'ts_sp_dp', 'ts_dn_dp', 'Tm_bb2', 'ts_b_c', 'ts_b_d', 'ts_u_sp', 'ts_u_w', 'ts_u_b', 'ts_u_c', 'ts_u_d', 'ts_a_c', 'am_on_sp', 'am_on_dn', 'am_on_dp', 'am_on_u', 'am_on_v', 'am_on_a', 'am_on_b', 'am_on_c', 'am_on_off', 'am_dn_sp', 'ar_on_dn__on_sp', 'ar_on_dp__on_sp', 'ar_dn_sp__on_sp', 'ar_dp_sp__on_sp', 'val_vpg_c', 'val_vpg_d', 'r_w_u', 'r_v_u', 'r_val_vpg_c__u', 'r_val_vpg_d__u', 'r_b_a', 'r_c_a', 'r_d_a', 'r_e_a', 'r_bcde_a', 'r_bcd_a', 'wa_on_off', 'wa_on_sp', 'wa_on_c', 'wa_on_dn', 'pa_on_sp_ppg', 'pa_u_sp_ppg', 'pa_sp_c_ppg', 'pa_sp_d_ppg', 'pa_on_sp_vpg', 'pa_u_sp_vpg', 'pa_sp_c_vpg', 'pa_sp_d_vpg', 'pa_on_sp_apg', 'pa_u_sp_apg', 'pa_sp_c_apg', 'pa_sp_d_apg', 'pa_on_off_ppg', 'pa_on_off_vpg', 'pa_on_off_apg', 'r_ts_on_a__ts_u_next_u', 'r_ts_on_u__ts_u_next_u', 'r_ts_on_b__ts_u_next_u', 'r_ts_on_sp__ts_u_next_u', 'r_ts_on_c__ts_u_next_u', 'r_ts_on_v__ts_u_next_u', 'r_ts_on_dn__ts_u_next_u', 'r_ts_u_w__ts_u_next_u', 'r_ts_sp_dp__ts_u_next_u', 'r_Tm_bb2_Tss', 'r_am_on_a__am_on_sp', 'r_am_on_u__am_on_sp', 'r_am_on_b__am_on_sp', 'r_am_on_c__am_on_sp', 'r_am_on_v__am_on_sp', 'r_am_on_off__am_on_sp', 'r_wa_dn_off__wa_on_dn', 'r_sp_on', 'r_wa_on_sp__wa_on_off', 'r_wa_on_c__wa_on_off', 'r_wa_on_dn__wa_on_off', 'r_pa_on_sp_ppg__pa_on_off_ppg', 'r_pa_u_sp_ppg__pa_on_off_ppg', 'r_pa_sp_c_ppg__pa_on_off_ppg', 'r_pa_sp_d_ppg__pa_on_off_ppg', 'r_pa_on_sp_vpg__pa_on_off_vpg', 'r_pa_u_sp_vpg__pa_on_off_vpg', 'r_pa_sp_c_vpg__pa_on_off_vpg', 'r_pa_sp_d_vpg__pa_on_off_vpg', 'r_pa_on_sp_apg__pa_on_off_apg', 'r_pa_u_sp_apg__pa_on_off_apg', 'r_pa_sp_c_apg__pa_on_off_apg', 'r_pa_sp_d_apg__pa_on_off_apg', 's_sp_c_ppg', 's_sp_d_ppg', 's_b_sp_ppg', 's_b_c_ppg', 's_b_d_ppg', 's_u_sp_ppg', 's_on_sp_ppg', 's_a_b_ppg', 's_a_b_apg', 's_b_sp_apg', 's_b_c_apg', 's_b_d_apg', 's_b_e_apg', 's_sp_c_apg', 's_u_sp_apg', 's_on_sp_apg'});
        feature_table = array2table(obj.feature.total, 'VariableNames', {'val_on', 'val_sp', 'val_dn', 'val_dp', 'val_off', 'val_u', ...
                         'val_x', 'val_v', 'val_w', 'val_a', 'val_b', 'val_c', 'val_d', 'val_e', 'val_f', 'time_on', ...
                         'time_sp', 'time_dn', 'time_dp', 'time_off', 'time_u', 'time_x', 'time_v', 'time_w', 'time_a', ...
                         'time_b', 'time_c', 'time_d', 'time_e', 'time_f', 'ts_on_sp', 'ts_on_dn', 'ts_on_dp', 'ts_on_u', ...
                         'ts_on_v', 'ts_on_a', 'ts_on_b', 'ts_on_c', 'ts_u_next_u', 'ts_sp_c', 'ts_sp_d', 'ts_sp_e', 'ts_sp_dp', ...
                         'ts_dn_dp', 'ts_b_c', 'ts_b_d', 'ts_u_sp', 'ts_u_w', 'ts_u_b', 'ts_u_c', 'ts_u_d', 'ts_a_c', ...
                         'am_on_sp', 'am_on_dn', 'am_on_dp', 'am_on_u', 'am_on_v', 'am_on_a', 'am_on_b', 'am_on_c', 'am_on_off', ...
                         'am_dn_sp', 'ar_on_dn__on_sp', 'ar_on_dp__on_sp', 'ar_dn_sp__on_sp', 'ar_dp_sp__on_sp', 'val_vpg_c', ...
                         'val_vpg_d', 'r_w_u', 'r_v_u', 'r_val_vpg_c__u', 'r_val_vpg_d__u', 'r_b_a', 'r_c_a', 'r_d_a', 'r_e_a', ...
                         'r_bcde_a', 'r_bcd_a', 'wa_on_off', 'wa_on_sp', 'wa_on_c', 'wa_on_dn', 'pa_on_sp_ppg', 'pa_u_sp_ppg', ...
                         'pa_sp_c_ppg', 'pa_sp_d_ppg', 'pa_on_sp_vpg', 'pa_u_sp_vpg', 'pa_sp_c_vpg', 'pa_sp_d_vpg', 'pa_on_sp_apg', ...
                         'pa_u_sp_apg', 'pa_sp_c_apg', 'pa_sp_d_apg', 'pa_on_off_ppg', 'pa_on_off_vpg', 'pa_on_off_apg', ...
                         'r_ts_on_a__ts_u_next_u', 'r_ts_on_u__ts_u_next_u', 'r_ts_on_b__ts_u_next_u', 'r_ts_on_sp__ts_u_next_u', ...
                         'r_ts_on_c__ts_u_next_u', 'r_ts_on_v__ts_u_next_u', 'r_ts_on_dn__ts_u_next_u', 'r_ts_u_w__ts_u_next_u', ...
                         'r_ts_sp_dp__ts_u_next_u', 'r_am_on_a__am_on_sp', 'r_am_on_u__am_on_sp', ...
                         'r_am_on_b__am_on_sp', 'r_am_on_c__am_on_sp', 'r_am_on_v__am_on_sp', 'r_am_on_off__am_on_sp', ...
                         'r_wa_dn_off__wa_on_dn', 'r_sp_on', 'r_wa_on_sp__wa_on_off', 'r_wa_on_c__wa_on_off', 'r_wa_on_dn__wa_on_off', ...
                         'r_pa_on_sp_ppg__pa_on_off_ppg', 'r_pa_u_sp_ppg__pa_on_off_ppg', 'r_pa_sp_c_ppg__pa_on_off_ppg', ...
                         'r_pa_sp_d_ppg__pa_on_off_ppg', 'r_pa_on_sp_vpg__pa_on_off_vpg', 'r_pa_u_sp_vpg__pa_on_off_vpg', ...
                         'r_pa_sp_c_vpg__pa_on_off_vpg', 'r_pa_sp_d_vpg__pa_on_off_vpg', 'r_pa_on_sp_apg__pa_on_off_apg', ...
                         'r_pa_u_sp_apg__pa_on_off_apg', 'r_pa_sp_c_apg__pa_on_off_apg', 'r_pa_sp_d_apg__pa_on_off_apg', ...
                         's_sp_c_ppg', 's_sp_d_ppg', 's_b_sp_ppg', 's_b_c_ppg', 's_b_d_ppg', 's_u_sp_ppg', 's_on_sp_ppg', ...
                         's_a_b_ppg', 's_a_b_apg', 's_b_sp_apg', 's_b_c_apg', 's_b_d_apg', 's_b_e_apg', 's_sp_c_apg', ...
                         's_u_sp_apg', 's_on_sp_apg'});
        

        % entry and segment index
        obj.entry_and_seg_id = resize(obj.entry_and_seg_id, [obj.total_seg_idx 2]);

        % location of min and max within segment
        obj.SEG_min_max = resize(obj.SEG_min_max, [obj.total_seg_idx 2]);

        % c and d presence in segment
        obj.c_d_APG = resize(obj.c_d_APG, [obj.total_seg_idx 1]);

        obj.quality_arr = resize(obj.quality_arr, [obj.total_seg_idx 6]);

        info_table = array2table([obj.entry_and_seg_id obj.SEG_min_max obj.c_d_APG obj.quality_arr], ...
                                 'VariableNames', {'entry_id', 'segment_id', 'min_1', 'min_2', 'c_d_pres', ...
                                                   'corr_qual', 'skew_qual', 'seg_corr_qual', 'seg_skew_qual', 'cycle_qual', 'seg_qual'});

        if if_mat
            save('PPG_Results_' + date_string, 'fiducial_table', 'feature_table', 'info_table')
        else
            writetable(fiducial_table, 'PPG_fiducials_' + date_string + '.xlsx', WriteMode='overwrite');
            writetable(feature_table, 'PPG_features_' + date_string + '.xlsx', WriteMode='overwrite')
            writetable(info_table, 'PPG_info_' + date_string + '.xlsx', WriteMode='overwrite');
        end

        
    end

    %calculate features for every entry in the loaded dataset
    function ProcessDataset(obj)
        working = true;
        progress = waitbar(0, sprintf('Entry: %d/%d Segment: %d', obj.entry_idx, obj.num_entries, obj.seg_idx), 'Name', 'Processing data');

        while working
            waitbar(obj.entry_idx/obj.num_entries, progress, ...
                    sprintf('Entry: %d/%d Segment: %d', obj.entry_idx, obj.num_entries, obj.seg_idx))

            if (obj.entry_idx == 105 && obj.seg_idx == 21)
                warning("breakpoint");
            end
            [min1, min2] = obj.FindBestCycle();
            obj.CalculateFiducial(min1, min2);
            working = obj.Next();
        end
        delete(progress);
    end
end

methods (Access = private)
    function LoadFile(obj)
        filepath = strjoin([obj.dir_folder '\' obj.dir_files(obj.entry_idx)], '');

        
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
        elseif (obj.is_seg == true && obj.total_seg_idx + (obj.num_entries - obj.entry_idx + 1) * obj.num_seg > obj.size_data(1))
            obj.size_data(1) = obj.total_seg_idx + (obj.num_entries - obj.entry_idx + 1) * obj.num_seg;
            obj.SEG_min_max = resize(obj.SEG_min_max, [obj.size_data(1) 2]);
            obj.PPG_filtered = resize(obj.PPG_filtered, [obj.size_data(1) obj.size_data(2)]);
            obj.PPG_SEG = resize(obj.PPG_SEG, [obj.size_data(1) obj.size_data(2)]);
            obj.APG_SEG = resize(obj.APG_SEG, [obj.size_data(1) obj.size_data(2)]);

            old_fiducial = obj.fiducial;
            old_fiducial_size = size(obj.fiducial.OnSpDnDpOff_value);

            obj.fiducial.OnSpDnDpOff_value = NaN(obj.size_data(1), 5);
            obj.fiducial.uxvw_value = NaN(obj.size_data(1), 4);
            obj.fiducial.abcdef_value = NaN(obj.size_data(1), 6);
            obj.fiducial.OnSpDnDpOff_index = NaN(obj.size_data(1), 5);
            obj.fiducial.uxvw_index = NaN(obj.size_data(1), 4);
            obj.fiducial.abcdef_index = NaN(obj.size_data(1), 6);
            
            obj.fiducial.OnSpDnDpOff_value(1:old_fiducial_size(1), :) = old_fiducial.OnSpDnDpOff_value;
            obj.fiducial.uxvw_value(1:old_fiducial_size(1), :) = old_fiducial.uxvw_value;
            obj.fiducial.abcdef_value(1:old_fiducial_size(1), :) = old_fiducial.abcdef_value;
            obj.fiducial.OnSpDnDpOff_index(1:old_fiducial_size(1), :) = old_fiducial.OnSpDnDpOff_index;
            obj.fiducial.uxvw_index(1:old_fiducial_size(1), :) = old_fiducial.uxvw_index;
            obj.fiducial.abcdef_index(1:old_fiducial_size(1), :) = old_fiducial.abcdef_index;

            old_feature = obj.feature;
            old_feature_size = size(obj.feature.total);

            obj.feature.total = NaN(obj.size_data(1), 145);
            obj.feature.fiducial_value = NaN(obj.size_data(1), 15);
            obj.feature.fiducial_time = NaN(obj.size_data(1), 15);
            obj.feature.timespan = NaN(obj.size_data(1), 22);
            obj.feature.amplitude = NaN(obj.size_data(1), 14);
            obj.feature.vpg_apg = NaN(obj.size_data(1), 12);
            obj.feature.waveform_area = NaN(obj.size_data(1), 4);
            obj.feature.power_area = NaN(obj.size_data(1), 15);
            obj.feature.ratio = NaN(obj.size_data(1), 32);
            obj.feature.slope = NaN(obj.size_data(1), 16);

            obj.feature.total(1:old_feature_size(1), :) = old_feature.total;
            obj.feature.fiducial_value(1:old_feature_size(1), :) = old_feature.fiducial_value;
            obj.feature.fiducial_time(1:old_feature_size(1), :) = old_feature.fiducial_time;
            obj.feature.timespan(1:old_feature_size(1), :) = old_feature.timespan;
            obj.feature.amplitude(1:old_feature_size(1), :) = old_feature.amplitude;
            obj.feature.vpg_apg(1:old_feature_size(1), :) = old_feature.vpg_apg;
            obj.feature.waveform_area(1:old_feature_size(1), :) = old_feature.waveform_area;
            obj.feature.power_area(1:old_feature_size(1), :) = old_feature.power_area;
            obj.feature.ratio(1:old_feature_size(1), :) = old_feature.ratio;
            obj.feature.slope(1:old_feature_size(1), :) = old_feature.slope;

            obj.c_d_APG = resize(obj.c_d_APG, [obj.size_data(1) 1]);
            obj.entry_and_seg_id = resize(obj.entry_and_seg_id, [obj.size_data(1) 2]);
            obj.quality_arr = resize(obj.quality_arr, [obj.size_data(1) 6]);
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

    function [start_index, end_index, peak_index] = FindCycles(obj)
        [index_PPG_max, index_PPG_min] = FindSegMaxMin(obj);

        %if obj.entry_idx == 12 && obj.seg_idx == 25
        %if obj.entry_idx == 67
        if 1 == 2
            warning('Break reached');
            temp = plot(obj.PPG_filtered(obj.total_seg_idx ,:));
            for i = 1:length(index_PPG_max)
                datatip(temp, index_PPG_max(i), obj.PPG_filtered(obj.total_seg_idx , index_PPG_max(i)),'FontSize',3)
            end
            for i = 1:length(index_PPG_min)
                datatip(temp, index_PPG_min(i), obj.PPG_filtered(obj.total_seg_idx , index_PPG_min(i)),'FontSize',3)
            end
            warning('Break done');
        end

        % Set start, end and peak index arrays to zero of size index_PPG_min
        [start_index, end_index, peak_index] = deal(zeros(size(index_PPG_min)));

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
                next_peak = index_PPG_max(find(index_PPG_max > index_PPG_min(i+1), 1));
                if ~isempty(maxima_between_minima) & ~isempty(next_peak)
                    cycle_index = cycle_index + 1;

                    start_index(cycle_index) = index_PPG_min(i);
                    end_index(cycle_index) = index_PPG_min(i+1);

                    peak_index(cycle_index) = maxima_between_minima(1);
                    peak_value = obj.PPG_filtered(maxima_between_minima(1));

                    % Find the highest peak in the cycle
                    if length(maxima_between_minima) > 1
                        for j = 2:length(maxima_between_minima)
                            if obj.PPG_filtered(maxima_between_minima(j)) > peak_value
                                peak_index(cycle_index) = maxima_between_minima(j);
                                peak_value = obj.PPG_filtered(maxima_between_minima(j));
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
end


end