classdef ppg_anal < handle

    properties (Access = public)
        loaded_data % Description  
        is_dir % wheter data is loaded from a directory or a file
        dir_files % to store files in directory
        seg %to store ppg segment
        PPG_SEG % for PPG segment
        data_idx % data index
        next % Description
        Ssqi % Description
        Sub_ID % Description
        PPG_filtered % Description
        size_data % Description
        SEG_min_max % Description
        VPG % VPG segment
        APG % APG segment
        filter_sos % filter 
        g % storing filter
        APG_SEG %to store APG segment
        c_d_not % to detect the presence of c and d in APG
        APG_maxima %store all maxima of APG
        APG_minima %store all minima of APG
        c_d_APG %to store result of c and d presence
        OFs %original sampling freq
        RFs %resampled sampling freq
        FL %filter low freq
        FH %filter high freq
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

        c_and_d_present % 1 if c and d points are present

        %Vectors for segment

        OnSpDnDp %for PPG
        uxvw %for VPG
        abcde %for APG
        feature %store feature table

        OnSpDnDp_time % dt values of PPG
        uxvw_time % dt values of VPG
        abcde_time % dt values of APG
        zero_c    % store zero crossing values
    end
  
methods
    function res = LoadPPG(obj, OFs, RFs, FL, FH)
        %data is loaded from file
        obj.is_dir = false;

        %start the counter for reading all PPG data
        obj.next = 1;
        obj.data_idx = 1;

        %set initial parameters for filter
        obj.OFs = OFs;
        obj.RFs = RFs;
        obj.FL = FL;
        obj.FH = FH;

        %REMEMBER TO CHANGE
        [file,path]= uigetfile('D:\Research\Examensarbete\Datasets\*.csv');  %read CVS file
        if (path == 0)
            res = false;
            return
        end
        obj.loaded_data = importdata([path file]); %load all RAW data

        %mimic = importdata([path file]);
        %obj.loaded_data = rot90(mimic.data(:, 2));

        obj.size_data = size(obj.loaded_data); %read size of loaded data
        obj.SEG_min_max = zeros(obj.size_data(1), 3); %variable to store location of min and max with data id
        filtered = zeros(obj.size_data); %Create matrix to store filtered PPG 
        obj.PPG_filtered = zeros(obj.size_data(1), ceil(obj.size_data(2) * (RFs / OFs))); %Create matrix to store resampled PPG
        
        obj.PPG_SEG = zeros(obj.size_data(1), ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (RFs / OFs))); %create matrix to store segment
        obj.APG_SEG = obj.PPG_SEG; %APG database

        %Create Filter 
        % Read all PPG RAW data and apply filters and store
        [A,B,C,F] = cheby2(4,20,[obj.FL obj.FH]/(obj.OFs/2));
        [obj.filter_sos , obj.g] = ss2sos(A,B,C,F);

        for i = 1:obj.size_data(1)
            %apply filter to all data
            filtered_data = filtfilt(obj.filter_sos, obj.g, obj.loaded_data(i,:));

            %normalizing the data
            m = mean(filtered_data);
            st = std(filtered_data);
            filtered(i, :) = (filtered_data - m)/st; %store filtered data
        end

        if (RFs ~= OFs)
            %resample the data from OFs to RFs
            for i = 1:obj.size_data(1)     
                obj.PPG_filtered(i, :) = resample(filtered(i, :), RFs, OFs);
            end
        else
            obj.PPG_filtered = filtered;
        end

        %initialize the variables
        obj.OnSpDnDp = zeros(1,4);
        obj.uxvw = zeros(1,4);
        obj.abcde = zeros(1,5);
        obj.feature = zeros(obj.size_data(1), 30);
        obj.c_d_APG = zeros(obj.size_data(1), 1);

        % Create empty Ssqi values
        A = ones(obj.size_data(1),1);
        B = 1:1:obj.size_data(1);
        B = B';
        obj.Ssqi = [A B];

        obj.Sub_ID = obj.Ssqi(obj.next,2);

        res = true;
    end

    function res = LoadDirectory(obj, OFs, RFs, FL, FH)
        %data is loaded from directory
        obj.is_dir = true;

        %start the counter for reading all PPG data
        obj.next = 1;
        obj.data_idx = 1;

        %set initial parameters for filter
        obj.OFs = OFs;
        obj.RFs = RFs;
        obj.FL = FL;
        obj.FH = FH;

        %REMEMBER TO CHANGE
        path = uigetdir('D:\Research\Examensarbete\Datasets');  %read CVS file
        if (path == 0)
            res = false;
            return
        end
        obj.dir_files = dir([path, '\*.csv']);

        obj.size_data = [size(obj.dir_files), 1]; %read number of entries
        obj.SEG_min_max = zeros(obj.size_data(1), 3); %variable to store location of min and max with data id

        obj.PPG_filtered = zeros(obj.size_data(1), 1); %Create matrix to store resampled PPG
        obj.PPG_SEG = obj.PPG_filtered;
        obj.APG_SEG = obj.PPG_filtered;


        %initialize the variables
        obj.OnSpDnDp = zeros(1,4);
        obj.uxvw = zeros(1,4);
        obj.abcde = zeros(1,5);
        obj.feature = zeros(obj.size_data(1), 30);
        obj.c_d_APG = zeros(obj.size_data(1), 1);

        % Create empty Ssqi values
        A = ones(1,1);
        B = 1:1:1;
        B = B';
        obj.Ssqi = [A B];

        obj.Sub_ID = obj.Ssqi(obj.next,2);

        obj.LoadFile();

        res = true;
    end

    function LoadSsqi(obj)
        %REMEMBER TO CHANGE
        [file1,path1]= uigetfile('D:\Research\Examensarbete\Datasets\*.csv');
        obj.Ssqi= importdata([path1 file1]);
    end

    function res = Next(obj)
        if (obj.next < obj.size_data(1))

            obj.next = obj.next + 1;

            if (obj.is_dir)
                obj.LoadFile();
            else
                obj.data_idx = obj.data_idx + 1;
            end

            obj.Sub_ID = obj.Ssqi(obj.next,2);

            res = true;
        else
            res = false;
        end
    end

    function [PPG_max, PPG_min] = FindSegMaxMin(obj)

        segment = obj.PPG_filtered(obj.data_idx ,:);
        %find maxima and minima of current segment
        PPG_max = islocalmax(segment ,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/5));
        PPG_min = islocalmin(segment ,"MinProminence",0.1,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/10));
    end

    function CalculateFiducial(obj, min1, min2)
        %select the segment from the filtered ppg
        obj.seg = obj.PPG_filtered(obj.data_idx, (min1 - ceil(obj.RFs/67)):(min2 + ceil(obj.RFs/67)));
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
        obj.SEG_min_max(obj.next,1) = obj.Sub_ID;
        obj.SEG_min_max(obj.next,2) = min1 - ceil(obj.RFs/67);
        obj.SEG_min_max(obj.next,3) = min2 + ceil(obj.RFs/67);

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
        if abs(index_min_apg(1) - index_max_apg(1)) <= 50
            disp('The first minimum and maximum are close to each other by 50 units.');
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
        obj.APG_c_d_test();

        obj.OnSpDnDp = [obj.seg(obj.on) obj.seg(obj.sp) obj.seg(obj.dn) obj.seg(obj.dp)];

        obj.u = index_max_vpg(1);
        obj.x = index_max(1);
        obj.v = index_min_vpg(1);
        obj.w = index_max_vpg(2);

        obj.uxvw = [obj.VPG(obj.u) obj.VPG(obj.x) obj.VPG(obj.v) obj.VPG(obj.w)];

        obj.a = index_max_apg(1);
        obj.b = index_min_apg(1);

        obj.abcde = [obj.APG(obj.a) obj.APG(obj.b) obj.APG(obj.c) obj.APG(obj.d) obj.APG(obj.e)];

        %time variables of all waveform
        obj.OnSpDnDp_time  = [obj.on obj.sp obj.dn obj.dp];
        obj.uxvw_time  = [obj.u obj.x obj.v obj.w];
        obj.abcde_time = [obj.a obj.b obj.c obj.d obj.e];

        %for f values of APG
        F_Value = obj.APG(obj.f);
        F_t = obj.f;
        %for On+1 values of PPG
        O_next_t = (min2 - min1);
        O_next = obj.seg(O_next_t);

        %feature table
        obj.feature(obj.next,:) = [obj.OnSpDnDp O_next obj.uxvw obj.abcde F_Value obj.OnSpDnDp_time O_next_t obj.uxvw_time obj.abcde_time F_t];

        seg_size = size(obj.PPG_SEG);
        %store segments by zero padding remaining values
        obj.PPG_SEG(obj.next,:) = [obj.seg zeros(1,(seg_size(2))-length(obj.seg))];
        obj.APG_SEG(obj.next,:) = [obj.APG zeros(1,(seg_size(2))-length(obj.APG))];
    end

    function GenerateOutput(obj)
        %save segmented data of PPG for zero padded values
        writematrix(obj.PPG_SEG, 'PPG_Segments.xlsx')
        PPG = obj.PPG_SEG;

        %save segmented data of APG for zero padded values
        writematrix(obj.APG_SEG, 'APG_Segments.xlsx')
        APG_s = obj.APG_SEG;

        %save location of min and max with segment
        writematrix(obj.SEG_min_max, 'ID_min1_min2.xlsx')
        SMM = obj.SEG_min_max;

        %save filtered data of the whole input PPG 219 x 2100
        writematrix(obj.PPG_filtered, 'PPG_Filtered_HighSQI.xlsx')
        PPG_fil = obj.PPG_filtered;

        %save feature table
        writematrix(obj.feature, 'PPG_features.xlsx')
        P_feat = obj.feature;

        %save c and d presence table
        writematrix(obj.c_d_APG, 'c_d_presence.xlsx')
        C_D = obj.c_d_APG;

        save 'Results30Jan' 'PPG' 'APG_s' 'SMM' 'P_feat' PPG_fil 'C_D','-mat';
    end

end

methods (Access = private)
    function LoadFile(obj)
        %load next file in folder
        mimic = importdata([obj.dir_files(obj.next).folder '\' obj.dir_files(obj.next).name]);
        obj.loaded_data = rot90(mimic.data(:, 2));

        size_file = size(obj.loaded_data);

        if (size_file(2) > obj.size_data(2))
            obj.size_data(2) = size_file(2);
            obj.PPG_filtered = resize(obj.PPG_filtered, [obj.size_data(1) ceil(obj.size_data(2) * (obj.RFs / obj.OFs))]);
            obj.PPG_SEG = resize(obj.PPG_SEG, [obj.size_data(1) ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (obj.RFs / obj.OFs))]);
            obj.APG_SEG = resize(obj.APG_SEG, [obj.size_data(1) ceil((obj.size_data(2) - (obj.size_data(2)/3)) * (obj.RFs / obj.OFs))]);
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
        filtered(1, :) = (filtered_data - m)/st; %store filtered data

        if (obj.RFs ~= obj.OFs)
            %resample the data from OFs to RFs    
            obj.PPG_filtered(1, :) = resample(filtered(1, :), obj.RFs, obj.OFs);
        else
            obj.PPG_filtered = filtered;
        end

        A = ones(obj.size_data(1),1);
        B = 1:1:obj.size_data(1);
        B = B';
        obj.Ssqi = [A B];
    end

    function APG_c_d_test(obj)
        obj.c_d_not = 0;

        if (obj.APG(obj.APG_maxima(2)) > 0)  %means c and d is NOT present in APG
            obj.c_d_not = 1;

            obj.c_d_APG(obj.next) = 0;
            obj.c_and_d_present = obj.c_d_APG(obj.next);

            cal_c_d(obj);
        end
        %updated on 15th march
        if ((obj.APG(obj.APG_maxima(2))  < 0)  && ((obj.APG_maxima(2)-obj.APG_minima(2))<100)) %means c and d is present in APG
            obj.c_d_not = 0; %to provide infor that c & d is present

            obj.dn = obj.APG_maxima(3);
            obj.dp = obj.APG_minima(3);
            obj.c = obj.APG_maxima(2);
            obj.d = obj.APG_minima(2);
            obj.e = obj.APG_maxima(3);
            obj.f = obj.APG_minima(3);
            obj.c_d_APG(obj.next) = 1;
            obj.c_and_d_present = obj.c_d_APG(obj.next);
        else %changed 15th march
            obj.c_d_APG(obj.next) = 1;
            obj.c_and_d_present = obj.c_d_APG(obj.next);
            cal_c_d(obj);
        end
    end

    %%created zero crossing function on 9th jan
    function loc = zerocrossing(obj,x)
        inew = 1;
        r = x;
        loc =0;
        for ch = 2:length(r)
            
            if (((r(ch-1)< 0) && (r(ch)> 0)) || ((r(ch-1)> 0) && (r(ch)< 0)))
                loc(inew) = ch;
                inew = inew + 1;
            end
        end
        obj.zero_c = loc;
    end

    %%created function on 9th jan to implement methods for finding c
        %%and d uding new method
    function cal_c_d(obj)
        z_apg = zerocrossing(obj,obj.APG);
        Jpg = diff(obj.APG)*1000;
        Jpg = smoothdata(Jpg,"movmean", ceil(obj.RFs/12)); %data smoothing using 85 ms window at 1000Hz

        max_index_jpg = islocalmax(Jpg,"MinProminence",40,"FlatSelection","all",...
            "MinSeparation", ceil(obj.RFs/20) ,"MaxNumExtrema",5);
        min_index_jpg = islocalmin(Jpg,"MinProminence",50,"FlatSelection","all",...
            "MinSeparation",ceil(obj.RFs/20) ,"MaxNumExtrema",5);

        max_value_jpg = find(max_index_jpg); %first point of jpg is not detected
        min_value_jpg  = find(min_index_jpg);
        z_jpg = zerocrossing(obj,Jpg);

        %%To find values for c nd d using new method

        while (Jpg(min_value_jpg(2)) < 0)
            if (max_value_jpg(1) > min_value_jpg(1))

                obj.c = max_value_jpg(1);
                obj.d = z_apg(2);

                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;

                obj.dn = obj.e;
                obj.dp = f;
                break
            end

            if (max_value_jpg(2) > min_value_jpg(1))

                obj.c = max_value_jpg(2);
                obj.d = z_apg(2);

                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;

                obj.dn = obj.e;
                obj.dp = f;
                break
            end

            if (max_value_jpg(3) > min_value_jpg(1))

                obj.c = max_value_jpg(3);
                obj.d = z_apg(2);

                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;

                obj.dn = obj.e;
                obj.dp = f;
                break
            end
        end


        while (Jpg(min_value_jpg(2)) >0)
            if (obj.APG_maxima(2) < 0)
                obj.c = obj.APG_maxima(2);
                obj.d = obj.APG_minima(2);

                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;

                obj.dn = obj.e;
                obj.dp = f;
                break
            end

            if (obj.APG_maxima(2) > 0)
                obj.c = min_value_jpg(2) - obj.T2_5;
                obj.d = min_value_jpg(2) + obj.T2_5;

                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;

                obj.dn = obj.e;
                obj.dp = f;
                break
            end

        end


    end
end


end