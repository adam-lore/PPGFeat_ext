classdef ppg_anal < handle

    properties (Access = public)
        loadedData % Description
        seg %to store ppg segment
        PPG_SEG % for PPG segment
        next % Description
        Ssqi % Description
        Sub_ID % Description
        PPG_filtered % Description
        Size_loaded_data % Description
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
        Fs %sampling freq
        FL %filter low freq
        FH %filter high freq
        T2_5 %2.5 % of total length of selected segment

        % For PPG Segments
        O % Description
        S % Description
        N % Description
        D % Description

        % For VPG Segments
        W % Description
        X % Description
        Y % Description
        Z % Description

        % For APG Segments
        a % Description
        b % Description
        c % Description
        d % Description
        e % Description
        f % Description

        c_and_d_present % 1 if c and d points are present

        %Vectors for segment

        OSND %for PPG
        WXYZ %for VPG
        abcde %for APG
        feature %store feature table

        OSND_time % dt values of PPG
        WXYZ_time % dt values of VPG
        abcde_time % dt values of APG
        zero_c    % store zero crossing values
    end
  
methods
    function  APG_c_d_test(obj)
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

            obj.N = obj.APG_maxima(3);
            obj.D = obj.APG_minima(3);
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
        Jpg = smoothdata(Jpg,"movmean",85);
        %plot(Jpg);
        max_index_jpg = islocalmax(Jpg,"MinProminence",40,"FlatSelection","all",...
            "MinSeparation",50,"MaxNumExtrema",5);
        min_index_jpg = islocalmin(Jpg,"MinProminence",50,"FlatSelection","all",...
            "MinSeparation",50,"MaxNumExtrema",5);

        max_value_jpg = find(max_index_jpg); %first point of jpg is not detected
        min_value_jpg  = find(min_index_jpg);
        z_jpg = zerocrossing(obj,Jpg);

        %%To find values for c nd d using new method

        while (Jpg(min_value_jpg(2)) < 0)
            if (max_value_jpg(1) > min_value_jpg(1))

                obj.c = max_value_jpg(1);
                obj.d = z_apg(2);
                A_APG_c = obj.APG(obj.c);
                A_APG_d = obj.APG(obj.d);
                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;
                A_APG_e = obj.APG(obj.e);
                A_APG_f = obj.APG(f);
                obj.N = obj.e;
                obj.D = f;
                A_PPG_N = obj.seg(obj.N);
                A_PPG_D = obj.seg(obj.D);
                %disp("One")
                break
            end

            if (max_value_jpg(2) > min_value_jpg(1))

                obj.c = max_value_jpg(2);
                obj.d = z_apg(2);
                A_APG_c = obj.APG(obj.c);
                A_APG_d = obj.APG(obj.d);
                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;
                A_APG_e = obj.APG(obj.e);
                A_APG_f = obj.APG(f);
                obj.N = obj.e;
                obj.D = f;
                A_PPG_N = obj.seg(obj.N);
                A_PPG_D = obj.seg(obj.D);
                %disp("Two")
                break
            end

            if (max_value_jpg(3) > min_value_jpg(1))

                obj.c = max_value_jpg(3);
                obj.d = z_apg(2);
                A_APG_c = obj.APG(obj.c);
                A_APG_d = obj.APG(obj.d);
                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;
                A_APG_e = obj.APG(obj.e);
                A_APG_f = obj.APG(f);
                obj.N = obj.e;
                obj.D = f;
                A_PPG_N = obj.seg(obj.N);
                A_PPG_D = obj.seg(obj.D);
                %disp("three")
                break
            end
        end


        while (Jpg(min_value_jpg(2)) >0)
            if (obj.APG_maxima(2) < 0)
                obj.c = obj.APG_maxima(2);
                obj.d = obj.APG_minima(2);
                A_APG_c = obj.APG(obj.c);
                A_APG_d = obj.APG(obj.d);
                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;
                A_APG_e = obj.APG(obj.e);
                A_APG_f = obj.APG(f);
                obj.N = obj.e;
                obj.D = f;
                A_PPG_N = obj.seg(obj.N);
                A_PPG_D = obj.seg(obj.D);
                %disp("c d founded in APG")
                break
            end

            if (obj.APG_maxima(2) > 0)
                obj.c = min_value_jpg(2)- obj.T2_5;
                obj.d = min_value_jpg(2)+ obj.T2_5;
                A_APG_c = obj.APG(obj.c);
                A_APG_d = obj.APG(obj.d);
                obj.e = z_jpg(3);
                f = z_jpg(4);
                obj.f = f;
                A_APG_e = obj.APG(obj.e);
                A_APG_f = obj.APG(f);
                obj.N = obj.e;
                obj.D = f;
                A_PPG_N = obj.seg(obj.N);
                A_PPG_D = obj.seg(obj.D);
                %disp("c d founded in JPG + -20")
                break
            end

        end


    end

    function LoadPPG(obj, Fs, FL, FH)
        %start the counter for reading all PPG data
        obj.next = 1;

        %set initial parameters for filter
        obj.Fs = Fs;
        obj.FL = FL;
        obj.FH = FH;

        %REMEMBER TO CHANGE
        [file,path]= uigetfile('D:\Research\Examensarbete\Datasets\*.csv');  %read CVS file
        obj.loadedData = importdata([path file]); %load all RAW data

        obj.Size_loaded_data = size(obj.loadedData); %read size of loaded data
        obj.SEG_min_max = zeros(obj.Size_loaded_data(1), 3); %variable to store location of min and max with data id
        obj.PPG_filtered = zeros(obj.Size_loaded_data); %Create matrix to store filtered PPG 
        
        obj.PPG_SEG = zeros(obj.Size_loaded_data(1), (obj.Size_loaded_data(2) - (obj.Size_loaded_data(2)/3))); %create matrix to store segment
        obj.APG_SEG = obj.PPG_SEG; %APG database

        %Create Filter 
        % Read all PPG RAW data and apply filters and store
        [A,B,C,F] = cheby2(4,20,[obj.FL obj.FH]/(obj.Fs/2));
        [obj.filter_sos , obj.g] = ss2sos(A,B,C,F);

        for i = 1:obj.Size_loaded_data(1)
            %apply filter to all data
            filtered_data = filtfilt(obj.filter_sos, obj.g, obj.loadedData(i,:));

            %normalizing the data
            m = mean(filtered_data);
            st = std(filtered_data);
            obj.PPG_filtered(i,:) = (filtered_data - m)/st; %store filtered data
        end

        %initialize the variables
        obj.OSND = zeros(1,4);
        obj.WXYZ = zeros(1,4);
        obj.abcde = zeros(1,5);
        obj.feature = zeros(219,30);
        obj.c_d_APG = zeros(219,1);
    end

    function LoadSsqi(obj)
        %REMEMBER TO CHANGE
        [file1,path1]= uigetfile('D:\Research\Examensarbete\Datasets\*.csv');
        obj.Ssqi= importdata([path1 file1]);
    end

    function GenerateOutput(obj)
        %save segmented data of PPG for zero padded values

        PPG_segment = 'PPG_Segments.xlsx';
        writematrix(obj.PPG_SEG,PPG_segment)
        PPG = obj.PPG_SEG;

        %save segmented data of APG for zero padded values

        APG_segment = 'APG_Segments.xlsx';
        writematrix(obj.APG_SEG,APG_segment)
        APG_s = obj.APG_SEG;

        %save location of min and max with segment
        ID_min_max = 'ID_min1_min2.xlsx';
        writematrix(obj.SEG_min_max, ID_min_max)
        SMM = obj.SEG_min_max;


        %save filtered data of the whole input PPG 219 x 2100
        PPG_filter = 'PPG_Filtered_HighSQI.xlsx';
        writematrix(obj.PPG_filtered,PPG_filter )
        PPG_fil = obj.PPG_filtered;

        %save feature table
        
        PPG_feature = 'PPG_features.xlsx';
        writematrix(obj.feature,PPG_feature )
        P_feat = obj.feature;

        %save c and d presence table
        c_and_d = 'c_d_presence.xlsx';
        writematrix(obj.c_d_APG,c_and_d)
        C_D = obj.c_d_APG;

        save 'Results30Jan' 'PPG' 'APG_s' 'SMM' 'P_feat' PPG_fil 'C_D','-mat';
    end

end


end