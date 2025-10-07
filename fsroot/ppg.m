classdef ppg < handle

  properties (Access = private)
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
                %disp(1)
                %obj.NEditField.Value = obj.APG_maxima(2);
                %obj.DEditField.Value = obj.APG_minima(2);

%                     hold(obj.UIAxes5,"on")
%                     ll = [obj.APG(obj.APG_minima(1)) obj.APG(obj.APG_maxima(2))];
%                     plot(obj.UIAxes5,[obj.APG_minima(1) obj.APG_maxima(2)],ll)
%                     hold(obj.UIAxes5,"off")

                %obj.cEditField.Value = (obj.APG(obj.APG_maxima(2))-obj.APG(obj.APG_minima(1)))/(obj.APG_maxima(2)-obj.APG_minima(1));
                %obj.dEditField.Value = (obj.APG_maxima(2)-obj.APG_minima(1))/1000;
                obj.c_d_APG(obj.next) = 0;
                obj.canddpresentEditField.Value = obj.c_d_APG(obj.next);
                %obj.eEditField.Value = obj.APG_maxima(2); 
                cal_c_d(obj);
            end
%updated on 15th march
            if ((obj.APG(obj.APG_maxima(2))  < 0)  && ((obj.APG_maxima(2)-obj.APG_minima(2))<100)) %means c and d is present in APG
                obj.c_d_not = 0; %to provide infor that c & d is present
                %disp(0)
                obj.NEditField.Value = obj.APG_maxima(3);
                obj.DEditField.Value = obj.APG_minima(3);
                obj.cEditField.Value = obj.APG_maxima(2);
                obj.dEditField.Value = obj.APG_minima(2);
                obj.eEditField.Value = obj.APG_maxima(3);
                obj.fEditField.Value = obj.APG_minima(3);
                obj.c_d_APG(obj.next) = 1;
                obj.canddpresentEditField.Value = obj.c_d_APG(obj.next);
            else %changed 15th march
                obj.c_d_APG(obj.next) = 1;
                obj.canddpresentEditField.Value = obj.c_d_APG(obj.next);
                cal_c_d(obj);
            end
        end
end

end