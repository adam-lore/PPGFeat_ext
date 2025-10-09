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
        %disp(1)
        %obj.N = obj.APG_maxima(2);
        %obj.D = obj.APG_minima(2);

%                     hold(obj.UIAxes5,"on")
%                     ll = [obj.APG(obj.APG_minima(1)) obj.APG(obj.APG_maxima(2))];
%                     plot(obj.UIAxes5,[obj.APG_minima(1) obj.APG_maxima(2)],ll)
%                     hold(obj.UIAxes5,"off")

        %obj.c = (obj.APG(obj.APG_maxima(2))-obj.APG(obj.APG_minima(1)))/(obj.APG_maxima(2)-obj.APG_minima(1));
        %obj.d = (obj.APG_maxima(2)-obj.APG_minima(1))/1000;
        obj.c_d_APG(obj.next) = 0;
        obj.c_and_d_present = obj.c_d_APG(obj.next);
        %obj.e = obj.APG_maxima(2); 
        cal_c_d(obj);
    end
%updated on 15th march
    if ((obj.APG(obj.APG_maxima(2))  < 0)  && ((obj.APG_maxima(2)-obj.APG_minima(2))<100)) %means c and d is present in APG
        obj.c_d_not = 0; %to provide infor that c & d is present
        %disp(0)
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

              app.c = max_value_jpg(2);
              app.d = z_apg(2);
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
              A_APG_e = obj.APG(app.e);
              A_APG_f = obj.APG(f);
              obj.N = obj.e;
              obj.D = f;
              A_PPG_N = obj.seg(obj.N);
              A_PPG_D = obj.seg(obj.D);
              %disp("c d founded in APG")
              break
          end

          if (obj.APG_maxima(2) > 0)
              %disp(app.T2_5)
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

end


end