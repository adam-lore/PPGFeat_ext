function [N,D,A_PPG_N,A_PPG_D,c,d,A_APG_c,A_APG_d,e,f,A_APG_e,A_APG_f] = cal_c_d(minapg,maxapg,minjpg,maxjpg,z_apg,z_jpg,APG,PPG,Jpg)
%comparison for JPG max
while (Jpg(minjpg(2)) < 0)
    if (maxjpg(1) > minjpg(1))
        c = maxjpg(1);
        d = z_apg(2);
        A_APG_c = APG(c);
        A_APG_d = APG(d);
        e = z_jpg(3);
        f = z_jpg(4);
        A_APG_e = APG(e);
        A_APG_f = APG(f);
        N = e;
        D = f;
        A_PPG_N = PPG(N);
        A_PPG_D = PPG(D);
        disp("One")
        break
    end
    
    if (maxjpg(2) > minjpg(1))
        c = maxjpg(2);
        d = z_apg(2);
        A_APG_c = APG(c);
        A_APG_d = APG(d);
        e = z_jpg(3);
        f = z_jpg(4);
        A_APG_e = APG(e);
        A_APG_f = APG(f);
        N = e;
        D = f;
        A_PPG_N = PPG(N);
        A_PPG_D = PPG(D);
        disp("Two")
        break
    end
    
    if (maxjpg(3) > minjpg(1))
        c = maxjpg(3);
        d = z_apg(2);
        A_APG_c = APG(c);
        A_APG_d = APG(d);
        e = z_jpg(3);
        f = z_jpg(4);
        A_APG_e = APG(e);
        A_APG_f = APG(f);
        N = e;
        D = f;
        A_PPG_N = PPG(N);
        A_PPG_D = PPG(D);
        disp("three")
        break
    end
end

%finding c d case 2 (not clearnly detected c and  d but present in the graph)
while (Jpg(minjpg(2)) >0)
    if (maxapg(2) < 0)
        c = maxapg(2);
        d = minapg(2);
        A_APG_c = APG(c);
        A_APG_d = APG(d);
        e = z_jpg(3);
        f = z_jpg(4);
        A_APG_e = APG(e);
        A_APG_f = APG(f);
        N = e;
        D = f;
        A_PPG_N = PPG(N);
        A_PPG_D = PPG(D);
        disp("c d founded in APG")
        break
    end

    if (maxapg(2) > 0)
        c = minjpg(2)-20;
        d = minjpg(2)+20;
        A_APG_c = APG(c);
        A_APG_d = APG(d);
        e = z_jpg(3);
        f = z_jpg(4);
        A_APG_e = APG(e);
        A_APG_f = APG(f);
        N = e;
        D = f;
        A_PPG_N = PPG(N);
        A_PPG_D = PPG(D);
        disp("c d founded in JPG + -20")
        break
    end

end
end