function loc = zerocrossing(x)


inew = 1;
r = x;
loc =[];
for ch = 2:length(r)
    
    if (((r(ch-1)< 0) && (r(ch)> 0)) || ((r(ch-1)> 0) && (r(ch)< 0)))
        loc(inew) = ch;
        inew = inew + 1;
    end
end
end