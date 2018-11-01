function [sttn, stpn] = section_update(stt, stp, nf)
% update finer stable sections
%   Jinghao Lu, 10/16/2018

    %%% augment stt first %%%
    sttn = stt;
    stpn = stp;
    
    if ~any(stt == 1)
        sttn = [1; stt];
        stpn = [1; stp];
    end
    
    if ~any(stp == nf)
        sttn = [sttn; nf];
        stpn = [stpn; nf];
    end
    
    lnon = sttn(2: end) - stpn(1: end - 1);
    sstep = 5;
    for i = 1: length(lnon)
        if lnon(i) > sstep
            nc = ceil(lnon(i) / sstep);
            sc = round(linspace(stpn(i), sttn(i + 1), nc + 1));
            sttn = [sttn; sc(2: end - 1)'];
            stpn = [stpn; sc(2: end - 1)'];
        end
    end
    sttn = sort(sttn);
    stpn = sort(stpn);
end