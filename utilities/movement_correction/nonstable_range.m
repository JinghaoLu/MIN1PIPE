function [bstt, bstp] = nonstable_range(sttn, stpn, nf)
% update nonstable section frame boundaries
%   Jinghao Lu, 02/02/2018

    bstt = stpn(1: end - 1) + 1;
    bstp = sttn(2: end) - 1;
    
    tmp = bstp - bstt;
    id = tmp >= 0;
    bstt = bstt(id);
    bstp = bstp(id);
    
    if sttn(1) ~= 1
        bstt = [1; bstt];
        bstp = [sttn(1) - 1; bstp];
    end
    
    if stpn(end) ~= nf
        bstt = [bstt; stpn(end) + 1];
        bstp = [bstp; nf];
    end
end