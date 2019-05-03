function [p0, winsize] = cen_gen(Y)
% generate center of the frames to be tracked by LK tracker
%   Jinghao Lu, 01/30/2017

    %%% find correct threshold to generate the candidate center patch that
    %%% contains sufficient information %%%
    szbd = [31, 31];
    Ymax = max(Y, [], 3);
    Ymax = Ymax(4: end - 3, 4: end - 3);
    ithres = 0.5;
    n = 1;
    ar = zeros(1, 4);
    count = 0;
    mstp = 11; %%% cut off the boundary %%%
    Yuset = normalize(imgaussfilt(Ymax, 2));
    Yuse = Yuset(mstp + 1: end - mstp, mstp + 1: end - mstp);
    while ~all(ar(3: 4) >= szbd) && ithres > 0
        mskt = Yuse > ithres;
        [l, n] = bwlabeln(mskt);
        count = count + 1;
        try
            ar = regionprops(double(mskt), 'BoundingBox');
            ar = ar.BoundingBox;
        catch
            ar = zeros(1, 4);
        end
        ithres = exp(-0.1 * count) - 0.5;
    end
    
    %%% compute entropy of the raw patch %%%
    entref = entropy(Ymax(round(ar(2)) + 1: round(ar(2)) + ar(4) - 1, round(ar(1)) + 1: round(ar(1)) + ar(3) - 1));
    
    %%% compute the entropy of all subpatches with one continuous component
    %%% deleted %%%
    entall = [];
    bbc = {};
    count = 1;
    if n > 1
        for i = 1: n
            imgt = xor(mskt, l == i);
            a = regionprops(double(imgt), 'Boundingbox');
            a = a.BoundingBox;
            if a(3) >= szbd(2) && a(4) >= szbd(1)
                entall(count) = entropy(Ymax(round(a(2)) + 1: round(a(2)) + a(4) - 1, round(a(1)) + 1: round(a(1)) + a(3) - 1));
                bbc{count} = a;
                count = count + 1;
            end
        end
    end
    
    %%% find the case with largest entropy and use that %%%
    [~, idbest] = max([entall, entref]);
    if idbest == count
        bb = ar;
    else
        bb = bbc{idbest};
    end
    p0 = mstp + flipud(round(bb(1: 2) + bb(3: 4) / 2)');
    winsize = flipud((2 * floor(bb(3: 4) / 2) - 1)');
end