function [scrout, spath] = TSP(smatrix)
% pseudo-travelling-salesman-problem solver
%   Jinghao Lu, 09/01/2017

    n = size(smatrix, 1);
    if n < 10
        allpaths = perms(1: n);
        allpathsi = allpaths(:, 1: end - 1);
        allpathsj = allpaths(:, 2: end);
        allinds = reshape(sub2ind([n, n], allpathsi(:), allpathsj(:)), [], n - 1);
        alldist = smatrix(allinds);
        [scrout, idpath] = min(sum(alldist, 2));
        spath = allpaths(idpath, :);
    else
        icur = 1;
        irem = setdiff(1: n, icur);
        spath = ones(n - 1, n);
        spath(:, 1) = icur;
        scr = zeros(1, n - 1);
        for i = 1: n - 1
            icurt = irem(i);
            spath(i, 2) = icurt;
            scr(i) = scr(i) + smatrix(icur, icurt);
            iremt = setdiff(irem, icurt);
            for j = 3: n
                snext = smatrix(icurt, iremt);
                [scrnext, idnext] = min(snext);
                inext = iremt(idnext);
                spath(i, j) = inext;
                scr(i) = scr(i) + scrnext;
                icurt = inext;
                iremt = setdiff(iremt, icurt);
            end
        end
        [scrout, idfn] = min(scr);
        spath = spath(idfn, :);
    end
end