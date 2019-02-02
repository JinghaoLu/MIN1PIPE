function [mxout, xfuse, lduse, iduse, smatrix, xfmatrix, ldmatrix] = logdemons_cluster(mxin, pixs, scl, sigma_x, sigma_f, sigma_d, maskc)
% full graph applying logdemons registration
%   Jinghao Lu, 06/11/2017

    if nargin < 2 || isempty(pixs)
        [pixh, pixw, ~] = size(mxin);
        pixs = min(pixh, pixw);
    end
    
    if nargin < 3 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 4 || isempty(sigma_x)
        defpar = default_parameters;
        sigma_x = defpar.mc_sigma_x;
    end

    if nargin < 5 || isempty(sigma_f)
        defpar = default_parameters;
        sigma_f = defpar.mc_sigma_f;
    end
    
    if nargin < 6 || isempty(sigma_d)
        defpar = default_parameters;
        sigma_d = defpar.mc_sigma_d;
    end
    
    if nargin < 7 || isempty(maskc)
        maskc = true(size(mxin, 1), size(mxin, 2));
    end
    
    %% initialization %%
    countfn = 1;
    idclustfn = {};
    imgallfn = {};
    xfallfn = {};
    ldallfn = {};
    mxint1 = mxin;
    [pixh, pixw, nf] = size(mxin);
    pixupper = 1; %%% hard value of single pixel; if resolution high, this restricts the upper bound, if not, this does not even affect %%%
    pixthres = scl * pixs;
    flagincrease = false;
    istep = 0.001;
    
    if nf > 1 && scl < 1
        %% main loop when there are still possibility of combination %%
        while 1
            
            if ~flagincrease
                %%% compute score and transform correspondance matrix %%%
                n = size(mxint1, 3);
                smatrix = ones(n);
                xfmatrix = cell(n);
                ldmatrix = cell(n);
                imgmatrix = cell(n);
                
                mxint2 = cell(1, n);
                for i = 1: n
                    mxint2{i} = mxint1(:, :, i + 1: end);
                end
                
                parfor i = 1: n
                    imref = mxint1(:, :, i);
                    ncur = size(mxint2{i}, 3);
                    temp = zeros(1, ncur);
                    ldtemp = cell(1, ncur);
                    xftemp = cell(1, ncur);
                    imgtemp = cell(1, ncur);
                    for j = 1: ncur
                        imcur = mxint2{i}(:, :, j);
                        [scrtc, imgt, xft] = klt_ref_track(imcur, imref, [], maskc);
                        scrto = get_trans_score_ref(imcur, imref, maskc);
                        if scrtc < scrto
                            imcur = imgt;
                            xftemp{j} = xft{1};
                        else
                            xftemp{j} = affine2d(diag(ones(1, 3)));
                        end
                        [imgo, sx, sy] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d, maskc);
                        temp(j) = get_trans_score_ref(imgo, imref, maskc);
                        ldtemp{j} = cellfun(@(x, y) cat(3, x, y), sx, sy, 'uniformoutput', false);
                        imgtemp{j} = imgo;
                    end
                    temp = [100 * ones(1, i), temp];
                    xftemp = [cell(1, i - 1), {affine2d(diag(ones(1, 3)))}, xftemp];
                    ldtemp = [cell(1, i - 1), {zeros(pixh, pixw, 2)}, ldtemp];
                    imgtemp = [cell(1, i - 1), {imref}, imgtemp];
                    smatrix(i, :) = temp;
                    xfmatrix(i, :) = xftemp;
                    ldmatrix(i, :) = ldtemp;
                    imgmatrix(i, :) = imgtemp;
%                     disp(num2str(i))
                end
            end
            
            flagincrease = false;
            
            %%% compute connected subgraph components %%%
            aa = smatrix < pixthres;
%             aa = smatrix < 0.8;
            sgraph = digraph(aa, 'OmitSelfLoops');
            names = get_node_name(1: n);
            sgraph.Nodes.Name = names;
            sgraph1 = graph(aa, 'upper');
            sgraph1.Nodes.Name = names;

            %%% decompose and combine subgraphs within the main cluster %%%
            idclust = {};
            imgall = {};
            xfall = cell(1, n);
            ldall = cell(1, n);
            for i = 1: n
                xfall{i} = {};
                ldall{i} = {};
            end
            
            %%%% stop if no edges %%%%
            if ~isempty(sgraph1.Edges)
                count = 1;
                flag = true;
                
                %%%% continue if there are still uncombined subgraph clusters %%%%
                while flag || ~all(cellfun(@(x) size(x, 3), imgall) == 1)
                    if flag
                        sgrapht = sgraph;
                        sgraph1t = sgraph1;
                    else %%%% find the largest subgraph in the remaining graph %%%%
                        idsub = find(cellfun(@(x) size(x, 3), imgall) ~= 1, 1);
                        sgrapht = subgraph(sgraph, idclust{idsub});
                        sgraph1t = subgraph(sgraph1, idclust{idsub});
                        ntemp = length(idclust);
                        idclust = idclust(setdiff(1: ntemp, idsub));
                        imgall = imgall(setdiff(1: ntemp, idsub));
                    end
                    
                    %%%% for the current subgraph cluster, combine till no edges %%%%
                    while ~isempty(sgraph1t.Nodes)
                        bb = conncomp(sgraph1t);
                        
                        %%%%% compute largest subgraph %%%%%
                        s = zeros(1, max(bb));
                        for i = 1: max(bb)
                            s(i) = sum(bb == i);
                        end
                        [~, connuse] = max(s);
                        idcu = setdiff(1: max(bb), connuse);
                        idcurall = cellfun(@str2double, sgraph1t.Nodes.Name);
                        for i = 1: length(idcu)
                            idt = str2double(sgraph1t.Nodes.Name(bb == idcu(i)));
                            idclust = [idclust; {idt(:)}];
                            imgall = [imgall; {mxint1(:, :, idt)}];
                        end
                        ids = str2double(sgraph1t.Nodes.Name(bb == connuse));
                        idrmv = setdiff(idcurall, ids);
                        namermv = get_node_name(idrmv);
                        sgrapht = rmnode(sgrapht, namermv);
                        sgraph1t = rmnode(sgraph1t, namermv);

                        %%%%% find the hub (with largest degree) %%%%%
                        idsource = ids(1);
                        
                        %%%%% combine the hub with all its neighbors %%%%%
                        nbs = successors(sgrapht, num2str(idsource));
                        idnbs = cellfun(@str2double, nbs);
                        imgtemp = reshape(cell2mat(imgmatrix(idsource, [idsource; idnbs(:)])), pixh, pixw, length(idnbs) + 1);
%                         imgall = [imgall; {max(imgtemp, [], 3)}];
                        imgall = [imgall; {mean(imgtemp, 3)}];
                        idclust = [idclust; {[idsource; idnbs]}];
                        for i = 1: length(idnbs)
                            xfall{idnbs(i)} = [xfall{idnbs(i)}, xfmatrix(idsource, idnbs(i))];
                            ldall{idnbs(i)} = [ldall{idnbs(i)}, ldmatrix(idsource, idnbs(i))];
                        end
                        
                        %%%%% update parameters within the current round of graph %%%
                        idt = [num2str(idsource); nbs];
                        
                        %%%%% update the current graph %%%%%
                        sgrapht = rmnode(sgrapht, idt);
                        sgraph1t = rmnode(sgraph1t, idt);
%                         disp(num2str(count))
                        count = count + 1;
                    end
                    flag = false;
                end
                
                %%%% resume the original order %%%%
                idmin = cellfun(@min, idclust);
                [~, idori] = sort(idmin);

                %%%% update the main parameters %%%%
                idclust = idclust(idori);
                imgall = imgall(idori);
                idclustfn{countfn} = idclust;
                imgallfn{countfn} = imgall;
                xfallfn{countfn} = xfall;
                ldallfn{countfn} = ldall;
                
                %%%% update inner loop variables %%%%
                mxint1 = zeros(pixh, pixw, length(imgall));
                for i = 1: length(imgall)
                    mxint1(:, :, i) = imgall{i};
                end
                countfn = countfn + 1;
            elseif n == 1
                break
            else
                scl = scl + istep;
                pixthres = scl * pixs;
                if pixthres > pixupper
                    %%%% loop exit: no need to further combine; use the current sub-optimal transform as the final ones %%%%
                    %%%%% first compute full matrix of the remaining images %%%%%
                    n = size(mxint1, 3);
                    smatrix = ones(n);
                    xfmatrix = cell(n);
                    ldmatrix = cell(n);
                    nmaxloop = 20;
%                     imgmatrix = cell(n);
                    mxint2 = mxint1;
                    parfor i = 1: n
                        imref = mxint1(:, :, i);
                        [scrtc, imgt, xft] = klt_ref_track(mxint2, imref, nmaxloop, maskc);
                        scrto = get_trans_score_ref(mxint2, imref, maskc);
                        mxtmp = mxint2;
                        xftmp = cell(1, n);
                        for j = 1: n
                            if scrtc(j) < scrto(j)
                                mxtmp(:, :, j) = imgt(:, :, j);
                                xftmp{j} = xft{j};
                            else
                                xftmp{j} = affine2d(diag(ones(1, 3)));
                            end
                        end
                        xfmatrix(i, :) = xftmp;
                        
                        stmp = zeros(1, n);
                        ldtmp = cell(1, n);
                        for j = 1: n
                            imcur = mxtmp(:, :, j);
                            [imgo, sx, sy] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d, maskc);
                            stmp(j) = get_trans_score_ref(imgo, imref, maskc);
                            ldtmp{j} = cellfun(@(x, y) cat(3, x, y), sx, sy, 'uniformoutput', false);
                        end
                        smatrix(i, :) = stmp;
                        ldmatrix(i, :) = ldtmp;
%                         disp(num2str(i))
                    end
                    
                    %%%%% solve TSP problem: <10 using greedy algorithm; >= 10 using NN approx %%%%%
                    [scrtsp, spath] = TSP(smatrix);
                    
                    %%%%% update final parameters %%%%%
                    xftmp = cell(1, n);
                    ldtmp = cell(1, n);
                    imgtmp = mxint1(:, :, spath);
                    idtmp = {1: n};
                    for i = 2: n
                        for j = i: -1: 2
                            xftmp{spath(i)} = [xftmp{spath(i)}, xfmatrix(spath(j - 1), spath(j))];
                            ldtmp{spath(i)} = [ldtmp{spath(i)}, ldmatrix(spath(j - 1), spath(j))];
                            imgtmp(:, :, i) = klt_warp(imgtmp(:, :, i), xftmp{spath(i)}{end});
                            for k = 1: length(ldtmp{spath(i)}{end})
                                imgtmp(:, :, i) = iminterpolate(imgtmp(:, :, i), ldtmp{spath(i)}{end}{k}(:, :, 1), ldtmp{spath(i)}{end}{k}(:, :, 2));
                            end
                        end
                    end
                    xfallfn{countfn} = xftmp;
                    ldallfn{countfn} = ldtmp;
                    idclustfn{countfn} = idtmp;
                    imgallfn{countfn} = {imgtmp};
                    break
                else
                    flagincrease = true;
                end
            end
        end
    
        %% compute transform matrix %%
        [xfuse, lduse, iduse] = logdemons_combine_layers(xfallfn, ldallfn, idclustfn);
        
        %%% final update the final sub-optimal transform %%%
        mxout = imgallfn{end}{1};
%         mxout = reshape(cell2mat(mxout(:)'), pixh, pixw, length(mxout));
    else
        mxout = mxin;
        xfuse{1} = {affine2d(diag(ones(1, 3)))};
        lduse{1} = {{zeros(pixh, pixw, 2)}};
        iduse = {1};
    end
end