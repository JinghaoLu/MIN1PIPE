function [mxout, xfuse, lduse, iduse] = logdemons_loop(mxin, pixs, scl, sigma_x, sigma_f, sigma_d, maskc)
% serial loop of logdemons to register similar neighboring frames
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
    
    countfn = 1;
    idclustfn = {};
    imgallfn = {};
    xfallfn = {};
    ldallfn = {};
    mxint1 = mxin;
    [pixh, pixw, nf] = size(mxin);
    pixthres = scl * pixs;
    idrun = 1: nf - 1;
    xfmatrix = cell(1, nf);
    ldmatrix = cell(1, nf);
    for i = 1: nf
        xfmatrix{i} = {};
        ldmatrix{i} = {};
    end
    imgmatrix = mxint1;
    scrs = zeros(1, nf - 1);
    
    if nf > 1
        %% main loop when there are still possibility of combination %%
        while 1
            
            %%% compute score and transform correspondance matrix %%%
            n = size(mxint1, 3);
            mxint2 = mxint1;
            mxint3 = mxint1;
            parfor ii = 1: n - 1
                if ismember(ii, idrun)
                    imref = mxint2(:, :, ii);
                    imcur = mxint3(:, :, ii + 1);
                    [scrtc, imgt, xft] = klt_ref_track(imcur, imref, [], maskc);
                    scrto = get_trans_score_ref(imcur, imref, maskc);
                    if scrtc < scrto
                        imcur = imgt;
                        xfmatrix{ii + 1} = xft{1};
                    else
                        xfmatrix{ii + 1} = affine2d(diag(ones(1, 3)));
                    end
                    [imgo, sx, sy] = logdemons_unit(imref, imcur, pixs, scl, sigma_x, sigma_f, sigma_d, maskc);
                    imgmatrix(:, :, ii + 1) = imgo;
                    scrs(ii) = get_trans_score_ref(imgo, imref, maskc);
                    temp = cellfun(@(x, y) cat(3, x, y), sx, sy, 'uniformoutput', false);
                    ldmatrix{ii + 1} = temp;
                end
%                 disp(num2str(ii))
            end
            
            %%% compute connected subgraph components %%%
            aa = scrs < pixthres;
%             aa = scrs < 1;
            aa = diag(aa, 1);
            sgraph = digraph(aa, 'OmitSelfLoops');
            names = get_node_name(1: n);
            sgraph.Nodes.Name = names;
            sgraph1 = graph(aa, 'upper');
            sgraph1.Nodes.Name = names;
            
            %%% decompose and combine subgraphs within the main cluster %%%
            idclust = {};
            imgall = {};
            idsources = [];
            idtargets = [];
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
                        imgtemp = cat(3, mxint1(:, :, idsource), imgmatrix(:, :, idnbs));
%                         imgall = [imgall; {max(imgtemp, [], 3)}];
                        imgall = [imgall; {mean(imgtemp, 3)}];
                        idclust = [idclust; {[idsource; idnbs]}];
                        if ~isempty(idnbs)
                            idsources = [idsources, idsource];
                            idtargets = [idtargets, idnbs];
                        end
                        for i = 1: length(idnbs)
                            ldall{idnbs(i)} = [ldall{idnbs(i)}, ldmatrix(idnbs(i))];
                            xfall{idnbs(i)} = [xfall{idnbs(i)}, xfmatrix(idnbs(i))];
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
                
                idrun = zeros(1, length(idsources));
                for i = 1: length(idsources)
                    idrun(i) = find(cellfun(@(x) any(x == idsources(i)), idclust));
                end
                idrun = union(idrun, idrun - 1);
                idrun = idrun(idrun >= 1 & idrun <= length(imgall));
                
                ldmatrix = ldmatrix(setdiff(1: n, idtargets));
                xfmatrix = xfmatrix(setdiff(1: n, idtargets));
                imgmatrix = imgmatrix(:, :, setdiff(1: n, idtargets));
                scrs = scrs(setdiff(1: n - 1, idsources));
                countfn = countfn + 1;
            else
                break
            end
        end
    
        %% compute transform matrix %%
        [xfuse, lduse, iduse] = logdemons_combine_layers(xfallfn, ldallfn, idclustfn);
                
        %%% generate mxout %%%
        if ~isempty(imgallfn)
            mxout = imgallfn{end};
            mxout = reshape(cell2mat(mxout(:)'), pixh, pixw, length(mxout));
        else
            mxout = mxin;
            for i = 1: nf
                iduse{i} = i;
                xfuse{i} = {affine2d(diag(ones(1, 3)))};
                lduse{i} = {{zeros(pixh, pixw, 2)}};
            end
        end
        
    else
        mxout = mxin;
        xfuse{1} = {affine2d(diag(ones(1, 3)))};
        lduse{1} = {{zeros(pixh, pixw, 2)}};
        iduse = {1};
    end
end