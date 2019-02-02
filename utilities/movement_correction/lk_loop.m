function [mxout, xfuse, iduse] = lk_loop(mxin, pixs, scl, maskc)
% Hierarchical LK method of the given image cluster
%   Jinghao Lu, 02/02/2018

    if nargin < 2 || isempty(pixs)
        [pixh, pixw, ~] = size(mxall);
        pixs = min(pixh, pixw);
    end
    
    if nargin < 3 || isempty(scl)
        defpar = default_parameters;
        scl = defpar.mc_scl;
    end
    
    if nargin < 4 || isempty(maskc)
        maskc = true(size(mxin, 1), size(mxin, 2));
    end
    
    %% LK track and combine loop %%
    countfn = 1;
    idclustfn = {};
    imgallfn = {};
    xfallfn = {};
    mxint1 = mxin;
    mxintt = mxin;
    [pixh, pixw, nf] = size(mxin);
    pixthres = scl * pixs;
    idrun = 1: nf - 1;
    xfmatrix = cell(1, nf);
    for i = 1: nf
        xfmatrix{i} = {};
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
                    
                    %%%% track current two neighboring "frames" %%%%
                    [PNt, imgt, scrto] = lk_ref_track(imcur, imref, maskc);
                    scrtc = get_trans_score_ref(imgt, imref, maskc);
                    
                    %%%% track real neighboring frames %%%%
                    if countfn > 1
                        [imreft, imcurt] = find_frame(mxintt, idclustfn, ii);
                        
                        %%%%% the other track %%%%%
                        [PNt1, imgt1, scrto1] = lk_ref_track(imcurt, imreft, maskc);
                        scrtc1 = get_trans_score_ref(imgt1, imreft, maskc);
                        
                        %%%%% get smallest score %%%%%
                        if scrto1 < scrto
                            scrto = scrto1;
                        end
                        
                        if scrtc1 < scrtc
                            PNt = PNt1;
                            imgt = imgt1;
                            scrtc = scrtc1;
                        end
                    end
                    
                    if scrto < scrtc
                        PNt = [0; 0];
                        imgt = imcur;
                        scrtc = scrto;
                    end
                    xft = eye(3);
                    xft(3, 1: 2) = -fliplr(squeeze(PNt)');
                    xfmatrix{ii + 1} = affine2d(xft);
                    scrs(ii) = scrtc;
                    imgmatrix(:, :, ii + 1) = imgt;
                end
%                 disp(num2str(ii))
            end
            
            %%% compute connected subgraph components %%%
            aa = scrs < pixthres;
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
            for i = 1: n
                xfall{i} = {};
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
                
                xfmatrix = xfmatrix(setdiff(1: n, idtargets));
                imgmatrix = imgmatrix(:, :, setdiff(1: n, idtargets));
                scrs = scrs(setdiff(1: n - 1, idsources));
                
                mxintt = temporary_warp(xfallfn, idclustfn, mxin);
                countfn = countfn + 1;
            else
                break
            end
        end
    
        %% compute transform matrix %%
        [xfuse, iduse] = lk_combine_layers(xfallfn, idclustfn);
                
        %%% generate mxout %%%
        if ~isempty(imgallfn)
            mxout = imgallfn{end};
            mxout = reshape(cell2mat(mxout(:)'), pixh, pixw, length(mxout));
        else
            mxout = mxin;
            for i = 1: nf
                iduse{i} = i;
                xfuse{i} = {affine2d(diag(ones(1, 3)))};
            end
        end
        
    else
        mxout = mxin;
        xfuse{1} = {affine2d(diag(ones(1, 3)))};
        iduse = {1};
    end
end