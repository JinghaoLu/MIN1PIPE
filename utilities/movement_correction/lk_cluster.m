function [mxout, xfuse, iduse, smatrix, xfmatrix] = lk_cluster(mxin, pixs, scl, maskc)
% Full cluster Hierarchical LK method
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
    
    %% initialization %%
    countfn = 1;
    idclustfn = {};
    imgallfn = {};
    xfallfn = {};
    mxint1 = mxin;
    [pixh, pixw, nf] = size(mxin);
    pixthres = scl * pixs;
    
    if nf > 1
        %% main loop when there are still possibility of combination %%
        while 1
            
            %%% compute score and transform correspondance matrix %%%
            n = size(mxint1, 3);
            smatrix = ones(n);
            xfmatrix = cell(n);
            imgmatrix = cell(n);
            
            mxint2 = mxint1;
            
            parfor i = 1: n
                imref = mxint1(:, :, i);
                ncur = n;
                temp = zeros(1, ncur);
                xftemp = cell(1, ncur);
                imgtemp = cell(1, ncur);
                
                [PNt, imgt, scrto] = lk_ref_track(mxint2, imref, maskc);
                scrtc = get_trans_score_ref(mxint2, imref, maskc);
                ids = scrto < scrtc;
                PNt(:, ids) = 0;
                imgt(:, :, ids) = mxint2(:, :, ids);
                scrtc(ids) = scrto(ids);
                for j = 1: ncur
                    xft = eye(3);
                    xft(3, 1: 2) = -fliplr(squeeze(PNt(:, j))');
                    xftemp{j} = affine2d(xft);
                    imgtemp{j} = imgt(:, :, j);
                    temp(j) = scrtc(j);
                end
                smatrix(i, :) = temp;
                xfmatrix(i, :) = xftemp;
                imgmatrix(i, :) = imgtemp;
%                 disp(num2str(i))
            end
            
            
            %%% compute connected subgraph components %%%
            aa = smatrix < pixthres;
%             aa = xor(aa, diag(ones(1, n)));
            aa(diag(ones(1, n)) > 0) = false;
%             aa = smatrix < 0.8;
            sgraph = graph(aa, 'upper');
            names = get_node_name(1: n);
            sgraph.Nodes.Name = names;

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
            if ~isempty(sgraph.Edges)
                count = 1;
                flag = true;
                
                %%%% continue if there are still uncombined subgraph clusters %%%%
                while flag || ~all(cellfun(@(x) size(x, 3), imgall) == 1)
                    if flag
                        sgrapht = sgraph;
                    else %%%% find the largest subgraph in the remaining graph %%%%
                        idsub = find(cellfun(@(x) size(x, 3), imgall) ~= 1, 1);
                        sgrapht = subgraph(sgraph, idclust{idsub});
                        ntemp = length(idclust);
                        idclust = idclust(setdiff(1: ntemp, idsub));
                        imgall = imgall(setdiff(1: ntemp, idsub));
                    end
                    
                    %%%% for the current subgraph cluster, combine till no edges %%%%
                    while ~isempty(sgrapht.Nodes)
                        bb = conncomp(sgrapht);
                        
                        %%%%% compute largest subgraph %%%%%
                        s = zeros(1, max(bb));
                        for i = 1: max(bb)
                            s(i) = sum(bb == i);
                        end
                        [~, connuse] = max(s);
                        idcu = setdiff(1: max(bb), connuse);
                        idcurall = cellfun(@str2double, sgrapht.Nodes.Name);
                        for i = 1: length(idcu)
                            idt = str2double(sgrapht.Nodes.Name(bb == idcu(i)));
                            idclust = [idclust; {idt(:)}];
                            imgall = [imgall; {mxint1(:, :, idt)}];
                        end
                        ids = str2double(sgrapht.Nodes.Name(bb == connuse));
                        idrmv = setdiff(idcurall, ids);
                        namermv = get_node_name(idrmv);
                        sgrapht = rmnode(sgrapht, namermv);

                        %%%%% find the hub (with largest degree) %%%%%
                        centra = centrality(sgrapht, 'degree');
                        source = find(centra == max(centra));
                        source = source(1);
                        idsource = ids(source);
                        
                        %%%%% combine the hub with all its neighbors %%%%%
                        nbs = neighbors(sgrapht, num2str(idsource));
                        idnbs = cellfun(@str2double, nbs);
                        imgtemp = reshape(cell2mat(imgmatrix(idsource, [idsource; idnbs(:)])), pixh, pixw, length(idnbs) + 1);
%                         imgall = [imgall; {max(imgtemp, [], 3)}];
                        imgall = [imgall; {mean(imgtemp, 3)}];
                        idclust = [idclust; {[idsource; idnbs]}];
                        for i = 1: length(idnbs)
                            xfall{idnbs(i)} = [xfall{idnbs(i)}, xfmatrix(idsource, idnbs(i))];
                        end
                        
                        %%%%% update parameters within the current round of graph %%%
                        idt = [num2str(idsource); nbs];
                        
                        %%%%% update the current graph %%%%%
                        sgrapht = rmnode(sgrapht, idt);
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
                countfn = countfn + 1;
            else
                break
            end
        end
    
        %% compute transform matrix %%
        [xfuse, iduse] = lk_combine_layers(xfallfn, idclustfn);
        
        %%% final update the final sub-optimal transform %%%
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