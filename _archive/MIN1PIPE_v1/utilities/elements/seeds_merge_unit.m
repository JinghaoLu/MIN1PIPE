function iduset = seeds_merge_unit(conmtx, imax, iduse)
% [roiout, sigout, seedsout, Pout] = merge_unit merge based on graph
%   extracted rois
%   Jinghao Lu 06/10/2016

    %% build image graph %%
    %%% initialize image graph %%%
    img = imax;
    img = imgaussfilt(img, 1);
    imgt = 1- normalize(img);
    [pixh, pixw] = size(img);
    g = imageGraph([pixh, pixw], 8);
    gedges = g.Edges.EndNodes;
    gedges = [gedges; fliplr(gedges)];
    gwts = g.Edges.Weight;
    gwts = [gwts; gwts];
    edgetable = table(gedges, gwts, 'VariableNames', {'EndNodes', 'Weight'});
    g = digraph(edgetable, g.Nodes);
    
    %%% update image graph weights %%%
    ndi = g.Edges.EndNodes(:, 1);
    ndj = g.Edges.EndNodes(:, 2);
    aint = imgt(ndi);
    bint = imgt(ndj);
    wts = graph_cost(aint, bint);
    g.Edges.Weight = wts;

    %% build graph and find clusters %%
    grf = graph(conmtx, 'upper');
    names = get_node_name(1: size(conmtx, 1));
    grf.Nodes.Name = names;
    bb = conncomp(grf); %%% connected subgraph %%%
    s = zeros(1, max(bb));
    for i = 1: max(bb)
        s(i) = sum(bb == i);
    end
    connuse = find(s > 1); %%% potential merging subgraphs %%%
    idgrf = [];
    for i = 1: length(connuse)
        idgrf = [idgrf, find(bb == connuse(i))];
    end
    sgrf = subgraph(grf, idgrf);
    iduset = setdiff(1: length(bb), idgrf);
    
    %% get merged profiles %%
    parfor i = 1: length(connuse)
        idcur = find(bb == connuse(i));
        nodeid = get_node_name(idcur);
        sgrft = subgraph(sgrf, nodeid);
        idusep = iduse;
        imgp = img;
        
        while ~isempty(sgrft.Nodes)
            %%%% get better seed id %%%%
            [~, id2use] = max(imgp(iduse(idcur)));
            
            spt = idcur(id2use);
            sp = idusep(spt);
            cds = setdiff(idcur, spt);
            idt = [];
            for j = 1: length(cds)
                tp = iduse(cds(j));
                trc = shortestpath(g, sp, tp);
                tmp = imgp(trc);
                dtmp = diff(tmp);
                scr = sum(dtmp > 0);
                if scr == 0
                    idt = [idt, j];
                end
            end
            idt = [spt, cds(idt)];
            iduset = [iduset, idcur(id2use)];
            sgrft = rmnode(sgrft, get_node_name(idt));
            idcur = setdiff(idcur, idt);
        end
%         disp(num2str(i))
    end
end