function iduse = merge_unit(conmtx, datasmthf)
% [roiout, sigout, seedsout, Pout] = merge_unit merge based on graph
%   extracted rois
%   Jinghao Lu 06/10/2016

    %% generate graph clusters %%
    grf = graph(conmtx, 'upper');
    names = get_node_name(1: size(conmtx, 1));
    grf.Nodes.Name = names;
    bb = conncomp(grf); %%% connected subgraph %%%
    s = zeros(1, max(bb));
    for i = 1: max(bb)
        s(i) = sum(bb == i);
    end
    connuse = find(s > 1); %%% potential merging subgraphs %%%
    sgrft = [];
    for i = 1: length(connuse)
        sgrft = [sgrft, find(bb == connuse(i))];
    end
    sgrf = subgraph(grf, sgrft);
    iduse = setdiff(1: length(bb), sgrft);
    
    %% get merged profiles %%
    for i = 1: length(connuse)
        idcur = find(bb == connuse(i));
        
        %%%% get better seed id %%%%
        [~, id2use] = max(median(datasmthf(idcur, :), 2));
        iduse = [iduse, idcur(id2use)];
    end
end