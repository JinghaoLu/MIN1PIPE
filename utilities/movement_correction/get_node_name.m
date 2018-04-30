function names = get_node_name(ids)
% get node name as string
%   Jinghao Lu, 09/01/2017

    n = length(ids);
    names = cell(n, 1);
    for i = 1: n
        names{i} = num2str(ids(i));
    end
end