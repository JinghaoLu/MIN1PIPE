function cost = cross_entropy_loss(P, lbin)
%%% Jinghao Lu 02/25/2017 %%%

    bsize = length(lbin);
    lblact = sub2ind(size(P), lbin' + 1, 1: bsize);
    cost = sum(-log(P(lblact))) / bsize;
end