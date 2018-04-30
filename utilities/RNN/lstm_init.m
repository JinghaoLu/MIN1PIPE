function [Wt, Ut, ht1, ct1] = lstm_init(pars)
%%% Jinghao Lu 02/28/2017 %%%

    %%% basic structural pars %%%
    layern = pars.layern;
    indim = pars.indim;
    hdim = pars.hdim;
    outdim = pars.outdim;
    bsize = pars.bsize;
    
    Ut = {initMatrix(outdim, hdim, pars.gpu)}; %%% key weight %%%
    Wt = cell(1, layern);
    ht1 = cell(1, layern);
    ct1 = cell(1, layern);
    hdimt = hdim;
    indimt = indim;
    for l = 1: layern
        Wt{l} = initMatrix(4 * hdimt, hdimt + indimt + 1, pars.gpu); %%% key weight including bias %%%
        ht1{l} = zeroMatrix([hdimt, bsize], pars.gpu);
        ct1{l} = zeroMatrix([hdimt, bsize], pars.gpu);
        indimt = hdimt;
        hdimt = hdim;
    end
end