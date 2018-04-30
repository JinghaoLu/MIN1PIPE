function [varargout] = optimizer_init(pars, optim)
%%% Jinghao Lu 03/02/2017 %%%

    layern = pars.layern;
    indim = pars.indim;
    hdim = pars.hdim;
    hdimt = hdim;
    outdim = pars.outdim;
    indimt = indim;

    if strcmp(optim, 'adam')
        moment.Um = {zeroMatrix([outdim, hdim], pars.gpu)};
        moment.Uv = {zeroMatrix([outdim, hdim], pars.gpu)};
        for l = 1: layern
            moment.Wm{l} = zeroMatrix([4 * hdimt, hdimt + indimt + 1], pars.gpu);
            moment.Wv{l} = zeroMatrix([4 * hdimt, hdimt + indimt + 1], pars.gpu);
            indimt = hdimt;
            hdimt = hdim;
        end
        varargout = {moment};
    elseif strcmp(optim, 'rmsprop')
        et.U = {zeroMatrix([outdim, hdim], pars.gpu)};
        for l = 1: layern
            et.W{l} = zeroMatrix([4 * hdimt, hdimt + indimt + 1], pars.gpu);
            indimt = hdimt;
            hdimt = hdim;
        end
        varargout = {et};
    end
end