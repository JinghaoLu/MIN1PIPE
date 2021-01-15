function [siguse, roiuse, seedsuse, buse, bfuse, datasmthuse, cutoffuse, pkcutoffuse] = iter_seeds_select(m_in, mask, Params, P, options, flag)
    %%% variables initialization %%%
    siguse = [];
    roiuse = [];
%     buse = [];
%     bfuse = [];
    seedsuse = [];
    datasmthuse = [];
    cutoffuse = [];
    pkcutoffuse = [];
    
    fnamenew = [m_in.Properties.Source(1: end - 4), '_res.mat'];
    status = copyfile(m_in.Properties.Source, fnamenew, 'f');
    m_tmp = matfile(fnamenew, 'writable', true);

    
    for i = 1: 1
        %%% select pixel %%%
        sz = Params.neuron_size;
        Fsi_new = Params.Fsi_new;
        
        if flag == 1
            sigthres = Params.pix_select_sigthres;
            corrthres = Params.pix_select_corrthres;
            [roi, sig, bg, bgf, seeds, datasmth0, cutoff0, pkcutoff0] = pix_select(m_tmp, mask, sz, Fsi_new, sigthres, corrthres);
        else
            [roi, sig, seeds, bg, bgf, datasmth0, cutoff0, pkcutoff0] = manual_seeds_select(m_tmp, Fsi_new, sz);
        end
                
        %%% refine roi %%%
        noise = P.sn;
        ispara = Params.refine_roi_ispara;
        [roirf, sigupdt, seedsupdt, datasmthf1, cutofff1, pkcutofff1] = refine_roi(m_tmp, sig, bgf, roi, seeds, noise, datasmth0, cutoff0, pkcutoff0, ispara);
        
        %%% refine sig %%%
        p = 0; %%% no ar model used %%%
        [sigrf, bgfrf, roirf, seedsupdt, datasmthf1, cutofff1, pkcutofff1, Puse] = refine_sig(m_tmp, roirf, bg(:), sigupdt, bgf, seedsupdt, datasmthf1, cutofff1, pkcutofff1, p, options);
        
        %%% update reg_post.mat %%%
        [d1, d2, T] = size(m_tmp, 'reg');
        nsize = d1 * d2 * T * 8 * 2; %%% size of double %%%
        nbatch = batch_compute(nsize);
        ebatch = ceil(T / nbatch);
        idbatch = [1: ebatch: T, T + 1];
        nbatch = length(idbatch) - 1;
%         YrA = zeros(T, K + 1);
        for ib = 1: nbatch
            idtmp = idbatch(ib): idbatch(ib + 1) - 1;
            ytmp = m_tmp.reg(1: d1, 1: d2, idtmp);
            ytmp = reshape(ytmp, d1 * d2, []);
            tominus = roirf * sigrf(:, idtmp);
            yrt = ytmp - tominus;
            m_tmp.reg(1: d1, 1: d2, idtmp) = reshape(yrt, d1, d2, []);
%             YrA(idtmp, :) = yrt' * roirf ./ sum(roirf, 1);
        end
        
        %%% variable collect %%%
        siguse = [siguse; sigrf];
        roiuse = [roiuse, roirf];
%         buse = bgr;
%         bfuse = bgrf;
        seedsuse = [seedsuse; seedsupdt];
        datasmthuse = [datasmthuse; datasmthf1];
        cutoffuse = [cutoffuse; cutofff1];
        pkcutoffuse = [pkcutoffuse; pkcutofff1];
        
        if flag ~= 1
            break
        end
    end
    
    %%% delete empty ROIs %%%
    nseed = size(siguse, 1);
    inddel = union(find(sum(roiuse, 1) == 0), find(sum(siguse, 2) == 0));
    roiuse = roiuse(:, setdiff(1: nseed, inddel));
    siguse = siguse(setdiff(1: nseed, inddel), :);
    datasmthuse = datasmthuse(setdiff(1: nseed, inddel), :);
    cutoffuse = cutoffuse(setdiff(1: nseed, inddel));
    pkcutoffuse = pkcutoffuse(setdiff(1: nseed, inddel));
    seedsuse = seedsuse(setdiff(1: nseed, inddel));
    
    %%% compute backgrounds %%%
    [buse, bfuse] = bg_update(m_in, roiuse, siguse);
    
    %%% delete temp files %%%
    delete(fnamenew)
end






