function [sig, spk] = pure_refine_sig(sigfn, options)
    [nn, nf] = size(sigfn);
    sig = sigfn;
    spk = sigfn;
    parfor i = 1: nn
        [cc, cb, c1, gn, sn, spkt] = constrained_foopsi(sigfn(i, :), [], [], [], [], options);
        gd = max(roots([1, -gn']));  % decay time constant for initial concentration
        gd_vec = gd .^ ((0: nf - 1));
        sig(i, :) = full(cc(:)' + cb + c1 * gd_vec);
        spk(i, :) = spkt(:)';
    end
end