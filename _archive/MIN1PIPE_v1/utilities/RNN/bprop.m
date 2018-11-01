function gradt = bprop(lstmc, Wt, Ut, P, lbin, pars)
%%% Jinghao Lu 02/25/2017 %%%

    %%% single step gradient init, through time %%%
    [layern, tlen] = size(lstmc);
    bsize = length(lbin);
    lblact = sub2ind(size(P), lbin' + 1, 1: bsize);
    P(lblact) = P(lblact) - 1; %%% predict error %%%
    
    gradt.dLdh = cell(1, layern);
    for l = 1: layern
        if l == layern
            gradt.dLdh{l} = Ut{1}' * P; % dimension=(hidden*class)*(class*N)=hidden*Num_of_example
        else
            gradt.dLdh{l} = zeroMatrix(size(lstmc{l, 1}.ht), pars.gpu);
        end
    end
    gradt.dLdU = {P * lstmc{end, end}.ht'}; %dimension=(class*N)*(N*hidden)=class*hidden;
    
    gradt.dLdW = cell(1, layern);
    gradt.dLdc = cell(1, layern);
    for l = 1: layern
        gradt.dLdW{l} = zeroMatrix(size(Wt{l}), pars.gpu);
        gradt.dLdc{l} = zeroMatrix(size(lstmc{l, 1}.ct), pars.gpu);
    end
    
    %%% propagation %%%
    for i = tlen: -1: 1
%         disp(num2str(i))
        for l = layern: -1: 1
            if i > 1
                ct1 = lstmc{l, i - 1}.ct;
            else
                ct1 = zeroMatrix(size(lstmc{l, i}.ct), pars.gpu);
            end
        
            %%% step t %%%
            gradt.dLdc{l} = arrayfun(@dTanhSum, gradt.dLdc{l}, lstmc{l, i}.tct, lstmc{l, i}.ogate, gradt.dLdh{l});
            dLdo = arrayfun(@dSigmf, lstmc{l, i}.ogate, lstmc{l, i}.tct, gradt.dLdh{l});
            dLdi = arrayfun(@dSigmf, lstmc{l, i}.igate, lstmc{l, i}.c_hat, gradt.dLdc{l});
            dLdf = arrayfun(@dSigmf, lstmc{l, i}.fgate, ct1, gradt.dLdc{l});
            dLdchat = arrayfun(@dTanh, lstmc{l, i}.c_hat, lstmc{l, i}.igate, gradt.dLdc{l});
            dLdz = [dLdi; dLdf; dLdo; dLdchat];
            dLdW = dLdz * lstmc{l, i}.incat';
            gradt.dLdW{l} = gradt.dLdW{l} + dLdW; %%% sum up dW %%%
            
            %%% update for step t - 1 %%%
            dLdI = Wt{l}' * dLdz;
            indim = size(Wt{l}, 2) - size(Wt{l}, 1) / 4 - 1;
            gradt.dLdh{l} = dLdI(indim + 1: end - 1, :);
            gradt.dLdc{l} = lstmc{l, i}.fgate .* gradt.dLdc{l};
            
            if l > 1
                gradt.dLdh{l - 1} = gradt.dLdh{l - 1} + dLdI(1: indim, :);
            end
        end
    end
    
    for l = 1: layern
        gradt.dLdW{l} = gradt.dLdW{l} / bsize;
    end
    gradt.dLdU{1} = gradt.dLdU{1} / bsize;
end