function [avec, gamma_plus, gamma_minus] = calcAvec(new, dQ, W, lambda, active_set, M, positive)

[r,c] = find(active_set);
Mm = -M(r,r);


Mm=(Mm + Mm')/2;

% verify that there is no numerical instability 
eigMm = eig(Mm);
if any(eigMm < 0)
    min(eigMm)
    %error('The matrix Mm has negative eigenvalues')  
    flag = 1;
end


b = sign(W);
if new
    b(new) = sign(dQ(new));
end
b = b(active_set == 1);

avec = Mm\b;

if positive 
    if new 
        in = sum(active_set(1:new));
        if avec(in) <0
            new;
            %error('new component of a is negative')
            flag = 1;
        end
    end
end

    

one_vec = ones(size(W));

dQa = zeros(size(W));
for j=1:length(r)
    dQa = dQa + avec(j)*M(:, r(j));
end

gamma_plus = (lambda - dQ)./(one_vec + dQa);
gamma_minus = (lambda + dQ)./(one_vec - dQa);

end