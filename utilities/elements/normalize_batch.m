function m = normalize_batch(filename, var, imx, imn, idbatch, ip)
    if nargin < 6 || isempty(ip)
        ip = 3;
    end
    
    m = matfile(filename, 'writable', true);
    nbatch = length(idbatch) - 1;
    for i = 1: nbatch
        switch ip
            case 2
                eval(['tmp = m.', var, '(:, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ', :);'])
                tmp = (tmp - imn) / (imx - imn);
                eval(['m.', var, '(:, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ', :) = tmp;'])
            case 3
                eval(['tmp = m.', var, '(:, :, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ');'])
                tmp = (tmp - imn) / (imx - imn);
                eval(['m.', var, '(:, :, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ') = tmp;'])
        end
    end
end