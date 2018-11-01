function m = normalize_batch(filename, var, imx, imn, idbatch)
    m = matfile(filename, 'writable', true);
    nbatch = length(idbatch) - 1;
    for i = 1: nbatch
        eval(['tmp = m.', var, '(:, :, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ');'])
        tmp = (tmp - imn) / (imx - imn);
        eval(['m.', var, '(:, :, ', num2str(idbatch(i)), ': ', num2str(idbatch(i + 1) - 1), ') = tmp;'])
    end
end