function ithres = intensity_filter(maxall)
    tmp1 = maxall;
    tmp1 = tmp1(tmp1 > 1/2^8); %%% above 0 in terms of single precision %%%
    if ~isempty(tmp1)
        %         try
        %             nbins = round(length(tmp1) / 10);
        %             [tmp1, ctrs] = histcounts(tmp1, nbins);
        %             tmp1 = smooth(tmp1)';
        %             [tm, ~] = max(tmp1);
        %             [~, id1] = min(tmp1);
        %             id = find(tmp1 < 0.1 * tm);
        %             id = id(find(id < id1, 1));
        %             f1 = fit(double(ctrs(id1: id))', double(tmp1(id1: id))', 'poly1');
        %             y1 = abs(f1.p2 / f1.p1);
        %             id = find(tmp1 < 0.01 * tm);
        %             id = id(find(id > id1, 1));
        %             f2 = fit(double(ctrs(id: min(id + 100, length(tmp1))))', double(tmp1(id: min(id + 100, length(tmp1))))', 'poly1');
        %             y2 = ctrs(find(tmp1 < f2.p2, 1));
        %             %         ithres = (y1 + y2) / 2;
        %             ithres = y1;
        %         catch
        tmp1 = sort(tmp1);
        tmp2 = linspace(tmp1(1), 1 * tmp1(end), length(tmp1));
        tmp = tmp1(:) - 0.2 * tmp2(:);
        x = 1: length(tmp);
        sker = 2 * round(length(tmp) / 100) + 1;
        xq = [1 - sker: 0, x, length(tmp) + 1: length(tmp) + sker];
        tmpt = interp1(x, tmp, xq, 'linear', 'extrap');
        tmpg = smooth(diff(smooth(tmpt, sker)), sker);
        tmpg = tmpg(x);
        idthres = find(tmpg >= 0, 1);
        ithres = min(prctile(tmp1, 90), tmp1(idthres));
        %
        %         imgmax = feature2_comp(maxall, 0, 40, 1 / ithres);
        %         imgmaxt = imgmax(imgmax > 0);
        %         tmp1 = sort(imgmaxt);
        %         tmp2 = linspace(tmp1(1), 1 * tmp1(end), length(tmp1));
        %         tmp = tmp1(:) - tmp2(:);
        %         x = 1: length(tmp);
        %         sker = 2 * round(length(tmp) / 100) + 1;
        %         xq = [1 - sker: 0, x, length(tmp) + 1: length(tmp) + sker];
        %         tmpt = interp1(x, tmp, xq, 'linear', 'extrap');
        %         tmpg = smooth(diff(smooth(tmpt, sker)), sker);
        %         tmpg = tmpg(x);
        %         idthres = find(tmpg >= 0, 1);
        %         ithres = tmp1(idthres);
        %         end
    else
        ithres = 0;
    end
end
