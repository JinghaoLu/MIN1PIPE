function acorr = get_trans_score_ref(Y, imref)
% get transform score relative to a reference frame
%   Jinghao Lu, 02/20/2018

    mq = 0.04;
    biderr = 2;
    [~, ~, nframes] = size(Y);
    acorr = zeros(1, nframes);
    if nframes == 1
        for i = 1: nframes
            img_old = imref;
            img = Y(:, :, i);
            d = klt2(img_old, img, biderr, mq);
            if ~isempty(d)
                temp = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
%                 acorr(i) = max(1, exp((-size(d, 1) + 10) / 2)) * temp;
                acorr(i) = temp;
            else
                acorr(i) = 10; %%% a large score %%%
            end
        end
    else
        parfor i = 1: nframes
            img_old = imref;
            img = Y(:, :, i);
            d = klt2(img_old, img, biderr, mq);
            if ~isempty(d)
                temp = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
%                 acorr(i) = max(1, exp((-size(d, 1) + 10) / 2)) * temp;
                acorr(i) = temp;
            else
                acorr(i) = 10; %%% a large score %%%
            end
        end
    end
end