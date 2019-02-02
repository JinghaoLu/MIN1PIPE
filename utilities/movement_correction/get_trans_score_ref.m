function acorr = get_trans_score_ref(Y, imref, maskc)
% get transform score relative to a reference frame
%   Jinghao Lu, 02/20/2018

    if nargin < 3 || isempty(maskc)
        maskc = true(size(imref));
    end
    
    mq = 0.01;
    biderr = 2;
    [~, ~, nframes] = size(Y);
    acorr = zeros(1, nframes);
    if nframes == 1
        for i = 1: nframes
            img_old = imref;
            img = Y(:, :, i);
            d = klt2(normalize(img_old), normalize(img), biderr, mq, [], [], maskc);
            if ~isempty(d)
                temp = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
%                 acorr(i) = max(1, exp((-size(d, 1) + 10) / 2)) * temp;
                acorr(i) = temp;
            else
                acorr(i) = 100; %%% a large score %%%
            end
        end
    else
        parfor i = 1: nframes
            img_old = imref;
            img = Y(:, :, i);
            d = klt2(normalize(img_old), normalize(img), biderr, mq, [], [], maskc);
            if ~isempty(d)
                temp = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
%                 acorr(i) = max(1, exp((-size(d, 1) + 10) / 2)) * temp;
                acorr(i) = temp;
            else
                acorr(i) = 100; %%% a large score %%%
            end
        end
    end
end