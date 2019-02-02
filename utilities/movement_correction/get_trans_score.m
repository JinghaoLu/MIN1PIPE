function [acorr, duse] = get_trans_score(Y, flag, ispara, isdisp, mq, maskc)
% compute translation score with KLT tracker: mean displacement of matched
% points
%   Jinghao Lu, 01/30/2018

    %% initialization %%
    %%% initialize parameters %%%
    if nargin < 2 || isempty(flag)
        flag = 1;
    end
    
    if nargin < 3 || isempty(ispara)
        ispara = 1;
    end
    
    if nargin < 4 || isempty(isdisp)
        isdisp = 0;
    end
    
    if nargin < 5 || isempty(mq)
        mq = 0.01;
    end
    
    if nargin < 6 || isempty(maskc)
        maskc = true(size(Y, 1), size(Y, 2));
    end
    
    %% compute score based on the flag signal %%
    regt = Y;
%     mq = 0.04;
%     mq = 0.001;
    biderr = 2;
    [~, ~, nframes] = size(Y);
    acorr = zeros(1, nframes - 1);
    duse = zeros(nframes - 1, 2);
    if nframes == 2 || ~ispara
        for i = 2: nframes
            try
                img_old = regt(:, :, i - 1);
                img = Y(:, :, i);
                d = klt2(normalize(img_old), normalize(img), biderr, mq, [], [], maskc);
                duse(i - 1, :) = mean(d, 1);
                if ~isempty(d)
                    if flag == 1
                        acorr(i - 1) = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
                    elseif flag == 2
                        acorr(i - 1) = prctile(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2), 95);
                    end
                else
                    acorr(i - 1) = 100; %%% a large score %%%
                end
            catch
                acorr(i - 1) = NaN;
            end
            if isdisp
                if mod(i, round(nframes / 10)) == 0
                    disp(['computed score of frame #', num2str(i), '/', num2str(nframes)])
                end
            end
        end
    else
        parfor i = 2: nframes
            try
                img_old = regt(:, :, i - 1);
                img = Y(:, :, i);
                d = klt2(normalize(img_old), normalize(img), biderr, mq, [], [], maskc);
                duse(i - 1, :) = mean(d, 1);
                if ~isempty(d)
                    if flag == 1
                        acorr(i - 1) = mean(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2));
                    else
                        acorr(i - 1) = prctile(sqrt(d(:, 1) .^ 2 + d(:, 2) .^ 2), 95);
                    end
                else
                    acorr(i - 1) = 100; %%% a large score %%%
                end
            catch
                acorr(i - 1) = NaN;
            end
            if isdisp
                if mod(i, round(nframes / 10)) == 0
                    disp(['computed score of frame #', num2str(i), '/', num2str(nframes)])
                end
            end
        end
    end
    
    if any(isnan(acorr))
        ids = find(~isnan(acorr));
        idss = find(isnan(acorr));
        acorr(idss) = interp1(ids, acorr(ids), idss, 'pchip', 'extrap');
    end
end