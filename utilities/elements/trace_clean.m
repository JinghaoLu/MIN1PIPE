function sigout = trace_clean(sigin, Fs, tflag)
    [n1, n2] = size(sigin);
    sigout = sigin;
    parfor i = 1: n1
        tmp = sigin(i, :);
        tmp(1: round(Fs / 4)) = linspace(prctile(tmp, 1), tmp(round(Fs / 4)), round(Fs / 4));
        
        %%% 1. get smoothed signal %%%
        tmpt = wdenoise(double(tmp), 5, 'Wavelet', 'sym8', 'DenoisingMethod', 'sure', 'noiseestimate', 'leveldependent');
        mn = median(tmp);
                
        if tflag == 1
            
            %%% 2.1 lower "envelope" %%%
            tt = find_valleys(tmpt, tflag);
        
            %%% 2.2 get smoothed baseline %%%
            kl = min(Fs * 200, length(tt) / 2);
            ttt = conv(tt, gausswin(kl) / sum(gausswin(kl)), 'valid');
            ttt = interp1(ceil(kl / 2): ceil(kl / 2) + length(ttt) - 1, ttt, 1: length(tt), 'linear', 'extrap');
            tmp1 = tmp - ttt;
        else
            %%% 2.1 RDP algorithm to get key turning points %%%
            mnt = tmp - tmpt;
            mnt = prctile(mnt, 95) - prctile(mnt, 5);
            tt = rdp([1: n2; tmpt], mnt);
            tt = interp1(tt(1, :), tt(2, :), 1: n2, 'linear', 'extrap');
            ttt = find_valleys(tt, tflag);
%             tmp1 = max(0, tmp - ttt);
            tmp1 = tmp - ttt;
        end
        
        %%% 3. correct baseline %%%
        mn2 = median(tmp1);
        offset = mn - mn2;
        
        %%% 4. suppress noise %%%
        sigout(i, :) = max(0, tmp1 + offset);
    end
    
%     %%% 1. calculate the baseline %%%
%     [n1, n2] = size(sigin);
%     mn = zeros(n1, 1);
%     parfor i = 1: n1
%         tmp = sigin(i, :);
%         edges = linspace(min(tmp), max(tmp), 101);
%         n = ksdensity(tmp, edges);
%         [~, id] = max(n);
%         mn(i) = edges(id);
%     end
% 
%     %%% 2. calculate trend by down- & upsample %%%
%     tstep = round(Fs / 2);
%     gkernel = round(Fs * 5 - 1);
%     datat = sigin(:, 1: tstep: end);
%     thresst = medfilt2(datat, [1, gkernel]);
%     thress = sigin;
%     parfor i = 1: n1
%         thress(i, :) = interp1(1: tstep: n2, thresst(i,:), 1: n2, 'linear', 'extrap');
%     end
%     data2 = sigin - thress;
% 
%     %%% 4. suppress noisy periods %%%
%     parfor i = 1: n1
%         tmp = data2(i, :);
%         tmpt = sgolayfilt(double(tmp), 1, round(Fs * 20 - 1));
%         data2(i, :) = sigmoid(tmp, 1000, tmpt) .* tmp;
%     end
% 
%     %%% 5. suppress noise %%%
%     parfor i = 1: n1
%         tmp = data2(i, :);
%         data2(i, :) = sgolayfilt(double(tmp), 1, round(Fs / 2 - 1));
%     end
%     sigout = data2 + mn;
end

function tt = find_valleys(tmpt, tflag)
    [pks, x] = findpeaks(-tmpt);

    mn = median(tmpt);
    if tflag == 1
        md = std(tmpt);
        id = pks >= -mn - 2 * md;
    else
        id = pks >= -mn;
    end

    x = x(id);
    x = [1, x, length(tmpt)];
    tt = interp1(x, tmpt(x), 1: length(tmpt), 'linear');
end