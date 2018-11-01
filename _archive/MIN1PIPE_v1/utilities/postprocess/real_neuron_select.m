function real_neuron_select
% post process: semi-auto false positive removal
%   Jinghao Lu, 10/02/2017

    %% initialization %%
    %%% select processed file %%%
    [file_name, path_name] = uigetfile('*.mat');
    file_id = file_name(1: end - 4);
    
    %%% load files %%%
    load([path_name, file_name])
    roi = roifn;
    sig = sigfn;
    seeds = seedsfn;
    bg = bgfn;
    bgf = bgffn;
    roir = roifnr;
    sigr = sigfnr;
    [pixh, pixw] = size(imax);

    %%% initialize parameters %%%
    nr = 10;
    nn = size(sig, 1);
    nc = ceil(nn / nr);
    sig = normalize(sig);
    for i = 1: nn
        sig(i, :) = normalize(sig(i, :));
    end
%     tp = reshape(max(roi * sig, [], 2), pixh, pixw);
    tp = imax;
    
    %%% generate cursor judge region with higher resolution %%%
    rt = full(roi);
    for i = 1: nn
        tt = rt(:, i);
        tt = imgaussfilt(reshape(tt, pixh, pixw), 3);
        rt(:, i) = tt(:) > max(tt(:)) .* 0.5; 
    end
    
    roit = reshape(max((1: size(roi, 2)) .* rt, [], 2), pixh, pixw);
    global id_mouse notuse;
    id_mouse = 1;
    notuse = [];
    flag = 0;
    
    %% draw each 10 cells %%
    for i = 1: nc
        figure(1)
        clf
        pause(0.1)
        suptitle(['Section #', num2str(i), '/', num2str(nc)])
%         set(gcf, 'units', 'normalized', 'outerposition', [0 0 0.5 1])

        %%% draw roi selection pool %%%
        hroi = subplot(nr, 6, [6 * (nr - 6) + 3, 6 * (nr - 2)], 'align');
        set(hroi, 'Units', 'Pixels')
        posroi = get(hroi, 'Position');
        imagesc(tp)
        
        %%% draw single trace demo pool %%%
        hroi_wave = subplot(nr, 6, [6 * (nr - 2) + 4, 6 * nr], 'align');
        set(hroi_wave, 'Units', 'Pixels')
        axis off
        
        %%% draw single roi demo pool %%%
        sroi = subplot(nr, 6, [6 * (nr - 2) + 3, 6 * (nr - 1) + 3], 'align');
        set(sroi, 'Units', 'Pixels')
        posroiwave = get(hroi_wave, 'Position');
        axis off
        
        %%% draw not used roi pool %%%
        hsroi = subplot(nr, 6, [3, 6 * (nr - 6)], 'align');
        set(sroi, 'Units', 'Pixels')
        possroi = get(hsroi, 'Position');
        imagesc(tp)   
        hold on
        for j = 1: length(notuse)
            [y, x] = ind2sub([pixh, pixw], seeds(notuse(j)));
            plot(x, y, '.r')
        end
        hold off
        
        %%% plot 10 traces %%%
        h = [];
        posf = [];
        for j = 1: min(nr, nn - (i - 1) * nr)
            h(j) = subplot(nr, 6, [6 * (j - 1) + 1, 6 * (j - 1) + 2], 'align');
            set(gca, 'Units', 'Pixels')
            posf(j, :) = get(h(j), 'Position');
%             plot(all_wave_lpf((i - 1) * nr + j, :), 'k')
            plot(sig((i - 1) * nr + j, :), 'color', [0.8,0.8,0.8])
            axis tight
            ylim([0, 1])
            axis off
        end
        
        %%% set mouse move function & mouse click function %%%
        set(gcf, 'WindowButtonMotionFcn', @(object, eventdata) mouseMove(object, eventdata, tp, sig, roit, roi, sroi, hroi, hroi_wave, h, posf, posroi, i, nn, nr, pixh, pixw));
        set(gcf, 'WindowButtonDownFcn', @(object, eventdata) mouseClick(object, eventdata, tp, seeds, roit, roi, sroi, hroi, hsroi, h, posf, posroi, possroi, i, nn, nr, pixh, pixw));
        
        %%% judge key press %%%
        key = get(gcf, 'CurrentKey');
        while ~strcmp(key, 'return')
            waitforbuttonpress
            key = get(gcf, 'CurrentKey');
            if strcmp(key, 'q')
                flag = 1; %%% quit the whole loop %%%
                break
            end
        end
        
        if flag == 1
            break
        end
%         set(gcf, 'WindowKeyReleaseFcn', @keyRelease);
        close
    end
    
    %% update output variables %%
    notuse = unique(notuse);
    neuron_use = setdiff(1: nn, notuse);
    
    %% save good neurons %%
    roifn = roi(:, neuron_use);
    sigfn = sig(neuron_use, :);
    seedsfn = seeds(neuron_use);
    roifnr = roir(:, neuron_use);
    sigfnr = sigr(neuron_use, :);
    bgfn = bg;
    bgffn = bgf;
    save([path_name, file_id, '_refined.mat'], 'roifn', 'sigfn', 'seedsfn', 'bgfn', 'bgffn', 'roifnr', 'sigfnr', 'imax', 'pixh', 'pixw', 'raw_score', 'corr_score', 'Params')
end


