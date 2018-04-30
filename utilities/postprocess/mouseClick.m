function mouseClick(~, ~, tp, seeds, all_neuron_roi, roi, sroi, hroi, hsroi, h, posf, posroi, possroi, i, nn, nr, pixh, pixw)       
% callback function for mouse clicking in the figure
%   Jinghao Lu, 05/01/2015

    global notuse
    C = get(gcf, 'CurrentPoint');
    x = C(1);
    y = C(2);
    flag = 0;
    
    %% first check if on the left 10 traces %%
    for j = 1: size(posf, 1)
        xmin = posf(j, 1);
        xmax = posf(j, 1) + posf(j, 3);
        ymin = posf(j, 2);
        ymax = posf(j, 2) + posf(j, 4);
        if inpolygon(x, y, [xmin, xmin, xmax, xmax], [ymin, ymax, ymax, ymin])
            notuse = [notuse, (i - 1) * nr + j];
            disp(num2str((i - 1) * nr + j))
            flag = 1;
            break
        end
    end
    
    %% not on the left, then check in the roi image %%
    if flag == 0
        xmin = posroi(1, 1);
        xmax = posroi(1, 1) + posroi(1, 3);
        ymin = posroi(1, 2);
        ymax = posroi(1, 2) + posroi(1, 4);
        img = all_neuron_roi;
        if inpolygon(x, y, [xmin, xmin, xmax, xmax], [ymin, ymax, ymax, ymin])
            roi_pos = get(hroi, 'CurrentPoint');
            roi_pos = [round(roi_pos(1, 1)), round(roi_pos(1, 2))];
            colorid = img(roi_pos(2), roi_pos(1));
            %%% if mouse is on some rois, then reimagesc the panel %%%
            if colorid ~= 0
                notuse = [notuse, colorid];
                disp(num2str(colorid))
                
                a = roi(:, colorid);
%                 axes(sroi);
                tmp = reshape(a, pixh, pixw);
%                 [xx, yy] = find(tmp == max(tmp(:)));
%                 xx = xx(1);
%                 yy = yy(1);
%                 tpt = zeros(51, 51);
%                 ttmp = tmp(max(1, xx - 25): min(pixh, xx + 25), max(1, yy - 25): min(pixw, yy + 25));
%                 tpt(max(1, 27 - xx): size(ttmp, 1), max(1, 27 - yy): size(ttmp, 2)) = ttmp;
%                 imagesc(tpt)
%                 axis off
%                 xlabel(['(X,Y) = (', num2str(C(1)), ', ',num2str(C(2)), ')']);
                
                tmpp = imgaussfilt(tmp, 3);
                lvl = max(tmpp(:)) * 0.5;
                ctrs = contours(tmpp, [lvl, lvl]);
                axes(hroi);
                cla
                imagesc(tp)
                hold on
                plot(ctrs(1, 2: end), ctrs(2, 2: end), 'r')
                hold off
            end
        end
    end
    
    %%% update in the roi image of previous selected neurons %%%
    axes(hsroi);
    imagesc(tp)
    hold on
    for j = 1: length(notuse)
        [y, x] = ind2sub([pixh, pixw], seeds(notuse(j)));
        plot(x, y, '.r')
    end
    hold off
end
