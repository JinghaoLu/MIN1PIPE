function mouseMove(~, ~, tp, all_wave_det, all_neuron_roi, roi, sroi, hroi, hroi_wave, h, posf, posroi, i, nn, nr, pixh, pixw)
% callback function for mouse position
%   Jinghao Lu, 05/06/2015

    global id_mouse
    C = get(gcf, 'CurrentPoint');
    x = C(1);
    y = C(2);
    flag = 0;
    
    %% judge if the mouse is on 10 traces %%
    for j = 1: size(posf, 1)
        xmin = posf(j, 1);
        xmax = posf(j, 1) + posf(j, 3);
        ymin = posf(j, 2);
        ymax = posf(j, 2) + posf(j, 4);
        if inpolygon(x, y, [xmin, xmin, xmax, xmax], [ymin, ymax, ymax, ymin])
            id_mouse = (i - 1) * nr + j;
            flag = 1;

            for jj = 1: min(nr, nn - (i - 1) * nr)
                axes(h(jj));
                if jj == j
                    plot(all_wave_det((i - 1) * nr + jj, :), 'color', 'k')
                else
                    plot(all_wave_det((i - 1) * nr + jj, :), 'color', [0.8,0.8,0.8])
                end
                axis tight
                ylim([0, 1])
                axis off
            end
        
            axes(hroi_wave)
            plot(all_wave_det(id_mouse, :), 'k')
            axis tight
%             ylim([-0.1, 0.2])
            axis off
            
            a = roi(:, id_mouse);
            axes(sroi);
            tmp = reshape(a, pixh, pixw);
            [xx, yy] = find(tmp == max(tmp(:)));
            xx = xx(1);
            yy = yy(1);
            tpt = zeros(51, 51);
            ttmp = tmp(max(1, xx - 25): min(pixh, xx + 25), max(1, yy - 25): min(pixw, yy + 25));
            tpt(max(1, 27 - xx): max(1, 27 - xx) + size(ttmp, 1) - 1, max(1, 27 - yy): max(1, 27 - yy) + size(ttmp, 2) - 1) = ttmp;
            imagesc(tpt)
            axis off
            xlabel(['(X,Y) = (', num2str(C(1)), ', ',num2str(C(2)), ')']);
            
            tmpp = imgaussfilt(tmp, 3);
            lvl = max(tmpp(:)) * 0.5;
            ctrs = contours(tmpp, [lvl, lvl]);
            axes(hroi);
            cla
            imagesc(tp)
            hold on
            plot(ctrs(1, 2: end), ctrs(2, 2: end), 'r')
            hold off
            break
        end
    end
    
    %% judge if the mouse is in roi selection pool %%
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
                id_mouse = colorid; %%% update id_mouse %%%
          
                %%% plot current roi wave %%%
                axes(hroi_wave)
                plot(all_wave_det(colorid, :), 'k')
                axis tight
%                 ylim([0, 0.2])
                axis off
                
                a = roi(:, colorid);
                axes(sroi);
                tmp = reshape(a, pixh, pixw);
                [xx, yy] = find(tmp == max(tmp(:)));
                xx = xx(1);
                yy = yy(1);
                tpt = zeros(51, 51);
                ttmp = tmp(max(1, xx - 25): min(pixh, xx + 25), max(1, yy - 25): min(pixw, yy + 25));
                tpt(max(1, 27 - xx): max(1, 27 - xx) + size(ttmp, 1) - 1, max(1, 27 - yy): max(1, 27 - yy) + size(ttmp, 2) - 1) = ttmp;
                imagesc(tpt)
                axis off
                xlabel(['(X,Y) = (', num2str(C(1)), ', ',num2str(C(2)), ')']);
            
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
end
