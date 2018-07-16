%%% demo of the full MIN1PIPE %%%

%% session-specific parameter initialization %% 
Fsi = 20;
Fsi_new = 20; %%% no temporal downsampling %%%
spatialr = 1; %%% no spatial downsampling %%%
se = 5; %%% structure element for background removal %%%
ismc = true; %%% run movement correction %%%
flag = 1; %%% use auto seeds selection; 2 if manual %%%
% isvis = true; %%% do visualize %%%
ifpost = false; %%% set true if want to see post-process %%%

%% main program %%
[fname, frname] = min1pipe(Fsi, Fsi_new, spatialr, se, ismc, flag);

%% plot some images %%
load(fname)
figure(1)
clf
%%% raw max %%%
subplot(2, 3, 1, 'align')
imagesc(imaxn)
axis square
title('Raw')

%%% neural enhanced before movement correction %%%
subplot(2, 3, 2, 'align')
imagesc(imaxy)
axis square
title('Before MC')

%%% neural enhanced after movement correction %%%
subplot(2, 3, 3, 'align')
imagesc(imax)
axis square
title('After MC')

%%% contour %%%
subplot(2, 3, 4, 'align')
plot_contour(roifn, sigfn, seedsfn, imax, pixh, pixw)
axis square

%%% movement measurement %%%
subplot(2, 3, 5, 'align')
plot(raw_score); hold on; plot(corr_score); hold off;
axis square
title('MC Scores')

%%% all identified traces %%%
subplot(2, 3, 6, 'align')
sigt = sigfn;
for i = 1: size(sigt, 1)
    sigt(i, :) = normalize(sigt(i, :));
end
plot((sigt + (1: size(sigt, 1))')')
axis tight
axis square
title('Traces')

%% make a movie %%
load(frname)
id = find(fname == filesep, 1, 'last');
fmovie = [fname(1: id), 'demo_vid.avi'];
v = VideoWriter(fmovie);
v.FrameRate = Fsi_new;
v.Quality = 100;
open(v)

figure(2)
clf
set(gcf, 'Units', 'normalized', 'position', [0.5, 0.1, 0.4, 0.2])
datar = reshape(roifn * sigfn, pixh, pixw, []);
for i = 1: size(datar, 3)
    subplot(1, 4, 1, 'align')
    imagesc(data.norm_frame(:, :, i), [0, 1])
    axis off
    axis square
    title('Raw')
    
    subplot(1, 4, 2, 'align')
    imagesc(data.Ydebg(:, :, i), [0, 1])
    axis off
    axis square
    title('Before MC')
    
    subplot(1, 4, 3, 'align')
    imagesc(data.reg(:, :, i), [0, 1])
    axis off
    axis square
    title('After MC')
    
    subplot(1, 4, 4, 'align')
    imagesc(datar(:, :, i), [0, 1])
    axis off
    axis square
    title('Processed')

    suptitle(['Frame #', num2str(i)])
    
    movtmp = getframe(gcf);
    writeVideo(v, movtmp);
end
close(v)

%% post-process %%
if ifpost
    real_neuron_select
end