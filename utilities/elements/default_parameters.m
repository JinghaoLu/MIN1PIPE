function Params = default_parameters 
% default parameters of MIN1PIPE; not important, just to make it robust
%   Jinghao Lu, 03/21/2018

    Params.Fsi = 20;
    Params.Fsi_new = 10;
    Params.spatialr = 0.5;
    Params.ttype = 'single'; %%% default data type %%%
    Params.neuron_size = 5; %%% half neuron size; 9 for Inscopix and 5 for UCLA, with 0.5 spatialr separately %%%
%     Params.data_cat_ispara = 0;                                         
    Params.anidenoise_iter = 4;                                         
    Params.anidenoise_dt = 1/7; %%% for 5-point formular %%%                                       
    Params.anidenoise_kappa = 0.5;                                      
    Params.anidenoise_opt = 1;                                          
    Params.anidenoise_ispara = 1;                                          
    Params.bg_remove_ispara = 1;  
    Params.mc_scl = 0.004;
    Params.mc_sigma_x = 5;
    Params.mc_sigma_f = 10;
    Params.mc_sigma_d = 1;
    Params.pix_select_sigthres = 0.8;                                   
    Params.pix_select_corrthres = 0.6;                                  
    Params.refine_roi_ispara = 1;                                       
    Params.merge_roi_corrthres = 0.9;                                    
end
