%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE D1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('SM_Community_Detection.mat');
figure;imagesc(R0);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_all);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_per);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_NICE);colormap jet;colorbar; caxis([-1, 1]); 
figure;imagesc(S_kmedoids);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_HCD);colormap jet;colorbar; caxis([-1, 1]); 
