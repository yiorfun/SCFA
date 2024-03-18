%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TCGA_Kidney_Gene_R_ALL_2448_1777.mat');
[A_TCGA, B_TCGA] = BEST_UNBIASED_ESTIMATOR(TCGA_Kidney_Gene_R_ALL_2448_1777, p_TCGA_Kidney_Gene);
TCGA_Kidney_Gene_R_ALL_EST_2448_1777 = BLOCK_HADAMARD_PRODUCT(A_TCGA, B_TCGA, p_TCGA_Kidney_Gene);
L = blkdiag(diag(ones(p_TCGA_Kidney_Gene(1))), diag(ones(p_TCGA_Kidney_Gene(2))), diag(ones(p_TCGA_Kidney_Gene(3))), diag(ones(p_TCGA_Kidney_Gene(4))), diag(ones(p_TCGA_Kidney_Gene(5))), diag(ones(p_TCGA_Kidney_Gene(6))));
Sigma_f = B_TCGA;
Sigma_u = blkdiag(A_TCGA(1, 1) * diag(diag(ones(p_TCGA_Kidney_Gene(1)))), A_TCGA(2, 2) * diag(diag(ones(p_TCGA_Kidney_Gene(2)))), A_TCGA(3, 3) * diag(diag(ones(p_TCGA_Kidney_Gene(3)))), A_TCGA(4, 4) * diag(diag(ones(p_TCGA_Kidney_Gene(4)))), A_TCGA(5, 5) * diag(diag(ones(p_TCGA_Kidney_Gene(5)))), A_TCGA(6, 6) * diag(diag(ones(p_TCGA_Kidney_Gene(6)))));

figure;imagesc(TCGA_Kidney_Gene_R_ALL_EST_2448_1777);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sigma
figure;imagesc(p_TCGA_Kidney_Gene);colormap jet;colorbar; 
%%% Community Membership (modifying the colors and width)
figure;imagesc(L);colormap jet;colorbar;caxis([-1, 1]); 
%%% L
figure;imagesc(Sigma_f);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sigma_f
figure;imagesc(L');colormap jet;colorbar;caxis([-1, 1]); 
%%% L'
figure;imagesc(Sigma_u);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sigma_u

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 3 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TCGA_Kidney_Gene_R_ALL_2448_1777.mat');
[A_TCGA, B_TCGA] = BEST_UNBIASED_ESTIMATOR(TCGA_Kidney_Gene_R_ALL_2448_1777, p_TCGA_Kidney_Gene);
TCGA_Kidney_Gene_R_ALL_EST_2448_1777 = BLOCK_HADAMARD_PRODUCT(A_TCGA, B_TCGA, p_TCGA_Kidney_Gene);
L_hat = blkdiag(diag(ones(p_TCGA_Kidney_Gene(1))), diag(ones(p_TCGA_Kidney_Gene(2))), diag(ones(p_TCGA_Kidney_Gene(3))), diag(ones(p_TCGA_Kidney_Gene(4))), diag(ones(p_TCGA_Kidney_Gene(5))), diag(ones(p_TCGA_Kidney_Gene(6))));
Sigma_f_hat = B_TCGA;
Sigma_u_hat = blkdiag(A_TCGA(1, 1) * diag(diag(ones(p_TCGA_Kidney_Gene(1)))), A_TCGA(2, 2) * diag(diag(ones(p_TCGA_Kidney_Gene(2)))), A_TCGA(3, 3) * diag(diag(ones(p_TCGA_Kidney_Gene(3)))), A_TCGA(4, 4) * diag(diag(ones(p_TCGA_Kidney_Gene(4)))), A_TCGA(5, 5) * diag(diag(ones(p_TCGA_Kidney_Gene(5)))), A_TCGA(6, 6) * diag(diag(ones(p_TCGA_Kidney_Gene(6)))));

figure;imagesc(TCGA_Kidney_Gene_R_RAW_2448_1777);colormap jet;colorbar;caxis([-1, 1]); 
%%% Raw Sample Covariance Matrix
figure;imagesc(TCGA_Kidney_Gene_R_ALL_2448_1777);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sample Covariance Matrix
figure;imagesc(TCGA_Kidney_Gene_R_ALL_EST_2448_1777);colormap jet;colorbar;caxis([-1, 1]); 
%%% (Estimated) Population Covariance Matrix
figure;imagesc(L_hat);colormap jet;colorbar;caxis([-1, 1]); 
%%% L_hat
figure;imagesc(Sigma_f_hat);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sigma_f_hat
figure;imagesc(Sigma_u_hat);colormap jet;colorbar;caxis([-1, 1]); 
%%% Sigma_u_hat



