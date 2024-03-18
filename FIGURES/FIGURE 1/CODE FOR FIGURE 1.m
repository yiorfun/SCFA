%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIGURE 1 A Brain Imaging Study %%%
load('EPSI_FOUR_MATS.mat');
[A_EPSI, B_EPSI] = BEST_UNBIASED_ESTIMATOR(EPSI_RAW, [89;89;89;89;89]);
EPSI_RAW_EST = BLOCK_HADAMARD_PRODUCT(A_EPSI, B_EPSI, [89;89;89;89;89]);
figure;imagesc(EPSI_RAW);colormap jet;colorbar; %%% in Figure 1 A Sample
figure;imagesc(EPSI_RAW_EST);colormap jet;colorbar; %%% in Figure 1 A Population

%%% FIGURE 1 B Gene Expression Study %%%
load('TCGA_Kidney_Gene_R_ALL_2448_1777.mat');
[A_TCGA, B_TCGA] = BEST_UNBIASED_ESTIMATOR(TCGA_Kidney_Gene_R_ALL_2448_1777, p_TCGA_Kidney_Gene);
TCGA_Kidney_Gene_R_ALL_EST_2448_1777 = BLOCK_HADAMARD_PRODUCT(A_TCGA, B_TCGA, p_TCGA_Kidney_Gene);
figure;imagesc(TCGA_Kidney_Gene_R_ALL_2448_1777);colormap jet;colorbar; %%% in Figure 1 B Sample
figure;imagesc(TCGA_Kidney_Gene_R_ALL_EST_2448_1777);colormap jet;colorbar; %%% in Figure 1 B Population

%%% FIGURE 1 C Multi-Omics Study %%%
load('Seed_FOUR_MATS.mat');
figure;imagesc(Seed_ALL);colormap jet;colorbar; %%% in Figure 1 C Sample
figure;imagesc(Seed_ALL_EST);colormap jet;colorbar; %%% in Figure 1 C Population

%%% FIGURE 1 D Plasma Metabolomics Study %%%
load('NMR_E_FOUR_MATS.mat');
figure;imagesc(NMR_E_170);colormap jet;colorbar; %%% in Figure 1 D Sample
figure;imagesc(NMR_E_170_5_EST);colormap jet;colorbar; %%% in Figure 1 D Population

%%% FIGURE 1 E Environmental Exposome Study %%%
load('Exposome_FOUR_MATS.mat');
figure;imagesc(Exposome_89);colormap jet;colorbar; %%% in Figure 1 E Sample
figure;imagesc(Exposome_89_EST);colormap jet;colorbar; %%% in Figure 1 E Population

%%% FIGURE 1 F Environmental Plasma Metabolomics Study %%%
load('Metabolites_FOUR_MATS.mat');
figure;imagesc(Metabolites_141);colormap jet;colorbar; %%% in Figure 1 F Sample
figure;imagesc(Metabolites_141_EST);colormap jet;colorbar; %%% in Figure 1 F Population



