#######################################################################
### Supplementary Material Section C: Additional Simulation Studies ###
#######################################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"HCD", 
	### HCD::HCD(), hierarchical community detection
	"mvnfast",		
	### mvnfast::rmvn(), fast generate multivariate norm
    "matlib", 
	### matlib::symMat(), create a symmetric matrix from a vector
	"matrixcalc",
	### matrixcalc::vech(), create a vector from a symmetric matrix
	"R.matlab"
	### convert between R and MATLAB files
)

CHECK_PACKAGES <- lapply(X = REQUIRED_PACKAGES,
					   FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    require(x, character.only = TRUE)
  }
)

#############################
### Loading the functions ###
#############################

BLOCK_HADAMARD_PRODUCT <- function(A, B, p_vec){
	K <- length(p_vec)
	COL_temp <- c()
	for(k in 1 : K){
		ROW_temp <- c()
		for(kp in 1 : K){
			if(k == kp){
				ROW_temp <- cbind(ROW_temp, A[k, k] * diag(p_vec[k]) + B[k, k] * matrix(1, p_vec[k], p_vec[k]))
			} else{
				ROW_temp <- cbind(ROW_temp, B[k, kp] * matrix(1, p_vec[k], p_vec[kp]))
			}
		}
		COL_temp <- rbind(COL_temp, ROW_temp)
	}
	M <- COL_temp
	return(M)
}

BEST_UNBIASED_ESTIMATOR <- function(S, p_vec){
	## this version can deal with NA values
	K <- length(p_vec)
	A_temp <- matrix(0, K, K)
	B_temp <- matrix(0, K, K)
	for(k in 1 : K){
		for(kp in 1 : K){
			SUB <- S[(sum(p_vec[0 : (k - 1)]) + 1) : sum(p_vec[0 : k]), 
						  (sum(p_vec[0 : (kp - 1)]) + 1) : sum(p_vec[0 : kp])]
			if(kp == k){
				SUB_ON <- mean(diag(SUB), na.rm = TRUE)
				SUB_OF <- (sum(SUB[upper.tri(SUB)], na.rm = TRUE) + sum(SUB[lower.tri(SUB)], na.rm = TRUE)) / (p_vec[k] ^ 2 - p_vec[k] - sum(is.na(SUB[upper.tri(SUB)])) - sum(is.na(SUB[lower.tri(SUB)])))
				A_temp[k, kp] <- SUB_ON - SUB_OF
				B_temp[k, kp] <- SUB_OF
			} else{
				B_temp[k, kp] <- mean(as.vector(SUB), na.rm = TRUE)
			}	
		}
	}
	A <- A_temp
	B <- (B_temp + t(B_temp)) / 2
	return(list(A = A, B = B))
}n <- 50

### sample size
p_vec <- c(30, 30, 40)
### the known partition-size vector
K <- length(p_vec) 				
p <- sum(p_vec)         		

mu0 <- rep(0, p)
theta_A0 <- c(0.1, 0.2, 0.5)
theta_B0 <- c(1.020878, - 0.730920, -0.153572, 3.130999, - 1.633474, 2.893128)
A0 <- diag(theta_A0)			
B0 <- matlib::symMat(theta_B0)			
Sigma0 <- BLOCK_HADAMARD_PRODUCT(A0, B0, p_vec)
R0 <- cov2cor(Sigma0)
set.seed(20)
obs_dat <- mvnfast::rmvn(n = n, mu = mu0, sigma = R0)
### data set without permutation
S_all <- cor(obs_dat)
### sample covariance matrix without permutation

var_names <- paste("Var", seq(p), sep = "_")
colnames(obs_dat) <- var_names
set.seed(2024)
names_permutation <- sample.int(n = p, size = p, replace = FALSE)
obs_per_dat <- obs_dat[, names_permutation]
### data set with permutation 
S_per <- cor(obs_per_dat)
### sample covariance matrix with permutation
### this is the version we start with

### community detection using HCD algorithm ###
### RES_HCD <- HCD::HCD(A = S_per, method = "SC", stopping = "Fix", reg = FALSE, n.min = 15, D = 2, notree = TRUE)
RES_HCD <- HCD::HCD(A = S_per, method = "SC", stopping = "NB", reg = FALSE, n.min = 15, D = NULL, notree = TRUE)
HCD_labels <- RES_HCD$labels
table(HCD_labels)
obs_per_HCD <- cbind(obs_per_dat[, HCD_labels == 1], 
					 obs_per_dat[, HCD_labels == 2], 
					 obs_per_dat[, HCD_labels == 3])
S_HCD <- cor(obs_per_HCD)
	
W <- S_per
W[upper.tri(S_per, diag = TRUE)] <- NA
W_vec <- vech(W)
W_vec <- W_vec[!is.na(W_vec)]

OUTPUT_ADDRESS <- "C:/Users/yxy1423/Documents/R Working Directory File/SCFA_Output"
DATA_OUTPUT_FILENAME <- paste(OUTPUT_ADDRESS, "SM_Community_Detection.mat", sep = "/")
	
R.matlab::writeMat(DATA_OUTPUT_FILENAME, 
			obs_dat = obs_dat, 
	    obs_per_dat = obs_per_dat, 
    		     R0 = R0, 
		      S_per = S_per, 
		      S_all = S_all,
			  S_HCD = S_HCD,			  
  names_permutation = names_permutation, 			
		      W_vec = W_vec)

### community detection using NICE and k-medoids clustering algorithms in Matlab ###

load('SM_Community_Detection.mat')

rng(0);
[CindxVICC, CIDVICC, ClistVICC] = NICE(transpose(W_vec), -0.3, 0, 20);
warning off;
S_NICE = S_per(ClistVICC, ClistVICC);

[idx, C, sumd, D, midx, info] = kmedoids(transpose(obs_per_dat), 3);
obs_per_kmedoids = [obs_per_dat(1:50, idx == 1), obs_per_dat(1:50, idx == 3), obs_per_dat(1:50, idx == 2)];
S_kmedoids = corrcoef(obs_per_kmedoids);

figure;imagesc(R0);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_all);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_per);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_NICE);colormap jet;colorbar; caxis([-1, 1]); 
figure;imagesc(S_kmedoids);colormap jet;colorbar; caxis([-1, 1]);
figure;imagesc(S_HCD);colormap jet;colorbar; caxis([-1, 1]); 

save("SM_Community_Detection.mat", "S_NICE", "S_kmedoids", "-append");




