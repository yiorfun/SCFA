#####################################################################
### Simulation Studies: Assessing the estimated factor scores and ###
### model parameters 											  ###
#####################################################################

##################################################################
# Goal 1: finite-sample performance of factor score estimators   #
# average of the Euclidean loss (mean)			      			 #
# standard deviation of the Euclidean loss (standard deviation)	 #
# computation time in seconds (time)							 #
#																 #
# Goal 2: finite-sample performance of model parameter estimator #	 
# average of estimation bias (bias)								 #
# Monte Carlo standard deviation (MCSD)							 #
# average of standard errors (ASE)								 #
# coverage probability (CP) of 95% Wald-type confidence interval #
#																 #
# Goal 3: misspecification analysis								 #	 
##################################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"MASS",		
	### MASS::ginv(), solve the generalized inverse
    "matlib", 
	### matlib::symMat(), create a symmetric matrix from a vector
	"matrixcalc",
	### matrixcalc::vech(), create a vector from a symmetric matrix
	"sem",
	### estimate a CFA model
	"OpenMx",
	### estimate a CFA model
	"lavaan"
	### estimate a CFA model
)

CHECK_PACKAGES <- lapply(X = REQUIRED_PACKAGES,
					   FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    require(x, character.only = TRUE)
  }
)

##################################
### Loading required functions ###
##################################

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
}

CALCULATE_VAR_A_B <- function(A, B, p_vec, n){
	
	K <- length(p_vec)
	VAR_A_MAT <- diag(K)
	VAR_B_MAT <- matrix(0, K, K)
	
	for(k in 1 : K){
		for(kp in 1 : K){
			if(kp == k){
				VAR_A_MAT[k, kp] <- 2 * A[k, k] ^ 2 / ((n - 1) * (p_vec[k] - 1))
				VAR_B_MAT[k, kp] <- 2 * ((A[k, k] + p_vec[k] * B[k, k]) ^ 2 - B[k, k] * (2 * A[k, k] + p_vec[k] * B[k, k])) / ((n - 1) * (p_vec[k] - 1) * p_vec[k])
			} else {
				VAR_A_MAT[k, kp] <- 0
				VAR_B_MAT[k, kp] <- (p_vec[k] * p_vec[kp] * (B[k, kp] ^ 2 + B[kp, k] ^ 2) + 2 * (A[k, k] + p_vec[k] * B[k, k]) * (A[kp, kp] + p_vec[kp] * B[kp, kp])) / (2 * (n - 1) * p_vec[k] * p_vec[kp])
			}
		}
	}
	return(list(VAR_A_MAT = VAR_A_MAT, VAR_B_MAT = VAR_B_MAT))
}

cfa_MODIFIED <- function(file="", text, covs=paste(factors, collapse=","), reference.indicators=TRUE, raw=FALSE, 
    subscript=c("name", "number"), ...){
    Lines <- if (!missing(text)) scan(text=text, what="", sep=";", strip.white=TRUE, comment.char="#", quiet = TRUE)
    else scan(file=file, what="", sep=";", strip.white=TRUE, comment.char="#", quiet = TRUE)
    lines <- character(0)
    current.line <- ""
    for (line in Lines){
        if (current.line != "") line <- paste(current.line, line)
        if (length(grep(",$", line)) > 0){
            current.line <- line
            next
        }
        current.line <- ""
        lines <- c(lines, line)
    }
    subscript <- match.arg(subscript)
    variance.lines <- grepl("^[Vv]ar.*:", lines)
    variances <- lines[variance.lines]
    lines <- lines[!variance.lines]
    nfactor <- length(lines)
    factors <- rep("", nfactor)
    all.obs.vars <- ram <- character(0)
    equality.sets <- list()
    for (i in 1:nfactor){
        par.number <- 0
        Line <- line <- lines[[i]]
        line <- gsub(" ", "", line)
        line <- strsplit(line, ":")[[1]]
        if (length(line) == 1){
            factors[i] <- paste("Factor.", i, sep="")
            variables <- strsplit(line, ",")[[1]]
            all.obs.vars <- c(all.obs.vars, variables)
        }
        else if (length(line) == 2){
            factors[i] <- line[1]
            variables <- strsplit(line[2], ",")[[1]]
            all.obs.vars <- c(all.obs.vars, unlist(strsplit(variables, "=")))
        }
        else stop("Parse error in ", Line)
        if (reference.indicators){
            if (!grepl("=", variables[1])){
                ram <- c(ram, paste(factors[i], " -> ", variables[1], ", NA, 1", sep=""))}
            else{
                vars <- strsplit(variables[1], "=")[[1]]
                equality.sets[[length(equality.sets) + 1]] <- vars
                for (var in vars){
                    ram <- c(ram, paste(factors[i], " -> ", var, ", NA, 1", sep=""))
                }
            }
            variables <- variables[-1]
        }
        for (variable in variables){
            if (length(grep("\\(", variable)) > 0){
                if (length(grep("\\)", variable)) == 0) stop ("Parse error in ", Line)
                variable <- sub("\\)", "", variable)   
                var.start <- strsplit(variable, "\\(")[[1]]
                if (length(var.start) != 2) stop("Parse error in ", Line)
                variable <- var.start[1]
                start <- var.start[2]
                if (not.number(start)) stop ("Bad start value ", start, " in ", Line)
            }
            else start <- "NA"
            if (!grepl("=", variable)){
                par.number <- par.number + 1
                par.name <- if (subscript == "name") variable else as.character(par.number)
                ram <- c(ram, paste(factors[i], " -> ", variable, ", lam[", par.name, ":", factors[i], "], ", start, sep=""))
            }
            else {
                vars <- strsplit(variable, "=")[[1]]
                equality.sets[[length(equality.sets) + 1]] <- vars
                par.number <- par.number + 1
                lam <- if (subscript == "name") paste(vars, collapse=".") else as.character(par.number)
                for (var in vars){
                    ram <- c(ram, paste(factors[i], " -> ", var, ", lam[", lam, ":", factors[i], "], ", start, sep=""))
                }
            }
        }
    }
    ram <- if (reference.indicators) {
        c(ram, sapply(factors, function(factor) paste(factor, " <-> ", factor, ", ", paste("V[", factor, "]", sep="") , ", NA", sep="")))
    }
    else{
        c(ram, sapply(factors, function(factor) paste(factor, " <-> ", factor, ", NA, 1", sep="")))
    }
    if (raw){
        all.obs.vars <- unique(all.obs.vars)
        if (length(equality.sets) == 0){
            int <- if (subscript == "name") all.obs.vars else as.character(seq(1, length(all.obs.vars)))
            names(int) <- all.obs.vars
            ram <- c(ram, sapply(all.obs.vars, function(var) paste("Intercept -> ", var, ", intercept(", int[var], "), NA", sep="")))
        }
        else{
            par.number <- 0
            for (set in equality.sets){
                par.number <- par.number + 1
                int <- if (subscript == "name") paste(set, collapse=".") else as.character(par.number)
                ram <- c(ram, sapply(set, function(var) 
                    paste("Intercept -> ", var, ", intercept(", int, "), NA", sep="")))
                all.obs.vars <- setdiff(all.obs.vars, set)
            }
            if (length(all.obs.vars) > 0) {
                int <- if (subscript == "name") all.obs.vars else as.character(seq(par.number + 1, par.number + length(all.obs.vars)))
                names(int) <- all.obs.vars
                ram <- c(ram, sapply(all.obs.vars, function(var) paste("Intercept -> ", var, ", intercept(", int[var], "), NA", sep="")))
            }
        }
        message('NOTE: specify fixed.x="Intercept" in call to sem')
    }
    if (length(variances) > 0){
        var.number <- 0
        variances <- sub("^[Vv]ar.*:", "", variances)
        variances <- gsub(" ", "", variances)
        variances <- strsplit(variances, ",")
        for (vars in variances){
            var <- strsplit(vars, "=")
            sub <- if (subscript == "name") sapply(var, function(x) paste(x, collapse=".")) else as.character(seq(var.number + 1, var.number + length(var) + 1))
            var.number <- var.number + length(var)
            for (i in 1:length(var)){
                ram <- c(ram, sapply(var[i], function(x) paste(x, " <-> ", x, ", V[", sub[i], "]", sep="")))
            }
        }
    }
    sem::specifyModel(text=ram, covs=covs, ..., quiet=TRUE)
}

####################################
### Setup for Simulation Studies ###
####################################

SET_NO <- 1

SET_UP <- matrix(c( 40,  2, 0, ### 1
				    40,  3, 0, ### 2
				    40,  4, 0, ### 3
				    80,  4, 0, ### 4
				    80,  8, 0, ### 5
				    80, 12, 0, ### 6 
				   120,  4, 0, ### 7
				   120, 12, 0, ### 8
				   120, 20, 0, ### 9
				   120, 20, 1
			), nrow = 10, ncol = 3, byrow = TRUE)
### sample size, number of observations, misspecified or not
### 1,2,4,7 -> sem() and cfa()
### all -> prop and mxRun()

set.seed(2023)
reps <- 1
repsMax <- 100
ALPHA <- 0.05
n <- SET_UP[SET_NO, 1]
multiple <- SET_UP[SET_NO, 2]
mis_or_not <- SET_UP[SET_NO, 3]
p_vec <- c(3, 3, 4) * multiple	
K <- length(p_vec) 				
### K = 3
p <- sum(p_vec)         		
### p = 10 * multiple
theta_A0 <- c(0.1, 0.2, 0.5)
theta_B0 <- c(2.020878, 0.730920, 1.153572, 3.130999, 1.633474, 3.693128)
theta0 <- c(theta_A0, theta_B0)
A0 <- diag(theta_A0)			
### A0 should be positive definite
B0 <- symMat(theta_B0)			
### B0 should be positive definite

MU_x <- rep(0, p)
MU_f <- rep(0, K)
MU_u <- rep(0, p)
SIGMA_f <- B0
SIGMA_u <- BLOCK_HADAMARD_PRODUCT(A = A0, B = matrix(0, K, K), p_vec = p_vec)
SIGMA_0 <- BLOCK_HADAMARD_PRODUCT(A = A0, B = B0, p_vec = p_vec)
E_sigma_1 <- rWishart(n = 1, df = p, Sigma = 0.01 * diag(p))[,,1]
E_sigma_3 <- rWishart(n = 1, df = p, Sigma = 0.03 * diag(p))[,,1]
E_sigma_5 <- rWishart(n = 1, df = p, Sigma = 0.05 * diag(p))[,,1]
L <- as.matrix(cbind(rep(c(1, 0, 0), p_vec), rep(c(0, 1, 0), p_vec), rep(c(0, 0, 1), p_vec)))
### p by K

theta_EST <- matrix(0, repsMax, length(theta0))
theta_ASE <- matrix(0, repsMax, length(theta0))
theta_WCP <- matrix(0, repsMax, length(theta0))
theta_EST_1 <- matrix(0, repsMax, length(theta0))
theta_ASE_1 <- matrix(0, repsMax, length(theta0))
theta_WCP_1 <- matrix(0, repsMax, length(theta0))
theta_EST_3 <- matrix(0, repsMax, length(theta0))
theta_ASE_3 <- matrix(0, repsMax, length(theta0))
theta_WCP_3 <- matrix(0, repsMax, length(theta0))
theta_EST_5 <- matrix(0, repsMax, length(theta0))
theta_ASE_5 <- matrix(0, repsMax, length(theta0))
theta_WCP_5 <- matrix(0, repsMax, length(theta0))

loading_VEC_EST_LAV <- loading_VEC_EST_SEM <- loading_VEC_EST_MX <- matrix(0, repsMax, p) 
### estimation of non-zero elements of loadings
loading_VEC_ASE_LAV <- loading_VEC_ASE_SEM <- loading_VEC_ASE_MX <- matrix(0, repsMax, p) 
### standard error of non-zero elements of loadings
error_VEC_EST_LAV <- error_VEC_EST_SEM <- error_VEC_EST_MX <- matrix(0, repsMax, p)   
### estimation of A0 * I(p)
error_VEC_ASE_LAV <- error_VEC_ASE_SEM <- error_VEC_ASE_MX <- matrix(0, repsMax, p)   
### standard error of A0 * I(p)
loading_VAR_EST_LAV <- loading_VAR_EST_SEM <- loading_VAR_EST_MX <- matrix(0, repsMax, K * (K + 1) / 2)
### estimation of B0
loading_VAR_ASE_LAV <- loading_VAR_ASE_SEM <- loading_VAR_ASE_MX <- matrix(0, repsMax, K * (K + 1) / 2)
### standard error of B0

TIME_UB <- TIME_LAV <- TIME_SEM <- TIME_MX <- 0
LOSS_UB <- LOSS_LAV <- LOSS_SEM <- LOSS_MX <- rep(0, repsMax)
LOSS_MIS_1 <- LOSS_MIS_3 <- LOSS_MIS_5 <- rep(0, repsMax)
RLOSS_MIS_1 <- RLOSS_MIS_3 <- RLOSS_MIS_5 <- rep(0, repsMax)

###########################
### Computing procedure ###
###########################

while(reps <= repsMax){
	
	tryCatch({
	
	if(mis_or_not == 0){
	
		#################################
		### data generation procedure ###
		#################################
	
		F <- t(mvrnorm(n = n, mu = MU_f, Sigma = SIGMA_f)) 
		### K by n
		U <- t(mvrnorm(n = n, mu = MU_u, Sigma = SIGMA_u)) 
		### p by n
		X <- matrix(rep(MU_x, n), p, n, byrow = FALSE) + L %*% F + U 
		### p by n
		DATA <- t(X)
		### n by p
		
		X_NAMES <- paste(rep('x', p_vec[1]), seq(p_vec[1]), sep = '')
		Y_NAMES <- paste(rep('y', p_vec[2]), seq(p_vec[2]), sep = '')
		Z_NAMES <- paste(rep('z', p_vec[3]), seq(p_vec[3]), sep = '')
		COL_NAMES <- c(X_NAMES, Y_NAMES, Z_NAMES)
		colnames(DATA) <- COL_NAMES
		ERR_NAMES <- paste(rep('e', p), seq(p), sep = '')
		LOD_NAMES <- paste(rep('l', p), seq(p), sep = '')
		MU_NAMES <- paste(rep('mean', p), COL_NAMES, sep = '')
		
		#########################################
		### the proposed estimation procedure ###
		#########################################
	
		S <- cov(DATA)
		TIME_TEMP <- Sys.time()
		RES_theta_hat <- BEST_UNBIASED_ESTIMATOR(S, p_vec)
		TIME_UB <- TIME_UB + (Sys.time() - TIME_TEMP)
		A_hat <- RES_theta_hat$A
		B_hat <- RES_theta_hat$B
		theta_hat_temp <- c(diag(A_hat), vech(B_hat))
		theta_EST[reps, ] <- theta_hat_temp
	
		RES_variance <- CALCULATE_VAR_A_B(A_hat, B_hat, p_vec, n)
		VAR_A_MAT <- RES_variance$VAR_A_MAT
		VAR_B_MAT <- RES_variance$VAR_B_MAT
		theta_VAR_temp <- c(diag(VAR_A_MAT), vech(VAR_B_MAT))
		theta_ASE[reps, ] <- sqrt(theta_VAR_temp)
	
		CI_LOWER_temp <- theta_hat_temp - qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp)
		CI_UPPER_temp <- theta_hat_temp + qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp)
		CI_temp <- rbind(CI_LOWER_temp, CI_UPPER_temp)
		theta_WCP_temp <- 1 * (CI_temp[1, ] < theta0 & CI_temp[2, ] > theta0)
		theta_WCP[reps, ] <- theta_WCP_temp
		
		F_HAT_UB <- solve(t(L) %*% L) %*% t(L) %*% X ### K by n
		R_UB <- as.matrix(F_HAT_UB - F)
		LOSS_UB[reps] <- sum(apply(R_UB, 2, norm, "2"))
		
		#########################################
		### sem() functions in package "sem"  ###
	    #########################################
		
		if(SET_NO %in% c(1, 2, 3, 4, 5, 7)){
		
		DATA_SEM <- data.frame(DATA)
		X_NAMES_SEM <- paste(X_NAMES, collapse = ', ')	
		Y_NAMES_SEM <- paste(Y_NAMES, collapse = ', ')
		Z_NAMES_SEM <- paste(Z_NAMES, collapse = ', ') 
		FX_SEM <- paste("FX_SEM: ", X_NAMES_SEM, sep = '')
		FY_SEM <- paste("FY_SEM: ", Y_NAMES_SEM, sep = '')
		FZ_SEM <- paste("FZ_SEM: ", Z_NAMES_SEM, sep = '')
		suppressMessages(
			MODEL_SEM <- cfa_MODIFIED(text = paste(FX_SEM, FY_SEM, FZ_SEM, sep = '\n'), reference.indicators = TRUE))
		
		TIME_TEMP <- Sys.time()
		#CFA_SEM <- sem::sem(MODEL_SEM, data = DATA_SEM)
		CFA_SEM <- sem::sem(MODEL_SEM, S = cor(DATA), N = n)
		TIME_SEM <- TIME_SEM + (Sys.time() - TIME_TEMP)
		EST_SEM <- CFA_SEM$coeff
		#ASE_SEM <- sqrt(diag(CFA_SEM$robust.vcov))
		ASE_SEM <- sqrt(diag(CFA_SEM$vcov))
		
		loading_VEC_EST_SEM[reps, ] <- c(1, EST_SEM[1 : (cumsum(p_vec - 1)[1])], 1, EST_SEM[((cumsum(p_vec - 1)[1]) + 1) : (cumsum(p_vec - 1)[2])], 1, EST_SEM[((cumsum(p_vec - 1)[2]) + 1) : (cumsum(p_vec - 1)[3])])	
		loading_VEC_ASE_SEM[reps, ] <- c(0, ASE_SEM[1 : (cumsum(p_vec - 1)[1])], 0, ASE_SEM[((cumsum(p_vec - 1)[1]) + 1) : (cumsum(p_vec - 1)[2])], 0, ASE_SEM[((cumsum(p_vec - 1)[2]) + 1) : (cumsum(p_vec - 1)[3])])	
		error_VEC_EST_SEM[reps, ] <- EST_SEM[- seq(cumsum(p_vec - 1)[3] + K * (K + 1) / 2)]
		error_VEC_ASE_SEM[reps, ] <- ASE_SEM[- seq(cumsum(p_vec - 1)[3] + K * (K + 1) / 2)]
		loading_VAR_EST_SEM[reps, ] <- EST_SEM[(cumsum(p_vec - 1)[3] + 1) : (cumsum(p_vec - 1)[3] + K * (K + 1) / 2)][c(1, 4, 5, 2, 6, 3)]
		loading_VAR_ASE_SEM[reps, ] <- ASE_SEM[(cumsum(p_vec - 1)[3] + 1) : (cumsum(p_vec - 1)[3] + K * (K + 1) / 2)][c(1, 4, 5, 2, 6, 3)]
		
		F_HAT_SEM <- sem::fscores(CFA_SEM, data = DATA_SEM) ### n by K
		R_SEM <- as.matrix(t(F_HAT_SEM) - F)
		LOSS_SEM[reps] <- sum(apply(R_SEM, 2, norm, "2"))
		
		}
		
		##############################################
		### mxRun() functions in package "OpenMx"  ###
		##############################################
		
		if(SET_NO %in% c(1, 2, 4, 7)){
		
		DATA_MX <- data.frame(DATA)
		DATA_MX_RAW <- mxData(observed = DATA_MX, type = "raw")
		resVars <- mxPath(from = COL_NAMES,
						arrows = 2,
                          free = TRUE, 
						values = rep(1, p),
                        labels = ERR_NAMES)
		# residual variances
		latVars <- mxPath(from = c("FX_MX", "FY_MX", "FZ_MX"), 
						arrows = 2, 
					   connect = "unique.pairs",
                          free = TRUE, 
						values = c(1, 0.5, 0.5, 1, 0.5, 1), 
						labels = c("varFX", "covXY", "covXZ", 
								            "varFY", "covYZ",
													 "varFZ"))
		# latent variances and covariance
		facLoadsX <- mxPath(from = "FX_MX", 
							  to = X_NAMES, 
						  arrows = 1,
                            free = c(FALSE, rep(TRUE, p_vec[1] - 1)), 
						  values = rep(1, p_vec[1]), 
						  labels = LOD_NAMES[1 : cumsum(p_vec)[1]])
		# factor loadings for x variables
		facLoadsY <- mxPath(from = "FY_MX", 
							  to = Y_NAMES, 
						  arrows = 1,
							free = c(FALSE, rep(TRUE, p_vec[2] - 1)),
						  values = rep(1, p_vec[2]), 
						  labels = LOD_NAMES[(cumsum(p_vec)[1] + 1) : cumsum(p_vec)[2]])
		# factor loadings for y variables
		facLoadsZ <- mxPath(from = "FZ_MX", 
							  to = Z_NAMES, 
						  arrows = 1,
							free = c(FALSE, rep(TRUE, p_vec[3] - 1)),
						  values = rep(1, p_vec[3]), 
						  labels = LOD_NAMES[(cumsum(p_vec)[2] + 1) : cumsum(p_vec)[3]])
		# factor loadings for z variables
		
		means <- mxPath(from = "one", 
						  to = c(COL_NAMES, "FX_MX", "FY_MX", "FZ_MX"), 
                      arrows = 1,
                        free = c(rep(TRUE, p), rep(FALSE, 3)), 
				      values = c(rep(1, p), rep(0, 3)),
                      labels = c(MU_NAMES, NA, NA, NA)) 
		# means
		MODEL_MX <- mxModel(model = "MODEL_MX", 
							  type = "RAM",
					  manifestVars = COL_NAMES, 
                        latentVars = c("FX_MX","FY_MX", "FZ_MX"),
                        DATA_MX_RAW, resVars, latVars, facLoadsX, facLoadsY, facLoadsZ, means)
		
		TIME_TEMP <- Sys.time()
		CFA_MX <- mxRun(MODEL_MX, silent = TRUE)
		TIME_MX <- TIME_MX + (Sys.time() - TIME_TEMP)
			
		EST_MX <- CFA_MX$output$estimate
		ASE_MX <- CFA_MX$output$standardErrors[,1]
		
		loading_VEC_EST_MX[reps, ] <- c(1, EST_MX[1 : (cumsum(p_vec - 1)[1])], 1, EST_MX[((cumsum(p_vec - 1)[1]) + 1) : (cumsum(p_vec - 1)[2])], 1, EST_MX[((cumsum(p_vec - 1)[2]) + 1) : (cumsum(p_vec - 1)[3])])	
		loading_VEC_ASE_MX[reps, ] <- c(1, ASE_MX[1 : (cumsum(p_vec - 1)[1])], 1, ASE_MX[((cumsum(p_vec - 1)[1]) + 1) : (cumsum(p_vec - 1)[2])], 1, ASE_MX[((cumsum(p_vec - 1)[2]) + 1) : (cumsum(p_vec - 1)[3])])	
		error_VEC_EST_MX[reps, ] <- EST_MX[(cumsum(p_vec - 1)[3] + 1) : (cumsum(p_vec - 1)[3] + p)]
		error_VEC_ASE_MX[reps, ] <- ASE_MX[(cumsum(p_vec - 1)[3] + 1) : (cumsum(p_vec - 1)[3] + p)]
		
		loading_VAR_EST_MX[reps, ] <- EST_MX[(cumsum(p_vec - 1)[3] + p + 1) : (cumsum(p_vec - 1)[3] + p + K * (K + 1) / 2)][c(1, 2, 4, 3, 5, 6)]
		loading_VAR_ASE_MX[reps, ] <- ASE_MX[(cumsum(p_vec - 1)[3] + p + 1) : (cumsum(p_vec - 1)[3] + p + K * (K + 1) / 2)][c(1, 2, 4, 3, 5, 6)]

		suppressMessages(
			F_HAT_MX <- mxFactorScores(CFA_MX)[,,1])
		R_MX <- as.matrix(t(F_HAT_MX) - F)
		LOSS_MX[reps] <- sum(apply(R_MX, 2, norm, "2"))
		
		}
		
		#########################################
		### cfa() function in package lavaan  ###
		#########################################
		
		if(SET_NO %in% c(1, 2, 4, 7)){
		
		X_NAME_LAV <- paste(COL_NAMES[1 : sum(p_vec[1 : 1])], collapse = ' + ')				
		Y_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 1]) + 1) : sum(p_vec[1 : 2])], collapse = ' + ')
		Z_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 2]) + 1) : sum(p_vec[1 : 3])], collapse = ' + ') 
		FX_LAV <- paste(' FX_LAV =~ ', X_NAME_LAV, sep = '')
		FY_LAV <- paste(' FY_LAV =~ ', Y_NAME_LAV, sep = '')
		FZ_LAV <- paste(' FZ_LAV =~ ', Z_NAME_LAV, sep = '')
		MODEL_LAV <- paste(FX_LAV, FY_LAV, FZ_LAV, sep = " \n ")
		
		TIME_TEMP <- Sys.time()
		CFA_LAV <- lavaan::cfa(MODEL_LAV, data = DATA)
		TIME_LAV <- TIME_LAV + (Sys.time() - TIME_TEMP)
		CFA_PARA_LAV <- lavaan::parameterEstimates(CFA_LAV)
	
		loading_VEC_EST_LAV[reps, ] <- CFA_PARA_LAV$est[1 : p]
		loading_VEC_ASE_LAV[reps, ] <- CFA_PARA_LAV$se[1 : p]
		error_VEC_EST_LAV[reps, ] <- CFA_PARA_LAV$est[(p + 1) : (2 * p)]
		error_VEC_ASE_LAV[reps, ] <- CFA_PARA_LAV$se[(p + 1) : (2 * p)]
		loading_VAR_EST_LAV_TEMP <- CFA_PARA_LAV$est[(2 * p + 1) : (2 * p + K * (K + 1) / 2)]
		loading_VAR_ASE_LAV_TEMP <- CFA_PARA_LAV$se[(2 * p + 1) : (2 * p + K * (K + 1) / 2)]
		loading_VAR_EST_LAV[reps, ] <- loading_VAR_EST_LAV_TEMP[c(1, 4, 5, 2, 6, 3)] 
		### change the order of elements
		loading_VAR_ASE_LAV[reps, ] <- loading_VAR_ASE_LAV_TEMP[c(1, 4, 5, 2, 6, 3)] 
		### change the order of elements
	
		F_HAT_LAV <- t(lavaan::lavPredict(CFA_LAV)) ### K by n
		R_LAV <- as.matrix(F_HAT_LAV - F)
		LOSS_LAV[reps] <- sum(apply(R_LAV, 2, norm, "2"))
		}
	
	} else {
		
		F <- t(mvrnorm(n = n, mu = MU_f, Sigma = SIGMA_f)) 
		### K by n
		U_1 <- t(mvrnorm(n = n, mu = MU_u, Sigma = SIGMA_u + E_sigma_1)) 
		### p by n
		U_3 <- t(mvrnorm(n = n, mu = MU_u, Sigma = SIGMA_u + E_sigma_3)) 
		### p by n
		U_5 <- t(mvrnorm(n = n, mu = MU_u, Sigma = SIGMA_u + E_sigma_5)) 
		### p by n
		X_1 <- matrix(rep(MU_x, n), p, n, byrow = FALSE) + L %*% F + U_1 
		### p by n
		X_3 <- matrix(rep(MU_x, n), p, n, byrow = FALSE) + L %*% F + U_3 
		### p by n
		X_5 <- matrix(rep(MU_x, n), p, n, byrow = FALSE) + L %*% F + U_5 
		### p by n
		DATA_1 <- t(X_1)
		### n by p
		DATA_3 <- t(X_3)
		### n by p
		DATA_5 <- t(X_5)
		### n by p
		
		F_HAT_MIS_1 <- solve(t(L) %*% L) %*% t(L) %*% X_1 ### K by n
		R_MIS_1 <- as.matrix(F_HAT_MIS_1 - F)
		LOSS_MIS_1[reps] <- sum(apply(R_MIS_1, 2, norm, "2"))
		RLOSS_MIS_1[reps] <- sum(apply(R_MIS_1, 2, norm, "2")) / sum(apply(F, 2, norm, "2"))
		
		F_HAT_MIS_3 <- solve(t(L) %*% L) %*% t(L) %*% X_3 ### K by n
		R_MIS_3 <- as.matrix(F_HAT_MIS_3 - F)
		LOSS_MIS_3[reps] <- sum(apply(R_MIS_3, 2, norm, "2"))
		RLOSS_MIS_3[reps] <- sum(apply(R_MIS_3, 2, norm, "2")) / sum(apply(F, 2, norm, "2"))
		
		F_HAT_MIS_5 <- solve(t(L) %*% L) %*% t(L) %*% X_5 ### K by n
		R_MIS_5 <- as.matrix(F_HAT_MIS_5 - F)
		LOSS_MIS_5[reps] <- sum(apply(R_MIS_5, 2, norm, "2"))
		RLOSS_MIS_5[reps] <- sum(apply(R_MIS_5, 2, norm, "2")) / sum(apply(F, 2, norm, "2"))
		
		S_1 <- cov(DATA_1)
		S_3 <- cov(DATA_3)
		S_5 <- cov(DATA_5)

		RES_theta_hat_1 <- BEST_UNBIASED_ESTIMATOR(S_1, p_vec)
		RES_theta_hat_3 <- BEST_UNBIASED_ESTIMATOR(S_3, p_vec)
		RES_theta_hat_5 <- BEST_UNBIASED_ESTIMATOR(S_5, p_vec)

		A_hat_1 <- RES_theta_hat_1$A
		B_hat_1 <- RES_theta_hat_1$B
		theta_hat_temp_1 <- c(diag(A_hat_1), vech(B_hat_1))
		theta_EST_1[reps, ] <- theta_hat_temp_1
		
		A_hat_3 <- RES_theta_hat_3$A
		B_hat_3 <- RES_theta_hat_3$B
		theta_hat_temp_3 <- c(diag(A_hat_3), vech(B_hat_3))
		theta_EST_3[reps, ] <- theta_hat_temp_3
		
		A_hat_5 <- RES_theta_hat_5$A
		B_hat_5 <- RES_theta_hat_5$B
		theta_hat_temp_5 <- c(diag(A_hat_5), vech(B_hat_5))
		theta_EST_5[reps, ] <- theta_hat_temp_5
	
		RES_variance_1 <- CALCULATE_VAR_A_B(A_hat_1, B_hat_1, p_vec, n)
		VAR_A_MAT_1 <- RES_variance_1$VAR_A_MAT
		VAR_B_MAT_1 <- RES_variance_1$VAR_B_MAT
		theta_VAR_temp_1 <- c(diag(VAR_A_MAT_1), vech(VAR_B_MAT_1))
		theta_ASE_1[reps, ] <- sqrt(theta_VAR_temp_1)
		
		RES_variance_3 <- CALCULATE_VAR_A_B(A_hat_3, B_hat_3, p_vec, n)
		VAR_A_MAT_3 <- RES_variance_3$VAR_A_MAT
		VAR_B_MAT_3 <- RES_variance_3$VAR_B_MAT
		theta_VAR_temp_3 <- c(diag(VAR_A_MAT_3), vech(VAR_B_MAT_3))
		theta_ASE_3[reps, ] <- sqrt(theta_VAR_temp_3)
		
		RES_variance_5 <- CALCULATE_VAR_A_B(A_hat_5, B_hat_5, p_vec, n)
		VAR_A_MAT_5 <- RES_variance_5$VAR_A_MAT
		VAR_B_MAT_5 <- RES_variance_5$VAR_B_MAT
		theta_VAR_temp_5 <- c(diag(VAR_A_MAT_5), vech(VAR_B_MAT_5))
		theta_ASE_5[reps, ] <- sqrt(theta_VAR_temp_5)
	
		CI_LOWER_temp_1 <- theta_hat_temp_1 - qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_1)
		CI_UPPER_temp_1 <- theta_hat_temp_1 + qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_1)
		CI_temp_1 <- rbind(CI_LOWER_temp_1, CI_UPPER_temp_1)
		theta_WCP_temp_1 <- 1 * (CI_temp_1[1, ] < theta0 & CI_temp_1[2, ] > theta0)
		theta_WCP_1[reps, ] <- theta_WCP_temp_1
		
		CI_LOWER_temp_3 <- theta_hat_temp_3 - qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_3)
		CI_UPPER_temp_3 <- theta_hat_temp_3 + qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_3)
		CI_temp_3 <- rbind(CI_LOWER_temp_3, CI_UPPER_temp_3)
		theta_WCP_temp_3 <- 1 * (CI_temp_3[1, ] < theta0 & CI_temp_3[2, ] > theta0)
		theta_WCP_3[reps, ] <- theta_WCP_temp_3
		
		CI_LOWER_temp_5 <- theta_hat_temp_5 - qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_5)
		CI_UPPER_temp_5 <- theta_hat_temp_5 + qnorm(1 - ALPHA / 2) * sqrt(theta_VAR_temp_5)
		CI_temp_5 <- rbind(CI_LOWER_temp_5, CI_UPPER_temp_5)
		theta_WCP_temp_5 <- 1 * (CI_temp_5[1, ] < theta0 & CI_temp_5[2, ] > theta0)
		theta_WCP_5[reps, ] <- theta_WCP_temp_5
		
	}
		
	cat(" iteration:  ", reps, "\r")
	reps <- reps + 1
	}, error = function(e){})

}

save.image(paste("SCFA_SET_NO", SET_NO, "repsMax_100.RData", sep = "_"))

if(SET_NO %in% c(1, 2, 4, 7)){

	### Computing results of Goal 1 for SET_NO: 1,2,4,7 ###
	
	round(c(mean(LOSS_UB), sd(LOSS_UB), TIME_UB, 
			mean(LOSS_SEM), sd(LOSS_SEM), TIME_SEM,
			mean(LOSS_MX), sd(LOSS_MX), TIME_MX, 
			mean(LOSS_LAV), sd(LOSS_LAV), TIME_LAV), 2)
	
	### Computing results of Goal 2 for SET_NO: 1,2,4,7 ###
	
	bias <- apply(theta_EST, 2, mean) - theta0
	MCSD <- apply(theta_EST, 2, sd)
	ASE <- apply(theta_ASE, 2, mean)
	WCP <- apply(theta_WCP, 2, mean)

	round(cbind(bias, MCSD, ASE, WCP) * 100, 1)

	L_MAT <- cbind(apply(loading_VEC_EST_SEM, 2, mean) - rep(1, p), 
				   apply(loading_VEC_EST_SEM, 2, sd),
				   apply(loading_VEC_ASE_SEM, 2, mean), 
			       apply(loading_VEC_EST_MX, 2, mean) - rep(1, p), 
			       apply(loading_VEC_EST_MX, 2, sd),
			       apply(loading_VEC_ASE_MX, 2, mean),
			       apply(loading_VEC_EST_LAV, 2, mean) - rep(1, p), 
			       apply(loading_VEC_EST_LAV, 2, sd),
			       apply(loading_VEC_ASE_LAV, 2, mean))

	round(L_MAT * 100, 2)

	SIGMA_U_MAT <- cbind(rep(apply(theta_EST[,1:K],2,mean)-theta_A0,p_vec),
						 rep(apply(theta_EST[,1:K],2,sd), p_vec),
					     rep(apply(theta_ASE[,1:K],2,mean), p_vec),
					     apply(error_VEC_EST_SEM,2,mean) - rep(theta_A0,p_vec), 
					     apply(error_VEC_EST_SEM,2,sd),
						 apply(error_VEC_ASE_SEM,2,mean),
					     apply(error_VEC_EST_MX,2,mean) - rep(theta_A0, p_vec), 
					     apply(error_VEC_EST_MX,2,sd),
					     apply(error_VEC_ASE_MX,2,mean),
					     apply(error_VEC_EST_LAV,2,mean) - rep(theta_A0,p_vec), 
					     apply(error_VEC_EST_LAV,2,sd),
					     apply(error_VEC_ASE_LAV,2,mean))

	round(SIGMA_U_MAT * 100, 2)

	SIGMA_F_MAT <- cbind(apply(theta_EST[,- (1 : K)], 2, mean) - theta_B0,
						 apply(theta_EST[,- (1 : K)], 2, sd),
						 apply(theta_ASE[,- (1 : K)], 2, mean),
						 apply(loading_VAR_EST_SEM, 2, mean) - theta_B0,
						 apply(loading_VAR_EST_SEM, 2, sd), 
						 apply(loading_VAR_ASE_SEM, 2, mean),
						 apply(loading_VAR_EST_MX, 2, mean) - theta_B0,
						 apply(loading_VAR_EST_MX, 2, sd), 
						 apply(loading_VAR_ASE_MX, 2, mean),
						 apply(loading_VAR_EST_LAV, 2, mean) - theta_B0,
						 apply(loading_VAR_EST_LAV, 2, sd), 
						 apply(loading_VAR_ASE_LAV, 2, mean))

	round(SIGMA_F_MAT * 100, 2)

} else if(SET_NO %in% c(3, 5)){

	### Computing results of Goal 1 for SET_NO: 3,5 ###

	round(c(mean(LOSS_UB), sd(LOSS_UB), TIME_UB, 
			NA, NA, NA,
			mean(LOSS_MX), sd(LOSS_MX), TIME_MX, 
			NA, NA, NA), 2)
	
	### Computing results of Goal 2 for SET_NO: 3,5 ###

	bias <- apply(theta_EST, 2, mean) - theta0
	MCSD <- apply(theta_EST, 2, sd)
	ASE <- apply(theta_ASE, 2, mean)
	WCP <- apply(theta_WCP, 2, mean)

	round(cbind(bias, MCSD, ASE, WCP) * 100, 1)

	L_MAT <- cbind(rep(NA, p), rep(NA, p), rep(NA, p), 			   
			       apply(loading_VEC_EST_MX, 2, mean) - rep(1, p), 
			       apply(loading_VEC_EST_MX, 2, sd), 
			       apply(loading_VEC_ASE_MX, 2, mean),
			       rep(NA, p), rep(NA, p), rep(NA, p))

	round(L_MAT * 100, 2)

	SIGMA_U_MAT <- cbind(rep(apply(theta_EST[,1:K],2,mean)-theta_A0, p_vec),
						 rep(apply(theta_EST[,1:K],2,sd), p_vec), 
						 rep(apply(theta_ASE[,1:K],2,mean), p_vec),
						 rep(NA, p), rep(NA, p), rep(NA, p),
						 apply(error_VEC_EST_MX,2,mean) - rep(theta_A0, p_vec), 
						 apply(error_VEC_EST_MX,2,sd),
						 apply(error_VEC_ASE_MX,2,mean),
						 rep(NA, p), rep(NA, p), rep(NA, p))

	round(SIGMA_U_MAT * 100, 2)

	SIGMA_F_MAT <- cbind(apply(theta_EST[,- (1 : K)], 2, mean) - theta_B0,
						 apply(theta_EST[,- (1 : K)], 2, sd),
						 apply(theta_ASE[,- (1 : K)], 2, mean),
						 rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2),
						 apply(loading_VAR_EST_MX, 2, mean) - theta_B0,
						 apply(loading_VAR_EST_MX, 2, sd), 
						 apply(loading_VAR_ASE_MX, 2, mean),
						 rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2), 
						 rep(NA, K * (K + 1) / 2))

	round(SIGMA_F_MAT * 100, 2)

} else if(SET_NO %in% c(6, 8, 9)){

	### Computing results of Goal 1 for SET_NO: 6, 8, 9 ###

	round(c(mean(LOSS_UB), sd(LOSS_UB), TIME_UB, 
			NA, NA, NA,
			NA, NA, NA,
			NA, NA, NA), 2)

	### Computing results of Goal 2 for SET_NO: 6, 8, 9 ###

	bias <- apply(theta_EST, 2, mean) - theta0
	MCSD <- apply(theta_EST, 2, sd)
	ASE <- apply(theta_ASE, 2, mean)
	WCP <- apply(theta_WCP, 2, mean)

	round(cbind(bias, MCSD, ASE, WCP) * 100, 1)

	L_MAT <- cbind(rep(NA, p), rep(NA, p), rep(NA, p), rep(NA, p),
				   rep(NA, p), rep(NA, p), rep(NA, p), rep(NA, p), 
			       rep(NA, p))

	round(L_MAT * 100, 2)

	SIGMA_U_MAT <- cbind(rep(apply(theta_EST[,1:K],2,mean)-theta_A0, p_vec),
						 rep(apply(theta_EST[,1:K],2,sd), p_vec), 
						 rep(apply(theta_ASE[,1:K],2,mean), p_vec),
						 rep(NA, p), rep(NA, p), rep(NA, p), rep(NA, p), 
						 rep(NA, p), rep(NA, p), rep(NA, p), rep(NA, p), 
						 rep(NA, p))

	round(SIGMA_U_MAT * 100, 2)

	SIGMA_F_MAT <- cbind(apply(theta_EST[,- (1 : K)], 2, mean) - theta_B0,
						 apply(theta_EST[,- (1 : K)], 2, sd),
						 apply(theta_ASE[,- (1 : K)], 2, mean),
						 rep(NA, K * (K + 1) / 2), rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2), rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2), rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2), rep(NA, K * (K + 1) / 2),
						 rep(NA, K * (K + 1) / 2))

	round(SIGMA_F_MAT * 100, 2)

} else if(SET_NO == 10){

	### Computing results of Goal 3 for SET_NO: 10 ###

	bias_1 <- apply(theta_EST_1, 2, mean) - theta0
	MCSD_1 <- apply(theta_EST_1, 2, sd)
	ASE_1 <- apply(theta_ASE_1, 2, mean)
	WCP_1 <- apply(theta_WCP_1, 2, mean)

	bias_3 <- apply(theta_EST_3, 2, mean) - theta0
	MCSD_3 <- apply(theta_EST_3, 2, sd)
	ASE_3 <- apply(theta_ASE_3, 2, mean)
	WCP_3 <- apply(theta_WCP_3, 2, mean)

	bias_5 <- apply(theta_EST_5, 2, mean) - theta0
	MCSD_5 <- apply(theta_EST_5, 2, sd)
	ASE_5 <- apply(theta_ASE_5, 2, mean)
	WCP_5 <- apply(theta_WCP_5, 2, mean)

	round(cbind(bias_1, MCSD_1, ASE_1, WCP_1) * 100, 1)
	round(cbind(bias_3, MCSD_3, ASE_3, WCP_3) * 100, 1)
	round(cbind(bias_5, MCSD_5, ASE_5, WCP_5) * 100, 1)

	round(c(mean(LOSS_MIS_1), sd(LOSS_MIS_1), 
		    mean(LOSS_MIS_3), sd(LOSS_MIS_3), 
		    mean(LOSS_MIS_5), sd(LOSS_MIS_5)), 2)
	
	round(c(mean(RLOSS_MIS_1), sd(RLOSS_MIS_1), 
		    mean(RLOSS_MIS_3), sd(RLOSS_MIS_3), 
		    mean(RLOSS_MIS_5), sd(RLOSS_MIS_5)), 2)	
}
