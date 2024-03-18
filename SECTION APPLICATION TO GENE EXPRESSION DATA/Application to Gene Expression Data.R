###########################################
### Application to Gene Expression Data ###
###########################################

############################
### Loading the packages ###
############################

PACKAGES <- c(
			 "R.matlab",		
			 ### convert between R and MATLAB files
			 "xlsx",
			 ### xlsx::write.xlsx() saves the results in a sheet
			 "matrixcalc", 
			 ### matrixcalc::vech() creates a vector from a symmetric matrix
			 "sem",       
			 ### estimates a CFA model
			 "lavaan", 	
			 ### estimates a CFA model
			 "OpenMx"		
			 ### estimates a CFA model
			 )

CHECK_PACKAGES <- lapply(X = PACKAGES,
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
}CALCULATE_VAR_A_B <- function(A, B, p_vec, n){
	
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

############################################
### Loading the real-world data examples ###
############################################

INPUT_ADDRESS <- "C:/Users/Yifan Yang/Desktop/SUB 3 ADAPTIVE ESTIMATION/SIMULATION AND REAL DATA/REAL DATA/VERSION 1.2.0"

OUTPUT_ADDRESS <- INPUT_ADDRESS
DATA_INPUT_FILENAME <- paste(INPUT_ADDRESS, "TCGA_Kidney_Gene_R_ALL_2448_1777.mat", sep = "/")
Data_Example <- R.matlab::readMat(DATA_INPUT_FILENAME)
Y_Data <- Data_Example$TCGA.KIDNEY.Y.2448.1777
### n by p
Working_Data <- Data_Example$TCGA.Kidney.Gene.R.ALL.2448.1777
### p by p
Variable_2448_1777_Name <- as.vector(unlist(Data_Example$TCGA.Kidney.Gene.R.ALL.2448.1777.NAMES))
Variable_Raw_Name <- as.vector(unlist(Data_Example$TCGA.Kidney.Gene.R.RAW.NAMES))
p_vec <- drop(Data_Example$p.TCGA.Kidney.Gene)
n <- 712 
p_cumsum <- cumsum(p_vec)
Members <- rep(0, length(Variable_2448_1777_Name))
for(k in 1 : length(p_vec)){
	Members[ifelse(length(p_cumsum[k - 1]) == 0, 1, p_cumsum[k - 1] + 1) : p_cumsum[k]] <- paste("Community", k, sep = "_")
}
Membership <- cbind(Variable_2448_1777_Name, Members)

#################################################
### Fitting the proposed model to the dataset ###
#################################################

K <- length(p_vec)
p <- sum(p_vec)
sig_level <- 0.05

### Estimate model parameters ###

RES_EST <- BEST_UNBIASED_ESTIMATOR(Working_Data, p_vec)
A_EST <- RES_EST$A
B_EST <- RES_EST$B
theta_EST <- c(diag(A_EST), drop(matrixcalc::vech(B_EST)))
RES_SE <- CALCULATE_VAR_A_B(A_EST, B_EST, p_vec, n)

A_SE_EST <- RES_SE$VAR_A_MAT
B_SE_EST <- RES_SE$VAR_B_MAT
theta_ASE <- sqrt(c(diag(A_SE_EST), drop(matrixcalc::vech(B_SE_EST))))
theta_LB <- theta_EST - qnorm(1 - sig_level / 2) * theta_ASE
theta_UB <- theta_EST + qnorm(1 - sig_level / 2) * theta_ASE
round(cbind(theta_EST, theta_ASE, theta_LB, theta_UB), 4)

xlsx::write.xlsx(
		  x = Membership,
	   file = paste(OUTPUT_ADDRESS, "Application to Gene Expression Data.xlsx", sep = "/"),
  sheetName = "TCGA_Membership",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

### Estimate the factor scores ###
L <- c()
for(k in 1 : K){
	e_temp <- rep(0, K)
	e_temp[k] <- 1
	L <- cbind(L, rep(e_temp, p_vec))
}
L <- as.matrix(L)
### p by K
YT_Data <- t(as.matrix(Y_Data))
### p by n
YT_mean <- apply(YT_Data, 1, mean)
### p by 1
YT_Data_Cent <- YT_Data - matrix(rep(YT_mean, n), p, n, byrow = FALSE)
F_HAT_UB <- solve(t(L) %*% L) %*% t(L) %*% YT_Data_Cent 
### K by n, each column refers to an estimated factor score
### number of factor scores = sample size

xlsx::write.xlsx(
		  x = F_HAT_UB,
	   file = paste(OUTPUT_ADDRESS, "Application to Gene Expression Data.xlsx", sep = "/"),
  sheetName = "Factor scores",
  col.names = TRUE,
  row.names = TRUE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

########################
### Pathway Analysis ###
########################

source(paste(INPUT_ADDRESS, "Path.fisher.R", sep = "/"))
load(paste(INPUT_ADDRESS, "pathways.rda", sep = "/"))

for(k in 1 : K){

	Path_k <- Path.fisher(
				x = Variable_2448_1777_Name[ifelse(length(p_cumsum[k - 1]) == 0, 1, p_cumsum[k - 1] + 1) : p_cumsum[k]], 
	   background = Variable_Raw_Name, 
          pathway = c(GOBP.genesets, GOMF.genesets, GOCC.genesets, KEGG.genesets, Reactome.genesets), #Oncogenic.genesets,
          top_num = 500)

	xlsx::write.xlsx(
				x = Path_k,
			 file = paste(OUTPUT_ADDRESS, "PathwaysAnalysis.xlsx", sep = "/"),
		sheetName = paste("Community", k, sep = "_"),
		col.names = TRUE,
		row.names = TRUE,
		   append = TRUE,
		   showNA = TRUE,
		 password = NULL
	)
}


######################################################
### Fitting the conventional models to the dataset ###
######################################################

DATA <- Y_Data
### n by p

X_NAMES <- paste(rep('x', p_vec[1]), seq(p_vec[1]), sep = '')
Y_NAMES <- paste(rep('y', p_vec[2]), seq(p_vec[2]), sep = '')
Z_NAMES <- paste(rep('z', p_vec[3]), seq(p_vec[3]), sep = '')
W_NAMES <- paste(rep('w', p_vec[4]), seq(p_vec[4]), sep = '')
S_NAMES <- paste(rep('s', p_vec[5]), seq(p_vec[5]), sep = '')
Q_NAMES <- paste(rep('q', p_vec[6]), seq(p_vec[6]), sep = '')

COL_NAMES <- c(X_NAMES, Y_NAMES, Z_NAMES, W_NAMES, S_NAMES, Q_NAMES)
colnames(DATA) <- COL_NAMES
ERR_NAMES <- paste(rep('e', p), seq(p), sep = '')
LOD_NAMES <- paste(rep('l', p), seq(p), sep = '')
MU_NAMES <- paste(rep('mean', p), COL_NAMES, sep = '')

### sem() functions in package "sem"  ###

DATA_SEM <- data.frame(DATA)
X_NAMES_SEM <- paste(X_NAMES, collapse = ', ')	
Y_NAMES_SEM <- paste(Y_NAMES, collapse = ', ')
Z_NAMES_SEM <- paste(Z_NAMES, collapse = ', ') 
W_NAMES_SEM <- paste(W_NAMES, collapse = ', ')	
S_NAMES_SEM <- paste(S_NAMES, collapse = ', ')
Q_NAMES_SEM <- paste(Q_NAMES, collapse = ', ') 
FX_SEM <- paste("FX_SEM: ", X_NAMES_SEM, sep = '')
FY_SEM <- paste("FY_SEM: ", Y_NAMES_SEM, sep = '')
FZ_SEM <- paste("FZ_SEM: ", Z_NAMES_SEM, sep = '')
FW_SEM <- paste("FW_SEM: ", W_NAMES_SEM, sep = '')
FS_SEM <- paste("FS_SEM: ", S_NAMES_SEM, sep = '')
FQ_SEM <- paste("FQ_SEM: ", Q_NAMES_SEM, sep = '')

suppressMessages(
	MODEL_SEM <- cfa_MODIFIED(text = paste(FX_SEM, FY_SEM, FZ_SEM, FW_SEM, FS_SEM, FQ_SEM, sep = '\n'), reference.indicators = TRUE))
		
### TIME_TEMP <- Sys.time()
CFA_SEM <- sem::sem(MODEL_SEM, data = DATA_SEM)
### Error in solve.default(S)
### CFA_SEM <- sem::sem(MODEL_SEM, S = cor(DATA, use = "pairwise.complete.obs"), N = n)
### Error in solve.default(S)
### TIME_SEM <- TIME_SEM + (Sys.time() - TIME_TEMP)

### mxRun() functions in package "OpenMx"  ###

DATA_MX <- data.frame(DATA)
DATA_MX_RAW <- mxData(observed = DATA_MX, type = "raw")
resVars <- OpenMx::mxPath(
				  from = COL_NAMES,
				arrows = 2,
                  free = TRUE, 
				values = rep(1, p),
                labels = ERR_NAMES)
latVars <- OpenMx::mxPath(
				  from = c("FX_MX", "FY_MX", "FZ_MX", 
						   "FW_MX", "FS_MX", "FQ_MX"), 
				arrows = 2, 
			   connect = "unique.pairs",
                  free = TRUE, 
				values = c(1, 0.5, 0.5, 0.5, 0.5, 0.5,
							    1, 0.5, 0.5, 0.5, 0.5,
									 1, 0.5, 0.5, 0.5,
									      1, 0.5, 0.5,
										       1, 0.5,
											        1), 
				labels = c("varFX","covXY","covXZ","covXW","covXS","covXQ",
						           "varFY","covYZ","covYW","covYS","covYQ",
										   "varFZ","covZW","covZS","covZQ", 
										           "varFW","covWS","covWQ",
												           "varFS","covWQ",
												                   "varFQ"))
facLoadsX <- OpenMx::mxPath(
					from = "FX_MX", 
					  to = X_NAMES, 
				  arrows = 1,
                    free = c(FALSE, rep(TRUE, p_vec[1] - 1)), 
				  values = rep(1, p_vec[1]), 
				  labels = LOD_NAMES[1 : cumsum(p_vec)[1]])
		
facLoadsY <- OpenMx::mxPath(
					from = "FY_MX", 
					  to = Y_NAMES, 
				  arrows = 1,
					free = c(FALSE, rep(TRUE, p_vec[2] - 1)),
				  values = rep(1, p_vec[2]), 
				  labels = LOD_NAMES[(cumsum(p_vec)[1] + 1) : cumsum(p_vec)[2]])

facLoadsZ <- OpenMx::mxPath(
					from = "FZ_MX", 
					  to = Z_NAMES, 
				  arrows = 1,
					free = c(FALSE, rep(TRUE, p_vec[3] - 1)),
				  values = rep(1, p_vec[3]), 
				  labels = LOD_NAMES[(cumsum(p_vec)[2] + 1) : cumsum(p_vec)[3]])

facLoadsW <- OpenMx::mxPath(
					from = "FW_MX", 
					  to = W_NAMES, 
				  arrows = 1,
					free = c(FALSE, rep(TRUE, p_vec[4] - 1)),
				  values = rep(1, p_vec[4]), 
				  labels = LOD_NAMES[(cumsum(p_vec)[3] + 1) : cumsum(p_vec)[4]])

facLoadsS <- OpenMx::mxPath(
					from = "FS_MX", 
					  to = S_NAMES, 
				  arrows = 1,
					free = c(FALSE, rep(TRUE, p_vec[5] - 1)),
				  values = rep(1, p_vec[5]), 
				  labels = LOD_NAMES[(cumsum(p_vec)[4] + 1) : cumsum(p_vec)[5]])

facLoadsQ <- OpenMx::mxPath(
					from = "FQ_MX", 
					  to = Q_NAMES, 
				  arrows = 1,
					free = c(FALSE, rep(TRUE, p_vec[6] - 1)),
				  values = rep(1, p_vec[6]), 
				  labels = LOD_NAMES[(cumsum(p_vec)[5] + 1) : cumsum(p_vec)[6]])

means <- OpenMx::mxPath(
				from = "one", 
				  to = c(COL_NAMES, "FX_MX", "FY_MX", "FZ_MX", "FW_MX", "FS_MX", "FQ_MX"), 
              arrows = 1,
                free = c(rep(TRUE, p), rep(FALSE, 6)), 
		      values = c(rep(1, p), rep(0, 6)),
              labels = c(MU_NAMES, NA, NA, NA, NA, NA, NA)) 

MODEL_MX <- OpenMx::mxModel(
					model = "MODEL_MX", 
				     type = "RAM",
			 manifestVars = COL_NAMES, 
               latentVars = c("FX_MX","FY_MX", "FZ_MX", "FW_MX","FS_MX", "FQ_MX"),
            DATA_MX_RAW, resVars, latVars, facLoadsX, facLoadsY, facLoadsZ, means)
		
TIME_TEMP <- Sys.time()
CFA_MX <- OpenMx::mxRun(MODEL_MX, silent = TRUE)
TIME_MX <- TIME_MX + (Sys.time() - TIME_TEMP)

### cfa() function in package lavaan  ###

X_NAME_LAV <- paste(COL_NAMES[1 : sum(p_vec[1 : 1])], collapse = ' + ')
Y_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 1]) + 1) : sum(p_vec[1 : 2])], collapse = ' + ')
Z_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 2]) + 1) : sum(p_vec[1 : 3])], collapse = ' + ') 
W_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 3]) + 1) : sum(p_vec[1 : 4])], collapse = ' + ')		
S_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 4]) + 1) : sum(p_vec[1 : 5])], collapse = ' + ')
Q_NAME_LAV <- paste(COL_NAMES[(sum(p_vec[1 : 5]) + 1) : sum(p_vec[1 : 6])], collapse = ' + ') 

FX_LAV <- paste(' FX_LAV =~ ', X_NAME_LAV, sep = '')
FY_LAV <- paste(' FY_LAV =~ ', Y_NAME_LAV, sep = '')
FZ_LAV <- paste(' FZ_LAV =~ ', Z_NAME_LAV, sep = '')
FW_LAV <- paste(' FW_LAV =~ ', W_NAME_LAV, sep = '')
FS_LAV <- paste(' FS_LAV =~ ', S_NAME_LAV, sep = '')
FQ_LAV <- paste(' FQ_LAV =~ ', Q_NAME_LAV, sep = '')

MODEL_LAV <- paste(FX_LAV, FY_LAV, FZ_LAV, FW_LAV, FS_LAV, FQ_LAV, sep = " \n ")
		
### TIME_TEMP <- Sys.time()
CFA_LAV <- lavaan::cfa(MODEL_LAV, data = DATA)
### lavaan error: sample covariance matrix is not positive-definite
### TIME_LAV <- TIME_LAV + (Sys.time() - TIME_TEMP)


