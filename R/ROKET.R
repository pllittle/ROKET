# Package Functions

# ----------
# Minor Functions
# ----------
setdirs = function(work_dir){
	
	# Make directories
	verbose = TRUE
	if( verbose ) cat(sprintf("%s: Make directories ...\n",date()))
	proj_dir 	= file.path(work_dir,"sim_ROKET"); 	smart_mkdir(proj_dir)
	sim_dir		= file.path(proj_dir,"sim"); 				smart_mkdir(sim_dir)
	reps_dir 	= file.path(sim_dir,"REPS"); 				smart_mkdir(reps_dir)
	rout_dir	= file.path(sim_dir,"REG");					smart_mkdir(rout_dir)
	sout_dir	= file.path(sim_dir,"OUT");					smart_mkdir(sout_dir)
	
	# Output
	if( verbose ) cat(sprintf("%s: Output ...\n",date()))
	list(proj_dir = proj_dir,sim_dir = sim_dir,reps_dir = reps_dir,
		rout_dir = rout_dir,sout_dir = sout_dir)
	
}


# ----------
# Optimal Transport Functions
# ----------

#' @title run_myOT
#' @description Runs balanced or unbalanced optimal transport
#'	on two input vectors
#' @param XX A numeric vector of positive masses
#' @param YY A numeric vector of positive masses
#' @param COST A numeric matrix of non-negative values
#'	representing the costs to transport masses between
#'	features of \code{XX} and \code{YY}. The rows of \code{COST}
#'	and features of \code{XX} need to be aligned.
#'	The columns of \code{COST} and features of \code{YY}
#'	need to be aligned.
#' @param EPS A positive numeric value representing the
#'	tuning parameter for entropic regularization.
#' @param LAMBDA1 A non-negative numeric value representing
#'	the tuning parameter penalizing the distance between \code{XX}
#'	and the row sums of the optimal transport matrix.
#' @param LAMBDA2 A non-negative numeric value representing
#'	the tuning parameter penalizing the distance between \code{YY}
#'	and the column sums of the optimal transport matrix.
#' @param balance Boolean set to \code{TRUE} to run balanced
#'	optimal transport. Otherwise run unbalanced optimal transport.
#' @param conv A positive numeric value to determine 
#'	algorithmic convergence. The default value is \code{1e-5}.
#' @param max_iter A positive integer denoting the maximum
#'	iterations to run the algorithm.
#' @param show Boolean value to display verbose algorithm output.
#' @param show_iter A positive integer to display iteration details
#'	at multiples of \code{show_iter} but only if \code{show = TRUE}.
#' 
#' @export
run_myOT = function(XX,YY,COST,EPS,LAMBDA1,LAMBDA2,
	balance = FALSE,conv = 1e-5,max_iter = 3e3,
	show = TRUE,show_iter = 50){
	
	nXX = names(XX)
	nYY = names(YY)
	nrCOST = rownames(COST)
	ncCOST = colnames(COST)
	
	# Check inputs
	if( is.null(nXX) ) stop("Set names(XX)")
	if( is.null(nYY) ) stop("Set names(YY)")
	if( is.null(ncCOST) ) stop("Set colnames(COST)")
	if( is.null(nrCOST) ) stop("Set rownames(COST)")
	if( length(XX) != nrow(COST) ) stop("length(XX) != nrow(COST)")
	if( length(YY) != ncol(COST) ) stop("length(YY) != ncol(COST)")
	if( !all(nXX == nrCOST) ) stop("row COST name mismatch")
	if( !all(nYY == ncCOST) ) stop("column COST name mismatch")
	
	# If there are zeros in XX or YY, subset
	XX2 = XX[XX > 0]
	YY2 = YY[YY > 0]
	COST_XY = COST[names(XX2),names(YY2)]
	
	if( sum(XX2) <= 0 ) stop("XX needs to have mass > 0")
	if( sum(YY2) <= 0 ) stop("YY needs to have mass > 0")
	
	# Run OT
	OT = Rcpp_run_OT(XX = XX2,YY = YY2,COST_XY = COST_XY,
		EPS = EPS,LAMBDA1 = LAMBDA1,LAMBDA2 = LAMBDA2,
		balance = balance,highLAM_lowMU = TRUE,
		conv = conv,max_iter = max_iter,show = show,
		show_iter = show_iter)
	
	# Add dimnames
	OT_fin = matrix(0,length(XX),length(YY))
	OT_fin = smart_names(MAT = OT_fin,
		ROW = nXX,COL = nYY)
	OT_fin[names(XX2),names(YY2)] = OT
	
	# Calculate OT metrics
	DIST_1 = sum(COST * OT_fin)
	sum_OT = sum(OT_fin)
	DIST_2 = DIST_1 / sum_OT
	
	out = list(OT = OT_fin,DIST_1 = DIST_1,
		sum_OT = sum_OT,DIST_2 = DIST_2)
	return(out)
	
}

#' @title run_myOTs
#' @param ZZ A numeric matrix of non-negative mass to transport.
#'	Rows correspond to features (e.g. genes) and columns
#'	correspond to samples or individuals. Each column must have
#'	strictly positive mass
#' @param COST A numeric square matrix of non-negative values
#'	representing the non-negative costs to transport 
#'	masses between pairs of features
#' @param EPS A positive numeric value representing the
#'	tuning parameter for entropic regularization.
#' @param LAMBDA1 A non-negative numeric value representing
#'	the tuning parameter penalizing the distance between \code{XX}
#'	and the row sums of the optimal transport matrix.
#' @param LAMBDA2 A non-negative numeric value representing
#'	the tuning parameter penalizing the distance between \code{YY}
#'	and the column sums of the optimal transport matrix.
#' @param balance Boolean set to \code{TRUE} to run balanced
#'	optimal transport. Otherwise run unbalanced optimal transport.
#' @param conv A positive numeric value to determine 
#'	algorithmic convergence. The default value is \code{1e-5}.
#' @param max_iter A positive integer denoting the maximum
#'	iterations to run the algorithm.
#' @param ncores A positive integer for the number of cores/threads
#'	to reduce computational runtime when running a for loop
#' @param show Boolean value to display verbose algorithm output.
#' @param show_iter A positive integer to display iteration details
#'	at multiples of \code{show_iter} but only if \code{show = TRUE}.
#' @export
run_myOTs = function(ZZ,COST,EPS,LAMBDA1,LAMBDA2,
	balance,conv = 1e-5,max_iter = 3e3,ncores = 1,
	show = TRUE,show_iter = 50){
	
	# Check inputs
	mass = apply(ZZ,2,function(xx){
		sum(xx[xx > 0])
	})
	if( any(mass <= 0) )
		stop("A subset of samples have non-positive mass")
	if( nrow(COST) != ncol(COST) ) stop("COST not square")
	if( !all(COST == t(COST)) ) stop("COST not symmetric")
	
	# Subset features of ZZ with any mass
	mass_rows = apply(ZZ,1,function(xx){
		sum(xx[xx > 0])
	})
	ZZ = ZZ[mass_rows > 0,,drop = FALSE]
	nZZ = rownames(ZZ)
	nCOST = rownames(COST)
	if( nrow(ZZ) > nrow(COST) )
		stop("More features in ZZ than features in COST")
	if( !all(nZZ %in% nCOST) )
		stop("Some features in ZZ not among features in COST")
	COST = COST[nZZ,nZZ,drop = FALSE]
	
	# Run OT across all individuals
	out_OT = Rcpp_run_full_OT(COST = COST,ZZ = ZZ,
		EPS = EPS,LAMBDA1 = LAMBDA1,LAMBDA2 = LAMBDA2,
		balance = balance,highLAM_lowMU = TRUE,
		conv = conv,max_iter = max_iter,ncores = ncores,
		show = show,show_iter = show_iter)
	
}


# ----------
# Optimal transport simulation
# ----------
OT_sim = function(){
	print("Add simulation code eventually")
}


# ----------
# Simulation functions
# ----------
kOT_sim_geneInfo = function(nGENE,nPATH,path_corr,prop_noise){
	
	# Group genes into pathways
	
	# Code if genes are exclusive to certain paths
	all_PATHS = seq(nPATH)
	PATH_prob = exp(seq(nPATH,1)*0.1)
	PATH_prob = PATH_prob / sum(PATH_prob)
	
	gene_names = paste0("G",seq(nGENE))
	gdat = smart_df(GENE = gene_names,
		PATH = sample(all_PATHS,nGENE,replace = TRUE,prob = PATH_prob))
	dim(gdat); gdat[1:5,]
	table(gdat$PATH)
	
	# Simulate gene similarities based on underlying pathways (range between 0 and 1)
	gene_sim = matrix(NA,nGENE,nGENE)
	gene_sim = smart_names(gene_sim,ROW = gene_names,COL = gene_names)
	# diag(gene_sim) = 1
	
	for(pp1 in seq(nPATH)){
	for(pp2 in seq(nPATH)){
		# pp1 = 1; pp2 = 2
		
		# Only fill in the upper triangle
		PATH1 = all_PATHS[pp1]; PATH1
		PATH2 = all_PATHS[pp2]; PATH2
		if( PATH1 > PATH2 ) next
		
		idx1 = which(gdat$PATH == PATH1)
		idx2 = which(gdat$PATH == PATH2)
		ngenes1 = length(idx1); ngenes1
		ngenes2 = length(idx2); ngenes2
		
		if( min(c(PATH1,PATH2)) > 0 ){
			if( PATH1 == PATH2 ){
				sim_range = sort(runif(2,path_corr,1))
				gene_sim[idx1,idx2] = runif(ngenes1*ngenes2,sim_range[1],sim_range[2])
			} else {
				# Pairs of genes not in same pathway
				sim_range = sort(runif(2,0,1-path_corr))
				gene_sim[idx1,idx2] = runif(ngenes1*ngenes2,sim_range[1],sim_range[2])
				gene_sim[idx2,idx1] = runif(ngenes1*ngenes2,sim_range[1],sim_range[2])
			}
		}
		
		rm(idx1,idx2)
	}}
	
	diag(gene_sim) = 1
	gene_sim[lower.tri(gene_sim)] = t(gene_sim)[lower.tri(gene_sim)]
	# any(is.na(gene_sim))
	
	# Add some noise
	gs_noise = matrix(runif(nGENE^2,0,1),nGENE,nGENE)
	gs_noise[lower.tri(gs_noise)] = t(gs_noise)[lower.tri(gs_noise)]
	diag(gs_noise) = 1
	# prop_noise = 0.35
	gene_sim2 = (1 - prop_noise) * gene_sim + prop_noise * gs_noise
	
	return(list(gene_sim = gene_sim2,gdat = gdat))
	
}
kOT_sim_baseCov = function(subj_names){
	
	verbose = TRUE
	if( verbose ) cat(sprintf("%s: Simulate baseline covariates ...\n",date()))
	NN = length(subj_names)
	XX = matrix(data = 1,nrow = NN,
		ncol = 1,dimnames = list(seq(NN),"Int"))
	XX = cbind(XX,make_dummy(x = sample(c("M","F"),NN,replace = TRUE)))
	XX = cbind(XX,make_dummy(x = sample(c("I","II","III","IV"),NN,replace = TRUE)))
	XX = cbind(XX,age = scale(sample(seq(20,80),NN,replace = TRUE))[,1])
	rownames(XX) = subj_names
	XX = as.matrix(XX)
	
	return(XX)
}
kOT_sim_sPATH = function(subj_names,path_names,SCEN){
	
	verbose = TRUE
	if( verbose ) cat(sprintf("%s: Simulate individual pathways mutated ...\n",date()))
	NN 				= length(subj_names)
	nPATH 		= length(path_names)
	vec_path 	= seq(nPATH)
	sPATH = matrix(NA,NN,nPATH,dimnames = list(subj_names,path_names))
	
	if( SCEN %in% c(1,2) ){
		# Approx constant PB mass
		PB = sample(seq(9,11),NN,replace = TRUE)
	} else if( SCEN %in% c(3,4) ){
		# Variable PB mass
		PB = sample(seq(5,20),NN,replace = TRUE)
	}
	PB[PB >= nPATH] = nPATH
	PB[PB <= 5] = 5
	
	for(ii in seq(NN)){
		# ii = 1
		
		subj_path_set = c(rep(0,nPATH - PB[ii]),rep(1,PB[ii]))
		
		# Shuffle paths
		sPATH[ii,] = sample(subj_path_set)
	}
	dim(sPATH); sPATH[1:5,]
	
	return(sPATH)
}
kOT_sim_sMUT = function(sPATH,geneInfo,SCEN){
	
	verbose = TRUE
	if( verbose ) cat(sprintf("%s: Simulate individual genes mutated ...\n",date()))
	NN = nrow(sPATH)
	gene_names = rownames(geneInfo$gene_sim)
	nGENE = length(gene_names)
	sMUT = matrix(0,nGENE,NN)
	sMUT = smart_names(MAT = sMUT,
		ROW = gene_names,
		COL = rownames(sPATH))
	
	for(ii in seq(NN)){
		# ii = 1
		
		# Determine all pathways mutated and sample from those genes
		# sPATH[ii,]
		subj_paths = which(sPATH[ii,] == 1)
		subj_paths = as.integer(gsub("P","",names(subj_paths)))
		subj_paths
		tmp_df = smart_df(PATH = subj_paths)
		
		prob_paths = rep(1,length(subj_paths))
		prob_paths = prob_paths / sum(prob_paths)
		tmp_df$PROB = prob_paths
		
		# Sample from individual's mutated and null paths
		if( SCEN %in% c(1,3) ){
			# TMB is about 2 or 3 times PMB
			subj_TMB = round(length(subj_paths) * 3)
		} else if( SCEN %in% c(2,4) ){
			# TMB is varied multiple of PMB
			subj_TMB = sample(seq(length(subj_paths) + 5,8e1,2),1)
		}
		samp_paths = sample(subj_paths,subj_TMB,
			replace = TRUE,prob = prob_paths)
		subj_tab = table(samp_paths); subj_tab
		nsubj_paths = length(subj_tab)
		
		tmp_df$CNT = sapply(tmp_df$PATH,function(xx)
			sum(samp_paths == xx))
		if( min(tmp_df$CNT) == 0 ){
			tmp_df$CNT[tmp_df$CNT == 0] = 1
		}
		tmp_df
		
		subj_genes = c()
		for(jj in seq(nrow(tmp_df))){
			# jj = 1
			tmp_genes = geneInfo$gdat$GENE[which(geneInfo$gdat$PATH == tmp_df$PATH[jj])]
			subj_genes = c(subj_genes,sample(tmp_genes,tmp_df$CNT[jj],replace = TRUE))
			rm(tmp_genes)
		}
		subj_genes = unique(subj_genes)
		sMUT[subj_genes,ii] = 1
		rm(tmp_df,subj_paths,subj_genes,samp_paths,subj_tab)
		
	}
	
	return(sMUT)
}
kOT_sim_GEN = function(geneInfo,NN,SCEN){
	
	nGENE = nrow(geneInfo$gene_sim)
	PATHs = sort(unique(geneInfo$gdat$PATH[geneInfo$gdat$PATH > 0])); PATHs
	subj_names = paste0("S",seq(NN))
	gene_names = paste0("G",seq(nGENE))
	path_names = paste0("P",PATHs)
	
	# Simulate pathways mutated
	sPATH = kOT_sim_sPATH(subj_names = subj_names,
		path_names = path_names,SCEN = SCEN)
	
	# Simulate genes mutated
	sMUT = kOT_sim_sMUT(sPATH = sPATH,
		geneInfo = geneInfo,SCEN = SCEN)
	
	# Output
	return(list(sMUT = sMUT,sPATH = sPATH))
	
}
kOT_sim_pXX = function(sPATH,pBETA,maxWAY){
	
	verbose = TRUE
	if( verbose ) cat(sprintf("%s: Construct pathway design matrix ...\n",date()))
	nPATH = ncol(sPATH)
	pXX = c()
	aa = 1
	
	for(ii in seq(maxWAY)){
		# ii = 2
		path_subsets = combn(x = seq(nPATH),m = ii)
		bb = aa + ncol(path_subsets) - 1; aa; bb
		idx_zero = which(pBETA[seq(aa,bb)] != 0)
		path_subsets = path_subsets[,idx_zero,drop = FALSE]
		# dim(path_subsets)
		# length(idx_zero)
	for(jj in seq(ncol(path_subsets))){
		# jj = 1
		path_group = path_subsets[,jj]; path_group
		tmp_mat = sPATH[,path_group,drop = FALSE]
		pXX = cbind(pXX,apply(tmp_mat,1,prod))
		colnames(pXX)[ncol(pXX)] = paste(colnames(tmp_mat),
			collapse = ".")
	}
		aa = bb + 1
		# print(dim(pXX)); pXX[1:4,]
	}
	
	dim(pXX); pXX[1:3,]
	
	return(pXX)
}
kOT_sim_hZZ = function(sPATH,pXX,TYPE,pBETA){
	
	NN = nrow(sPATH)
	nPATH = ncol(sPATH)
	path_names = colnames(sPATH)
	
	if( TYPE == "CAT" ){
		hZZ = c(pXX %*% pBETA[colnames(pXX)])
	} else {
		stop("No code yet")
	}
	
	return(hZZ)
	
}
kOT_sim_inputs = function(NN,nGENE,nPATH,path_corr,
	prop_noise,maxWAY,pBETA_sig){
	
	# Make COST
	geneInfo = kOT_sim_geneInfo(nGENE = nGENE,
		nPATH = nPATH,path_corr = path_corr,
		prop_noise = prop_noise)
	
	# Sim pBETA
	pBETA = c()
	fin_names = c()
	prop_nz = 0.5; # prop 25 pathways with non-zero 1-way effect
	prop_nz_inc = 0.5 # % of n-way effects are non-zero
	for(ii in seq(maxWAY)){ # maximum-way interactions
		# num_coef = choose(n = nPATH,k = ii)
		tmp_combs = combn(x = nPATH,m = ii)
		tmp_names = apply(tmp_combs,2,function(xx)
			paste(paste0("P",xx),collapse = "."))
		num_coef = ncol(tmp_combs)
		
		# Ensure prop_nz of 1-way, 2-way, etc interactions have non-zero effects
		tmp_pBETA = rnorm(n = num_coef,mean = 0,sd = pBETA_sig)
		binBETA = sample(c(0,1),num_coef,
			prob = c(1 - prop_nz,prop_nz),
			replace = TRUE)
		tmp_pBETA = tmp_pBETA * binBETA
		
		pBETA = c(pBETA,tmp_pBETA)
		fin_names = c(fin_names,tmp_names)
		prop_nz = prop_nz * prop_nz_inc
	}
	pBETA = round(pBETA,3)
	names(pBETA) = fin_names
	
	# Simulate baseline covariates (gender, tumor stage, age)
	subj_names = paste0("S",seq(NN))
	XX = kOT_sim_baseCov(subj_names = subj_names)
	xBETA = c(-5,0.5,seq(3)/3,-1)
	names(xBETA) = colnames(XX)
	
	# Output
	inputs = list(nGENE = nGENE,nPATH = nPATH,
		path_corr = path_corr,prop_noise = prop_noise,
		geneInfo = geneInfo,pBETA = pBETA,XX = XX,
		xBETA = xBETA)
	
	return(inputs)
	
}

#' @title kOT_sim_make
#' @description Generates simulation files
#' @param work_dir A full path to create "sim_ROKET" and subdirectories
#' @param NN A positive integer for sample size
#' @param nGENE A positive integer for number of genes to simulate
#' @param nPATH A positive integer for number of pathways to simulate
#' @param RR A positive integer for number of replicates to simulate
#' @export
kOT_sim_make = function(work_dir,NN = 200,
	nGENE = 500,nPATH = 25,RR = 200){
	
	my_dirs = setdirs(work_dir = work_dir)
	if( NN != round(NN) && NN > 0 )
		stop("NN isn't a positive integer")
	if( nGENE != round(nGENE) && nGENE > 0 )
		stop("nGENE isn't a positive integer")
	if( nPATH != round(nPATH) && nPATH > 0 )
		stop("nPATH isn't a positive integer")
	if( RR != round(RR) && RR > 0 )
		stop("RR isn't a positive integer")
	
	# Simulation parameters
	maxWAY 			= 4			# maximum number of n-way interactions
	path_corr 	= 0.50 	# lower bound on correlation among genes of same pathway
	prop_noise 	= 0.50 	# adds static noise to correlation between genes
	pBETA_sig		= 0.40 	# increase variance on pathway coefficients
	
	# Gene relationships, covariates
	inputs_fn = file.path(my_dirs$sim_dir,
		sprintf("inputs_nS.%s_nG.%s_nP.%s.rds",
		NN,nGENE,nPATH))
	if( !file.exists(inputs_fn) ){
		inputs = kOT_sim_inputs(NN = NN,nGENE = nGENE,
			nPATH = nPATH,path_corr = path_corr,
			prop_noise = prop_noise,maxWAY = maxWAY,
			pBETA_sig = pBETA_sig)
		saveRDS(inputs,inputs_fn)
	}
	inputs = readRDS(inputs_fn)
	
	# Get sim results for a SCENARIO and original replicate
	hBETAs = seq(0,1,0.1)
	prob_events = c(0.25,0.5,0.75,1)
	
	for(SCEN in seq(4)){
		# SCEN = 1
		cat(sprintf("%s: SCEN = %s ...\n",date(),SCEN))
		
		# Import scenario data
		scen_fn = file.path(my_dirs$reps_dir,
			sprintf("sim_SCEN%s.rds",SCEN))
		
		if( !file.exists(scen_fn) ){
			scen = kOT_sim_GEN(geneInfo = inputs$geneInfo,
				NN = NN,SCEN = SCEN)
			names(scen)
			saveRDS(scen,scen_fn)
		}
		
		scen = readRDS(scen_fn)
		pXX = kOT_sim_pXX(sPATH = scen$sPATH,
			pBETA = inputs$pBETA,maxWAY = maxWAY)
		
		for(rr in seq(RR)){
			smart_progress(ii = rr,nn = RR,iter = 5,iter2 = 1e2)
			out_fn = file.path(my_dirs$reps_dir,
				sprintf("sim_SCEN%s_rr%s.rds",SCEN,rr))
			if( file.exists(out_fn) ) next
			
			OUT = list()
			for(hBETA in hBETAs){
				# hBETA = hBETAs[1]
				
				# Simulate h(ZZ)
				if( hBETA == 0 ){
					hZZ = rep(0,NN)
				} else {
					hZZ = kOT_sim_hZZ(sPATH = scen$sPATH,pXX = pXX,
						TYPE = "CAT",pBETA = inputs$pBETA * hBETA)
				}
				
				# Expected mean
				regMU = c(inputs$XX %*% inputs$xBETA + hZZ )
				
				# Linear
				YY = regMU + rnorm(n = NN,0,sd = 1.5)
				
				# Survival
				RHO = 2.0
				UU = runif(n = NN,0,1)
				TT = ( -log(UU) / exp(regMU) )^(1/RHO)
				SD = smart_df(TT = TT)
				
				for(prob_event in prob_events){
					
					if( prob_event < 1 ){
						max_delta = max(TT); max_delta
						
						while(TRUE){
							# max_delta
							CC = runif(n = NN,0,max_delta)
							prop_delta = mean(TT <= CC); # print(prop_delta)
							prop_diff = prop_delta - prob_event; # prop_diff
							if( abs(prop_diff) < 0.02 ){
								break
							} else if( prop_diff > 0 ){
								max_delta = max_delta * 0.9
							} else if( prop_diff < 0 ){
								max_delta = max_delta * 1.1
							}
						}
						
					} else {
						CC = rep(max(TT) + 1,NN)
						
					}
					
					tmp_SD = smart_df(SS = pmin(TT,CC),DD = 1*(TT <= CC))
					prob_event2 = round(prob_event * 100)
					names(tmp_SD) = paste0(names(tmp_SD),"_",prob_event2)
					SD = cbind(SD,tmp_SD); rm(tmp_SD)
				
				}
				
				tmp_name = sprintf("hBETA_%s",hBETA)
				OUT[[tmp_name]] = smart_df(hZZ = hZZ,YY = YY,SD)
				rm(hZZ,regMU,YY,SD,UU,TT,CC)
				
			}
			
			saveRDS(OUT,out_fn)
			rm(OUT)
			
		}
		
		rm(scen)
	}
	
	return(NULL)
	
}


# ----------
# Simulation Analysis functions
# ----------

#' @title kOT_sim_OT
#' @param work_dir A full path to create "sim_ROKET" and subdirectories
#' @param NN A positive integer for sample size
#' @param nGENE A positive integer for number of genes to simulate
#' @param nPATH A positive integer for number of pathways to simulate
#' @param SCEN An integer taking values 1, 2, 3, or 4
#' @param ncores A positive integer specifying the number of
#'	cores/threads to use for optimal transport calculations
#' @export
kOT_sim_OT = function(work_dir,NN,nGENE,nPATH,SCEN,ncores = 1){
	
	my_dirs = setdirs(work_dir = work_dir)
	
	# Import geneInfo
	inputs_fn = file.path(my_dirs$sim_dir,
		sprintf("inputs_nS.%s_nG.%s_nP.%s.rds",
		NN,nGENE,nPATH))
	inputs = readRDS(inputs_fn)
	
	# Import scenario's genomics
	scen_fn = file.path(my_dirs$reps_dir,
		sprintf("sim_SCEN%s.rds",SCEN))
	scen = readRDS(scen_fn)
	
	# OT params
	EPS 				= 1e-3
	LAMs 				= c(0.5,1,5)
	BALs 				= c(TRUE,FALSE)
	conv 				= 1e-5
	max_iter 		= 3e3
	show_iter		= 10
	
	# Run NN by NN OT
	OT_fn = file.path(my_dirs$reps_dir,
		sprintf("OT_SCEN%s.rds",SCEN))
	
	if( !file.exists(OT_fn) ){
		OT = list()
		for(LAM in LAMs){
		for(BAL in BALs){
			# LAM = LAMs[1]; BAL = BALs[1]
			
			if( BAL == TRUE ){
				if( LAM != LAMs[1] ) next # aka we only need to run BAL = true once
			}
			cat(sprintf("%s: BAL = %s, LAM = %s ...\n",date(),BAL,LAM))
			
			tmp_name = sprintf("BAL.%s_LAM.%s",BAL,LAM)
			tmp_OT_fn = file.path(my_dirs$reps_dir,
				sprintf("tmp_OT_SCEN%s_%s.rds",SCEN,tmp_name))
			
			if( !file.exists(tmp_OT_fn) ){
			
				COST = 1 - inputs$geneInfo$gene_sim
				sMUT = scen$sMUT; sMUT = sMUT[rowSums(sMUT) >= 1,]
				int_genes = intersect(rownames(COST),rownames(sMUT))
				COST = COST[int_genes,int_genes]
				sMUT = sMUT[int_genes,]
				dim(COST); dim(sMUT)
				
				if( !all(colnames(COST) == rownames(sMUT)) ) stop("rownames mismatch")
				out = Rcpp_run_full_OT(COST = COST,ZZ = sMUT,EPS = EPS,
					LAMBDA1 = LAM,LAMBDA2 = LAM,balance = BAL,
					highLAM_lowMU = TRUE,conv = conv,max_iter = max_iter,
					ncores = ncores,show = FALSE,show_iter = show_iter)
				
				NN = ncol(sMUT)
				out$DIST = smart_names(out$DIST,
					ROW = paste0("S",seq(NN)),
					COL = paste0("S",seq(NN)))
				out$DIST_2 = out$DIST / out$sum_OT
				names(out)
				saveRDS(out,tmp_OT_fn)
				
			} else {
				out = readRDS(tmp_OT_fn)
				
			}
			
			OT[[tmp_name]] = out
			
		}}
		saveRDS(OT,OT_fn)
	}
	
	tmp_fns = list.files(my_dirs$reps_dir)
	tmp_fns = tmp_fns[grepl(sprintf("tmp_OT_SCEN%s_",SCEN),tmp_fns)]
	if( length(tmp_fns) > 0 )
		unlink(file.path(my_dirs$reps_dir,tmp_fns))
	
	return(NULL)
	
}

kOT_sim_GENE = function(sim,out = "OLS",hBETAs = NULL,nPERM,samp_thres){
	
	num_basecov = ncol(sim$XX); num_basecov
	vec_genes = rownames(sim$sGEN$sMUT)
	num_genes = nrow(sim$sGEN$sMUT); num_genes
	
	if( is.null(hBETAs) ){
		hBETAs = names(sim$OUT); hBETAs
		hBETAs = hBETAs[as.numeric(gsub("hBETA_","",hBETAs)) >= 0]; hBETAs
	}
	
	OUTs = out
	if( out == "SURV" ){
		OUTs = names(sim$OUT$hBETA_0)
		OUTs = OUTs[grepl("SS_",OUTs)]
		OUTs = gsub("SS_","",OUTs)
		OUTs
	}
	# GMs = c(0,5,10,20,30)
	log10_TMB = log10(colSums(sim$sGEN$sMUT))
	XX = smart_df(cbind(sim$XX[,-1],log10_TMB))
	NN = nrow(XX)
	
	res = c()
	nom_PVAL = matrix(NA,num_genes,length(OUTs) *length(hBETAs))
	rownames(nom_PVAL) = vec_genes
	colnames(nom_PVAL) = sprintf("C%s",seq(ncol(nom_PVAL)))
	iter = 1
	
	for(OUT in OUTs){
	for(hBETA in hBETAs){
		# OUT = OUTs[1]; hBETA = hBETAs[1]
		cat(sprintf("%s: OUT = %s; hBETA = %s ...\n",
			date(),OUT,hBETA))
		
		colnames(nom_PVAL)[iter] = sprintf("%s_%s",OUT,hBETA)
		
		if( out == "SURV" ){
			tmp_surv = sim$OUT[[hBETA]]
			SS = tmp_surv[,paste0("SS_",OUT)]
			DD = tmp_surv[,paste0("DD_",OUT)]
		}
		
		aPVAL = matrix(NA,nPERM + 1,num_genes)
		colnames(aPVAL) = vec_genes
		idx_SUBJ = matrix(NA,nPERM + 1,NN)
		idx_SUBJ[1,] = seq(NN)
		idx_SUBJ[-1,] = t(sapply(seq(nPERM),function(xx) sample(NN)))
		cnt = 1
		
	for(gene in vec_genes){
		# gene = vec_genes[79]; gene
		smart_progress(ii = cnt,nn = num_genes,iter = 5,iter2 = 1e2)
		tab1 = table(sim$sGEN$sMUT[gene,]); tab1
		if( length(tab1) != 2 || min(tab1) < samp_thres ){
			aPVAL[,gene] = 1
			cnt = cnt + 1
			next
		}
		
	for(PERM in seq(nrow(aPVAL))){
		# PERM = 1
		if( out == "OLS" ){
			lm_out = lm(sim$OUT[[hBETA]]$YY ~ 
				sim$sGEN$sMUT[gene,idx_SUBJ[PERM,]] + .,data = XX)
			aPVAL[PERM,gene] = summary(lm_out)$coefficients[2,4]
		}
		if( out == "SURV" ){
			
			cx_out = tryCatch(coxph(Surv(SS,DD) ~ 
				sim$sGEN$sMUT[gene,idx_SUBJ[PERM,]] + .,data = XX),
				warning = function(ww){NULL},error = function(ee){NULL})
			
			if( is.null(cx_out) ){
				aPVAL[PERM,gene] = 1
				next
			}
			
			aPVAL[PERM,gene] = summary(cx_out)$coefficients[1,5]
			
		}
		
	}
		cnt = cnt + 1
	}
		
		keep_genes = colnames(aPVAL)[aPVAL[1,] != 1]; length(keep_genes)
		
		if( length(keep_genes) == 0 ){
			PVAL_perm = 1
			
		} else {
			vec_mPVAL = apply(aPVAL[,keep_genes,drop = FALSE],1,min)
			PVAL_perm = mean(vec_mPVAL <= vec_mPVAL[1]); PVAL_perm
			
		}
		
		hBETA_2 = as.numeric(gsub("hBETA_","",hBETA)); hBETA_2
		OUT2 = OUT
		if( out == "SURV" ) OUT2 = paste0("SURV_nCENS.",OUT2)
		OUT2
		
		res = rbind(res,smart_df(OUT = OUT2,hBETA = hBETA_2,
			PVAL_perm = PVAL_perm,nTEST_genes = length(keep_genes)))
		nom_PVAL[,iter] = aPVAL[1,]
		iter = iter + 1
	}}
	
	print(res)
	return(list(RES = res,nom_PVAL = nom_PVAL))
	
}
kOT_sim_KERN = function(sim,OT,nPERM,hBETAs = NULL){
	
	OT_EUC = OT
	OT_EUC$EUC = as.matrix(dist(t(sim$sGEN$sMUT)))
	all_out = c("OLS","SURV")
	all_LABS = c("EUC","OT")
	
	if( is.null(hBETAs) ){
		hBETAs = names(sim$OUT); hBETAs
	}
	
	# Set null model covariates
	log10_TMB = log10(colSums(sim$sGEN$sMUT))
	XX = cbind(sim$XX[,-1],log10_TMB)
	
	reg_out = c()
	for(out in all_out){
	for(hBETA in hBETAs){
	for(LAB in all_LABS){
		if(FALSE){
			out = all_out[2]; hBETA = hBETAs[1]; 
			LAB = all_LABS[2]
		}
		
		cat(sprintf("%s: OUT = %s; hBETA = %s; LAB = %s ...\n",
			date(),out,hBETA,LAB))
		hBETA_2 = as.numeric(gsub("hBETA_","",hBETA))
		
		# Make Kernel list
		KK = list()
		if( LAB == "EUC" ){
			dd = OT_EUC$EUC
			KK$EUC = D2K(D = dd)
		} else {
			for(vv in names(OT)){
				dd = OT[[vv]]$DIST
				vv2 = paste(strsplit(vv,"_")[[1]][1:2],collapse = "_")
				KK[[vv2]] = D2K(D = dd)
			}
		}
		rm(dd)
		
		# Run regression
		if( out == "OLS" ){
			fit = MiRKAT(y = sim$OUT[[hBETA]]$YY,X = XX,
				Ks = KK,out_type = "C",method = "permutation",
				omnibus = "permutation",nperm = nPERM,
				returnKRV = FALSE,returnR2 = FALSE)
			fit
			
			if( LAB == "EUC" ){
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "MiRKAT",LAB = LAB,PVAL_perm = as.numeric(fit$p_values)))
			} else {
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "MiRKAT",LAB = "omni_OT",PVAL_perm = fit$omnibus_p))
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "MiRKAT",LAB = names(fit$p_values),
					PVAL_perm = as.numeric(fit$p_values)))
			}
			
			lm_out = lm(sim$OUT[[hBETA]]$YY ~ .,data = smart_df(XX))
			RESI = as.numeric(lm_out$residuals)
			fit2 = Rcpp_KernTest(RESI = RESI,KK = KK,nPERMS = nPERM,
				iter1 = 50,iter2 = 1000,verbose = FALSE); rm(RESI)
			names(fit2$PVALs) = names(KK)
			fit2
			
			if( LAB == "EUC" ){
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "Rcpp",LAB = LAB,PVAL_perm = as.numeric(fit2$PVALs)))
			} else {
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "Rcpp",LAB = "omni_OT",PVAL_perm = fit2$omni_PVAL))
				reg_out = rbind(reg_out,smart_df(OUT = out,hBETA = hBETA_2,
					SOFT = "Rcpp",LAB = names(fit2$PVALs),
					PVAL_perm = as.numeric(fit2$PVALs)))
			}
			
		}
		if( out == "SURV" ){
			OUTs = names(sim$OUT[[hBETA]])
			OUTs = OUTs[grepl("SS_",OUTs)]
			OUTs = gsub("SS_","",OUTs)
			OUTs
			
			for(OUT in OUTs){
				# OUT = OUTs[1]
				
				OUT2 = paste0("SURV_nCENS.",OUT)
				tmp_surv = sim$OUT[[hBETA]]
				SS = tmp_surv[,paste0("SS_",OUT)]
				DD = tmp_surv[,paste0("DD_",OUT)]
				
				fit = suppressWarnings(MiRKATS(obstime = SS,delta = DD,X = XX,
					Ks = KK,perm = TRUE,omnibus = "permutation",nperm = nPERM))
				fit
				
				if( LAB == "EUC" ){
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "MiRKAT",LAB = LAB,PVAL_perm = as.numeric(fit$p_values)))
				} else {
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "MiRKAT",LAB = "omni_OT",PVAL_perm = fit$omnibus_p))
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "MiRKAT",LAB = names(fit$p_values),
						PVAL_perm = as.numeric(fit$p_values)))
				}
			
				cout = coxph(Surv(SS,DD) ~ .,data = smart_df(XX))
				RESI = as.numeric(cout$residuals)
				fit2 = Rcpp_KernTest(RESI = RESI,KK = KK,nPERMS = nPERM,
					iter1 = 50,iter2 = 1000,verbose = FALSE); rm(RESI)
				names(fit2$PVALs) = names(KK)
				fit2
				
				if( LAB == "EUC" ){
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "Rcpp",LAB = LAB,PVAL_perm = as.numeric(fit2$PVALs)))
				} else {
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "Rcpp",LAB = "omni_OT",PVAL_perm = fit2$omni_PVAL))
					reg_out = rbind(reg_out,smart_df(OUT = OUT2,hBETA = hBETA_2,
						SOFT = "Rcpp",LAB = names(fit2$PVALs),
						PVAL_perm = as.numeric(fit2$PVALs)))
				}
			
			}
			
		}
		
	}}}
	print(reg_out)
	
	return(reg_out)
	
}

#' @title kOT_sim_REG
#' @param work_dir A full path to create "sim_ROKET" and subdirectories
#' @param NN A positive integer for sample size
#' @param nGENE A positive integer for number of genes to simulate
#' @param nPATH A positive integer for number of pathways to simulate
#' @param SCEN An integer taking values 1, 2, 3, or 4
#' @param rr A positive integer indexing a replicate
#' @export
kOT_sim_REG = function(work_dir,NN,nGENE,nPATH,SCEN,rr){
	
	my_dirs = setdirs(work_dir = work_dir)
	
	inputs_fn = file.path(my_dirs$sim_dir,
		sprintf("inputs_nS.%s_nG.%s_nP.%s.rds",
		NN,nGENE,nPATH))
	inputs = readRDS(inputs_fn)
	
	rds_fn = file.path(my_dirs$rout_dir,
		sprintf("sim_rr_%s_SCEN%s.rds",rr,SCEN))
	if( file.exists(rds_fn) ) unlink(rds_fn)
	
	# Import mutations
	sim_fn = file.path(my_dirs$reps_dir,
		sprintf("sim_SCEN%s.rds",SCEN))
	sim = readRDS(sim_fn)
	sim$XX = inputs$XX
	
	# Import OT
	OT_fn = file.path(my_dirs$reps_dir,
		sprintf("OT_SCEN%s.rds",SCEN))
	if( !file.exists(OT_fn) ) stop(sprintf("OT %s not done!",SCEN))
	OT = readRDS(OT_fn)
	
	# Import outcomes
	out_fn = file.path(my_dirs$reps_dir,
		sprintf("sim_SCEN%s_rr%s.rds",SCEN,rr))
	sim$OUT = readRDS(out_fn)
	
	# Run gene regressions
	reg_gene = c()
	OUTs = c("OLS","SURV")
	nom_PVAL = list()
	for(out in OUTs){
		# out = OUTs[2]
		tmp_out = kOT_sim_GENE(sim = sim,
			out = out,hBETAs = NULL,nPERM = 1e2,
			samp_thres = 20)
		reg_gene = rbind(reg_gene,tmp_out$RES)
		nom_PVAL[[out]] = tmp_out$nom_PVAL
	}
	
	# Run kernel regressions
	kern_rds_fn = file.path(my_dirs$rout_dir,
		sprintf("sim_rr_%s_SCEN%s_kern.rds",rr,SCEN))
	if( !file.exists(kern_rds_fn) ){
		reg_kern = kOT_sim_KERN(sim = sim,OT = OT,
			nPERM = 2e3,hBETAs = NULL)
		saveRDS(reg_kern,kern_rds_fn)
	}
	reg_kern = readRDS(kern_rds_fn)
	
	saveRDS(list(GENE = reg_gene,
		nom_PVAL = nom_PVAL,KERN = reg_kern),rds_fn)
	return(NULL)
	
}

kOT_sim_AGG = function(work_dir){
	
	my_dirs = setdirs(work_dir = work_dir)
	res_fn 	= file.path(my_dirs$sout_dir,"res.rds")
	ures_fn = file.path(my_dirs$sout_dir,"ures.rds")
	
	# Aggregate replicates
	if( !file.exists(ures_fn) ){
		RR = list.files(my_dirs$rout_dir,pattern = "sim")
		RR = RR[grepl("_rr_",RR)]
		RR = RR[grepl(".rds$",RR)]
		RR = RR[!grepl("kern",RR)]
		RR = RR[grepl("SCEN1",RR)]
		RR = length(RR); RR
		
		res = c()
		fin_nom_PVAL = list()
		for(SCEN in seq(4)){
			# SCEN = 1
			cat(sprintf("%s: SCEN = %s ...\n",date(),SCEN))
		for(rr in seq(RR)){
			# rr = 1
			smart_progress(ii = rr,nn = RR,iter2 = 1e2)
			rr_fn = file.path(my_dirs$rout_dir,
				sprintf("sim_rr_%s_SCEN%s.rds",rr,SCEN))
			if( !file.exists(rr_fn) ){
				cat(sprintf("%s miss ",rr))
				next
			}
			tmp_rds = readRDS(rr_fn)
			
			if( rr == 1 ){
				nom_PVAL = tmp_rds$nom_PVAL
				nom_PVAL$OLS = ifelse(nom_PVAL$OLS < 0.05,1,0)
				nom_PVAL$SURV = ifelse(nom_PVAL$SURV < 0.05,1,0)
				
			} else {
				nom_PVAL$OLS = nom_PVAL$OLS + ifelse(tmp_rds$nom_PVAL$OLS < 0.05,1,0)
				nom_PVAL$SURV = nom_PVAL$SURV + ifelse(tmp_rds$nom_PVAL$SURV < 0.05,1,0)
				
			}
			
			tmp_gene = tmp_rds$GENE
			tmp_gene$RR = rr
			tmp_gene$SCEN = SCEN
			res_gene = c()
			for(VAR in c("PVAL_perm")){
				tmp_df = tmp_gene
				tmp_df$THRES = tmp_df[,VAR]
				VAR2 = VAR
				if( VAR == "PVAL_perm" ) 	VAR2 = "minP_gene_perm1"
				tmp_df$DIST = VAR2
				tmp_df$SOFT = "none"
				tmp_df = tmp_df[,c("RR","OUT","DIST","hBETA",
					"SOFT","SCEN","THRES","nTEST_genes")]
				res_gene = rbind(res_gene,tmp_df); rm(tmp_df)
			}
			
			tmp_kern = tmp_rds$KERN[which(tmp_rds$KERN$hBETA >= 0),]
			tmp_kern$RR = rr
			tmp_kern$SCEN = SCEN
			tmp_kern$THRES = tmp_kern$PVAL_perm
			tmp_kern$DIST = tmp_kern$LAB
			tmp_kern$nTEST_genes = 500
			tmp_kern = tmp_kern[,c("RR","OUT","DIST","hBETA",
				"SOFT","SCEN","THRES","nTEST_genes")]
			tmp_kern[1:5,]
			
			res = rbind(res,res_gene,tmp_kern)
			rm(tmp_rds,tmp_gene,tmp_kern)
			
		}
			nom_PVAL$OLS = nom_PVAL$OLS / RR
			nom_PVAL$SURV = nom_PVAL$SURV / RR
			fin_nom_PVAL[[sprintf("SCEN%s",SCEN)]] = list(nom_PVAL = nom_PVAL)
			rm(nom_PVAL,mt_PVAL)
		}
		saveRDS(list(RES = res,fin_nom_PVAL = fin_nom_PVAL),res_fn)
		
		# Calculate hypothesis rejection
		rds = readRDS(res_fn)
		res = rds$RES
		ures = unique(res[,c("OUT","DIST","SOFT","SCEN","hBETA")])
		ures$REJECT = NA
		dim(ures); ures[1:5,]
		tot = nrow(ures); tot
		cat(sprintf("%s: Calculate %s rejections ...\n",date(),tot))
		for(ii in seq(tot)){
			# ii = 1
			smart_progress(ii = ii,nn = tot,iter2 = 1e2)
			idx = which(res$OUT == ures$OUT[ii]
				& res$DIST == ures$DIST[ii]
				& res$SOFT == ures$SOFT[ii]
				& res$SCEN == ures$SCEN[ii]
				& res$hBETA == ures$hBETA[ii])
			length(idx)
			
			# ures$REJECT[ii] = calc_POWER_2(PVAL = res$THRES[idx])
			ures$REJECT[ii] = mean(res$THRES[idx] < 0.05)
			
			rm(idx)
		}
		
		mDIST = c("minP_gene_perm1")
		tmp_df = ures[which(ures$DIST %in% mDIST),]
		tmp_df$SOFT = "Rcpp"
		ures$SOFT[ures$DIST %in% mDIST] = "MiRKAT"
		ures = rbind(ures,tmp_df); rm(tmp_df)
		ures$SOFT[grepl("Rcpp",ures$SOFT)] = "ROKET"
		saveRDS(ures,ures_fn)
	}
	ures = readRDS(ures_fn)
	
	if( TRUE ){ # Polish
		ures[1:3,]
		smart_table(ures$DIST)
		tmp_lev = sort(unique(ures$DIST)); tmp_lev
		tmp_lev = tmp_lev[c(6,7,5,4,3,2,1,8)]; tmp_lev
		# tmp_lev = tmp_lev[c(6,7,8,9,10,11,12,13,5,4,3,2,1,14)]; tmp_lev
		# tmp_lev = tmp_lev[c(6,7,8,9,10,5,4,3,2,1,11)]; tmp_lev
		ures$DIST2 = factor(ures$DIST,levels = tmp_lev,
			labels = c("Gene-Based (Bonf)","Gene-Based (Perm)",
			# sprintf("Path-Based (Bonf%s)",c(1,2)),
			# sprintf("Path-Based (Perm%s)",c(1,2)),
			"Euclidean",sprintf("OT (\u03BB = %s)",c("\u221E","5.0","1.0","0.5")),
			"OT (omnibus)"))
		
		ures$SCEN2 = factor(ures$SCEN,
			levels = sort(unique(ures$SCEN)),
			labels = paste0("Scenario ",sort(unique(ures$SCEN))))
	}
	# dim(ures); ures[1:3,]
	
	my_theme = theme(legend.position = c("none","bottom","right")[2],
		text = element_text(size = 40),
		axis.text.x = element_text(size = 25),
		# axis.text.y = element_text(size = 25),
		# strip.text.y = element_text(size = 18),
		panel.background = element_blank(),
		panel.grid.major = element_line(colour = "grey50",
			size = 1,linetype = "dotted"),
		panel.spacing.x = unit(2,"lines"),
		panel.spacing.y = unit(2,"lines"),
		legend.key.width = unit(2,"line"),
		# legend.title = element_text(size = 34),
		# legend.text = element_text(size = 34),
		plot.title = element_text(hjust = 0.5))##,size = 38))
	
	# New plot for paper
	hBETA2 = DIST2 = NULL
	if( TRUE ){
		ures2 = ures[which(ures$SOFT == "ROKET" & ures$hBETA >= 0 
			& ures$OUT %in% c("OLS",sprintf("SURV_nCENS.%s",c(25,50,75,100)))
			),]
		
		# dim(ures2); ures2[1:5,]
		lev_OUT = sort(unique(ures2$OUT)); lev_OUT
		all_OUT = c("OLS",sprintf("SURV_nCENS.%s",c(100,75,50,25))); all_OUT
		lev_OUT = all_OUT[which(all_OUT %in% lev_OUT)]; lev_OUT
		ures2$OUT = factor(ures2$OUT,levels = lev_OUT,
			labels = c("OLS",sprintf("CENS = %s%%",c(0,25,50,75))))
		# table(ures2$OUT)
		# ures2$DIST
		
		my_theme = theme(legend.position = c("none","bottom","right")[2],
			text = element_text(size = 45),
			axis.text.x = element_text(size = 30),
			# axis.text.y = element_text(size = 42),
			# strip.text.y = element_text(size = 18),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				size = 1,linetype = "dotted"),
			panel.spacing.x = unit(4,"lines"),
			panel.spacing.y = unit(2,"lines"),
			legend.key.width = unit(4,"line"),
			# legend.title = element_text(size = 34),
			# legend.text = element_text(size = 34),
			plot.title = element_text(hjust = 0.5))##,size = 38))
		
		gg = ggplot(data = ures2,
			mapping = aes(x = hBETA,y = REJECT,color = DIST2)) +
			facet_grid(SCEN2 ~ OUT) + geom_line(size = 3) + 
			geom_hline(aes(yintercept = 0.05),
				size = 2,linetype = 2,color = "black") +
			ylim(c(0,1)) +
			xlab(expression(gamma)) + ylab("Power / Type I Error") +
			labs(color = "Approach") + my_theme
		# gg
		
		png_fn = file.path(my_dirs$sout_dir,"final_sim.png")
		ggsave(png_fn,plot = gg,device = "png",
			width = 25,height = 25,units = "in",dpi = 75)
		rm(gg)
		
	}
	
	if( TRUE ){ # Check distribution of gene power
		rds = readRDS(res_fn)
		names(rds)
		nPVAL = rds$fin_nom_PVAL
		names(nPVAL)
		
		wPVAL = c("nom_PVAL")[1]
		out = c()
		for(SCEN in seq(4)){
			# SCEN = 1
			SCEN2 = sprintf("SCEN%s",SCEN)
		for(OUT in names(nPVAL[[SCEN2]][[wPVAL]])){
			# OUT = c("OLS","SURV")[2]
			cnames = colnames(nPVAL[[SCEN2]][[wPVAL]][[OUT]])
		for(VAR in cnames){
			if( OUT == "OLS" ){
				hBETA = as.numeric(gsub("OLS_hBETA_","",VAR))
				OUT2 = OUT
			} else if( OUT == "SURV" ){
				nCENS = as.integer(strsplit(VAR,"_")[[1]][1])
				OUT2 = sprintf("SURV_nCENS%s",nCENS)
				hBETA = as.numeric(strsplit(VAR,"_")[[1]][3])
			}
			
			REJECT = nPVAL[[SCEN2]][[wPVAL]][[OUT]][,VAR]
			REJECT = as.numeric(REJECT[REJECT > 0])
			out = rbind(out,smart_df(SCEN = SCEN,
				OUT = OUT2,hBETA = hBETA,REJECT = REJECT))
			rm(REJECT)
		}}}
		
		out$SCEN2 = sprintf("Scenario = %s",out$SCEN)
		tmp_lev = sort(unique(out$OUT)); tmp_lev
		tmp_lev = tmp_lev[c(1,2,5,4,3)]; tmp_lev
		out$OUT2 = factor(out$OUT,levels = tmp_lev,
			labels = c("OLS",sprintf("CENS = %s%%",c(0,25,50,75))))
		out$hBETA2 = factor(out$hBETA,levels = sort(unique(out$hBETA)),
			labels = smart_digits(sort(unique(out$hBETA)),digits = 1))
		dim(out); out[1:20,]
		
		my_theme = theme(legend.position = c("none","bottom","right")[1],
			text = element_text(size = 45),
			axis.text.x = element_text(size = 30),
			# axis.text.y = element_text(size = 42),
			# strip.text.y = element_text(size = 18),
			panel.background = element_blank(),
			panel.grid.major = element_line(colour = "grey50",
				size = 1,linetype = "dotted"),
			panel.spacing.x = unit(4,"lines"),
			panel.spacing.y = unit(2,"lines"),
			legend.key.width = unit(2,"line"),
			# legend.title = element_text(size = 34),
			# legend.text = element_text(size = 34),
			plot.title = element_text(hjust = 0.5))
		
		gg = ggplot(data = out,mapping = aes(x = hBETA2,y = REJECT)) +
			geom_boxplot(fill = "#56B4E9") + facet_grid(SCEN2 ~ OUT2) +
			geom_hline(aes(yintercept = 0.05),size = 2,linetype = 2,color = "black") +
			xlab(expression(gamma)) + ylab("Power / Type I Error") +
			my_theme
		png_fn = file.path(my_dirs$sout_dir,
			sprintf("final_sim_%s.png",wPVAL))
		ggsave(png_fn,plot = gg,device = "png",
			width = 45,height = 25,units = "in",dpi = 75)
		rm(gg)
	}
	
}


#' @importFrom stats dist lm rnorm runif
#' @importFrom smartr smart_df smart_progress 
#'	smart_mkdir smart_table smart_digits smart_names
#'	make_dummy
#' @importFrom survival Surv coxph
#' @importFrom utils combn
#' @importFrom MiRKAT D2K MiRKAT MiRKATS
#' @importFrom ggplot2 ggsave theme element_text
#'	unit ggplot aes facet_grid geom_line geom_hline
#'	ylim xlab ylab labs geom_boxplot element_blank
#'	element_line
#' @importFrom Rcpp sourceCpp
#' @useDynLib ROKET
NULL

# Steps to create/check/install package from directory
# bb = strsplit(getwd(),"/")[[1]]; pack_dir = paste(bb[-length(bb)],collapse = "/")
# pack = strsplit(pack_dir,"/")[[1]]; pack = pack[length(pack)]
# if( pack %in% installed.packages()[,1] ) remove.packages(pack)
# Rcpp::compileAttributes(pkgdir = pack_dir)
# devtools::document(pkg = pack_dir); usethis::use_gpl3_license()
# devtools::check(pkg = pack_dir,manual = TRUE,cran = TRUE,error_on = "note")
# devtools::install(pack_dir)


###

