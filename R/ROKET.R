# Package Functions

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
#'	optimal transport regardless of LAMDA1 and LAMBDA2. 
#'	Otherwise run unbalanced optimal transport.
#' @param conv A positive numeric value to determine 
#'	algorithmic convergence. The default value is \code{1e-5}.
#' @param max_iter A positive integer denoting the maximum
#'	iterations to run the algorithm.
#' @param verbose Boolean value to display verbose function output.
#' @param show_iter A positive integer to display iteration details
#'	at multiples of \code{show_iter} but only if \code{verbose = TRUE}.
#' 
#' @export
run_myOT = function(XX,YY,COST,EPS,LAMBDA1,LAMBDA2,
	balance = FALSE,conv = 1e-5,max_iter = 3e3,
	verbose = TRUE,show_iter = 50){
	
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
		conv = conv,max_iter = max_iter,show = verbose,
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
#' @inheritParams run_myOT
#' @param ZZ A numeric matrix of non-negative mass to transport.
#'	Rows correspond to features (e.g. genes) and columns
#'	correspond to samples or individuals. Each column must have
#'	strictly positive mass
#' @param COST A numeric square matrix of non-negative values
#'	representing the non-negative costs to transport 
#'	masses between pairs of features
#' @param ncores A positive integer for the number of cores/threads
#'	to reduce computational runtime when running for loops
#' @param show_iter A positive integer to display iteration details
#'	at multiples of \code{show_iter} but only if \code{verbose = TRUE}.
#' @export
run_myOTs = function(ZZ,COST,EPS,LAMBDA1,LAMBDA2,
	balance,conv = 1e-5,max_iter = 3e3,ncores = 1,
	verbose = TRUE,show_iter = 50){
	
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
	
	if( !all(colnames(COST) == rownames(ZZ)) )
		stop("rownames mismatch")
	
	# Run OT across all individuals
	out_OT = Rcpp_run_full_OT(COST = COST,ZZ = ZZ,
		EPS = EPS,LAMBDA1 = LAMBDA1,LAMBDA2 = LAMBDA2,
		balance = balance,highLAM_lowMU = TRUE,
		conv = conv,max_iter = max_iter,ncores = ncores,
		show = verbose,show_iter = show_iter)
	
	out_OT$DIST = smart_names(out_OT$DIST,
		ROW = colnames(ZZ),COL = colnames(ZZ))
	out_OT$sum_OT = smart_names(out_OT$sum_OT,
		ROW = colnames(ZZ),COL = colnames(ZZ))
	out_OT$DIST_2 = out_OT$DIST / out_OT$sum_OT
	
	return(out_OT)
	
}


# ----------
# Kernel Regression Hypothesis Testing
# ----------

#' @title kernTEST
#' @inheritParams run_myOTs
#' @param RESI A numeric vector of null model residuals
#'	\code{names(RESI)} must be set to maintain sample ordering
#' @param KK An array containing double-centered positive semi-definite
#'	kernel matrices. Refer to \code{MiRKAT::D2K()} for transforming 
#'	distance matrices to kernel matrices. The \code{dimnames(KK)[[1]]} and 
#'	\code{dimnames(KK)[[2]]} must match \code{names(RESI)}.
#'	Also set dimnames(KK)[[3]] to keep track of each kernel matrix.
#' @param OMNI A matrix of zeros and ones. Each column corresponds to a
#'	distance matrix while each row corresponds to an omnibus test. Set
#'	\code{rownames(OMNI)} for labeling outputted p-values and 
#'	\code{colnames(OMNI)} which should match \code{dimnames(KK)[[3]]}.
#' @param nPERMS A positive integer to specify the number of
#'	permutation-based p-value calculation
#' @export
kernTEST = function(RESI,KK,OMNI,nPERMS = 1e5,ncores = 1){
	
	samp_names = names(RESI)
	
	if( is.null(samp_names) )
		stop("Specify names(RESI)")
	if( is.null(dimnames(KK)[[3]]) || is.null(dimnames(KK)[[1]]) 
		|| is.null(dimnames(KK)[[2]]) )
		stop("Specify dimnames(KK) to track sample order and kernel matrices")
	
	if( !all(samp_names == dimnames(KK[[1]])) )
		stop("names(RESI) != dimnames(KK[[1]]")
	if( !all(samp_names == dimnames(KK[[2]])) )
		stop("names(RESI) != dimnames(KK[[2]]")
	
	# Check OMNI object
	if( !all(colnames(OMNI) == dimnames(KK)[[3]]) )
		stop("colnames(OMNI) != dimnames(KK[[3]]")
	if( is.null(rownames(OMNI)) )
		stop("Specify rownames(OMNI)")
	if( !all(c(OMNI) %in% c(0,1)) )
		stop("OMNI should take values 0 or 1")
	if( any(rowSums(OMNI) == 0) )
		stop("Each row of OMNI should contain at least one non-zero element")
	
	out_test = Rcpp_KernTest(RESI = RESI,cKK = KK,
		OMNI = OMNI,nPERMS = nPERMS,ncores = ncores)
	names(out_test$PVALs) = dimnames(KK)[[3]]
	names(out_test$omni_PVALs) = rownames(OMNI)
	
	return(out_test)
	
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
# if( pack %in% installed.packages()[,1] ){ remove.packages(pack); q("no")}
# Rcpp::compileAttributes(pkgdir = pack_dir)
# devtools::document(pkg = pack_dir); usethis::use_gpl3_license()
# Sys.setenv("RSTUDIO_PANDOC" = "C:/Program Files/RStudio/bin/pandoc")
# check_pandoc = rmarkdown::pandoc_available(); check_pandoc
#### usethis::use_vignette(name = "test",title = "Testing")
# make_vign = check_pandoc && !TRUE; make_vign
# devtools::check(pkg = pack_dir,manual = TRUE,cran = FALSE,error_on = c("warning","note")[1],vignettes = make_vign)
# devtools::install(pack_dir,build_vignettes = make_vign)


###

