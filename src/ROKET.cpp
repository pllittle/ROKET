// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// --------------------
// Minor functions
// --------------------
double Rcpp_logSumExp(const arma::vec& log_x){
	if( log_x.n_elem == 1 ){
		return log_x.at(0);
	} else {
		double max_val = arma::max(log_x);
		arma::vec log_x_2 = log_x - max_val;
		return std::log(arma::sum(arma::exp(log_x_2))) + max_val;
	}
}

arma::vec Rcpp_calc_rank(const arma::vec& aa){
	arma::vec out = arma::zeros<arma::vec>(aa.n_elem);
	out = arma::conv_to<arma::vec>::from(arma::sort_index(aa));
	out = arma::conv_to<arma::vec>::from(arma::sort_index(out)) + 1;
	
	// Check ties
	arma::vec uniq_aa = arma::unique(aa);
	if( uniq_aa.n_elem != aa.n_elem ){
		arma::uword ii;
		double tmp_avg;
		for(ii = 0; ii < uniq_aa.n_elem; ii++){
			arma::uvec idx = arma::find(aa == uniq_aa.at(ii));
			if( idx.n_elem == 1 ) continue;
			tmp_avg = arma::sum(out.elem(idx)) / idx.n_elem;
			out.elem(idx).fill(tmp_avg);
		}
	}
	
	return out;
	
}


// --------------------
// Optimal Transport Functions
// --------------------
arma::mat run_OT_OPT(const arma::vec& XX,const arma::vec& YY,
	const arma::mat& COST_XY,const double& EPS,const double& LAMBDA1,
	const double& LAMBDA2,const double& conv,const arma::uword& max_iter,
	const bool& show,const arma::uword& show_iter){
	
	double cnst_u = 1.0,cnst_v = 1.0,diff_uu = 0.0,diff_vv = 0.0,
		max_diff = 0.0;
	if( LAMBDA1 > 0.0 ) cnst_u = LAMBDA1 / (LAMBDA1 + EPS);
	if( LAMBDA2 > 0.0 ) cnst_v = LAMBDA2 / (LAMBDA2 + EPS);
	arma::uword ii,jj,kk,ll,n1 = XX.n_elem,n2 = YY.n_elem,iter = 0;
	arma::mat log_KK = -1.0 / EPS * COST_XY;
	arma::vec log_XX = arma::log(XX),log_YY = arma::log(YY),
		log_uu = arma::zeros<arma::vec>(n1),log_uu2 = log_uu,
		log_vv = arma::zeros<arma::vec>(n2),log_vv2 = log_vv,
		log_tt = log_uu,log_ss = log_vv;
	
	while(iter < max_iter){
		
		// Update log_tt = log( KK * vv )
		for(kk = 0; kk < n1; kk++){
			log_tt.at(kk) = Rcpp_logSumExp(log_KK.row(kk).t() + log_vv);
		}
		
		// Update log_uu = cnst_u * ( log(XX) - log(tt) )
		log_uu2 = cnst_u * ( log_XX - log_tt );
		if( iter > 0 ) diff_uu = arma::max(arma::abs(log_uu - log_uu2));
		log_uu = log_uu2;
		
		// Update log_ss = log( KK.t() * uu )
		for(ll = 0; ll < n2; ll++){
			log_ss.at(ll) = Rcpp_logSumExp(log_KK.col(ll) + log_uu);
		}
		
		// Update log_vv = cnst_v * ( log(YY) - log(ss) )
		log_vv2 = cnst_v * ( log_YY - log_ss );
		if( iter > 0 ) diff_vv = arma::max(arma::abs(log_vv - log_vv2));
		log_vv = log_vv2;
		
		// Check convergence
		if( iter > 0 ){
			max_diff = diff_uu;
			if( diff_vv > diff_uu ) max_diff = diff_vv;
			
			if( show ){
				if( (iter + 1) % show_iter == 0 ){
					Rcpp::Rcout << "   Iter = " << iter+1 
						<< ";max_diff = " << max_diff << "\n";
					if( n1 < 10 && n2 < 10 ){
						Rcpp::Rcout << "log_uu = " << log_uu.t();
						Rcpp::Rcout << "log_vv = " << log_vv.t();
					}
				}
			}
			if( max_diff < conv ){
				if( show ){
					Rcpp::Rcout << "   Iter = " << iter+1 
						<< ";max_diff = " << max_diff << "\n";
					if( n1 < 10 && n2 < 10 ){
						Rcpp::Rcout << "log_uu = " << log_uu.t();
						Rcpp::Rcout << "log_vv = " << log_vv.t();
					}
				}
				break;
			}
		}
		
		iter++;
	}
	
	// first calculate log_OT(ii,jj), then exponentiate
	arma::mat OT = arma::zeros<arma::mat>(n1,n2);
	for(ii = 0; ii < n1; ii++){
	for(jj = 0; jj < n2; jj++){
		OT.at(ii,jj) = log_uu.at(ii) + log_KK.at(ii,jj) + log_vv.at(jj);
	}}
	OT = arma::exp(OT);
	
	if( show ){
		Rcpp::Rcout << "   OT_DIST = " << arma::accu(OT % COST_XY) << "\n";
	}
	
	return OT;
}

// [[Rcpp::export(Rcpp_run_OT)]]
arma::mat Rcpp_run_OT(const arma::vec& XX,const arma::vec& YY,
	const arma::mat& COST_XY,const double& EPS,const double& LAMBDA1,
	const double& LAMBDA2,const bool& balance,const bool& highLAM_lowMU,
	const double& conv,const arma::uword& max_iter,const bool& show,
	const arma::uword& show_iter){
	
	// Check dim of XX and YY match COST_XY
	if( XX.n_elem != COST_XY.n_rows || YY.n_elem != COST_XY.n_cols ){
		Rcpp::stop("Dimension mismatch!");
	}
	
	arma::vec XX_2 = XX, YY_2 = YY;
	double sum_XX = arma::sum(XX), sum_YY = arma::sum(YY);
	if( balance == true ){
		XX_2 /= sum_XX;
		YY_2 /= sum_YY;
	}
	
	if( show ){
		Rcpp::Rcout << "sum(XX) = " << arma::sum(XX_2)
			<< "; sum(YY) = " << arma::sum(YY_2) << "\n";
		Rcpp::Rcout << "sum(XX) - sum(YY) -> " 
			<< (arma::sum(XX_2) - arma::sum(YY_2)) << "\n";
	}
	
	if( balance == true ){
		if( show ) Rcpp::Rcout << "Run balanced OT ...\n";
		return run_OT_OPT(XX_2,YY_2,COST_XY,EPS,0.0,0.0,
			conv,max_iter,show,show_iter);
	} else {
		if( show ) Rcpp::Rcout << "Run unbalanced OT ...\n";
		arma::vec LAMBDAs = {LAMBDA1,LAMBDA2};
		double LAMBDA1_final, LAMBDA2_final;
		if( highLAM_lowMU ){
			// Match subject with lower MU(mass) with higher LAMBDA penalty
			if( sum_XX < sum_YY ){
				LAMBDA1_final = arma::max(LAMBDAs);
				LAMBDA2_final = arma::min(LAMBDAs);
			} else { // sum_YY < sum_XX
				LAMBDA1_final = arma::min(LAMBDAs);
				LAMBDA2_final = arma::max(LAMBDAs);
			}
		} else { // Otherwise
			if( sum_XX < sum_YY ){
				LAMBDA1_final = arma::min(LAMBDAs);
				LAMBDA2_final = arma::max(LAMBDAs);
			} else { // sum_YY < sum_XX
				LAMBDA1_final = arma::max(LAMBDAs);
				LAMBDA2_final = arma::min(LAMBDAs);
			}
		}
		return run_OT_OPT(XX_2,YY_2,COST_XY,EPS,
			LAMBDA1_final,LAMBDA2_final,conv,max_iter,
			show,show_iter);
	}
	
}

// [[Rcpp::export(Rcpp_run_full_OT)]]
Rcpp::List Rcpp_run_full_OT(const arma::mat& COST,
	const arma::mat& ZZ,const double& EPS,
	const double& LAMBDA1,const double& LAMBDA2,
	const bool& balance = false,const bool& highLAM_lowMU = true,
	const double& conv = 1e-5,const arma::uword& max_iter = 3e3,
	const int& ncores = 1,const bool& show = true,
	const arma::uword& show_iter = 50){
	
	arma::uword NN = ZZ.n_cols;
	arma::mat DIST = arma::zeros<arma::mat>(NN,NN),
		sum_OT = DIST;
	bool show2 = show && ncores == 1;
	
	#ifdef _OPENMP
	# pragma omp parallel for collapse(2) schedule(dynamic) \
		num_threads(ncores) \
		shared(NN,COST,ZZ,EPS,LAMBDA1,LAMBDA2,\
			balance,highLAM_lowMU,conv,max_iter,show2,\
			show_iter,DIST,sum_OT)
	#endif
	for(arma::uword ii = 0; ii < NN; ii++){
	for(arma::uword jj = 0; jj < NN; jj++){
		
		if( show2 ){
			if( jj == 0 ) Rcpp::Rcout << "ii = " << ii + 1 << ": ";
			if( (jj + 1) % 5 == 0 ) Rcpp::Rcout << ".";
			if( (jj + 1) % 100 == 0 || (jj + 1) == NN )
				Rcpp::Rcout << (jj + 1) << " out of " << NN << "\n";
		}
		
		if( ii < jj ){
			// if LAM1 == LAM2 or balanced OT, we don't need 
			//		to calculate the upper triangle
			continue;
		}
		
		// Get XX, YY, COST_XY
		arma::vec XX = ZZ.col(ii), YY = ZZ.col(jj);
		arma::mat COST_XY = COST.submat(arma::find(XX > 0.0),
			arma::find(YY > 0.0));
		
		// Run OT
		arma::mat OT = Rcpp_run_OT(XX.elem(arma::find(XX > 0.0)),
			YY.elem(arma::find(YY > 0.0)),COST_XY,EPS,LAMBDA1,LAMBDA2,
			balance,highLAM_lowMU,conv,max_iter,false,show_iter);
		
		// Calculate DIST and sum_OT
		DIST.at(ii,jj) 		= arma::accu(OT % COST_XY);
		sum_OT.at(ii,jj) 	= arma::accu(OT);
		
		if( ii != jj ){
			DIST.at(jj,ii) = DIST.at(ii,jj);
			sum_OT.at(jj,ii) = sum_OT.at(ii,jj);
		}
		
	}}
	
	// Output
	return Rcpp::List::create(Rcpp::Named("DIST",DIST),
		Rcpp::Named("sum_OT",sum_OT));
	
}


// --------------------
// Kernel Matrices and Omnibus Hypothesis Testing
// --------------------

// [[Rcpp::export(Rcpp_KernTest)]]
Rcpp::List Rcpp_KernTest(const arma::vec& RESI,
	const arma::cube& cKK,const arma::umat& OMNI,
	const arma::uword& nPERMS = 2e3){
	
	arma::uword kk, pp, NN = RESI.n_elem,
		oo, nKK = cKK.n_slices;
	arma::vec PVALs = arma::zeros<arma::vec>(nKK);
	arma::mat pRESI = arma::zeros<arma::mat>(nPERMS,NN),
		pSTAT = arma::zeros<arma::mat>(nPERMS + 1,nKK),
		tmp_mat = arma::zeros<arma::mat>(nPERMS,nPERMS);
	arma::uvec tmp_vec = arma::zeros<arma::uvec>(nPERMS + 1);
	
	// Store permuted residuals
	for(pp = 0; pp < nPERMS; pp++){
		pRESI.row(pp) = arma::shuffle(RESI).t();
	}
	
	// Calculate test-statistics
	for(kk = 0; kk < nKK; kk++){
		
		// Calculate unpermuted test-statistics
		pSTAT.at(0,kk) = arma::dot(RESI,cKK.slice(kk) * RESI);
		
		// Calculate permuted test-statistics
		tmp_mat = pRESI * cKK.slice(kk) * pRESI.t();
		pSTAT(arma::span(1,nPERMS),kk) = tmp_mat.diag();
		
		// Transform statistics to empirical distribution 
		//	b/c kernels can be on different scales
		pSTAT.col(kk) = 1.0 - ( Rcpp_calc_rank(pSTAT.col(kk)) 
			- 1.0 ) / (nPERMS + 1.0);
		
		// Permutation p-value per kernel
		PVALs.at(kk) = pSTAT.at(0,kk);
		
	}
	
	// Calculate omnibus p-values
	arma::vec omni_PVALs = arma::zeros<arma::vec>(OMNI.n_rows);
	
	for(oo = 0; oo < OMNI.n_rows; oo++){
		arma::uvec tmp_cols = arma::find(OMNI.row(oo).t() == 1);
		arma::mat tmp_mat = pSTAT.cols(tmp_cols);
		omni_PVALs.at(oo) = arma::sum( arma::min(tmp_mat,1) 
			<= arma::min(tmp_mat.row(0).t()) ) * 1.0 / (nPERMS + 1.0);
	}
	
	return Rcpp::List::create(
		Rcpp::Named("PVALs",
			Rcpp::NumericVector(PVALs.begin(),PVALs.end())),
		Rcpp::Named("omni_PVALs",
			Rcpp::NumericVector(omni_PVALs.begin(),omni_PVALs.end()))
		);
	
}



