#include <RcppArmadillo.h>
#include <omp.h>
#include <smartr.h>

// [[Rcpp::depends(RcppArmadillo)]]

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
			log_tt.at(kk) = smartr::Rcpp_logSumExp(log_KK.row(kk).t() + log_vv);
		}
		
		// Update log_uu = cnst_u * ( log(XX) - log(tt) )
		log_uu2 = cnst_u * ( log_XX - log_tt );
		if( iter > 0 ) diff_uu = arma::max(arma::abs(log_uu - log_uu2));
		log_uu = log_uu2;
		
		// Update log_ss = log( KK.t() * uu )
		for(ll = 0; ll < n2; ll++){
			log_ss.at(ll) = smartr::Rcpp_logSumExp(log_KK.col(ll) + log_uu);
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
	
	omp_set_num_threads(ncores);
	# pragma omp parallel for collapse(2) schedule(static) \
		shared(NN,COST,ZZ,EPS,LAMBDA1,LAMBDA2,balance,highLAM_lowMU,\
		conv,max_iter,show2,show_iter,DIST,sum_OT)
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
		arma::mat COST_XY = COST.submat(arma::find(XX == 1.0),
			arma::find(YY == 1.0));
		
		// Run OT
		arma::mat OT = Rcpp_run_OT(XX.elem(arma::find(XX == 1.0)),
			YY.elem(arma::find(YY == 1.0)),COST_XY,EPS,LAMBDA1,LAMBDA2,
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
	const Rcpp::List& KK,const arma::uword& nPERMS = 2e3,
	const arma::uword& iter1 = 50,const arma::uword& iter2 = 1e3,
	const bool& verbose = false){
	
	arma::uword kk, pp, NN = RESI.n_elem,
		nKK = KK.length();
	arma::vec PVALs = arma::zeros<arma::vec>(nKK);
	arma::mat pRESI = arma::zeros<arma::mat>(nPERMS,NN),
		pSTAT = arma::zeros<arma::mat>(nPERMS + 1,nKK);
	arma::uvec tmp_vec = arma::zeros<arma::uvec>(nPERMS + 1);
	
	// Store KK as cube
	arma::cube cKK = arma::zeros<arma::cube>(NN,NN,nKK);
	for(kk = 0; kk < nKK; kk++){
		cKK.slice(kk) = Rcpp::as<arma::mat>(KK[kk]);
	}
	
	// Store permuted residuals
	for(pp = 0; pp < nPERMS; pp++){
		pRESI.row(pp) = arma::shuffle(RESI).t();
	}
	
	// Calculate test-statistics
	for(kk = 0; kk < nKK; kk++){
		if( verbose ) Rcpp::Rcout << "kk = " << kk + 1 << ": ";
		
		// Unpermuted test-statistics
		pSTAT.at(0,kk) = arma::dot(RESI,cKK.slice(kk) * RESI);
		
		// Permuted test-statistics
		for(pp = 0; pp < nPERMS; pp++){
			if( verbose ){
				if( (pp+1) % iter1 == 0 ) Rcpp::Rcout << ".";
				if( (pp+1) % iter2 == 0 || (pp+1) == nPERMS )
					Rcpp::Rcout << (pp+1) << " out of " << nPERMS << "\n";
			}
			pSTAT.at(pp + 1,kk) = arma::dot(pRESI.row(pp) * cKK.slice(kk),pRESI.row(pp).t());
		}
		
	}
	
	// Calculate p-values
	for(kk = 0; kk < nKK; kk++){
		// Permutation p-value per kernel
		PVALs.at(kk) = arma::sum(pSTAT.col(kk) >= pSTAT.at(0,kk)) * 1.0 / (nPERMS + 1.0);
		
		// Transform statistics to empirical distribution b/c kernels can be on different scales
		tmp_vec = arma::sort_index(arma::sort_index(pSTAT.col(kk)));
		pSTAT.col(kk) = 1.0 - arma::conv_to<arma::vec>::from(tmp_vec) / (nPERMS + 1.0);
	}
	
	// Calculate ominibus p-value
	double omni_PVAL = arma::sum(arma::min(pSTAT,1) 
		<= arma::min(pSTAT.row(0).t())) * 1.0 / (nPERMS + 1.0);
	
	return Rcpp::List::create(
		Rcpp::Named("PVALs",Rcpp::NumericVector(PVALs.begin(),PVALs.end())),
		Rcpp::Named("omni_PVAL",omni_PVAL));
	
}







