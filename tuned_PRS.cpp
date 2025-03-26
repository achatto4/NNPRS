#include <RcppArmadillo.h>
#include <cmath> // For std::isnan()
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gradient_descent_h_g_ADAM(
    arma::vec y0, arma::mat X0, 
    arma::vec beta_summ, arma::vec h_hat, arma::vec g_hat,
    double eta, double beta1, double beta2, double epsilon, int max_iter) {
  
  int p = X0.n_cols;
  
  // Initialize h and g from given inputs
  arma::vec h = h_hat;
  arma::vec g = g_hat;
  
  // Initialize ADAM momentum terms
  arma::vec m_h = arma::zeros<arma::vec>(p), v_h = arma::zeros<arma::vec>(p);
  arma::vec m_g = arma::zeros<arma::vec>(p), v_g = arma::zeros<arma::vec>(p);
  
  // Optimization loop for h and g
  for (int m = 1; m <= max_iter; ++m) {
    arma::vec residual0 = y0 - X0 * (beta_summ + arma::square(h) - arma::square(g));
    
    arma::vec grad_h = -2 * X0.t() * (residual0 % h);
    arma::vec grad_g = 2 * X0.t() * (residual0 % g);
    
    // ADAM update
    m_h = beta1 * m_h + (1 - beta1) * grad_h;
    v_h = beta2 * v_h + (1 - beta2) * arma::square(grad_h);
    m_g = beta1 * m_g + (1 - beta1) * grad_g;
    v_g = beta2 * v_g + (1 - beta2) * arma::square(grad_g);
    
    arma::vec m_h_hat = m_h / (1 - std::pow(beta1, m));
    arma::vec v_h_hat = v_h / (1 - std::pow(beta2, m));
    arma::vec m_g_hat = m_g / (1 - std::pow(beta1, m));
    arma::vec v_g_hat = v_g / (1 - std::pow(beta2, m));
    
    h -= (eta / X0.n_rows) * m_h_hat / (arma::sqrt(v_h_hat) + epsilon);
    g -= (eta / X0.n_rows) * m_g_hat / (arma::sqrt(v_g_hat) + epsilon);
  }
  
  // Final estimates
  arma::vec h_final = h;
  arma::vec g_final = g;
  arma::vec beta_hat = beta_summ + arma::square(h_final) - arma::square(g_final);
  
  return Rcpp::List::create(
    Rcpp::Named("hat_h") = h_final,
    Rcpp::Named("hat_g") = g_final,
    Rcpp::Named("hat_beta") = beta_hat
  );
}


// // [[Rcpp::depends(RcppArmadillo)]]
// 
// // Gradient descent function for transfer learning
// //' Compute beta joint
//  //'
//  //' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
//  //' @param r0: summary statistics for base ancestry (set to 0 for ethnic-specific SNPs)
//  //' @param n0: sample size for base ancestry 
//  //' @param rk_list: list of summary statistics for auxillary ancestry (set to 0 for ethnic-specific SNPs)
//  //' @param R0: correlation matrix of all SNPs for base ancestry (set to 0 for ethnic-specific SNPs)
//  //' @param nk_list: list of sample size of auxillary ancestries
//  //' @param Rk_list: correlation matrix of all SNPs for auxillary ancestry (set to 0 for ethnic-specific SNPs)
//  //' @param indx: indicator of whether a SNP in a certain ethnic group
//  //' @param M: number of ethnic groups
//  //' @param alpha1-4: initialization parameters
//  //' @return beta vector
//  // [[Rcpp::export]]
//  Rcpp::List gradient_descent_transfer_learning_rcpp_PRS(
//      double n0,
//      arma::vec r0,  // r0 is now a vector
//      arma::mat R0,
//      std::vector<double> nk_list,
//      std::vector<arma::vec> rk_list,
//      std::vector<arma::mat> Rk_list,
//      double alpha1,
//      double alpha2,
//      double alpha3,
//      double alpha4,
//      double eta_l,
//      double eta_m,
//      int max_iter
//  ) {
//    
//    int p = R0.n_rows;  // Assuming square matrices, p is the number of SNPs (rows in R0)
//    
//    // Initialize variables
//    arma::vec u_l = alpha1 * arma::ones<arma::vec>(p);
//    arma::vec v_l = alpha2 * arma::ones<arma::vec>(p);
//    arma::vec h_m = alpha3 * arma::ones<arma::vec>(p);
//    arma::vec g_m = alpha4 * arma::ones<arma::vec>(p);
//    
//    // Gradient descent for auxiliary data
//    for (int l = 0; l <= max_iter; ++l) {
//      arma::vec grad_u = arma::zeros<arma::vec>(p);
//      arma::vec grad_v = arma::zeros<arma::vec>(p);
//      
//      for (size_t k = 0; k < nk_list.size(); ++k) {
//        double nk = nk_list[k];
//        arma::vec rk = rk_list[k];
//        arma::mat Rk = Rk_list[k];
//        
//        arma::vec u_l_sq = arma::square(u_l);
//        arma::vec v_l_sq = arma::square(v_l);
//        arma::vec diff_sq = u_l_sq - v_l_sq;
//        
//        grad_u += (-4 * nk * rk % u_l + 4 * nk * (Rk * diff_sq) % u_l);
//        grad_v += (4 * nk * rk % v_l - 4 * nk * (Rk * diff_sq) % v_l);
//      }
//      
//      double total_n = std::accumulate(nk_list.begin(), nk_list.end(), 0.0) + n0;
//      u_l -= (eta_l / total_n) * grad_u;
//      v_l -= (eta_l / total_n) * grad_v;
//    }
//    
//    arma::vec hat_u = u_l;
//    arma::vec hat_v = v_l;
//    
//    // Gradient descent for main data
//    for (int m = 0; m <= max_iter; ++m) {
//      arma::vec u_hat_sq = arma::square(hat_u);
//      arma::vec v_hat_sq = arma::square(hat_v);
//      arma::vec grad_h = arma::zeros<arma::vec>(p);
//      arma::vec grad_g = arma::zeros<arma::vec>(p);
//      
//      arma::vec h_m_sq = arma::square(h_m);
//      arma::vec g_m_sq = arma::square(g_m);
//      arma::vec diff_sq = u_hat_sq - v_hat_sq + h_m_sq - g_m_sq;
//      
//      grad_h = (-4 * n0 * r0 % h_m + 4 * n0 * (R0 * diff_sq) % h_m);
//      grad_g = (4 * n0 * r0 % g_m - 4 * n0 * (R0 * diff_sq) % g_m);
//      
//      h_m -= (eta_m / n0) * grad_h;
//      g_m -= (eta_m / n0) * grad_g;
//      
//    }
//    
//    arma::vec hat_h = h_m;
//    arma::vec hat_g = g_m;
//    
//    // Define the output
//    arma::vec hat_beta = arma::square(hat_u) - arma::square(hat_v) + arma::square(hat_h) - arma::square(hat_g);
//    
//    return Rcpp::List::create(
//      Rcpp::Named("hat_u") = hat_u,
//      Rcpp::Named("hat_v") = hat_v,
//      Rcpp::Named("hat_h") = hat_h,
//      Rcpp::Named("hat_g") = hat_g,
//      Rcpp::Named("hat_beta") = hat_beta
//    );
//  }
//  
//  // [[Rcpp::export]]
//  void replace_nan_with_zero(arma::mat &mat) {
//    mat.elem(arma::find_nonfinite(mat)).zeros(); // Replace NaN and Inf with 0
//  }
//  
//  
//  // Function to apply gradient descent transfer learning across all blocks
//  // [[Rcpp::export]]
//  Rcpp::List gradient_descent_transfer_learning_all_blocks(
//      Rcpp::List summ_list,
//      Rcpp::List LD_list,
//      int M,
//      Rcpp::List indx,
//      Rcpp::IntegerVector indx_block,
//      double n0,
//      std::vector<double> nk_list,
//      double alpha1,
//      double alpha2,
//      double alpha3,
//      double alpha4,
//      double eta_l,
//      double eta_m,
//      int max_iter,
//      bool adaptive,
//      bool alpha_adaptive,
//      double eta,
//      double beta1,
//      double beta2,
//      double epsilon
//  ) {
//    int num_blocks = indx_block.size();
//    Rcpp::List beta_results(num_blocks);
//    
//    for (int bl = 0; bl < num_blocks; ++bl) {
//      if (indx_block[bl] == 0) {
//        beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = R_NilValue);
//        continue;
//      }
//      
//      // Extract summary statistics and LD matrices
//      arma::mat summ = Rcpp::as<arma::mat>(summ_list[bl]);
//      replace_nan_with_zero(summ); // Ensure no NaNs
//      arma::mat indx_mat = Rcpp::as<arma::mat>(indx[bl]);
//      
//      
//      Rcpp::List R = LD_list[bl];
//      // Prepare rk_list: all columns except the first from summ, only if M != 1
//      std::vector<arma::vec> rk_list;
//      if (M != 1) {
//        for (size_t i = 1; i < summ.n_cols; ++i) {
//          arma::vec col = summ.col(i);
//          col.elem(arma::find_nonfinite(col)).zeros(); // Replace NaNs in column
//          rk_list.push_back(col);
//        }
//      }
//      
//      // Prepare Rk_list: all LD matrices except the first (target population), only if M != 1
//      std::vector<arma::mat> Rk_list;
//      if (M != 1) {
//        for (size_t i = 1; i < R.size(); ++i) {
//          arma::mat Rk = Rcpp::as<arma::mat>(R[i]);
//          replace_nan_with_zero(Rk); // Ensure no NaNs in LD matrices
//          Rk_list.push_back(Rk);
//        }
//      }
//      double inv_num_rows = 1.0 / std::sqrt(indx_mat.n_rows);
//      // Call gradient descent function with correct inputs
//      Rcpp::List beta_block;
//         beta_block = gradient_descent_transfer_learning_rcpp_PRS(
//            n0,
//            summ.col(0),
//            Rcpp::as<arma::mat>(R[0]),
//            nk_list,
//            rk_list,
//            Rk_list,
//            alpha1, alpha2, alpha3, alpha4,
//            eta_l, eta_m, max_iter
//          );
//      
//      // Extract beta vector for the target population
//      arma::vec beta_vec = as<arma::vec>(beta_block["hat_beta"]);
//      arma::vec indx_binary = arma::conv_to<arma::vec>::from(indx_mat.col(0) != 0);
//      
//      // Element-wise multiplication to apply indexing
//      arma::vec beta_final = beta_vec % indx_binary;
//      
//      // Store the final beta vector in the results
//      beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = beta_final);
//    }
//    
//    return beta_results;
//  }
 