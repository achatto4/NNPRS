#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Gradient descent function for transfer learning
//' Compute beta joint
 //'
 //' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
 //' @param r0: summary statistics for base ancestry (set to 0 for ethnic-specific SNPs)
 //' @param n0: sample size for base ancestry 
 //' @param rk_list: list of summary statistics for auxillary ancestry (set to 0 for ethnic-specific SNPs)
 //' @param R0: correlation matrix of all SNPs for base ancestry (set to 0 for ethnic-specific SNPs)
 //' @param nk_list: list of sample size of auxillary ancestries
 //' @param Rk_list: correlation matrix of all SNPs for auxillary ancestry (set to 0 for ethnic-specific SNPs)
 //' @param indx: indicator of whether a SNP in a certain ethnic group
 //' @param M: number of ethnic groups
 //' @param alpha1-4: initialization parameters
 //' @return beta vector
// [[Rcpp::export]]
Rcpp::List gradient_descent_transfer_learning_rcpp_PRS(
    double n0,
    arma::vec r0,  // r0 is now a vector
    arma::mat R0,
    std::vector<double> nk_list,
    std::vector<arma::vec> rk_list,
    std::vector<arma::mat> Rk_list,
    double alpha1,
    double alpha2,
    double alpha3,
    double alpha4,
    double eta_l,
    double eta_m,
    int max_iter
) {
  
  int p = R0.n_rows;  // Assuming square matrices, p is the number of SNPs (rows in R0)
  
  // Initialize variables
  arma::vec u_l = alpha1 * arma::ones<arma::vec>(p);
  arma::vec v_l = alpha2 * arma::ones<arma::vec>(p);
  arma::vec h_m = alpha3 * arma::ones<arma::vec>(p);
  arma::vec g_m = alpha4 * arma::ones<arma::vec>(p);
  
  // Gradient descent for auxiliary data
  for (int l = 0; l <= max_iter; ++l) {
    arma::vec grad_u = arma::zeros<arma::vec>(p);
    arma::vec grad_v = arma::zeros<arma::vec>(p);
    
    for (size_t k = 0; k < nk_list.size(); ++k) {
      double nk = nk_list[k];
      arma::vec rk = rk_list[k];
      arma::mat Rk = Rk_list[k];
      
      arma::vec u_l_sq = arma::square(u_l);
      arma::vec v_l_sq = arma::square(v_l);
      arma::vec diff_sq = u_l_sq - v_l_sq;
      
      grad_u += (-4 * nk * rk % u_l + 4 * nk * (Rk * diff_sq) % u_l);
      grad_v += (4 * nk * rk % v_l - 4 * nk * (Rk * diff_sq) % v_l);
    }
    
    double total_n = std::accumulate(nk_list.begin(), nk_list.end(), 0.0) + n0;
    u_l -= (eta_l / total_n) * grad_u;
    v_l -= (eta_l / total_n) * grad_v;
  }
  
  arma::vec hat_u = u_l;
  arma::vec hat_v = v_l;
  
  // Gradient descent for main data
  for (int m = 0; m <= max_iter; ++m) {
    arma::vec u_hat_sq = arma::square(hat_u);
    arma::vec v_hat_sq = arma::square(hat_v);
    arma::vec grad_h = arma::zeros<arma::vec>(p);
    arma::vec grad_g = arma::zeros<arma::vec>(p);
    
    arma::vec h_m_sq = arma::square(h_m);
    arma::vec g_m_sq = arma::square(g_m);
    arma::vec diff_sq = u_hat_sq - v_hat_sq + h_m_sq - g_m_sq;
    
    grad_h = (-4 * n0 * r0 % h_m + 4 * n0 * (R0 * diff_sq) % h_m);
    grad_g = (4 * n0 * r0 % g_m - 4 * n0 * (R0 * diff_sq) % g_m);
    
    h_m -= (eta_m / n0) * grad_h;
    g_m -= (eta_m / n0) * grad_g;
  }
  
  arma::vec hat_h = h_m;
  arma::vec hat_g = g_m;
  
  // Define the output
  arma::vec hat_beta = arma::square(hat_u) - arma::square(hat_v) + arma::square(hat_h) - arma::square(hat_g);
  
  return Rcpp::List::create(
    Rcpp::Named("hat_u") = hat_u,
    Rcpp::Named("hat_v") = hat_v,
    Rcpp::Named("hat_h") = hat_h,
    Rcpp::Named("hat_g") = hat_g,
    Rcpp::Named("hat_beta") = hat_beta
  );
}
 
 
 // Function to replace NaN values with 0 in an arma::vec or arma::mat
 void replace_nan_with_zero(arma::mat& mat) {
   for (size_t i = 0; i < mat.n_elem; ++i) {
     if (std::isnan(mat[i])) {
       mat[i] = 0.0;
     }
   }
 }
 
 void replace_nan_with_zero(arma::vec& vec) {
   for (size_t i = 0; i < vec.n_elem; ++i) {
     if (std::isnan(vec[i])) {
       vec[i] = 0.0;
     }
   }
 }
 
 // Function to apply gradient descent transfer learning across all blocks
 // [[Rcpp::export]]
 Rcpp::List gradient_descent_transfer_learning_all_blocks(
     Rcpp::List summ_list,
     Rcpp::List LD_list,
     int M,
     Rcpp::List indx,
     Rcpp::IntegerVector indx_block,
     double n0,
     std::vector<double> nk_list,
     double alpha1,
     double alpha2,
     double alpha3,
     double alpha4,
     double eta_l,
     double eta_m,
     int max_iter
 ) {
   int num_blocks = indx_block.size();
   Rcpp::List beta_results(num_blocks);
   
   for (int bl = 0; bl < num_blocks; ++bl) {
     if (indx_block[bl] == 0) {
       beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = R_NilValue);
       continue;
     }
     
     // Extract summary statistics and LD matrices
     arma::mat summ = Rcpp::as<arma::mat>(summ_list[bl]);
     Rcpp::List R = LD_list[bl];
     arma::mat indx_mat = Rcpp::as<arma::mat>(indx[bl]);
     
     // Replace NaN values with 0 in summ, R0, and indx_mat
     replace_nan_with_zero(summ);
     arma::mat R0 = Rcpp::as<arma::mat>(R[0]);
     replace_nan_with_zero(R0);
     replace_nan_with_zero(indx_mat);
     
     // Prepare rk_list: all columns except the first from summ
     std::vector<arma::vec> rk_list;
     for (size_t i = 1; i < summ.n_cols; ++i) {
       arma::vec rk = summ.col(i);
       replace_nan_with_zero(rk);  // Replace NaN in rk
       rk_list.push_back(rk);
     }
     
     // Prepare Rk_list: all LD matrices except the first (target population)
     std::vector<arma::mat> Rk_list;
     for (size_t i = 1; i < R.size(); ++i) {
       arma::mat Rk = Rcpp::as<arma::mat>(R[i]);
       replace_nan_with_zero(Rk);  // Replace NaN in Rk
       Rk_list.push_back(Rk);
     }
     
     // Call gradient descent function with corrected inputs
     Rcpp::List beta_block = gradient_descent_transfer_learning_rcpp_PRS(
       n0, 
       summ.col(0), // r0 is the first column
       R0, // R0 is the first LD matrix
       nk_list, 
       rk_list, 
       Rk_list,
       alpha1, alpha2, alpha3, alpha4, 
       eta_l, eta_m, max_iter
     );
     
     // Extract beta vector for the target population
     arma::vec beta_vec = as<arma::vec>(beta_block["hat_beta"]);
     
     arma::vec indx_binary = arma::conv_to<arma::vec>::from(indx_mat.col(0) != 0);
     
     // Element-wise multiplication to apply indexing
     arma::vec beta_final = beta_vec % indx_binary;
     
     // Store the final beta vector in the results
     beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = beta_final);
   }
   
   return beta_results;
 }