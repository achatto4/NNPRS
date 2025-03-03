#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Gradient descent function for transfer learning in PRS
// [[Rcpp::export]]
Rcpp::List gradient_descent_transfer_learning_rcpp(
    double n0,
    arma::vec r0,  
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

// Gradient descent function for PRS without auxiliary dataset
// [[Rcpp::export]]
Rcpp::List gradient_descent_main_only(
    double n0,
    arma::vec r0,  
    arma::mat R0,
    double alpha3,
    double alpha4,
    double eta_m,
    int max_iter
) {
  
  int p = R0.n_rows;  // Assuming square matrices, p is the number of SNPs (rows in R0)
  
  // Initialize variables
  arma::vec h_m = alpha3 * arma::ones<arma::vec>(p);
  arma::vec g_m = alpha4 * arma::ones<arma::vec>(p);
  
  // Gradient descent for main data
  for (int m = 0; m <= max_iter; ++m) {
    arma::vec h_m_sq = arma::square(h_m);
    arma::vec g_m_sq = arma::square(g_m);
    arma::vec diff_sq = h_m_sq - g_m_sq;
    
    arma::vec grad_h = (-4 * n0 * r0 % h_m + 4 * n0 * (R0 * diff_sq) % h_m);
    arma::vec grad_g = (4 * n0 * r0 % g_m - 4 * n0 * (R0 * diff_sq) % g_m);
    
    h_m -= (eta_m / n0) * grad_h;
    g_m -= (eta_m / n0) * grad_g;
  }
  
  arma::vec hat_h = h_m;
  arma::vec hat_g = g_m;
  
  // Compute the final beta estimate
  arma::vec hat_beta = arma::square(hat_h) - arma::square(hat_g);
  
  return Rcpp::List::create(
    Rcpp::Named("hat_h") = hat_h,
    Rcpp::Named("hat_g") = hat_g,
    Rcpp::Named("hat_beta") = hat_beta
  );
}


// [[Rcpp::depends(RcppArmadillo)]]

// Gradient descent function for transfer learning (Main Data Only)
// [[Rcpp::export]]
Rcpp::List gradient_descent_main_P4(
    double n0,
    arma::vec r0,
    arma::mat R0,
    double alpha3,
    double alpha4,
    double eta_m,
    int max_iter
) {
  int p = R0.n_rows; // Number of SNPs
  
  // Initialize variables
  arma::vec h_m = alpha3 * arma::ones<arma::vec>(p);
  arma::vec g_m = alpha4 * arma::ones<arma::vec>(p);
  
  // Gradient descent for main data
  for (int m = 0; m <= max_iter; ++m) {
    arma::vec h_m_4th = arma::pow(h_m, 4);
    arma::vec g_m_4th = arma::pow(g_m, 4);
    arma::vec grad_h = arma::zeros<arma::vec>(p);
    arma::vec grad_g = arma::zeros<arma::vec>(p);
    
    arma::vec diff_4th = h_m_4th - g_m_4th;
    
    grad_h = (-8 * n0 * r0 % arma::pow(h_m, 3) + 8 * n0 * (R0 * diff_4th) % arma::pow(h_m, 3));
    grad_g = (8 * n0 * r0 % arma::pow(g_m, 3) - 8 * n0 * (R0 * diff_4th) % arma::pow(g_m, 3));
    
    h_m -= (eta_m / n0) * grad_h;
    g_m -= (eta_m / n0) * grad_g;
  }
  
  arma::vec hat_h = h_m;
  arma::vec hat_g = g_m;
  
  // Compute the final beta estimate with 4th power terms
  arma::vec hat_beta = arma::pow(hat_h, 4) - arma::pow(hat_g, 4);
  
  return Rcpp::List::create(
    Rcpp::Named("hat_h") = hat_h,
    Rcpp::Named("hat_g") = hat_g,
    Rcpp::Named("hat_beta") = hat_beta
  );
}
