// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}


// [[Rcpp::export]]
arma::sp_mat mult_sp_sp_to_sp(const arma::sp_mat& a, const arma::sp_mat& b) {
    // sparse x sparse -> sparse
    arma::sp_mat result(a * b);
    return result;
}

// [[Rcpp::export]]
arma::sp_mat mult_sp_den_to_sp(const arma::sp_mat& a, const arma::mat& b) {
    // sparse x dense -> sparse
    arma::sp_mat result(a * b);
    return result;
}

// [[Rcpp::export]]
arma::sp_mat mult_den_sp_to_sp(const arma::mat& a, const arma::sp_mat& b) {
    // dense x sparse -> sparse
    arma::sp_mat result(a * b);
    return result;
}