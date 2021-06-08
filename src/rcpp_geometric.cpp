#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>       /* pow */
using namespace Rcpp;
#include <limits>
#include "geometric.h"
//' lgeo_cpp
//'
//' Cpp version of lgeo, computing likelihood of geometric distribution
//' @param Y (n x 1) observations
//' @param Theta (n x 1) param
//' @param logl (bool) use log likelihood or likelihood
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double lgeo_cpp( Eigen::VectorXd&  Y, Eigen::VectorXd&  Theta, bool logl = true) {

  double lik = 0.;
  for(int i = 0; i < Y.size(); i++)
    lik += (Y[i] - 1) * log(1 - 1. /  Theta[i]) - log(Theta[i]);

  if(logl)
    return lik;

  return exp(lik);
}

//' dd_lgeo_cpp
//'
//' Cpp version of dd_lgeo, computing likelihood, gradient of the likehood, and Hessian
//' @param Y (n x 1) observations
//' @param Theta (n x 1) param
//' @param calc_ddlog (bool) caculate second
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List dd_lgeo_cpp( Eigen::VectorXd&  Y,
                  Eigen::VectorXd&  Theta){
  return(dd_lgeo_cpp_internal(Y, Theta));
}

List dd_lgeo_cpp_internal( const Eigen::VectorXd&  Y,
                           const Eigen::VectorXd&  Theta,
                           bool calc_ddlog) {

  double lik = 0.;
  Eigen::VectorXd dlog  = Eigen::VectorXd::Zero(Y.size());
  Eigen::VectorXd ddlog = Eigen::VectorXd::Zero(Y.size());

  for(int i = 0; i < Y.size(); i++){

    double p = 1. /  Theta[i];
    double p_m1 =1 /(Theta[i] - 1);
    lik += (Y[i] - 1) * log(1 - p) + log(p);
    double Y_p_p_m1 = (Y[i] - 1) * p * p_m1;
    dlog[i]  = - p + Y_p_p_m1;
    if(calc_ddlog == true)
      ddlog[i] = - p * dlog[i] - Y_p_p_m1 * p_m1;

  }
  List result;
  if(std::isnan(lik) )
    lik = -std::numeric_limits<double>::max();
  result["lik"]   = lik;
  result["dlogl"]  = dlog;
  result["ddlogl"] = ddlog;
  return result;
}




//' dd_lgeom_cpp
//'
//' Cpp version of dd_lgeom, computing likelihood, gradient of the likehood, and Hessian
//' For geometric distribution with covariates and mean like
//' E[Y_i|X_i] = 1 + exp(X_i%*%beta)
//' @param beta (k x 1) coeffients
//' @param Y (n x 1) observations
//' @param X (n x k) param
//' @param offset (n x 1) could be zero length then not used
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List dd_lgeom_cpp(const Eigen::VectorXd&  beta,
                  Eigen::VectorXd&  Y,
                  Eigen::MatrixXd&  X,
                  Eigen::VectorXd&  offset,
                  bool              calc_hess=true) {

  return(dd_lgeom(beta, Y, X, offset, calc_hess));
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List dd_lgeomc_cpp(const    Eigen::VectorXd& beta ,
                   const Eigen::VectorXd&    Y,
                   const Eigen::MatrixXd&    X,
                   const Eigen::VectorXd&   offset,
                   const Eigen::VectorXd&    K,
                   bool                     calc_hess= true) {

  return(dd_lgeomc(beta, Y, X, offset, K, calc_hess));
}

Rcpp::List dd_lgeomc(const    Eigen::VectorXd& beta ,
                     const Eigen::VectorXd&    Y,
                     const Eigen::MatrixXd&    X,
                     const Eigen::VectorXd&   offset,
                     const Eigen::VectorXd&    K,
                     bool                     calc_hess
){

  Eigen::VectorXd theta_m1;
  theta_m1 = X * beta;
  if(offset.size() > 0 )
    theta_m1 += offset;
  theta_m1.array()      = theta_m1.array().exp();
  Eigen::VectorXd theta = theta_m1;
  theta.array()         += 1;
  List likelihood_list  =  ddlgeo_const_internal(Y, theta, K, calc_hess);
  Eigen::VectorXd grad_ = Rcpp::as<Eigen::VectorXd>(likelihood_list["dlogl"] );
  grad_.array() *= theta_m1.array();
  Eigen::VectorXd dlogl = X.transpose()  * grad_;


  List result;
  if(calc_hess == true){
    Eigen::VectorXd hess_ = Rcpp::as<Eigen::VectorXd>(likelihood_list["ddlogl"] );
    hess_.array() *= theta_m1.array().square();
    Eigen::MatrixXd temp   = grad_.asDiagonal();
    temp += hess_.asDiagonal();
    Eigen::MatrixXd ddlogl = X.transpose() * temp * X;
    result["hessian"]  = ddlogl;
  }



  result["lik"]      = Rcpp::as<double>(likelihood_list["lik"]);
  result["gradient"] = dlogl;
  return result;

}

//' dd_lgeom
//'
List dd_lgeom(  const  Eigen::VectorXd&  beta,
                const  Eigen::VectorXd&  Y,
                const  Eigen::MatrixXd&  X,
                const  Eigen::VectorXd&  offset,
                  bool              calc_hess) {

  Eigen::VectorXd theta_m1;
  theta_m1 = X * beta;
  if(offset.size() > 0 )
    theta_m1 += offset;
  theta_m1.array()      = theta_m1.array().exp();
  Eigen::VectorXd theta = theta_m1;
  theta.array()         += 1;

  List likelihood_list  = dd_lgeo_cpp_internal(Y, theta, calc_hess);
  Eigen::VectorXd grad_ = Rcpp::as<Eigen::VectorXd>(likelihood_list["dlogl"] );
  grad_.array() *= theta_m1.array();
  Eigen::VectorXd dlogl = X.transpose()  * grad_;


  List result;
  if(calc_hess == true){
    Eigen::VectorXd hess_ = Rcpp::as<Eigen::VectorXd>(likelihood_list["ddlogl"] );
    hess_.array() *= theta_m1.array().square();
    Eigen::MatrixXd temp   = grad_.asDiagonal();
    temp += hess_.asDiagonal();
    Eigen::MatrixXd ddlogl = X.transpose() * temp * X;
    result["hessian"]  = ddlogl;
  }

  result["lik"]      = Rcpp::as<double>(likelihood_list["lik"]);
  result["gradient"] = dlogl;
  return result;
}

//' Computing the constraint geometric likelihood the observations
//' are constrained with K
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double lgeo_const_cpp( Eigen::VectorXd&  Y,
                       Eigen::VectorXd&  Theta,
                       Eigen::VectorXd&  K,
                       bool logl = true) {

  double lik = 0.;

  double p, p_m;
  for(int i = 0; i < Y.size(); i++){
    p = 1. / Theta[i];

    p_m = 1. - p;
    lik += (Y[i] - 1) * log(p_m) + log(p);
    lik -= log(1 - pow(p_m, K[i]));
  }

  if(logl)
    return lik;

  return exp(lik);
}

//' Computing the constraint geometric likelihood the observations
//' are constrained with K
//'
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List ddlgeo_const_cpp( const Eigen::VectorXd&  Y,
                       const Eigen::VectorXd&  Theta,
                       const Eigen::VectorXd&  K,
                       bool calc_ddlog = true){
  return(ddlgeo_const_internal(Y, Theta, K, calc_ddlog));
}
List ddlgeo_const_internal( const Eigen::VectorXd&  Y,
                           const Eigen::VectorXd&  Theta,
                           const Eigen::VectorXd&  K,
                           bool calc_ddlog) {
  double lik = 0.;
  Eigen::VectorXd dlog  = Eigen::VectorXd::Zero(Y.size());
  Eigen::VectorXd ddlog = Eigen::VectorXd::Zero(Y.size());

  for(int i = 0; i < Y.size(); i++){
    double theta_i  = std::min( Theta[i], 1e16);
    double p              = 1. /  theta_i;
    double one_mp         = 1 - p;
    double theta_inv_m1   = 1. /(theta_i - 1);
    double base_const     = 1 - pow(one_mp, K[i]);
    double Y_p_t_m1       = (Y[i] - 1) * p * theta_inv_m1;

    if(log(base_const) == -std::numeric_limits<double>::infinity()){
      lik = std::numeric_limits<double>::infinity();
      dlog[i] = 0;
      ddlog[i] = 0;
      break;
    }
    if(log(one_mp) == -std::numeric_limits<double>::infinity()){
      lik = -std::numeric_limits<double>::infinity();
      dlog[i] = 0;
      ddlog[i] = 0;
      break;
    }

    lik += (Y[i] - 1) * log(one_mp) + log(p);

    lik -= log(base_const); // constraint correction

    double Y_p_p_m1 = (Y[i] - 1) * p * theta_inv_m1;

    double p_m1_powK1 = pow(one_mp, K[i] - 1);
    dlog[i]  = - p + Y_p_p_m1;
    double dlog_c  = - K[i] * pow(p, 2) * p_m1_powK1 / base_const;

    if(calc_ddlog == true){
      double ddlog_c = (- 2 * p + pow(p, 2) * (K[i] - 1) / one_mp - dlog_c) * dlog_c;
      ddlog[i] = - p * dlog[i] - Y_p_t_m1 * theta_inv_m1;
      ddlog[i] -= ddlog_c;
    }
    dlog[i]  -= dlog_c; // constraint correction


  }
  List result;
  result["lik"]   = lik;
  result["dlogl"]  = dlog;
  if(calc_ddlog == true)
    result["ddlogl"] = ddlog;

  return result;
}
