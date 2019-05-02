#ifndef GEOMETRIC_h
#define GEOMETRIC_h

#include <Rcpp.h>
#include <RcppEigen.h>

Rcpp::List dd_lgeo_cpp_internal(const Eigen::VectorXd&,
                  const Eigen::VectorXd&,
                  bool = true);
Rcpp::List ddlgeo_const_internal(const Eigen::VectorXd&, //Y
                                 const Eigen::VectorXd&, //Theta
                                 const Eigen::VectorXd&, //K
                                 bool = true             // caclc Hessian
                                );
Rcpp::List dd_lgeom(const    Eigen::VectorXd&  , //beta
                  const Eigen::VectorXd&  ,      //Y
                  const Eigen::MatrixXd&   ,     //X
                  const Eigen::VectorXd&  ,      //offset
                  bool              = true      // calcl Hessian
                      );
Rcpp::List dd_lgeomc(const    Eigen::VectorXd&  , //beta
                    const Eigen::VectorXd&  , //Y
                    const Eigen::MatrixXd&   , //X
                    const Eigen::VectorXd&  ,  // offset
                    const Eigen::VectorXd&  ,  // K
                    bool              = true   // calcl Hessian
                       );
#endif
