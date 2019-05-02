#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
#include "geometric.h"


// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]

#include <RcppNumerical.h>

typedef Eigen::MatrixXd MapMat;
typedef Eigen::VectorXd MapVec;

class Geomm: public Numer::MFuncGrad
{

  private:
    const MapMat X;
    const MapVec Y;
    const MapVec offset;
  public:
    Geomm(const MapMat x_, const MapVec y_,const MapVec offset_) : X(x_), Y(y_), offset(offset_) {}

    double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
    {
      Rcpp::List res = dd_lgeom(beta,
                                Y,
                                X,
                                offset,
                                false);
      const double f = -Rcpp::as<double>(res["lik"]) ;


      grad.noalias() = -Rcpp::as<Eigen::VectorXd>(res["gradient"]);
      if(isnan(f))
        Rcpp::stop(" f is na");
      return f;
    }

};
//' Function for optimize likelihood with Geometric distribution,
//' using the parameterization E[Y] = 1 + \exp(X\beta + offset)
//' @param beta0  - inital guess
//' @param Y      - (n x 1) the observations >1
//' @param X      - (n x K) covariates
//' @param offset - (n x 1) the offset
//[[Rcpp::export]]
Rcpp::List optim_geomm(  Eigen::VectorXd&  beta0,
                         const  Eigen::VectorXd&  Y,
                         const  Eigen::MatrixXd&  X,
                         const  Eigen::VectorXd&  offset){

  Geomm geomm(X, Y, offset);
  double fopt;
  int status = optim_lbfgs(geomm, beta0, fopt);

  if(status < 0)
    Rcpp::stop("Warning: geomm failed to converge");

  Rcpp::List result;
  result["beta"] = beta0;
  return(result);
}

class Geommc: public Numer::MFuncGrad
{

private:
  const MapMat X;
  const MapVec Y;
  const MapVec K;
  const MapVec offset;
public:
  Geommc(const MapMat x_, const MapVec y_,const MapVec k_,const MapVec offset_) : X(x_), Y(y_), K(k_), offset(offset_) {}

  double f_grad(Numer::Constvec& beta, Numer::Refvec grad)
  {
    Rcpp::List res = dd_lgeomc(beta,
                              Y,
                              X,
                              offset,
                              K,
                              false);
    const double f = -Rcpp::as<double>(res["lik"]) ;


    grad.noalias() = -Rcpp::as<Eigen::VectorXd>(res["gradient"]);
    return f;
  }

};


//' Function for optimize likelihood with constrainted Geometric distribution,
//' using the parameterization E[Y] = 1 + \exp(X\beta + offset) (if unconstrained)
//' However now the observations are constrainted to Y[i] \leq K[i]
//' @param beta0  - inital guess
//' @param Y      - (n x 1) the observations >1
//' @param X      - (n x K) covariates
//' @param K      - (n x 1) upper constraint on Y
//' @param offset - (n x 1) the offset
//[[Rcpp::export]]
Rcpp::List optim_geommc(  Eigen::VectorXd&  beta0,
                         const  Eigen::VectorXd&  Y,
                         const  Eigen::MatrixXd&  X,
                         const  Eigen::VectorXd&  K,
                         const  Eigen::VectorXd&  offset){

  Geommc geommc(X, Y, K, offset);
  double fopt;
  int status = optim_lbfgs(geommc, beta0, fopt, 200);

  if(status < 0)
    Rcpp::stop("Warning: geomm failed to converge");

  Rcpp::List result;
  result["beta"] = beta0;
  return(result);
}

