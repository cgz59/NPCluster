/*
 * @author Marc A. Suchard
 */

#include "Rcpp.h"
#include "RcppEigen.h"

#include "AbstractEngine.h"

typedef Rcpp::XPtr<np_cluster::AbstractEngine> EnginePtr;

EnginePtr parsePtr(SEXP sexp) {
	EnginePtr ptr(sexp);
	if (!ptr) {
		Rcpp::stop("External pointer is uninitialized");
	}
	return ptr;
}

 
// [[Rcpp::export(.createEngine)]]
Rcpp::List createEngine(bool sort) {

	EnginePtr engine(
		new np_cluster::CPUEngine<double>(sort)
	);

	Rcpp::List list = Rcpp::List::create(
		Rcpp::Named("engine") = engine
	);

	return list;
}

// [[Rcpp::export(.accessEngine)]]
int accessEngine(SEXP sexp) {
  EnginePtr engine = parsePtr(sexp);
  return 1;
}
// [[Rcpp::export(.computePmfAndNeighborhoods)]]
Rcpp::List computePmfAndNeighborhoods(SEXP sexp,
                               int n0, const Rcpp::IntegerVector& nVec, double epsilon, double epsilon2,
                               int K, int N,
                               const Rcpp::NumericVector& Y, const Rcpp::NumericVector& Xsd, const Rcpp::IntegerVector& rowSubsetI,
                               const Rcpp::IntegerVector& CmVec, const int n2,
                               const Rcpp::NumericVector& phiV, const double tau, const double tau0, const int maxNeighborhoodSize,
                               const double cutOff) {
  EnginePtr engine = parsePtr(sexp);
  return engine->computePmfAndNeighborhoods(n0, nVec, epsilon, epsilon2,
                                     K, N,
                                     Y, Xsd, rowSubsetI,
                                     CmVec, n2,
                                     phiV, tau, tau0,
                                     maxNeighborhoodSize, cutOff);
}

// [[Rcpp::export(.computePdpLogLikelihood)]]
Rcpp::List computePdpLogLikelihood(SEXP sexp,
							const int k, const Rcpp::NumericMatrix& X,
							const Rcpp::NumericMatrix& A, const Rcpp::IntegerMatrix& S,
							const int G, const int N,
							const double tau, const double tau0, const double tauInt, bool colSums) {
	EnginePtr engine = parsePtr(sexp);

	return engine->computePdpLogLikelihood(k, X, A, S, G, N, tau, tau0, tauInt, colSums);
}

// [[Rcpp::export(.fastTabulate)]]
Rcpp::IntegerVector fastTabulate(const Rcpp::IntegerMatrix& mat, const int K) {
	using namespace Rcpp;
	IntegerVector count(K + 1);  // C uses 0-indexing

	for (auto it = std::begin(mat); it != std::end(mat); ++it) {
		count[*it]++;
	}

	return IntegerVector(std::begin(count) + 1, std::end(count)); // R uses 1-indexing
}

// [[Rcpp::export(.fastXtX)]]
Rcpp::NumericMatrix fastXtX(const Rcpp::NumericMatrix& rX) {
	using namespace Eigen;
  using namespace Rcpp;
	const Map<MatrixXd> X(as<Map<MatrixXd> >(rX));

	return wrap(X.transpose() * X);
}
