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

// [[Rcpp::export(.fastTabulateVector)]]
Rcpp::IntegerVector fastTabulateVector(const Rcpp::IntegerVector& vec, const int K,
                                       const bool includeZero) {
  using namespace Rcpp;
  IntegerVector count(K + 1);  // C uses 0-indexing

  for (auto it = std::begin(vec); it != std::end(vec); ++it) {
    count[*it]++;
  }

  if (includeZero) {
    return count;
  } else {
    return IntegerVector(std::begin(count) + 1, std::end(count)); // R uses 1-indexing
  }
}

// [[Rcpp::export(.fastXtX)]]
Rcpp::NumericMatrix fastXtX(const Rcpp::NumericMatrix& rX) {
	using namespace Eigen;
  using namespace Rcpp;
	const Map<MatrixXd> X(as<Map<MatrixXd> >(rX));

	return wrap(X.transpose() * X);
}


// for (gg in 0:parm$clust$K)
// {flag.gg <- (new.s.k==gg)
//   count.gg <- sum(flag.gg)
//
//   if (gg > 0)
//   {parm$clust$n.vec[gg] <- parm$clust$n.vec.k.comp[gg] + count.gg
//   }
//
//   if (gg == 0)
//   {parm$clust$n0 <- parm$clust$n0.k.comp + count.gg
//   }
//
//   if (count.gg > 0)
//   {rho.prop <- rho.prop + log(parm$clust$post.k[gg+1])*count.gg
//   }
// }


// [[Rcpp::export(.fastPrior)]]
Rcpp::List fastPrior(const Rcpp::IntegerVector& Sk,
                     const Rcpp::IntegerVector& Nk,
                     const int N0,
                     const Rcpp::NumericVector& Pk,
                     const int K) {
  using namespace Rcpp;

  IntegerVector count(K + 1);  // C uses 0-indexing

  double rho = 0.0;

  return List::create(
    _["count"] = count,
    _["rho"] = rho
  );
}

// [[Rcpp::export(.fastSumSafeLog)]]
double fastSumSafeLog(const Rcpp::NumericVector& prob,
                      const Rcpp::IntegerVector& count) {
  double total = 0.0;
  for (int i = 0; i < prob.length(); ++i) {
    if (count[i] > 0) {
      total += std::log(prob[i]) * count[i];
    }
  }
  return total;
}

