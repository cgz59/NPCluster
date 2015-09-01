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

// [[Rcpp::export(.vectorizedElementFnLogLik)]]
Rcpp::NumericVector vectorizedElementFnLogLik(SEXP sexp,
							const Rcpp::NumericVector& phi, const double sd, const int num, const double Y, const double Xsd) {
	EnginePtr engine = parsePtr(sexp);
	return engine->vectorizedElementFnLogLik(phi, sd, num, Y, Xsd);
}

// [[Rcpp::export(.computeMarginalLikelihood)]]
Rcpp::List computeMarginalLikelihood(SEXP sexp,
                                     const Rcpp::NumericMatrix& X,
                                     const Rcpp::NumericVector& phi,
                                     const Rcpp::NumericVector& Paux,
                                     const double tau, const double tau0,
                                     const bool sample,
                                     const bool exactBitStream) {
  EnginePtr engine = parsePtr(sexp);
  return engine->computeMarginalLikelihood(X, phi, Paux, tau, tau0, sample, exactBitStream);
}

// [[Rcpp::export(.computeDPAcceptanceRatio)]]
Rcpp::List computeDPAcceptanceRatio(SEXP sexp,
							const Rcpp::NumericVector& Y, const Rcpp::NumericVector& X,
							const Rcpp::IntegerVector& I, const Rcpp::IntegerVector& C,
							const Rcpp::NumericVector& phi,
							const Rcpp::IntegerVector& newS, const Rcpp::IntegerVector& oldS,
							const double tau, const double tau0,
							const int N) {
	EnginePtr engine = parsePtr(sexp);
	return engine->computeDPAcceptanceRatio(Y, X, I, C, phi, newS, oldS, tau, tau0, N);
}

// [[Rcpp::export(.computeColumnPmfAndNeighborhoods)]]
Rcpp::List computeColumnsPmfAndNeighborhoods(SEXP sexp,
                               int n0, const Rcpp::IntegerVector& nVec, double epsilon, double epsilon2,
                               int K, int N,
						  	   const Rcpp::NumericVector& Y, const Rcpp::NumericMatrix& X,
						  	   const Rcpp::NumericMatrix& A, const Rcpp::IntegerMatrix& S,
                               const Rcpp::IntegerVector& rowSubsetI,
                               const Rcpp::IntegerVector& CmVec, const int n2,
                               const Rcpp::NumericVector& phiV, const double tau, const double tau0, const double tauInt,
                               const int maxNeighborhoodSize,
                               const double cutOff, const bool collectMax) {
  EnginePtr engine = parsePtr(sexp);
  return engine->computeColumnPmfAndNeighborhoods(n0, nVec, epsilon, epsilon2,
                                     K, N,
                                     Y, X, A, S, rowSubsetI,
                                     CmVec, n2,
                                     phiV, tau, tau0, tauInt,
                                     maxNeighborhoodSize, cutOff, collectMax);
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

// #define PRINT_DIM
#ifdef PRINT_DIM
  std::cerr << "G = " << G << " and N = " << N << std::endl;
  std::cerr << X.rows() << ":" << X.cols() << std::endl;
  std::cerr << A.rows() << ":" << A.cols() << std::endl;
  std::cerr << S.rows() << ":" << S.cols() << std::endl << std::endl;
#endif // PRINT_DIM

	return engine->computePdpLogLikelihood(k, X, A, S, G, N, tau, tau0, tauInt, colSums);
}

// [[Rcpp::export(.fastIndexedSetCopy)]]
Rcpp::IntegerVector fastIndexedSetCopy(const Rcpp::IntegerVector inY,
					const Rcpp::IntegerVector& indices,
					const Rcpp::IntegerVector& x) {
	using namespace Rcpp;

	IntegerVector y(inY.size());
	auto index = std::begin(indices);
	auto value = std::begin(x);
	const auto end = std::end(indices);

	std::copy(std::begin(inY), std::end(inY), std::begin(y));

	for (; index != end; ++index, ++value) {
		y[*index - 1] = *value;
	}

	return y;
}


// [[Rcpp::export(.fastIndexedSetNoCopy)]]
void fastIndexedSetNoCopy(Rcpp::IntegerVector& y,
					const Rcpp::IntegerVector& indices,
					const Rcpp::IntegerVector& x) {
	using namespace Rcpp;

//	IntegerVector y(inY.size());
	auto index = std::begin(indices);
	auto value = std::begin(x);
	const auto end = std::end(indices);

	for (; index != end; ++index, ++value) {
		y[*index - 1] = *value;
	}

	//return y;
}


// [[Rcpp::export(.fastTabulate)]]
Rcpp::IntegerVector fastTabulate(const Rcpp::IntegerMatrix& mat, const int K,
                                 bool includeZero = false) {
	using namespace Rcpp;
	IntegerVector count(K + 1);  // C uses 0-indexing

	for (auto it = std::begin(mat); it != std::end(mat); ++it) {
		++count[*it];
	}

	if (includeZero) {
	  return count;
	} else {
	  return IntegerVector(std::begin(count) + 1, std::end(count)); // R uses 1-indexing
	}
}

// [[Rcpp::export(.fastTabulateVector)]]
Rcpp::IntegerVector fastTabulateVector(const Rcpp::IntegerVector& vec, const int K,
                                       const bool includeZero) {
  using namespace Rcpp;
  IntegerVector count(K + 1);  // C uses 0-indexing

  for (auto it = std::begin(vec); it != std::end(vec); ++it) {
    ++count[*it];
  }

  if (includeZero) {
    return count;
  } else {
    return IntegerVector(std::begin(count) + 1, std::end(count)); // R uses 1-indexing
  }
}

// [[Rcpp::export(.fastTabulateExcludeEmptiedIndices)]]
Rcpp::IntegerVector fastTabulateExcludeEmptiedIndices(const Rcpp::IntegerMatrix& mat, const Rcpp::IntegerVector& empty,
								                       const int K, bool includeZero = false) {
	using namespace Rcpp;
	IntegerVector count(K + 1);  // C uses 0-indexing

	for (auto it = std::begin(mat); it != std::end(mat); ++it) {
		++count[*it];
	}

	for (auto it = std::begin(empty); it != std::end(empty); ++it) {
		const auto& column = mat(_, *it - 1);
		for (auto itEntry = std::begin(column); itEntry != std::end(column); ++itEntry) {
			--count[*itEntry];
		}
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

// [[Rcpp::export(.fastSumSafeLog)]]
double fastSumSafeLog(const Rcpp::NumericVector& prob,
                      const Rcpp::IntegerVector& count,
                      const int length) {
  double total = 0.0;

  for (int i = 0; i < length; ++i) {
    if (count[i] > 0) {
      total += std::log(prob[i]) * count[i];
    }
  }
  return total;
}

