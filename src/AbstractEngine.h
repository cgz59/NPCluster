#ifndef _ABSTRACTENGINE_HPP
#define _ABSTRACTENGINE_HPP

#include <iostream>
#include <sstream>

#include "Rcpp.h"

#define MATCH_SG

#ifdef NDEBUG
# define verify(EX)
#else
# define verify(EX) (void)((EX) || (__verify #EX, __FILE__, __LINE__),0))
#endif

void __verify(const char *msg, const char *file, int line) {
    char buffer[100];
    snprintf(buffer, 100, "Assert Failure: %s at %s line #%d", msg, file, line);
    throw std::invalid_argument(buffer);
}

namespace np_cluster {

using namespace Rcpp;


class AbstractEngine {
public:

	typedef std::vector<int> StdIntVector;


  AbstractEngine(bool sort = true) : extraSorting(sort) { }
  virtual ~AbstractEngine() { }

  // virtual void computePmfAndNeighborhoods() = 0; // pure virtual


  Rcpp::List computePdpLogLikelihood(const int k, const Rcpp::NumericMatrix& X,
							const Rcpp::NumericMatrix& A, const Rcpp::IntegerMatrix& S,
							const int G, const int N,
							const double tau, const double tau0, const double tauInt, bool colSums) {
	using namespace Rcpp;

	if (colSums) {
		stop("colSums = TRUE is not yet implemented.");
	}

	// TODO X may be constant, so may consider making single, padded copy for vectorization

	NumericVector logLikelihood(G + 1);
	internalComputePdpLogLikelihood(&logLikelihood[0], k, X, A, S, G, N, tau, tau0, tauInt, colSums);

	return Rcpp::List::create(
		Rcpp::Named("logLikelihood") = logLikelihood
	);
  }

  template <typename OutputPtr>
  void internalComputePdpLogLikelihood(OutputPtr logLikelihood, const int k, const Rcpp::NumericMatrix& X,
							const Rcpp::NumericMatrix& A, const Rcpp::IntegerMatrix& S,
							const int G, const int N,
							const double tau, const double tau0, const double tauInt, bool colSums) {

	// TODO X may be constant, so may consider making single, padded copy for vectorization
	const auto Xmt = std::begin(X) + (k - 1) * N;

	// gg == 0 case

	const auto tmp = std::accumulate(
		Xmt, Xmt + N, 0.0,
		[](double sum, double x) {
			return (sum + (x - 1.0) * (x - 1.0));
		}
	);
	logLikelihood[0] = -0.5 * tmp / (tauInt * tauInt) - N * 0.5 * std::log(2.0 * M_PI)
											          - N * std::log(tauInt);
	// gg > 0 cases
	for (int gg = 1; gg <= G; ++gg) { // TODO Parallelize
		auto Xmtgg = Xmt;
		auto Agg = std::begin(A) + (gg - 1) * N;
		auto Sgg = std::begin(S) + (gg - 1) * N;

		const auto   occupied = -0.5 * std::log(2.0 * M_PI) - std::log(tau);
		const auto unoccupied = -0.5 * std::log(2.0 * M_PI) - std::log(tau0);

		auto tmp = 0.0;
		for (int i = 0; i < N; ++i) { // TODO Could use std::accumulate with zip(Agg, Sgg)
			if (*Sgg == 0) {
				tmp += -0.5 * (*Xmtgg) * (*Xmtgg) / (tau0 * tau0) + unoccupied;
			} else {
				tmp += -0.5 * (*Xmtgg - *Agg) * (*Xmtgg - *Agg) / (tau * tau) + occupied;
			}
			++Xmtgg;
			++Agg;
			++Sgg;
		}
		logLikelihood[gg] = tmp;
	}
  }

  double logLikelihood(double mean, double sd, int num, double Y, double Xsd) {
  	return -num / 2.0 * std::log(2.0 * M_PI * sd * sd)
  		- 0.5 * num * (Y - mean) * (Y - mean) / (sd * sd)
  		- 0.5 * (num - 1) * Xsd * Xsd / (sd * sd);
  	// TODO Write in terms of precisions
  	// TODO Ignore log(2pi) (?)
  }

  Rcpp::List computePmfAndNeighborhoods(int n0, const IntegerVector& nVec,
  		const double epsilon, const double epsilon2,
  		int K, int N,
  		const NumericVector& Y, const NumericVector& Xsd, const IntegerVector& rowSubsetI,
  		const IntegerVector& CmVec, const int n2,
  		const NumericVector& phiV, const double tau, const double tau0,
  		const int maxNeighborhoodSize,
  		const double cutOff) {

	// TODO n2 is const throughout run

	verify(K == nVec.size());
	verify(K == phiV.size());

    std::vector<double> priorProbV(K + 1); // TODO Reuse memory
    priorProbV[0] = n0;
    std::copy(std::begin(nVec), std::end(nVec), std::begin(priorProbV) + 1);
    std::transform(std::begin(priorProbV), std::end(priorProbV), std::begin(priorProbV),
    	[epsilon](double x) {
    		if (x < epsilon) {
    			x = epsilon;
    		}
    		return std::log(x);
    	}
    );
    std::vector<double> phiVKp1(K + 1);
    phiVKp1[0] = 0.0;
    std::copy(std::begin(phiV), std::end(phiV), std::begin(phiVKp1) + 1);

    std::vector<double> tauV(K + 1, tau); // TODO Remove, without switch in for-loop below
    tauV[0] = tau0;

	if ((K + 1) * N > logSsMt.size()) { // TODO Handle alignment, minor dimension == K
		logSsMt.resize((K + 1) * N);
	}

	// TODO Parallelize
	std::for_each(std::begin(rowSubsetI), std::end(rowSubsetI),
		[&phiVKp1,&tauV,&CmVec,&Y,&Xsd,&priorProbV,tau,tau0,n2,K,epsilon2,this](const int x) {
			const auto row = x - 1; //rowSubsetI[i] - 1;
			const auto gv = row / n2;

			// perform transformation and accumulate max entry
			auto max = -std::numeric_limits<double>::max();
			for (int j = 0; j < K + 1; ++j) {
				const auto entry =
					logLikelihood(phiVKp1[j], tauV[j], CmVec[gv], Y[row], Xsd[row])
						+ priorProbV[j];
				if (entry > max) {
					max = entry;
				}
				// Store
				logSsMt[row * (K + 1) + j] = entry;
			}

			// subtract off max, exponentiate and accumulate colSum1
			auto colSum1 = 0.0;
			for (int j = 0; j < K + 1; ++j) {
				const auto entry = std::exp(logSsMt[row * (K + 1) + j] - max);
				colSum1 += entry;
				logSsMt[row * (K + 1) + j] = entry;
			}

			// first normalize, replace zeros by "small" and accumulate colSum2
			auto colSum2 = 0.0;
			for (int j = 0; j < K + 1; ++j) {
				auto entry = logSsMt[row * (K + 1) + j] / colSum1;
				if (entry < epsilon2) {
					entry = epsilon2;
				}
				colSum2 += entry;
				logSsMt[row * (K + 1) + j] = entry;
			}

			// again normalize TODO Is second normalization really necessary?
			for (int j = 0; j < K + 1; ++j) {
				logSsMt[row * (K + 1) + j] /= colSum2;
			}
		}
	);

// 	for (int j = 0; j < K + 1; ++j) {
// 		std::cerr << " " << logSsMt[(rowSubsetI[0] - 1) * (K + 1) + j];
// 	}
// 	std::cerr << std::endl;
//
// 	double total = std::accumulate(std::begin(logSsMt), std::end(logSsMt), 0.0);
// 	std::cerr << total << std::endl;

	// now compute the delta-neighborhoods

	StdIntVector neighborhoodList;
	StdIntVector neighborhoodOffset;
	StdIntVector neighborhoodIndex;

	computeDeltaNeighborhoods(rowSubsetI,
		neighborhoodList, neighborhoodOffset, neighborhoodIndex,
		cutOff, maxNeighborhoodSize, K);

	return Rcpp::List::create(
		Rcpp::Named("neighbor") = neighborhoodList,
		Rcpp::Named("offset") = neighborhoodOffset,
		Rcpp::Named("index") = neighborhoodIndex
	);

  }

  Rcpp::List computeColumnPmfAndNeighborhoods(int n0, const IntegerVector& nVec,
  		const double epsilon, const double epsilon2,
  		int K, int N,  // NB: K -> G for column-space
  		const NumericVector& Y, const NumericMatrix& X,
  		const NumericMatrix& A, const IntegerMatrix& S,
  		const IntegerVector& rowSubsetI,
  		const IntegerVector& CmVec, const int P,
  		const NumericVector& phiV, const double tau, const double tau0, const double tauInt,
  		const int maxNeighborhoodSize,
  		const double cutOff) {

	// TODO n2 is const throughout run

	verify(K == nVec.size());
	verify(K == phiV.size());

    std::vector<double> priorProbV(K + 1); // TODO Reuse memory

    priorProbV[0] = n0;
    std::copy(std::begin(nVec), std::end(nVec), std::begin(priorProbV) + 1);
    std::transform(std::begin(priorProbV), std::end(priorProbV), std::begin(priorProbV),
    	[epsilon](double x) {
    		if (x < epsilon) {
    			x = epsilon;
    		}
    		return std::log(x);
    	}
    );

// 	for (int i = 0; i < K + 1; ++i) {
// 		std::cerr << " " << priorProbV[i];
// 	}
// 	std::cerr << std::endl;

//     std::vector<double> phiVKp1(K + 1);
//     phiVKp1[0] = 0.0;
//     std::copy(std::begin(phiV), std::end(phiV), std::begin(phiVKp1) + 1);
//
//     std::vector<double> tauV(K + 1, tau); // TODO Remove, without switch in for-loop below
//     tauV[0] = tau0;
//
	if ((K + 1) * N > logSsMt.size()) { // TODO Handle alignment, minor dimension == K
		logSsMt.resize((K + 1) * N);
	}

	// TODO Parallelize
	std::for_each(std::begin(rowSubsetI), std::end(rowSubsetI),
		[&CmVec,&Y,&X,&A,&S,&priorProbV,&rowSubsetI,tau,tau0,tauInt,P,K,epsilon2,N,this](const int x) {
			const auto col = x - 1; //rowSubsetI[i] - 1;
			// const auto gv = col / n2;

			internalComputePdpLogLikelihood(&logSsMt[col * (K + 1)], col + 1, X, A, S, K, N, tau, tau0, tauInt, false);

			// add prior and accumulate max entry
			auto max = -std::numeric_limits<double>::max();
			for (int j = 0; j < K + 1; ++j) {
				const auto entry =
					logSsMt[col * (K + 1) + j] // TODO Could compute here instead of above, one fewer write/read pair
						+ priorProbV[j];
				if (entry > max) {
					max = entry;
				}
				// Store
				logSsMt[col * (K + 1) + j] = entry;
			}

			// subtract off max, exponentiate and accumulate colSum1
			auto colSum1 = 0.0;
			for (int j = 0; j < K + 1; ++j) {
				const auto entry = std::exp(logSsMt[col * (K + 1) + j] - max);
				colSum1 += entry;
				logSsMt[col * (K + 1) + j] = entry;
			}

			// first normalize, replace zeros by "small" and accumulate colSum2
			auto colSum2 = 0.0;
			for (int j = 0; j < K + 1; ++j) {
				auto entry = logSsMt[col * (K + 1) + j] / colSum1;
				if (entry < epsilon2) {
					entry = epsilon2;
				}
				colSum2 += entry;
				logSsMt[col * (K + 1) + j] = entry;
			}

			// again normalize TODO Is second normalization really necessary?
			for (int j = 0; j < K + 1; ++j) {
				logSsMt[col * (K + 1) + j] /= colSum2;
			}

// 			if (col == 1) {
// 				for (int j = 0; j < K + 1; ++j) {
// 					std::cerr << " " << logSsMt[col * (K + 1) + j];
// 				}
// 				std::cerr << std::endl;
// 			}
		}
	);

	// now compute the delta-neighborhoods
	StdIntVector neighborhoodList;
	StdIntVector neighborhoodOffset;
	StdIntVector neighborhoodIndex;

	computeDeltaNeighborhoods(rowSubsetI,
		neighborhoodList, neighborhoodOffset, neighborhoodIndex,
		cutOff, maxNeighborhoodSize, K);

	return Rcpp::List::create(
		Rcpp::Named("neighbor") = neighborhoodList,
		Rcpp::Named("offset") = neighborhoodOffset,
		Rcpp::Named("index") = neighborhoodIndex
	);

  }


private:

	typedef std::pair<double, int> Score;
	typedef std::vector<Score> ScoreList;
	typedef typename ScoreList::iterator ScoreIterator;

	typedef std::set<int> NeighborhoodContainer;
	typedef ScoreList FastContainer;
	typedef NeighborhoodContainer::iterator NeighborhoodIterator;
	typedef typename FastContainer::iterator FastIterator;




	void computeDeltaNeighborhoods(const IntegerVector& rowSubsetI,
		StdIntVector& list, StdIntVector& offset, StdIntVector& index,
		const double cutOff,
		const int maxSize, const int K) {

	  	FastContainer I; // TODO Better to reuse?
	  	fillInitialList(I, rowSubsetI);

		list.clear(); offset.clear(); index.clear();

		FastIterator beginI = std::begin(I);
		const FastIterator endI = std::end(I);

	  	while (beginI != endI) {
			drawNextNeighborhood(beginI, endI, list, offset, index,
				cutOff, maxSize, K);
	  	}
 	  	offset.push_back(list.size() + 1); // Inclusive counting
	}

	void drawNextNeighborhood(FastIterator& begin, const FastIterator& end,
							  StdIntVector& list, StdIntVector& offset, StdIntVector& index,
							  double cutOff,
							  const int maxSize, const int K) {
		offset.push_back(list.size() + 1); // Add next offset to a neighborhood
		if (std::distance(begin, end) == 1) {
			// Singleton
			list.push_back(begin->second);
			index.push_back(begin->second);
			++begin; // advance
		} else {
			// 2 or more left
			auto k = sampleUniform(begin, end);

			// Score each remaining entry in [begin, end) against k
			std::transform(begin, end, begin, // TODO Parallelize
				[k,K,this](Score score) {
					const auto measure = distributionDistance(k.second - 1, score.second - 1, K);
					return std::make_pair(measure, score.second);
				}
			);
// 			std::cerr << "k = " << k.second << std::endl;

//			return;

// 			std::cerr << "Current length: " << std::distance(begin,end) << std::endl;

			// sort the first maxNeighborhoodSize elements in increasing measure
			auto lastSort = (std::distance(begin, end) > maxSize) ? begin + maxSize : end;

// 			std::cerr << "Partial length: " << std::distance(begin,lastSort) << std::endl;

// 			return;

			std::partial_sort(begin, lastSort, end);

			// find cut-off point
			auto lastCut = std::find_if(begin, lastSort, [cutOff](Score score) {
				return score.first > cutOff;
			});

// 			std::for_each(begin,lastCut, [](Score score) {
// 				std::cerr << " (" << score.first << ":" << score.second << ")";
// 			});
// 			std::cerr << std::endl;

// 			std::cerr << "Cut length: " << std::distance(begin,lastCut) << std::endl;
// 			std::cerr << "cutOff = " << cutOff << std::endl;

// 			return;


			// elements [begin, lastCut) <= cutOff .. and that range always has size >= 1, (k matches with k)
			if (begin == lastCut) {
				Rcpp::stop("Bad neighborhood");
			}

			auto cutLength = std::distance(begin, lastCut);
			for( ; begin != lastCut; ++begin) {
				list.push_back(begin->second);
			}
			index.push_back(k.second);

			if (extraSorting) {
				std::sort(std::end(list) - cutLength, std::end(list)); // TODO Ask SG: Do these need to be sorted?

				// TODO Ask SG: Why sort I.k again? Seems completely unnecessary
				std::sort(begin, end, [](const Score x, const Score y) {
					return x.second < y.second;
				});
			}
		}

// 		std::for_each(begin, begin + 10, [](Score score) {
// 			std::cerr << " (" << score.first << ":" << score.second << ")";
// 		});
// 		std::cerr << std::endl;
	}


	double distributionDistance(const int i, const int j, const int K) {
		return 2.0 * (1.0 -
				std::inner_product(
					std::begin(logSsMt) + i * (K + 1),
					std::begin(logSsMt) + (i + 1) * (K + 1),
					std::begin(logSsMt) + j * (K + 1), 0.0,
					std::plus<double>(),
					//std::product<double>()
					[](double x, double y) {
						return std::sqrt(x * y);
					}
					));
// 		return std::sqrt(
// 			std::inner_product(
//
// 			)
// 		);
	}

	void fillInitialList(FastContainer& container, const IntegerVector& rowSubsetI) {
		container.reserve(rowSubsetI.size());
		for (const auto row : rowSubsetI) {
			container.push_back(std::make_pair(0.0, row));
		}
	}

	template <typename Itr>
	auto sampleUniform(Itr begin, Itr end) -> decltype(*begin) {
		return *(begin + std::distance(begin, end) * unif_rand());
	}


protected:

	std::vector<double> logSsMt;

	bool extraSorting;

};

template <typename RealType>
class CPUEngine : public AbstractEngine {
public:
  CPUEngine(bool sort) : AbstractEngine(sort) { }
  virtual ~CPUEngine() { }


//   void computePmfAndNeighborhoods() {
//     std::cout << "Here" << std::endl;
//   }

private:

  // computePmfAndNeighborhoods

  // element_fn.post.prob.and.delta


	std::vector<RealType> t;
};

} // namespace np_cluster


#endif // _ABSTRACTENGINE_HPP
