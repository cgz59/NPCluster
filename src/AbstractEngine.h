#ifndef _ABSTRACTENGINE_HPP
#define _ABSTRACTENGINE_HPP

#include <iostream>
#include <sstream>

#include "Rmath.h"
#include "Rcpp.h"

#define USE_NT2

#ifdef USE_NT2
  #include "RcppNT2.h"
#endif

#include "tbb/parallel_for.h"

#include "boost/format.hpp"

#include "aligned_allocator.hpp"

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

class Tick {
public:
  Tick() : start(std::chrono::steady_clock::now()) { }

  virtual ~Tick() { }

  long operator()() {
    auto end = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return diff;
  }

private:
  std::chrono::time_point<std::chrono::steady_clock> start;
};


static void ProbSampleReplace(int n, double *p, int *perm, int nans, int *ans)
{
  double rU;
  int i, j;
  int nm1 = n - 1;

  /* record element identities */
  for (i = 0; i < n; i++)
    perm[i] = i + 1;

  /* sort the probabilities into descending order */
  revsort(p, perm, n);

  /* compute cumulative probabilities */
  for (i = 1 ; i < n; i++)
    p[i] += p[i - 1];

  /* compute the sample */
  for (i = 0; i < nans; i++) {
    rU = unif_rand();
    for (j = 0; j < nm1; j++) {
      if (rU <= p[j])
        break;
    }
    ans[i] = perm[j];
  }
}

namespace np_cluster {

using namespace Rcpp;


enum class SpecialComputeMode {
  TBB = 2,
  SSE = 4,
  AVX = 8,
};

class AbstractEngine {
public:

  typedef std::vector<int> StdIntVector;
  typedef std::vector<double> StdNumericVector;

  AbstractEngine(bool sort = true, int specialMode = 0) :
    extraSorting(sort),
    specialMode(specialMode),
    useTBB(specialMode & static_cast<long>(SpecialComputeMode::TBB)),
    useSSE(specialMode & static_cast<long>(SpecialComputeMode::SSE)) {

    Rcpp::Rcout << "AbstractEngine:\n";
    Rcpp::Rcout << "\tTBB: " << (useTBB ? "true" : "false") << "\n";
    Rcpp::Rcout << "\tSSE: " << (useSSE ? "true" : "false") << std::endl;
  }

  virtual ~AbstractEngine() { }

  void printTiming() {
    boost::format format("%1.3e");
    for (auto& d : duration) {
      Rcpp::Rcout << d.first << " " << (format % static_cast<double>(d.second)) << std::endl;
    }
  }

  // virtual void computePmfAndNeighborhoods() = 0; // pure virtual


  Rcpp::NumericVector vectorizedElementFnLogLik(const Rcpp::NumericVector& phi,
                                                const double sd, const int num, const double Y, const double Xsd) {
    using namespace Rcpp;

    NumericVector logLike(phi.size());
    std::transform(std::begin(phi), std::end(phi), std::begin(logLike),
                   [this,sd,num,Y,Xsd](const double mean) {
                     return logLikelihood(mean, sd, num, Y, Xsd);
                   });

    return logLike;
  }


  Rcpp::List computeDPAcceptanceRatio(const Rcpp::NumericVector& Y, const Rcpp::NumericVector& X,
                                      const Rcpp::IntegerVector& I, const Rcpp::IntegerVector& C,
                                      const Rcpp::NumericVector& phiR,
                                      const Rcpp::IntegerVector& newS, const Rcpp::IntegerVector& oldS,
                                      const double tau, const double tau0,
                                      const int N) {

    // 	double ratio = 0.0;

    // TODO Make local copy to facilitate vectorization
    // 	std::vector<double> phi(phiR.length() + 1); // filled with zero
    // 	std::copy(std::begin(phiR). std::end(phiR), std::begin(phi) + 1);
    // TODO Make local copy of tau to facilitate vectorization

    // 	ratio += std::accumulate(
    // 		std::begin(newS), std::end(newS),
    // 		[](const int sk) {
    // 			--sk; // R is 1-indexed
    // 			return logLikelihood(
    // 				(sk == -1) ? 0.0 : phiR[sk],
    // 				(sk == -1) ? tau0 : tau,
    //
    // 			);
    // 		}
    //   	);

    double newRatio = 0.0;
    double oldRatio = 0.0;

    for (int k = 0; k < I.length(); ++k) {
      const auto Ik = I[k] - 1; // R is 1-indexed
      const auto newSk = newS[k];
      const auto oldSk = oldS[k];
      const auto gk = Ik / N;
      const auto numK = C[gk];

      newRatio += logLikelihood(
        (newSk == 0) ? 0.0 : phiR[newSk - 1],
                                 (newSk == 0) ? tau0 : tau,
                                 numK,
                                 Y[Ik],
                                  X[Ik]
      );

      oldRatio += logLikelihood(
        (oldSk == 0) ? 0.0 : phiR[oldSk - 1],
                                 (oldSk == 0) ? tau0 : tau,
                                 numK,
                                 Y[Ik],
                                  X[Ik]
      );
    }

    double ratio = newRatio - oldRatio;



    return Rcpp::List::create(
      Rcpp::Named("ratio") = ratio,
      Rcpp::Named("new") = newRatio,
      Rcpp::Named("old") = oldRatio
    );

  }


  Rcpp::List computeMarginalLikelihood(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericVector& phi,
                                       const Rcpp::NumericVector& Paux,
                                       const double tau, const double tau0,
                                       const bool doSample,
                                       const bool exactBitStream) {
    using namespace Rcpp;

    // TODO if doSample == false, then no need to save to posterior[]

    const int rows = X.rows();
    //     const int cols = X.cols();
    const int len = phi.length();


    std::vector<int> sample(rows); // TODO Only allocate if doSample == true
    //     std::vector<double> logLikelihood(rows);

    std::vector<int> perm(len + 1);
    std::vector<double> posterior(len + 1);

    double logMarginalLikelihood = 0.0;
    for (int i = 0; i < rows; ++i) {

      const auto xMtTt = X[i];

      // Handle 0 entry // TODO For vectorization, could turn tau, phi into len + 1 vectors
      double sum = posterior[0] = Rf_dnorm4(xMtTt, 0.0, tau0, 0) * Paux[0];

      // Handle remaining entries
      for (int j = 1; j <= len; ++j) {
        sum += posterior[j] = Rf_dnorm4(xMtTt, phi[j - 1], tau, 0) * Paux[j];
      }


      if (doSample) {

#if 1
        for (int j = 0; j <=len; ++j) {
          posterior[j] /= sum;
        }

        int index;
        ProbSampleReplace(len + 1, &posterior[0], &perm[0], 1, &index); // explicit R code
        --index;
#else
        // Normalize posterior and sample
        double draw = unif_rand() * sum;
        int index = 0;
        while (draw >= posterior[index] && index < len) {
          draw -= posterior[index];
          ++index;
        }
#endif

        // Store results
        sample[i] = index;
      }

      // Store results
      logMarginalLikelihood += std::log(sum);
    }


    // 	double logMarginalLikelihood = std::accumulate(
    // 		boost::make_zip(std::begin(X), std::begin(sample)),
    // 		boost::make_zip(std::end(X), std::end(sample)),
    // 		0.0,
    // 		[](std::pair<double,int>&
    // 	);

    return Rcpp::List::create(
      Rcpp::Named("sVk") = sample,
      Rcpp::Named("logMarginalLikelihood") = logMarginalLikelihood,
      Rcpp::Named("posterior") = posterior
    );
  }


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

    auto f = [Xmt,&A,&S,N,tau,tau0,&logLikelihood](const int gg) {
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
    };

    Tick iCPL;

    if (useTBB) {
      tbb::parallel_for(tbb::blocked_range<int>(1, G + 1),
                        [=](const tbb::blocked_range<int>& r) {
                          const auto end = r.end();
                          for (auto i = r.begin(); i != end; ++i) {
                            f(i);
                          }
                        });
    } else {
      for (int gg = 1; gg <= G; ++gg) {
        f(gg);
      }
    }

    // gg > 0 cases
    // for (int gg = 1; gg <= G; ++gg) { // TODO Parallelize
      // auto Xmtgg = Xmt;
      // auto Agg = std::begin(A) + (gg - 1) * N;
      // auto Sgg = std::begin(S) + (gg - 1) * N;
      //
      // const auto   occupied = -0.5 * std::log(2.0 * M_PI) - std::log(tau);
      // const auto unoccupied = -0.5 * std::log(2.0 * M_PI) - std::log(tau0);
      //
      // auto tmp = 0.0;
      // for (int i = 0; i < N; ++i) { // TODO Could use std::accumulate with zip(Agg, Sgg)
      //   if (*Sgg == 0) {
      //     tmp += -0.5 * (*Xmtgg) * (*Xmtgg) / (tau0 * tau0) + unoccupied;
      //   } else {
      //     tmp += -0.5 * (*Xmtgg - *Agg) * (*Xmtgg - *Agg) / (tau * tau) + occupied;
      //   }
      //   ++Xmtgg;
      //   ++Agg;
      //   ++Sgg;
      // }
      // logLikelihood[gg] = tmp;
    // }

    duration["iCPL____"] += iCPL();
  }

  //   static inline double dnorm(double x, double mean, double sd) {
  // 	  return -0.5 * std::log(2.0 * M_PI * sd * sd)
  // 	        - 0.5 * (x - mean) * (x - mean) / (sd * sd);
  //   }

  static inline double logLikelihood(double mean, double sd, int num, double Y, double Xsd) {
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

    auto stride = getAlignedStride(K + 1);

    if (stride * N > logSsMt.size()) { // TODO Handle alignment, minor dimension == K
      logSsMt.resize(stride * N);
    }

    auto f = [&phiVKp1,&tauV,&CmVec,&Y,&Xsd,&priorProbV,tau,tau0,n2,K,stride,epsilon2,this](const int x) {
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
        logSsMt[row * stride + j] = entry;
      }

      // subtract off max, exponentiate and accumulate colSum1
      auto colSum1 = 0.0;
      for (int j = 0; j < K + 1; ++j) {
        const auto entry = std::exp(logSsMt[row * stride + j] - max);
        colSum1 += entry;
        logSsMt[row * stride + j] = entry;
      }

      // first normalize, replace zeros by "small" and accumulate colSum2
      auto colSum2 = 0.0;
      for (int j = 0; j < K + 1; ++j) {
        auto entry = logSsMt[row * stride + j] / colSum1;
        if (entry < epsilon2) {
          entry = epsilon2;
        }
        colSum2 += entry;
        logSsMt[row * stride + j] = entry;
      }

      // again normalize TODO Is second normalization really necessary?
      for (int j = 0; j < K + 1; ++j) {
        logSsMt[row * stride + j] /= colSum2;
      }
    };

    Tick cPAN;

    if (useTBB) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, rowSubsetI.size()),
                        [=](const tbb::blocked_range<size_t>& r) {
                          std::for_each(std::begin(rowSubsetI) + r.begin(), std::begin(rowSubsetI) + r.end(), f);
                        });
    } else {
      std::for_each(std::begin(rowSubsetI), std::end(rowSubsetI), f);
    }

    duration["cPAN____"] += cPAN();

    // now compute the delta-neighborhoods

    StdIntVector neighborhoodList;
    StdIntVector neighborhoodOffset;
    StdIntVector neighborhoodIndex;

    computeDeltaNeighborhoods(rowSubsetI,
                              neighborhoodList, neighborhoodOffset, neighborhoodIndex,
                              cutOff, maxNeighborhoodSize, K, stride, false);

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
                                              const double cutOff,
                                              const bool collectMax) {

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

    auto stride = getAlignedStride(K + 1);

    if (stride * N > logSsMt.size()) {
      logSsMt.resize(stride * N);
    }

    Tick cCPAN;

    // TODO Parallelize
    std::for_each(std::begin(rowSubsetI), std::end(rowSubsetI),
                  [&CmVec,&Y,&X,&A,&S,&priorProbV,&rowSubsetI,tau,tau0,tauInt,P,K,stride,epsilon2,N,this](const int x) {
                    const auto col = x - 1; //rowSubsetI[i] - 1;
                    // const auto gv = col / n2;

                    internalComputePdpLogLikelihood(&logSsMt[col * stride], col + 1, X, A, S, K, N, tau, tau0, tauInt, false);

                    // add prior and accumulate max entry
                    auto max = -std::numeric_limits<double>::max();
                    for (int j = 0; j < K + 1; ++j) {
                      const auto entry =
                        logSsMt[col * stride + j] // TODO Could compute here instead of above, one fewer write/read pair
                      + priorProbV[j];
                      if (entry > max) {
                        max = entry;
                      }
                      // Store
                      logSsMt[col * stride + j] = entry;
                    }

                    // subtract off max, exponentiate and accumulate colSum1
                    auto colSum1 = 0.0;
                    for (int j = 0; j < K + 1; ++j) {
                      const auto entry = std::exp(logSsMt[col * stride + j] - max);
                      colSum1 += entry;
                      logSsMt[col * stride + j] = entry;
                    }

                    // first normalize, replace zeros by "small" and accumulate colSum2
                    auto colSum2 = 0.0;
                    for (int j = 0; j < K + 1; ++j) {
                      auto entry = logSsMt[col * stride + j] / colSum1;
                      if (entry < epsilon2) {
                        entry = epsilon2;
                      }
                      colSum2 += entry;
                      logSsMt[col * stride + j] = entry;
                    }

                    // again normalize TODO Is second normalization really necessary?
                    for (int j = 0; j < K + 1; ++j) {
                      logSsMt[col * stride + j] /= colSum2;
                    }
                  }
    );

    duration["cCPAN___"] += cCPAN();

    // now compute the delta-neighborhoods
    StdIntVector neighborhoodList;
    StdIntVector neighborhoodOffset;
    StdIntVector neighborhoodIndex;

    const auto max = computeDeltaNeighborhoods(rowSubsetI,
                                               neighborhoodList, neighborhoodOffset, neighborhoodIndex,
                                               cutOff, maxNeighborhoodSize, K, stride, collectMax);

    return Rcpp::List::create(
      Rcpp::Named("neighbor") = neighborhoodList,
      Rcpp::Named("offset") = neighborhoodOffset,
      Rcpp::Named("index") = neighborhoodIndex,
      Rcpp::Named("neighborhoodMax") = max
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

  size_t getAlignedStride(size_t length) {
    const int multiple = (useSSE ? 2 : 1);

    assert(multiple && ((multiple & (multiple -1)) == 0));
    return (length + multiple - 1) & ~(multiple - 1);
  }

  double computeDeltaNeighborhoods(const IntegerVector& rowSubsetI,
                                   StdIntVector& list, StdIntVector& offset, StdIntVector& index,
                                   const double cutOff,
                                   const int maxSize, const int K, const int stride,
                                   const bool returnMax) {

    FastContainer I; // TODO Better to reuse?
    fillInitialList(I, rowSubsetI);

    list.clear(); offset.clear(); index.clear();

    StdNumericVector maxList; // TODO Better to reuse?

    FastIterator beginI = std::begin(I);
    const FastIterator endI = std::end(I);

    while (beginI != endI) {
      drawNextNeighborhood(beginI, endI, list, offset, index,
                           cutOff, maxSize, K, stride, (returnMax ? &maxList : nullptr));
    }
    offset.push_back(list.size() + 1); // Inclusive counting

    if (returnMax) {
      return *std::max_element(std::begin(maxList), std::end(maxList));
    } else {
      return 0.0;
    }
  }

  void drawNextNeighborhood(FastIterator& begin, const FastIterator& end,
                            StdIntVector& list, StdIntVector& offset, StdIntVector& index,
                            double cutOff,
                            const int maxSize, const int K, const int stride,
                            StdNumericVector* collectionMax) {
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

      auto f = [k,K,stride,this](Score score) {
        const auto measure = distributionDistance(k.second - 1, score.second - 1, K, stride);
        return std::make_pair(measure, score.second);
      };

      Tick dNN;

      if (useTBB) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, std::distance(begin, end)),
                          [=](const tbb::blocked_range<size_t>& r) {
                            std::transform(begin + r.begin(), begin + r.end(), begin + r.begin(), f);
                          });
      } else {
        std::transform(begin, end, begin, f);
      }

      duration["dNN_____"] += dNN();

      Tick remainder;

      // sort the first maxNeighborhoodSize elements in increasing measure
      auto lastSort = (std::distance(begin, end) > maxSize) ? begin + maxSize : end;

      // 			std::cerr << "Partial length: " << std::distance(begin,lastSort) << std::endl;

      // 			return;

      std::partial_sort(begin, lastSort, end);

      // find cut-off point
      auto lastCut = std::find_if(begin, lastSort, [cutOff](Score score) {
        return score.first > cutOff; // TODO Could replace with binary search via std::lower_bound
      });

      // 			if (collectionMax) {
      // 				std::for_each(begin,lastCut, [](Score score) {
      // 					std::cerr << " (" << score.first << ":" << score.second << ")";
      // 				});
      // 				std::cerr << std::endl;
      // 			}

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

      if (collectionMax) {
        // 				std::cerr << "k = " << k.second << " : " << (lastCut - 1)->second << " with " << (lastCut - 1)->first << std::endl;
        collectionMax->push_back((lastCut - 1)->first);
      }

      if (extraSorting) {
        std::sort(std::end(list) - cutLength, std::end(list)); // TODO Ask SG: Do these need to be sorted?

        // TODO Ask SG: Why sort I.k again? Seems completely unnecessary
        std::sort(begin, end, [](const Score x, const Score y) {
          return x.second < y.second;
        });
      }

      duration["rem_____"] += remainder();
    }

    // 		std::for_each(begin, begin + 10, [](Score score) {
    // 			std::cerr << " (" << score.first << ":" << score.second << ")";
    // 		});
    // 		std::cerr << std::endl;
  }

#ifdef USE_NT2
  template <typename V>
  class DDFunctor : public RcppNT2::PlusReducer<V> {
  public:
    template <typename T>
    T map(const T& lhs, const T& rhs) {
      return sqrt(lhs * rhs);
    }
  };

  template <typename MapReducer, typename U, typename T, typename... Ts>
  U simdMapReduce(MapReducer&& f, U init, const T* it, const T* end, const Ts*... ts)
  {

    using namespace RcppNT2::variadic;

    typedef boost::simd::pack<T> vT; // SIMD vector of T
    typedef boost::simd::pack<U> vU; // SIMD vector of U

#define NEW

#ifdef NEW
    // Buffer for the SIMD mapping operations
    vU buffer = boost::simd::splat<vU>(init);

    static const std::size_t N = vT::static_size;

    // Aligned, SIMD operations
    for (; it != end; increment<N>(it, ts...))
      buffer = f.combine(
        f.map(boost::simd::aligned_load<vT>(it), boost::simd::aligned_load<vT>(ts)...),
        buffer
      );
    init = f.reduce(buffer);
#else
    static const std::size_t N = vT::static_size;
    const T* aligned_begin = std::min(boost::simd::align_on(it, N * sizeof(T)), end);
    const T* aligned_end   = aligned_begin + (end - aligned_begin) / N * N;

    // Buffer for the SIMD mapping operations
    vU buffer = boost::simd::splat<vU>(init);

    // Scalar operations for the initial unaligned region
    for (; it != aligned_begin; increment<1>(it, ts...))
      init = f.combine(f.map(*it, *ts...), init);

    // Aligned, SIMD operations
    for (; it != aligned_end; increment<N>(it, ts...))
      buffer = f.combine(
        f.map(boost::simd::aligned_load<vT>(it), boost::simd::aligned_load<vT>(ts)...),
        buffer
      );

    // Reduce the buffer, joining it into the scalar vale
    init = f.combine(f.reduce(buffer), init);

    // Leftover unaligned region
    for (; it != end; increment<1>(it, ts...))
      init = f.combine(f.map(*it, *ts...), init);
#endif

    return init;
  }
#endif // USE_NT2

  double distributionDistance(const int i, const int j, const int K, const int stride) {

    auto innerProduct =
      (useSSE ?
      simdMapReduce(DDFunctor<double>(), 0.0,
                                       &logSsMt[i * stride],
                                       &logSsMt[(i + 1) * stride],
                                       &logSsMt[j * stride])
      :
                  std::inner_product(
                    std::begin(logSsMt) + i * stride,
                    std::begin(logSsMt) + (i + 1) * stride,
                    std::begin(logSsMt) + j * stride, 0.0,
                    std::plus<double>(),
                    //std::product<double>()
                    [](double x, double y) {
                      return std::sqrt(x * y);
                    }
                  )
    );

    return 2.0 * (1.0 - innerProduct);
    // Unroll template:  http://cpplove.blogspot.com/2012/07/a-generic-loop-unroller-based-on.html
    // Unroll by hand + SSE: http://fastcpp.blogspot.com/2011/04/how-to-unroll-loop-in-c.html
    // RcppNT2
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

  std::vector<double,
              util::aligned_allocator<double, 16> // SSE alignment
              > logSsMt;

  bool extraSorting;
  const int specialMode;
  const bool useTBB;
  const bool useSSE;
  std::map<std::string,long long> duration;
};

template <typename RealType>
class CPUEngine : public AbstractEngine {
public:
  CPUEngine(bool sort, int specialMode) : AbstractEngine(sort, specialMode) { }

  virtual ~CPUEngine() { }

private:
  std::vector<RealType> t;
};

} // namespace np_cluster


#endif // _ABSTRACTENGINE_HPP
