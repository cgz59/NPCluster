
#include <chrono>
#include <random>
#include <iostream>

//#include "Rcpp.h"
#include <RInside.h>
#include "AbstractEngine.h"

// template <typename T, typename PRNG, typename D>
// void generateLocation(T& locations, D& d, PRNG& prng) {
// 	for (auto& location : locations) {
// 		location = d(prng);
// 	}
// }

typedef std::shared_ptr<np_cluster::AbstractEngine> EnginePtr;

__attribute__((constructor))
static void initialize(void) {
	RInside R(0, nullptr); // This is a horrible hack to allow Rcpp memory objects as global variables
	// TODO Real solution is to move RInside() into main() and *not* use global variables
}

int main(int argc, char* argv[]) {

    using namespace Rcpp;
    
	auto seed = 666L;
    
	// Set the R RNG seed
	Rcpp::Environment baseEnv("package:base");
	Rcpp::Function setSeed = baseEnv["set.seed"];
	setSeed(seed);    
		
	auto N = 100000;
	auto P = 125;
	auto K = 50;
	
	auto tau = 0.5;
	auto tau0 = 2.0;
	auto tauInt = 0.5;
	
	auto extraSort = true;
	auto colSum = false;
		
	auto prng = std::mt19937(seed);
	auto normal = std::normal_distribution<double>(0.0, 1.0);
	auto uniform = std::uniform_int_distribution<int>(0, K);
	auto binomial = std::bernoulli_distribution(0.75);
	auto normalData = std::normal_distribution<double>(0.0, 1.0);
	
	EnginePtr engine = EnginePtr(
	    new np_cluster::AbstractEngine(extraSort));
	
	std::cout << "Loading data" << std::endl;
	
	NumericMatrix X(N, P);
	NumericMatrix A(N, K);
	IntegerMatrix S(N, K);
	
	IntegerVector nVec(K);
	NumericVector Y(N);
	NumericVector Xsd(N);
	IntegerVector rowSubsetI(N);
	IntegerVector cmVec(1000);
	NumericVector phiV(K);
	
	for (auto& x : X) {
	    x = normal(prng);
	}
	
	for (auto it = std::begin(X); it != std::end(X); ++it) {
	    *it = normal(prng);
	}
	
	for (auto it = std::begin(A); it != std::end(A); ++it) {
	    *it = normal(prng);
	}
	
	for (auto it = std::begin(S); it != std::end(S); ++it) {
	    *it = uniform(prng) + 1;
	}
	
	for (auto& n : nVec) {
		n = uniform(prng) + 1;
	}
	
	for (auto& y : Y) {
		y = normal(prng);
	}
	
	for (auto& x : Xsd) {
		x = normal(prng);
	}
	
	for (auto& c : cmVec) {
		c = uniform(prng) + 1;
	}
	
	for (auto& phi : phiV) {
		phi = normal(prng);
	}
	
	for (int i = 0; i < N; ++i) {
		rowSubsetI[i] = i + 1;
	}

	auto n0 = 0;
	auto n2 = 100;
	auto epsilon = 0.001;
	auto epsilon2 = 0.00001;
	auto maxNeighborhoodSize = 28;
	auto cutOff = 0.1;
	
    std::cout << "Running benchmark ..." << std::endl;
    auto startTime = std::chrono::steady_clock::now();
    
    // Run
    
    engine->computePmfAndNeighborhoods(
    	n0, nVec, epsilon, epsilon2, K, N, 
    	Y, Xsd, rowSubsetI, cmVec, n2, phiV, 
    	tau, tau0, maxNeighborhoodSize, cutOff);
    
    
//   Rcpp::List computePmfAndNeighborhoods(int n0, const IntegerVector& nVec,
//                                               const double epsilon, const double epsilon2,
//                                               int K, int N,  // NB: K -> G for column-space
//                                               const NumericVector& Y, const NumericMatrix& X,
//                                               const NumericMatrix& A, const IntegerMatrix& S,
//                                               const IntegerVector& rowSubsetI,
//                                               const IntegerVector& CmVec, const int P,
//                                               const NumericVector& phiV, const double tau, const double tau0, const double tauInt,
//                                               const int maxNeighborhoodSize,
//                                               const double cutOff,
//                                               const bool collectMax)     

	auto endTime = std::chrono::steady_clock::now();
	auto duration = endTime - startTime;
	std::cout << "End benchmark" << std::endl;
	std::cout << std::chrono::duration<double, std::milli> (duration).count() << " ms "
			  << std::endl;

}
