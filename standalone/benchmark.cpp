
#include <chrono>
#include <random>
#include <iostream>

#include "Rcpp.h"
#include "AbstractEngine.h"

// template <typename T, typename PRNG, typename D>
// void generateLocation(T& locations, D& d, PRNG& prng) {
// 	for (auto& location : locations) {
// 		location = d(prng);
// 	}
// }

typedef std::shared_ptr<np_cluster::AbstractEngine> EnginePtr;


int main(int argc, char* argv[]) {

    using namespace Rcpp;

	auto seed = 666L;
	
	auto N = 25;
	auto P = 125;
	auto K = 70;
	
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

    std::cout << "Running benchmark ..." << std::endl;
    auto startTime = std::chrono::steady_clock::now();
    
    // Run

	auto endTime = std::chrono::steady_clock::now();
	auto duration = endTime - startTime;
	std::cout << "End benchmark" << std::endl;
	std::cout << std::chrono::duration<double, std::milli> (duration).count() << " ms "
			  << std::endl;

}
