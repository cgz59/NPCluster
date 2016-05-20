#ifndef _NEIGHBORHOODGRAPH_HPP
#define _NEIGHBORHOODGRAPH_HPP

#include <iostream>
#include <sstream>
#include <algorithm>

#include "Rcpp.h"

namespace np_cluster {

using namespace Rcpp;

// struct Tuple {
//   int i;
//   int j;
//
// };

class NeighborhoodGraph {
public:

  typedef size_t KeyType;
  typedef int CountType;

  NeighborhoodGraph(const int size) : size(size) {
    // Rcpp::Rcout << "ctor NG" << std::endl;
  }

//   NeighborhoodGraph(const IntegerVector& i,
//                     const IntegerVector& j,
//                     const int size) : size(size) {
//     // TODO Ensure i.size() == j.size()
//     for (int idx = 0; idx < i.size(); ++idx) {
//       graph[getIndex(i[idx] - 1, j[idx] - 1)] = 1;
//     }
//   }

  virtual ~NeighborhoodGraph() {
    // Rcpp::Rcout << "dtor NG" << std::endl;
  }

  void incrementGraph(const IntegerVector& indices) {
    const int length = indices.size();
    for (int i = 0; i < length; ++i) {
      for (int j = 0 /* i */; j < length; ++j) { // TODO Can just save upper-triangular
        graph[getIndex(indices[i] - 1, indices[j] - 1)] += static_cast<CountType>(1);
      }
    }
    // TODO Could -1 from list once at beginning
  }

  double taxiDistance(const NeighborhoodGraph& rhs) {
    // TODO Ensure this.size == rhs.size

    // Get key sets
    std::vector<KeyType> lhsKeys;
    for (const auto e : graph) {
      lhsKeys.push_back(e.first);
    }

    std::vector<KeyType> rhsKeys;
    for (const auto e : rhs.graph) {
      rhsKeys.push_back(e.first);
    }

    std::vector<KeyType> difference;
    std::set_symmetric_difference(std::begin(lhsKeys), std::end(lhsKeys),
                                  std::begin(rhsKeys), std::end(rhsKeys),
                                  std::back_inserter(difference));

    double result = difference.size();

    return result / (size * size);
  }

  List getGraph() {
    std::vector<int> i;
    std::vector<int> j;
    std::vector<CountType> x;

    for (const auto entry : graph) {
      auto ij = getIJ(entry.first);
      i.push_back(ij.first + 1);
      j.push_back(ij.second + 1);
      x.push_back(entry.second);
    }
    return List::create(
      Named("i") = i,
      Named("j") = j,
      Named("x") = x
    );

  }

private:

  KeyType getIndex(int i, int j) {
    return j * size + i;
  }

  std::pair<int,int> getIJ(KeyType key) {
    int i = key / size;
    int j = key % size;
    return std::make_pair(i, j);
  }

  const int size;
  std::map<KeyType, CountType> graph;
};

} // namespace np_cluster

#endif // _NEIGHBORHOODGRAPH_HPP
