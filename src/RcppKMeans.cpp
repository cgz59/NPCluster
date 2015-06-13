// #include "RcppMLPACK.h"
// 
// using namespace mlpack::kmeans;
// using namespace Rcpp;
// 
// // [[Rcpp::export(.fastKMeans)]]
// List fastKMeans(const arma::mat& data, const int& clusters) {
// 
// 	arma::Col<size_t> assignments;
// 	arma::mat centroids;
// 
// 	KMeans<> k;
// 
// 	k.Cluster(data, clusters, assignments, centroids);
// 
// 	return List::create(
// 		_["cluster"] = assignments,
// 		_["centers"] = centroids,
// 		_["size"] = clusters);
// }

