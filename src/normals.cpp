#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// [[Rcpp::export]]
Rcpp::NumericMatrix jet_normals_cpp(Rcpp::NumericMatrix pts,
                                    unsigned nb_neighbors) {
  const size_t npoints = pts.nrow();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)),
                               Vector3(0.0, 0.0, 0.0));
  }

  CGAL::jet_estimate_normals<Concurrency_tag>(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  CGAL::mst_orient_normals(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  Rcpp::NumericMatrix normals(npoints, 3);
  for(size_t i = 0; i < npoints; i++) {
    const Vector3 normal = points[i].second;
    normals(i, 0) = normal.x();
    normals(i, 1) = normal.y();
    normals(i, 2) = normal.z();
  }

  return normals;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pca_normals_cpp(Rcpp::NumericMatrix pts,
                                    unsigned nb_neighbors) {
  const size_t npoints = pts.nrow();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)),
                               Vector3(0.0, 0.0, 0.0));
  }

  CGAL::pca_estimate_normals<Concurrency_tag>(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  CGAL::mst_orient_normals(
      points, nb_neighbors,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
          .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));

  Rcpp::NumericMatrix normals(npoints, 3);
  for(size_t i = 0; i < npoints; i++) {
    const Vector3 normal = points[i].second;
    normals(i, 0) = normal.x();
    normals(i, 1) = normal.y();
    normals(i, 2) = normal.z();
  }

  return normals;
}
