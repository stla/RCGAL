#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// [[Rcpp::export]]
Rcpp::List PolyMesh(const Rcpp::NumericMatrix points,
                    const Rcpp::IntegerMatrix faces) {
  Polyhedron mesh = makePolyMesh(points, faces);
  return RPolyMesh(mesh);
}

// [[Rcpp::export]]
Rcpp::List SurfMesh(const Rcpp::NumericMatrix points,
                    const Rcpp::IntegerMatrix faces) {
  Polyhedron poly = makePolyMesh(points, faces);
  Mesh3 mesh = Poly2Mesh3(poly);
  return RSurfMesh(mesh);
}
