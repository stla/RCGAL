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
Rcpp::List SurfMesh(const Rcpp::NumericMatrix points, const Rcpp::List faces, const bool merge) {
  // Polyhedron poly = makePolyMesh(points, faces);
  // Mesh3 mesh = Poly2Mesh3(poly);
  Mesh3 mesh = makeSurfMesh(points, faces, merge);
  return RSurfMesh(mesh);
}

// [[Rcpp::export]]
Rcpp::List SurfMeshWithNormals(const Rcpp::NumericMatrix points,
                               const Rcpp::List faces,
                               const bool merge) {
  // Polyhedron poly = makePolyMesh(points, faces);
  // Mesh3 mesh;
  // CGAL::copy_face_graph(poly, mesh);
  Mesh3 mesh = makeSurfMesh(points, faces, merge);
  auto vnormals = mesh.add_property_map<boost_vertex_descriptor, Vector3>(
                          "v:normals", CGAL::NULL_VECTOR)
                      .first;
  auto fnormals = mesh.add_property_map<boost_face_descriptor, Vector3>(
                          "f:normals", CGAL::NULL_VECTOR)
                      .first;
  CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals);
  Rcpp::List surfmesh = RSurfMesh(mesh);
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  {
    size_t i = 0;
    for(boost_vertex_descriptor vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const Vector3 normal = vnormals[vd];
      col_i(0) = normal.x();
      col_i(1) = normal.y();
      col_i(2) = normal.z();
      Normals(Rcpp::_, i) = col_i;
      i++;
    }
  }
  surfmesh["normals"] = Normals;
  return surfmesh;
}

// [[Rcpp::export]]
Rcpp::List SurfTMesh(const Rcpp::NumericMatrix points, const Rcpp::List faces, const bool merge) {
  Mesh3 mesh = makeSurfMesh(points, faces, merge);
  bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
  if(!success){
    Rcpp::stop("Triangulation has failed.");
  }
  return RSurfMesh(mesh);
}

// [[Rcpp::export]]
Rcpp::List Intersection(const Rcpp::List rmeshes, const bool merge, const bool triangulate) {
  const size_t nmeshes = rmeshes.size();
  Rcpp::List rmesh = rmeshes(0);
  Rcpp::NumericMatrix points = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  Rcpp::List faces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  Mesh3 mesh = makeSurfMesh(points, faces, merge);
  Nef NP(mesh);
  for(size_t i=1; i < nmeshes; i++){
    rmesh = rmeshes(i);
    points = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
    faces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
    mesh = makeSurfMesh(points, faces, merge);
    Nef np(mesh);
    NP = NP * np;
  }
  Mesh3 outmesh;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(NP, outmesh, true);
  if(triangulate){
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(outmesh);
    if(!success){
      Rcpp::stop("Triangulation has failed.");
    }
  }
  return RSurfMesh(outmesh);
}
