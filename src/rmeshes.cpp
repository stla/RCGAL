#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// // [[Rcpp::export]]
// Rcpp::List PolyMesh(const Rcpp::NumericMatrix points,
//                     const Rcpp::IntegerMatrix faces) {
//   Polyhedron mesh = makePolyMesh(points, faces);
//   return RPolyMesh(mesh);
// }

// [[Rcpp::export]]
Rcpp::List SurfMesh(const Rcpp::NumericMatrix points,
                    const Rcpp::List faces,
                    const bool isTriangle,
                    const bool triangulate,
                    const bool merge,
                    const bool normals,
                    const double epsilon) {
  // Polyhedron poly = makePolyMesh(points, faces);
  // Mesh3 mesh = Poly2Mesh3(poly);
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(points, faces, merge);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::IntegerMatrix Edges0;
  if(really_triangulate) {
    Edges0 = getEdges<Mesh3>(mesh);
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh = RSurfMesh<K, Mesh3, Point3>(
      mesh, isTriangle || triangulate, epsilon, false);
  if(really_triangulate) {
    routmesh["edges0"] = Edges0;
  }
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  if(normals) {
    auto vnormals = mesh.add_property_map<boost_vertex_descriptor, Vector3>(
                            "v:normals", CGAL::NULL_VECTOR)
                        .first;
    auto fnormals = mesh.add_property_map<boost_face_descriptor, Vector3>(
                            "f:normals", CGAL::NULL_VECTOR)
                        .first;
    CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals);
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
    routmesh["normals"] = Normals;
  }
  return routmesh;
}

/*
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
Rcpp::List SurfTMesh(const Rcpp::NumericMatrix points, const Rcpp::List faces,
const bool merge) { Mesh3 mesh = makeSurfMesh(points, faces, merge); bool
success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh); if(!success){
    Rcpp::stop("Triangulation has failed.");
  }
  return RSurfMesh(mesh);
}
*/

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::List Intersection(const Rcpp::List rmeshes,
                        const bool triangulate,
                        const bool merge,
                        const bool normals,
                        const bool exact) {
  const size_t nmeshes = rmeshes.size();
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Rcpp::NumericMatrix points = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  Rcpp::List faces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  MeshT mesh = makeSurfMesh<MeshT, PointT>(points, faces, merge);
  CGAL::Nef_polyhedron_3<KernelT> NP(mesh);
  Rcpp::Rcout << "NP defined - nfacets: " << NP.number_of_facets() << ".\n";
  for(size_t i = 1; i < nmeshes; i++) {
    {
      Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
      Rcpp::NumericMatrix points_i =
          Rcpp::as<Rcpp::NumericMatrix>(rmesh_i["vertices"]);
      Rcpp::List faces_i = Rcpp::as<Rcpp::List>(rmesh_i["faces"]);
      MeshT mesh_i = makeSurfMesh<MeshT, PointT>(points_i, faces_i, merge);
      CGAL::Nef_polyhedron_3<KernelT> NP_i(mesh_i);
      Rcpp::Rcout << "NP_i defined. - nfacets: " << NP_i.number_of_facets()
                  << ".\n";
      NP = NP_i * NP;
      Rcpp::Rcout << "intersection done - nfacets: " << NP.number_of_facets()
                  << ".\n";
    }
  }
  MeshT outmesh;
  // CGAL::convert_nef_polyhedron_to_polygon_mesh(NP, outmesh, false);
  std::vector<PointT> verts;
  std::vector<std::vector<size_t>> indices;
  CGAL::convert_nef_polyhedron_to_polygon_soup(NP, verts, indices, false);
  Rcpp::Rcout << "conversion to soup done.\n";
  bool success =
      CGAL::Polygon_mesh_processing::orient_polygon_soup(verts, indices);
  if(!success) {
    Rcpp::stop("Polygon orientation failed XXX.");
  }
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(verts, indices,
                                                              outmesh);
  Rcpp::Rcout << "conversion to mesh done - nvertices: "
              << outmesh.number_of_vertices() << ".\n";
  Rcpp::IntegerMatrix Edges0;
  if(triangulate) {
    Edges0 = getEdges<MeshT>(outmesh);
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(outmesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh =
      RSurfMesh<KernelT, MeshT, PointT>(outmesh, triangulate, 0, exact);
  if(triangulate) {
    routmesh["edges0"] = Edges0;
  }
  // const size_t nvertices = outmesh.number_of_vertices();
  // Rcpp::NumericMatrix Normals(3, nvertices);
  // if(normals) {
  //   if(exact) {
  //     auto vnormals = outmesh
  //                         .add_property_map<boost_vertex_descriptor, EVector3>(
  //                             "v:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     auto fnormals = outmesh
  //                         .add_property_map<boost_face_descriptor, EVector3>(
  //                             "f:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     CGAL::Polygon_mesh_processing::compute_normals(outmesh, vnormals,
  //                                                    fnormals);
  //     {
  //       size_t i = 0;
  //       for(boost_vertex_descriptor vd : vertices(outmesh)) {
  //         Rcpp::NumericVector col_i(3);
  //         const EVector3 normal = vnormals[vd];
  //         col_i(0) = CGAL::to_double(normal.x());
  //         col_i(1) = CGAL::to_double(normal.y());
  //         col_i(2) = CGAL::to_double(normal.z());
  //         Normals(Rcpp::_, i) = col_i;
  //         i++;
  //       }
  //     }
  //   } else {
  //     auto vnormals = outmesh
  //                         .add_property_map<boost_vertex_descriptor, Vector3>(
  //                             "v:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     auto fnormals = outmesh
  //                         .add_property_map<boost_face_descriptor, Vector3>(
  //                             "f:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     CGAL::Polygon_mesh_processing::compute_normals(outmesh, vnormals,
  //                                                    fnormals);
  //     {
  //       size_t i = 0;
  //       for(boost_vertex_descriptor vd : vertices(outmesh)) {
  //         Rcpp::NumericVector col_i(3);
  //         const Vector3 normal = vnormals[vd];
  //         col_i(0) = normal.x();
  //         col_i(1) = normal.y();
  //         col_i(2) = normal.z();
  //         Normals(Rcpp::_, i) = col_i;
  //         i++;
  //       }
  //     }
  //   }
  //   routmesh["normals"] = Normals;
  // }
  return routmesh;
}

// [[Rcpp::export]]
Rcpp::List Intersection_K(const Rcpp::List rmeshes,
                          const bool triangulate,
                          const bool merge,
                          const bool normals) {
  return Intersection<K, Mesh3, Point3>(rmeshes, triangulate, merge, normals,
                                        false);
}

// [[Rcpp::export]]
Rcpp::List Intersection_EK(const Rcpp::List rmeshes,
                           const bool triangulate,
                           const bool merge,
                           const bool normals) {
  return Intersection<EK, EMesh3, EPoint3>(rmeshes, triangulate, merge, normals,
                                           true);
}

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::List Intersection2(const Rcpp::List rmeshes,  // must be triangles
                         const bool merge,
                         const bool normals,
                         const bool exact) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  Rcpp::NumericMatrix points = Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  Rcpp::List faces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  meshes[0] = makeSurfMesh<MeshT, PointT>(points, faces, merge);
  for(size_t i = 1; i < nmeshes; i++) {
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    Rcpp::NumericMatrix points_i =
        Rcpp::as<Rcpp::NumericMatrix>(rmesh_i["vertices"]);
    Rcpp::List faces_i = Rcpp::as<Rcpp::List>(rmesh_i["faces"]);
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(points_i, faces_i, merge);
    bool ok = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(
        meshes[i - 1], mesh_i, meshes[i]);
    Rcpp::Rcout << "intersection: " << ok << "\n";
  }
  MeshT mesh = meshes[nmeshes - 1];
  Rcpp::List routmesh = RSurfMesh<KernelT, MeshT, PointT>(mesh, true, 0, exact);
  // const size_t nvertices = mesh.number_of_vertices();
  // Rcpp::NumericMatrix Normals(3, nvertices);
  // if(normals) {
  //   if(exact) {
      // auto vnormals = mesh.add_property_map<typename MeshT::Vertex_index, typename KernelT::Vector_3>(
      //                         "v:normals", CGAL::NULL_VECTOR)
      //                     .first;
      // auto fnormals = mesh.add_property_map<typename MeshT::Face_index, typename KernelT::Vector_3>(
      //                         "f:normals", CGAL::NULL_VECTOR)
      //                     .first;
      // CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals);
      // {
      //   size_t i = 0;
      //   for(typename MeshT::Vertex_index vd : vertices(mesh)) {
      //     Rcpp::NumericVector col_i(3);
      //     const typename KernelT::Vector_3 normal = vnormals[vd];
      //     col_i(0) = CGAL::to_double(normal.x());
      //     col_i(1) = CGAL::to_double(normal.y());
      //     col_i(2) = CGAL::to_double(normal.z());
      //     Normals(Rcpp::_, i) = col_i;
      //     i++;
      //   }
      // }
  //   } else {
  //     auto vnormals = mesh.add_property_map<boost_vertex_descriptor, Vector3>(
  //                             "v:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     auto fnormals = mesh.add_property_map<boost_face_descriptor, Vector3>(
  //                             "f:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals, fnormals);
  //     {
  //       size_t i = 0;
  //       for(boost_vertex_descriptor vd : vertices(mesh)) {
  //         Rcpp::NumericVector col_i(3);
  //         const Vector3 normal = vnormals[vd];
  //         col_i(0) = normal.x();
  //         col_i(1) = normal.y();
  //         col_i(2) = normal.z();
  //         Normals(Rcpp::_, i) = col_i;
  //         i++;
  //       }
  //     }
  //   }
  //  routmesh["normals"] = Normals;
  //}
  return routmesh;
}

// [[Rcpp::export]]
Rcpp::List Intersection2_K(const Rcpp::List rmeshes,
                           const bool merge,
                           const bool normals) {
  return Intersection2<K, Mesh3, Point3>(rmeshes, merge, normals, false);
}

// [[Rcpp::export]]
Rcpp::List Intersection2_EK(const Rcpp::List rmeshes,
                            const bool merge,
                            const bool normals) {
  return Intersection2<EK, EMesh3, EPoint3>(rmeshes, merge, normals, true);
}
