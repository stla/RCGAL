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
Rcpp::List SurfMesh(const Rcpp::List rmesh,
                    const bool isTriangle,
                    const bool triangulate,
                    const bool merge,
                    const bool normals,
                    const double epsilon) {
  // Polyhedron poly = makePolyMesh(points, faces);
  // Mesh3 mesh = Poly2Mesh3(poly);
  Mesh3 mesh = makeSurfMesh<Mesh3, Point3>(rmesh, merge);
  const bool really_triangulate = !isTriangle && triangulate;
  Rcpp::IntegerMatrix Edges0;
  if(really_triangulate) {
    Edges0 = getEdges1<Mesh3>(mesh);
    bool success = CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  Rcpp::List routmesh = RSurfKMesh(mesh, isTriangle || triangulate, epsilon);
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
    PMP::compute_normals(mesh, vnormals, fnormals);
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

// template <typename KernelT, typename MeshT, typename PointT>
// Rcpp::List Intersection(const Rcpp::List rmeshes,
//                         const bool triangulate,
//                         const bool merge,
//                         const bool normals,
//                         const bool exact) {
//   const size_t nmeshes = rmeshes.size();
//   Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
//   Rcpp::NumericMatrix points =
//   Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]); Rcpp::List faces =
//   Rcpp::as<Rcpp::List>(rmesh["faces"]); MeshT mesh = makeSurfMesh<MeshT,
//   PointT>(points, faces, merge); CGAL::Nef_polyhedron_3<KernelT> NP(mesh);
//   Rcpp::Rcout << "NP defined - nfacets: " << NP.number_of_facets() << ".\n";
//   for(size_t i = 1; i < nmeshes; i++) {
//     {
//       Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
//       Rcpp::NumericMatrix points_i =
//           Rcpp::as<Rcpp::NumericMatrix>(rmesh_i["vertices"]);
//       Rcpp::List faces_i = Rcpp::as<Rcpp::List>(rmesh_i["faces"]);
//       MeshT mesh_i = makeSurfMesh<MeshT, PointT>(points_i, faces_i, merge);
//       CGAL::Nef_polyhedron_3<KernelT> NP_i(mesh_i);
//       Rcpp::Rcout << "NP_i defined. - nfacets: " << NP_i.number_of_facets()
//                   << ".\n";
//       NP = NP_i * NP;
//       Rcpp::Rcout << "intersection done - nfacets: " << NP.number_of_facets()
//                   << ".\n";
//     }
//   }
//   MeshT outmesh;
//   // CGAL::convert_nef_polyhedron_to_polygon_mesh(NP, outmesh, false);
//   std::vector<PointT> verts;
//   std::vector<std::vector<size_t>> indices;
//   CGAL::convert_nef_polyhedron_to_polygon_soup(NP, verts, indices, false);
//   Rcpp::Rcout << "conversion to soup done.\n";
//   bool success =
//       CGAL::Polygon_mesh_processing::orient_polygon_soup(verts, indices);
//   if(!success) {
//     Rcpp::stop("Polygon orientation failed XXX.");
//   }
//   CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(verts, indices,
//                                                               outmesh);
//   Rcpp::Rcout << "conversion to mesh done - nvertices: "
//               << outmesh.number_of_vertices() << ".\n";
//   Rcpp::IntegerMatrix Edges0;
//   if(triangulate) {
//     Edges0 = getEdges<MeshT>(outmesh);
//     bool success = CGAL::Polygon_mesh_processing::triangulate_faces(outmesh);
//     if(!success) {
//       Rcpp::stop("Triangulation has failed.");
//     }
//   }
//   Rcpp::List routmesh =
//       RSurfMesh<KernelT, MeshT, PointT>(outmesh, triangulate, 0, exact);
//   if(triangulate) {
//     routmesh["edges0"] = Edges0;
//   }
//   // const size_t nvertices = outmesh.number_of_vertices();
//   // Rcpp::NumericMatrix Normals(3, nvertices);
//   // if(normals) {
//   //   if(exact) {
//   //     auto vnormals = outmesh
//   //                         .add_property_map<boost_vertex_descriptor,
//   EVector3>(
//   //                             "v:normals", CGAL::NULL_VECTOR)
//   //                         .first;
//   //     auto fnormals = outmesh
//   //                         .add_property_map<boost_face_descriptor,
//   EVector3>(
//   //                             "f:normals", CGAL::NULL_VECTOR)
//   //                         .first;
//   //     CGAL::Polygon_mesh_processing::compute_normals(outmesh, vnormals,
//   //                                                    fnormals);
//   //     {
//   //       size_t i = 0;
//   //       for(boost_vertex_descriptor vd : vertices(outmesh)) {
//   //         Rcpp::NumericVector col_i(3);
//   //         const EVector3 normal = vnormals[vd];
//   //         col_i(0) = CGAL::to_double(normal.x());
//   //         col_i(1) = CGAL::to_double(normal.y());
//   //         col_i(2) = CGAL::to_double(normal.z());
//   //         Normals(Rcpp::_, i) = col_i;
//   //         i++;
//   //       }
//   //     }
//   //   } else {
//   //     auto vnormals = outmesh
//   //                         .add_property_map<boost_vertex_descriptor,
//   Vector3>(
//   //                             "v:normals", CGAL::NULL_VECTOR)
//   //                         .first;
//   //     auto fnormals = outmesh
//   //                         .add_property_map<boost_face_descriptor,
//   Vector3>(
//   //                             "f:normals", CGAL::NULL_VECTOR)
//   //                         .first;
//   //     CGAL::Polygon_mesh_processing::compute_normals(outmesh, vnormals,
//   //                                                    fnormals);
//   //     {
//   //       size_t i = 0;
//   //       for(boost_vertex_descriptor vd : vertices(outmesh)) {
//   //         Rcpp::NumericVector col_i(3);
//   //         const Vector3 normal = vnormals[vd];
//   //         col_i(0) = normal.x();
//   //         col_i(1) = normal.y();
//   //         col_i(2) = normal.z();
//   //         Normals(Rcpp::_, i) = col_i;
//   //         i++;
//   //       }
//   //     }
//   //   }
//   //   routmesh["normals"] = Normals;
//   // }
//   return routmesh;
// }

// // [[Rcpp::export]]
// Rcpp::List Intersection_K(const Rcpp::List rmeshes,
//                           const bool triangulate,
//                           const bool merge,
//                           const bool normals) {
//   return Intersection<K, Mesh3, Point3>(rmeshes, triangulate, merge, normals,
//                                         false);
// }

// // [[Rcpp::export]]
// Rcpp::List Intersection_EK(const Rcpp::List rmeshes,
//                            const bool triangulate,
//                            const bool merge,
//                            const bool normals) {
//   return Intersection<EK, EMesh3, EPoint3>(rmeshes, triangulate, merge,
//   normals,
//                                            true);
// }

template <typename MeshT>
void checkMesh(MeshT mesh, size_t i) {
  const bool si = PMP::does_self_intersect(mesh);
  if(si) {
    std::string msg = "Mesh n\u00b0" + std::to_string(i) + " self-intersects.";
    Rcpp::stop(msg);
  }
  const bool bv = PMP::does_bound_a_volume(mesh);
  if(!bv) {
    std::string msg =
        "Mesh n\u00b0" + std::to_string(i) + " does not bound a volume.";
    Rcpp::stop(msg);
  }
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Intersection2(const Rcpp::List rmeshes,  // must be triangles
                    const bool merge,
                    const bool normals,
                    const bool exact) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  meshes[0] = makeSurfMesh<MeshT, PointT>(rmesh, merge);
  if(exact) {
    checkMesh<MeshT>(meshes[0], 1);
  }
  for(size_t i = 1; i < nmeshes; i++) {
    if(!exact) {
      checkMesh<MeshT>(meshes[i - 1], i);
    }
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, merge);
    checkMesh<MeshT>(mesh_i, i + 1);
    bool ok = PMP::corefine_and_compute_intersection(meshes[i - 1], mesh_i,
                                                     meshes[i]);
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
  // const size_t nvertices = mesh.number_of_vertices();
  // Rcpp::NumericMatrix Normals(3, nvertices);
  // if(normals) {
  //   if(exact) {
  // auto vnormals = mesh.add_property_map<typename MeshT::Vertex_index,
  // typename KernelT::Vector_3>(
  //                         "v:normals", CGAL::NULL_VECTOR)
  //                     .first;
  // auto fnormals = mesh.add_property_map<typename MeshT::Face_index, typename
  // KernelT::Vector_3>(
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
  //     auto vnormals = mesh.add_property_map<boost_vertex_descriptor,
  //     Vector3>(
  //                             "v:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     auto fnormals = mesh.add_property_map<boost_face_descriptor, Vector3>(
  //                             "f:normals", CGAL::NULL_VECTOR)
  //                         .first;
  //     CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals,
  //     fnormals);
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
  // return routmesh;
}

// [[Rcpp::export]]
Rcpp::List Intersection2_K(const Rcpp::List rmeshes,
                           const bool merge,
                           const bool normals) {
  Mesh3 mesh = Intersection2<K, Mesh3, Point3>(rmeshes, merge, normals, false);
  return RSurfKMesh(mesh, true, 0);
}

// [[Rcpp::export]]
Rcpp::List Intersection2_EK(const Rcpp::List rmeshes,
                            const bool merge,
                            const bool normals) {
  EMesh3 mesh = Intersection2<EK, EMesh3, EPoint3>(rmeshes, merge, normals, true);
  return RSurfEKMesh(mesh, true, 0);
}

// [[Rcpp::export]]
Rcpp::List Intersection_Q(const Rcpp::List rmeshes,  // must be triangles
                          const bool merge,
                          const bool normals) {
  const size_t nmeshes = rmeshes.size();
  std::vector<QMesh3> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  meshes[0] = makeSurfQMesh(rmesh, merge);
  checkMesh<QMesh3>(meshes[0], 0);
  for(size_t i = 1; i < nmeshes; i++) {
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    QMesh3 mesh_i = makeSurfQMesh(rmesh_i, merge);
    checkMesh<QMesh3>(mesh_i, i);
    bool ok = PMP::corefine_and_compute_intersection(meshes[i - 1], mesh_i,
                                                     meshes[i]);
    if(!ok) {
      Rcpp::stop("Intersection computation has failed.");
    }
  }
  return RSurfQMesh(meshes[nmeshes - 1], merge, 0);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Difference(const Rcpp::List rmesh1,  // must be triangles
                 const Rcpp::List rmesh2,
                 const bool merge,
                 const bool normals) {
  MeshT smesh1 = makeSurfMesh<MeshT, PointT>(rmesh1, merge);
  checkMesh<MeshT>(smesh1, 1);
  MeshT smesh2 = makeSurfMesh<MeshT, PointT>(rmesh2, merge);
  checkMesh<MeshT>(smesh2, 2);
  MeshT outmesh;
  bool ok = PMP::corefine_and_compute_difference(smesh1, smesh2, outmesh);
  if(!ok) {
    Rcpp::stop("Difference computation has failed.");
  }
  return outmesh;
}

// [[Rcpp::export]]
Rcpp::List Difference_K(const Rcpp::List rmesh1,
                        const Rcpp::List rmesh2,
                        const bool merge,
                        const bool normals) {
  Mesh3 mesh = Difference<K, Mesh3, Point3>(rmesh1, rmesh2, merge, normals);
  return RSurfKMesh(mesh, true, 0);
}

// [[Rcpp::export]]
Rcpp::List Difference_EK(const Rcpp::List rmesh1,
                         const Rcpp::List rmesh2,
                         const bool merge,
                         const bool normals) {
  EMesh3 mesh = Difference<EK, EMesh3, EPoint3>(rmesh1, rmesh2, merge, normals);
  return RSurfEKMesh(mesh, true, 0);
}

template <typename KernelT, typename MeshT, typename PointT>
MeshT Union(const Rcpp::List rmeshes,  // must be triangles
            const bool merge,
            const bool normals) {
  const size_t nmeshes = rmeshes.size();
  std::vector<MeshT> meshes(nmeshes);
  Rcpp::List rmesh = Rcpp::as<Rcpp::List>(rmeshes(0));
  meshes[0] = makeSurfMesh<MeshT, PointT>(rmesh, merge);
  checkMesh<MeshT>(meshes[0], 1);
  for(size_t i = 1; i < nmeshes; i++) {
    Rcpp::List rmesh_i = Rcpp::as<Rcpp::List>(rmeshes(i));
    MeshT mesh_i = makeSurfMesh<MeshT, PointT>(rmesh_i, merge);
    checkMesh<MeshT>(mesh_i, i + 1);
    bool ok = PMP::corefine_and_compute_union(meshes[i - 1], mesh_i, meshes[i]);
    if(!ok) {
      Rcpp::stop("Union computation has failed.");
    }
  }
  return meshes[nmeshes - 1];
}

// [[Rcpp::export]]
Rcpp::List Union_K(const Rcpp::List rmeshes,
                   const bool merge,
                   const bool normals) {
  Mesh3 mesh = Union<K, Mesh3, Point3>(rmeshes, merge, normals);
  return RSurfKMesh(mesh, true, 0);
}

// [[Rcpp::export]]
Rcpp::List Union_EK(const Rcpp::List rmeshes,
                    const bool merge,
                    const bool normals) {
  EMesh3 mesh = Union<EK, EMesh3, EPoint3>(rmeshes, merge, normals);
  return RSurfEKMesh(mesh, true, 0);
}
