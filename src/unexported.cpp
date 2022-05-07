#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// double sepsilon = sqrt(std::numeric_limits<double>::epsilon());

bool approxEqual(double x, double y, double epsilon) {
  return fabs(x - y) <= epsilon;
}

bool approxEqualVectors(Vector3 v, Vector3 w, double epsilon) {
  return approxEqual(v.x(), w.x(), epsilon) &&
         approxEqual(v.y(), w.y(), epsilon) &&
         approxEqual(v.z(), w.z(), epsilon);
}

Rcpp::String stringPair(const size_t i, const size_t j) {
  const size_t i0 = std::min(i, j);
  const size_t i1 = std::max(i, j);
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
      std::to_string(i0), "-", std::to_string(i1));
  return Rcpp::collapse(stringids);
}

std::array<Rcpp::String, 3> triangleEdges(const size_t i0,
                                          const size_t i1,
                                          const size_t i2) {
  return {stringPair(i0, i1), stringPair(i0, i2), stringPair(i1, i2)};
}

Rcpp::String stringTriple(const size_t i, const size_t j, const size_t k) {
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
      std::to_string(i), "-", std::to_string(j), "-", std::to_string(k));
  return Rcpp::collapse(stringids);
}

double volume_under_triangle(Dvector v0, Dvector v1, Dvector v2) {
  double x0 = v0(0);
  double y0 = v0(1);
  double z0 = v0(2);
  double x1 = v1(0);
  double y1 = v1(1);
  double z1 = v1(2);
  double x2 = v2(0);
  double y2 = v2(1);
  double z2 = v2(2);
  return (z0 + z1 + z2) *
         (x0 * y1 - x1 * y0 + x1 * y2 - x2 * y1 + x2 * y0 - x0 * y2) / 6.0;
}

void mark_domains0(CDT& ct,
                   CDFace_handle start,
                   int index,
                   std::list<CDT::Edge>& border) {
  if(start->info().nesting_level != -1) {
    return;
  }
  std::list<CDFace_handle> queue;
  queue.push_back(start);
  while(!queue.empty()) {
    CDFace_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++) {
        CDT::Edge e(fh, i);
        CDFace_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1) {
          if(ct.is_constrained(e))
            border.push_back(e);
          else
            queue.push_back(n);
        }
      }
    }
  }
}

void mark_domains(CDT& cdt) {
  for(CDFace_handle f : cdt.all_face_handles()) {
    f->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains0(cdt, cdt.infinite_face(), 0, border);
  while(!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDFace_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1) {
      mark_domains0(cdt, n, e.first->info().nesting_level + 1, border);
    }
  }
}

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<PointT> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<Point3> matrix_to_points3<Point3>(Rcpp::NumericMatrix);
template std::vector<EPoint3> matrix_to_points3<EPoint3>(Rcpp::NumericMatrix);

// std::vector<std::vector<size_t>> matrix_to_faces(const Rcpp::IntegerMatrix I)
// {
//   const size_t nfaces = I.ncol();
//   std::vector<std::vector<size_t>> faces;
//   faces.reserve(nfaces);
//   for(size_t i = 0; i < nfaces; i++) {
//     const Rcpp::IntegerVector face_rcpp = I(Rcpp::_, i);
//     std::vector<size_t> face(face_rcpp.begin(), face_rcpp.end());
//     faces.emplace_back(face);
//   }
//   return faces;
// }

std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List L) {
  const size_t nfaces = L.size();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    // face_rcpp = face_rcpp - 1;
    std::vector<size_t> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}

// Polyhedron makePolyMesh(const Rcpp::NumericMatrix M,
//                         const Rcpp::IntegerMatrix I) {
//   Points3 points = matrix_to_points3(M);
//   std::vector<std::vector<size_t>> faces = matrix_to_faces(I);
//   bool success =
//       CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
//   if(!success) {
//     Rcpp::stop("Polygon orientation failed.");
//   }
//   Polyhedron mesh;
//   CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces,
//                                                               mesh);
//   return mesh;
// }

// Rcpp::List RPolyMesh(Polyhedron mesh) {
//   int id = 1;
//   for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
//       vit != mesh.vertices_end(); ++vit) {
//     vit->id() = id;
//     id++;
//   }

//   const size_t nfacets = mesh.size_of_facets();
//   const size_t nvertices = mesh.size_of_vertices();

//   // const size_t nvertices = num_vertices(mesh);  // or number_of_vertices
//   // // const size_t nedges = num_edges(mesh);        // or number_of_edges
//   // const size_t nfaces = num_faces(mesh);
//   Rcpp::NumericMatrix vertices(3, nvertices);
//   {
//     size_t i = 0;
//     for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
//         vit != mesh.vertices_end(); vit++) {
//       Rcpp::NumericVector col_i(3);
//       col_i(0) = vit->point().x();
//       col_i(1) = vit->point().y();
//       col_i(2) = vit->point().z();
//       vertices(Rcpp::_, i) = col_i;
//       i++;
//     }
//   }

//   // Rcpp::NumericMatrix vertices(3, nvertices);
//   // Rcpp::IntegerVector ids(nvertices);
//   // {
//   //   size_t i = 0;
//   //   for(vertex_descriptor vd : mesh.vertices()) {
//   //     const IPoint3 ivertex = mesh.point(vd);
//   //     Rcpp::NumericVector col_i(3);
//   //     col_i(0) = ivertex.first.x();
//   //     col_i(1) = ivertex.first.y();
//   //     col_i(2) = ivertex.first.z();
//   //     vertices(Rcpp::_, i) = col_i;
//   //     const unsigned id = ivertex.second;
//   //     ids(i) = id;
//   //     i++;
//   //   }
//   // }

//   Rcpp::IntegerMatrix facets(3, nfacets);
//   {
//     size_t i = 0;
//     for(Polyhedron::Facet_iterator fit = mesh.facets_begin();
//         fit != mesh.facets_end(); fit++) {
//       Rcpp::IntegerVector col_i(3);
//       col_i(0) = fit->halfedge()->vertex()->id();
//       col_i(1) = fit->halfedge()->next()->vertex()->id();
//       col_i(2) = fit->halfedge()->opposite()->vertex()->id();
//       facets(Rcpp::_, i) = col_i;
//       i++;
//     }
//   }

//   // Rcpp::IntegerMatrix faces(nfaces, 3);
//   // {
//   //   size_t i = 0;
//   //   for(face_descriptor fa : mesh.faces()) {
//   //     size_t j = 0;
//   //     Rcpp::IntegerVector col_i(3);
//   //     for(vertex_descriptor vd :
//   //         vertices_around_face(mesh.halfedge(fa), mesh)) {
//   //       const IPoint3 ivertex = mesh.point(vd);
//   //       col_i(j) = ivertex.second;
//   //       j++;
//   //     }
//   //     faces(Rcpp::_, i) = col_i;
//   //     i++;
//   //   }
//   // }

//   return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
//                             // Rcpp::Named("ids") = ids,
//                             Rcpp::Named("faces") = facets);
// }

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(const Rcpp::NumericMatrix M,
                   const Rcpp::List L,
                   const bool merge) {
  std::vector<PointT> points = matrix_to_points3<PointT>(M);
  std::vector<std::vector<size_t>> faces = list_to_faces(L);
  bool success =
      CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces);
  if(!success) {
    Rcpp::stop("Polygon orientation failed.");
  }
  if(merge) {
    size_t nremoved =
        CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(
            points, faces);
    Rcpp::Rcout << "Number of points removed: " << nremoved << ".\n";
  }
  MeshT mesh;
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces,
                                                              mesh);
  return mesh;
}

template Mesh3 makeSurfMesh<Mesh3, Point3>(const Rcpp::NumericMatrix, const Rcpp::List, const bool);
template EMesh3 makeSurfMesh<EMesh3, EPoint3>(const Rcpp::NumericMatrix, const Rcpp::List, const bool);

template <typename MeshT>
Rcpp::IntegerMatrix getEdges(MeshT mesh) {
  const size_t nedges = mesh.number_of_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  {
    size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      Rcpp::IntegerVector col_i(2);
      col_i(0) = source(ed, mesh) + 1;
      col_i(1) = target(ed, mesh) + 1;
      Edges(Rcpp::_, i) = col_i;
      i++;
    }
  }
  return Edges;
}

template Rcpp::IntegerMatrix getEdges<Mesh3>(Mesh3);
template Rcpp::IntegerMatrix getEdges<EMesh3>(EMesh3);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::IntegerMatrix getEdges2(MeshT mesh, const double epsilon) {
  const size_t nedges = mesh.number_of_edges();
  Rcpp::IntegerMatrix Edges(3, nedges);
  {
    size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      typename MeshT::Vertex_index s = source(ed, mesh);
      typename MeshT::Vertex_index t = target(ed, mesh);
      Rcpp::IntegerVector col_i(3);
      col_i(0) = (int)s + 1;
      col_i(1) = (int)t + 1;
      std::vector<PointT> points(4);
      points[0] = mesh.point(s);
      points[1] = mesh.point(t);
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      points[2] = mesh.point(mesh.target(mesh.next(h0)));
      typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
      points[3] = mesh.point(mesh.target(mesh.next(h1)));
      bool exterior;
      if(epsilon == 0) {
        exterior = !CGAL::coplanar(points[0], points[1], points[2], points[3]);
      } else {
        typename KernelT::FT vol =
            CGAL::volume(points[0], points[1], points[2], points[3]);
        exterior = CGAL::abs(vol) > epsilon;
        // K::FT svol = CGAL::square(CGAL::volume(points[0], points[1],
        // points[2], points[3])); CGAL::Triangle_3<K> tr1(points[0], points[1],
        // points[2]); K::FT sarea1 = tr1.squared_area(); if(svol <
        // epsilon*sarea1){
        //   CGAL::Triangle_3<K> tr2(points[0], points[1], points[3]);
        //   K::FT sarea2 = tr2.squared_area();
        //   if(svol < epsilon*sarea2){
        //     CGAL::Triangle_3<K> tr3(points[0], points[2], points[3]);
        //     K::FT sarea3 = tr3.squared_area();
        //     if(svol < epsilon*sarea3){
        //       CGAL::Triangle_3<K> tr4(points[1], points[2], points[3]);
        //       K::FT sarea4 = tr4.squared_area();
        //       exterior = svol > epsilon*sarea4;
        //     }else{
        //       exterior = true;
        //     }
        //   }else{
        //     exterior = true;
        //   }
        // }else{
        //   exterior = true;
        // }
      }
      col_i(2) = (int)exterior;
      Edges(Rcpp::_, i) = col_i;
      i++;
    }
  }
  Rcpp::CharacterVector rowNames =
      Rcpp::CharacterVector::create("i1", "i2", "exterior");
  Rcpp::rownames(Edges) = rowNames;
  return Edges;
}

template Rcpp::IntegerMatrix getEdges2<K, Mesh3, Point3>(Mesh3, const double);
template Rcpp::IntegerMatrix getEdges2<EK, EMesh3, EPoint3>(EMesh3, const double);


template <typename KernelT, typename MeshT, typename PointT>
Rcpp::List RSurfMesh(MeshT mesh,
                     const bool isTriangle,
                     const double epsilon,
                     const bool exact) {
  Rcpp::IntegerMatrix Edges;
  if(isTriangle) {
    Edges = getEdges2<KernelT, MeshT, PointT>(mesh, epsilon);
  } else {
    Edges = getEdges<MeshT>(mesh);
  }
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nvertices);
  {
    size_t i = 0;
//    if(exact) {
      for(typename MeshT::Vertex_index vd : mesh.vertices()) {
        Rcpp::NumericVector col_i(3);
        const PointT vertex = mesh.point(vd);
        col_i(0) = CGAL::to_double(vertex.x());
        col_i(1) = CGAL::to_double(vertex.y());
        col_i(2) = CGAL::to_double(vertex.z());
        Vertices(Rcpp::_, i) = col_i;
        i++;
      }
    // } else {
    //   for(Mesh3::Vertex_index vd : mesh.vertices()) {
    //     Rcpp::NumericVector col_i(3);
    //     const Point3 vertex = mesh.point(vd);
    //     col_i(0) = vertex.x();
    //     col_i(1) = vertex.y();
    //     col_i(2) = vertex.z();
    //     Vertices(Rcpp::_, i) = col_i;
    //     i++;
    //   }
    // }
  }
  const size_t nfaces = mesh.number_of_faces();
  Rcpp::List Faces(nfaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(i) = col_i;
      i++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                            Rcpp::Named("edges") = Edges,
                            Rcpp::Named("faces") = Faces);
}

template Rcpp::List RSurfMesh<K, Mesh3, Point3>(Mesh3, const bool, const double, const bool);
template Rcpp::List RSurfMesh<EK, EMesh3, EPoint3>(EMesh3, const bool, const double, const bool);

// Mesh3 Poly2Mesh3(Polyhedron poly) {
//   Mesh3 mesh;
//   CGAL::copy_face_graph(poly, mesh);
//   // auto vnormals = mesh.add_property_map<boost_vertex_descriptor,
//   // Vector3>("v:normals", CGAL::NULL_VECTOR).first; auto fnormals =
//   // mesh.add_property_map<boost_face_descriptor, Vector3>("f:normals",
//   // CGAL::NULL_VECTOR).first;
//   // CGAL::Polygon_mesh_processing::compute_normals(mesh, vnormals,
//   fnormals); return mesh;
// }
