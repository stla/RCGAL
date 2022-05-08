#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// double sepsilon = sqrt(std::numeric_limits<double>::epsilon());

// [[Rcpp::export]]
Rcpp::List cxhull2d_cpp(Rcpp::NumericMatrix pts) {
  Points2 points;
  for(int i = 0; i < pts.nrow(); i++) {
    points.push_back(Point2(pts(i, 0), pts(i, 1)));
  }

  std::vector<std::size_t> indices(points.size()), out;
  std::iota(indices.begin(), indices.end(), 0);
  CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                      CHT2(CGAL::make_property_map(points)));

  const size_t npoints = out.size();

  Rcpp::IntegerVector ids(npoints);
  Rcpp::IntegerMatrix edges(npoints, 2);
  Rcpp::NumericMatrix chull(npoints, 2);
  // Rcpp::NumericMatrix midpoints(npoints, 2);
  Rcpp::NumericMatrix normals(npoints, 2);
  // Rcpp::NumericVector lengths(npoints);
  double perimeter = 0.0;
  double surface = 0.0;
  size_t id0 = out[npoints - 1];
  Point2 pt0 = points[id0];
  for(size_t i = 0; i < npoints; i++) {
    const size_t id1 = out[i];
    ids[i] = id1 + 1;
    edges(i, 0) = id0 + 1;
    edges(i, 1) = id1 + 1;
    Point2 pt1 = points[id1];
    const double l = sqrt(CGAL::squared_distance(pt0, pt1));
    // lengths(i) = l;
    perimeter += l;
    const Point2 middle = CGAL::midpoint(pt0, pt1);
    const std::array<double, 2> center = {middle.x(), middle.y()};
    // midpoints(i, 0) = middle.x();
    // midpoints(i, 1) = middle.y();
    const CGAL::Vector_2<K> vedge = pt0 - pt1;
    const Rcpp::NumericVector normal =
      Rcpp::NumericVector::create(-vedge.y(), vedge.x());
    normals(i, Rcpp::_) = normal / l;
    chull(i, 0) = pt1.x();
    chull(i, 1) = pt1.y();
    pt0 = pt1;
    id0 = id1;
    surface +=
      std::inner_product(center.begin(), center.end(), normal.begin(), 0.0);
  }
  surface /= 2.0;
  edges.attr("normals") = normals;

  return Rcpp::List::create(
    Rcpp::Named("verticesIds") = ids, Rcpp::Named("vertices") = chull,
    Rcpp::Named("edges") = edges, Rcpp::Named("surface") = surface,
    Rcpp::Named("perimeter") = perimeter);
}

// [[Rcpp::export]]
Rcpp::List cxhull3d_cpp(Rcpp::NumericMatrix pts, double epsilon) {
  const size_t npoints = pts.nrow();
  std::vector<IPoint3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)), i + 1);
  }

  // define surface mesh to hold convex hull
  Mesh mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh, CHT(Pmap()));

  const size_t nvertices = num_vertices(mesh);  // or number_of_vertices
  const size_t nedges = num_edges(mesh);        // or number_of_edges
  const size_t nfaces = num_faces(mesh);

  Rcpp::List vertices(nvertices);
  std::vector<unsigned> ids(nvertices);
  {
    size_t i = 0;
    for(vertex_descriptor vd : mesh.vertices()) {
      const IPoint3 ivertex = mesh.point(vd);
      Rcpp::NumericVector vertex(3);
      vertex[0] = ivertex.first.x();
      vertex[1] = ivertex.first.y();
      vertex[2] = ivertex.first.z();
      const unsigned id = ivertex.second;
      ids[i] = id;
      const Rcpp::List L = Rcpp::List::create(Rcpp::Named("id") = id,
                                              Rcpp::Named("point") = vertex);
      vertices[i] = L;
      i++;
    }
  }

  Rcpp::IntegerMatrix edges(3, nedges);
  {
    size_t i = 0;
    for(edge_descriptor ed : mesh.edges()) {
      Rcpp::IntegerVector edge_i(3);
      const vertex_descriptor v1 = source(ed, mesh);
      const vertex_descriptor v2 = target(ed, mesh);
      edge_i(0) = ids[v1];
      edge_i(1) = ids[v2];

      const Mesh::Halfedge_index h0 = mesh.halfedge(ed, 0);
      // const Mesh::Face_index face0 = mesh.face(h0);
      // std::array<unsigned, 3> vids0;
      // std::array<Point3, 3> vpoints0;
      // unsigned counter0 = 0;
      // for(vertex_descriptor vd :
      //       vertices_around_face(mesh.halfedge(face0), mesh)) {
      //   const IPoint3 ivertex = mesh.point(vd);
      //   vids0[counter0] = ivertex.second;
      //   vpoints0[counter0] = ivertex.first;
      //   counter0++;
      // }
      // const CGAL::Vector_3<K> normal0 =
      //   CGAL::unit_normal(vpoints0[0], vpoints0[1], vpoints0[2]);

      // const Mesh::Halfedge_index h1 = mesh.halfedge(ed, 1);
      // const Mesh::Face_index face1 = mesh.face(h1);
      // std::array<unsigned, 3> vids1;
      // std::array<Point3, 3> vpoints1;
      // unsigned counter1 = 0;
      // for(vertex_descriptor vd :
      //       vertices_around_face(mesh.halfedge(face1), mesh)) {
      //   const IPoint3 ivertex = mesh.point(vd);
      //   vids1[counter1] = ivertex.second;
      //   vpoints1[counter1] = ivertex.first;
      //   counter1++;
      // }
      // const CGAL::Vector_3<K> normal1 =
      //   CGAL::unit_normal(vpoints1[0], vpoints1[1], vpoints1[2]);

      // edges(i, 2) = approxEqualVectors(normal0, normal1, epsilon) ? 0 : 1;
      const Mesh::Halfedge_index h1 = mesh.halfedge(ed, 1);
      const vertex_descriptor v3 = mesh.target(mesh.next(h0));
      const vertex_descriptor v4 = mesh.target(mesh.next(h1));
      const Point3 point1 = mesh.point(v1).first;
      const Point3 point2 = mesh.point(v2).first;
      const Point3 point3 = mesh.point(v3).first;
      const Point3 point4 = mesh.point(v4).first;
      bool exterior;
      if(epsilon == 0) {
        exterior = !CGAL::coplanar(point1, point2, point3, point4);
      } else {
        K::FT vol = CGAL::volume(point1, point2, point3, point4);
        exterior = CGAL::abs(vol) > epsilon;
      }
      edge_i(2) = (int)exterior;
      edges(Rcpp::_, i) = edge_i;

      i++;
    }
    Rcpp::CharacterVector rowNames =
      Rcpp::CharacterVector::create("i1", "i2", "exterior");
    Rcpp::rownames(edges) = rowNames;
    // edges.attr("info") = "The `border` column indicates border edges.";
  }

  Rcpp::NumericMatrix circumcenters(3, nfaces);
  Rcpp::NumericMatrix normals(3, nfaces);
  Rcpp::NumericVector areas(nfaces);
  double totalArea = 0;
  double volume = 0;
  Rcpp::IntegerMatrix faces(3, nfaces);
  Rcpp::IntegerMatrix tfaces;
  {
    size_t i = 0;
    Rcpp::IntegerVector face_i(3);
    for(face_descriptor fd : mesh.faces()) {
      size_t j = 0;
      std::array<Point3, 3> fa_vertices;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(fd), mesh)) {
        const IPoint3 ivertex = mesh.point(vd);
        face_i(j) = ivertex.second;
        fa_vertices[j] = ivertex.first;
        j++;
      }
      faces(Rcpp::_, i) = face_i;
      const double area = sqrt(
        CGAL::squared_area(fa_vertices[0], fa_vertices[1], fa_vertices[2]));
      totalArea += area;
      areas(i) = area;
      Rcpp::NumericVector normal_i(3);
      const Vector3 normal =
        CGAL::unit_normal(fa_vertices[0], fa_vertices[1], fa_vertices[2]);
      normal_i(0) = normal.x();
      normal_i(1) = normal.y();
      normal_i(2) = normal.z();
      normals(Rcpp::_, i) = normal_i;
      Rcpp::NumericVector circumcenter_i(3);
      const Point3 circumcenter =
        CGAL::circumcenter(fa_vertices[0], fa_vertices[1], fa_vertices[2]);
      circumcenter_i(0) = circumcenter.x();
      circumcenter_i(1) = circumcenter.y();
      circumcenter_i(2) = circumcenter.z();
      circumcenters(Rcpp::_, i) = circumcenter_i;

      i++;
    }
    tfaces = Rcpp::transpose(faces);
    tfaces.attr("areas") = areas;
    tfaces.attr("normals") = Rcpp::transpose(normals);
    tfaces.attr("circumcenters") = Rcpp::transpose(circumcenters);
    for(size_t k = 0; k < nfaces; ++k) {
      const Rcpp::NumericVector center_k = circumcenters(Rcpp::_, k);
      volume += areas(k) * std::inner_product(center_k.begin(), center_k.end(),
                      normals(Rcpp::_, k).begin(), 0.0);
    }
    volume /= 3;
  }

  return Rcpp::List::create(
    Rcpp::Named("vertices") = vertices, 
    Rcpp::Named("edges") = Rcpp::transpose(edges),
    Rcpp::Named("faces") = tfaces, Rcpp::Named("surface") = totalArea,
    Rcpp::Named("volume") = volume);
}

