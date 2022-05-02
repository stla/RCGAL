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

  Rcpp::IntegerMatrix edges(nedges, 3);
  {
    size_t i = 0;
    for(edge_descriptor ed : mesh.edges()) {
      const vertex_descriptor s = source(ed, mesh);
      const vertex_descriptor t = target(ed, mesh);
      edges(i, 0) = ids[s];
      edges(i, 1) = ids[t];

      const Mesh::Halfedge_index h0 = mesh.halfedge(ed, 0);
      const Mesh::Face_index face0 = mesh.face(h0);
      std::array<unsigned, 3> vids0;
      std::array<Point3, 3> vpoints0;
      unsigned counter0 = 0;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(face0), mesh)) {
        const IPoint3 ivertex = mesh.point(vd);
        vids0[counter0] = ivertex.second;
        vpoints0[counter0] = ivertex.first;
        counter0++;
      }
      const CGAL::Vector_3<K> normal0 =
        CGAL::unit_normal(vpoints0[0], vpoints0[1], vpoints0[2]);

      const Mesh::Halfedge_index h1 = mesh.halfedge(ed, 1);
      const Mesh::Face_index face1 = mesh.face(h1);
      std::array<unsigned, 3> vids1;
      std::array<Point3, 3> vpoints1;
      unsigned counter1 = 0;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(face1), mesh)) {
        const IPoint3 ivertex = mesh.point(vd);
        vids1[counter1] = ivertex.second;
        vpoints1[counter1] = ivertex.first;
        counter1++;
      }
      const CGAL::Vector_3<K> normal1 =
        CGAL::unit_normal(vpoints1[0], vpoints1[1], vpoints1[2]);

      edges(i, 2) = approxEqualVectors(normal0, normal1, epsilon) ? 0 : 1;

      i++;
    }
    Rcpp::CharacterVector columnNames =
      Rcpp::CharacterVector::create("i1", "i2", "border");
    Rcpp::colnames(edges) = columnNames;
    edges.attr("info") = "The `border` column indicates border edges.";
  }

  Rcpp::NumericMatrix circumcenters(nfaces, 3);
  Rcpp::NumericMatrix normals(nfaces, 3);
  Rcpp::NumericVector areas(nfaces);
  double totalArea = 0;
  double volume = 0;
  Rcpp::IntegerMatrix faces(nfaces, 3);
  {
    size_t i = 0;
    for(face_descriptor fa : mesh.faces()) {
      size_t j = 0;
      std::array<Point3, 3> fa_vertices;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(fa), mesh)) {
        const IPoint3 ivertex = mesh.point(vd);
        faces(i, j) = ivertex.second;
        fa_vertices[j] = ivertex.first;
        j++;
      }
      const double area = sqrt(
        CGAL::squared_area(fa_vertices[0], fa_vertices[1], fa_vertices[2]));
      totalArea += area;
      areas(i) = area;
      const CGAL::Vector_3<K> normal =
        CGAL::unit_normal(fa_vertices[0], fa_vertices[1], fa_vertices[2]);
      normals(i, 0) = normal.x();
      normals(i, 1) = normal.y();
      normals(i, 2) = normal.z();
      const Point3 circumcenter =
        CGAL::circumcenter(fa_vertices[0], fa_vertices[1], fa_vertices[2]);
      circumcenters(i, 0) = circumcenter.x();
      circumcenters(i, 1) = circumcenter.y();
      circumcenters(i, 2) = circumcenter.z();

      i++;
    }
    faces.attr("areas") = areas;
    faces.attr("normals") = normals;
    faces.attr("circumcenters") = circumcenters;
    for(size_t k = 0; k < nfaces; ++k) {
      const Rcpp::NumericVector center_k = circumcenters(k, Rcpp::_);
      volume += areas(k) * std::inner_product(center_k.begin(), center_k.end(),
                      normals(k, Rcpp::_).begin(), 0.0);
    }
    volume /= 3;
  }

  return Rcpp::List::create(
    Rcpp::Named("vertices") = vertices, Rcpp::Named("edges") = edges,
    Rcpp::Named("faces") = faces, Rcpp::Named("surface") = totalArea,
    Rcpp::Named("volume") = volume);
}

