// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/Surface_mesh.h>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron3;
typedef Polyhedron3::Vertex_iterator Vertex_iterator;
typedef K::Point_2 Point2;
typedef K::Point_3 Point3;
typedef std::vector<Point2> Points2;

typedef std::pair<Point3, unsigned> IPoint3;
typedef CGAL::First_of_pair_property_map<IPoint3> Pmap;
typedef CGAL::Extreme_points_traits_adapter_3<Pmap,
                                              CGAL::Convex_hull_traits_3<K> >
    CHT;
typedef CGAL::Surface_mesh<IPoint3> Mesh;
// typedef Surface_mesh::Vertex_range Vertex_range;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Edge_index edge_descriptor;
typedef Mesh::Face_index face_descriptor;

// [[Rcpp::export]]
Rcpp::NumericMatrix test() {
  Points2 points, result;
  points.push_back(Point2(0, 0));
  points.push_back(Point2(10, 0));
  points.push_back(Point2(10, 10));
  points.push_back(Point2(6, 5));
  points.push_back(Point2(4, 1));
  CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result));

  size_t npoints = result.size();

  Rcpp::NumericMatrix chull(npoints, 2);
  for(size_t i = 0; i < npoints; i++) {
    chull(i, 0) = result[i].x();
    chull(i, 1) = result[i].y();
  }
  return chull;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cxhull2d(Rcpp::NumericMatrix pts) {
  Points2 points, result;

  for(int i = 0; i < pts.nrow(); i++) {
    points.push_back(Point2(pts(i, 0), pts(i, 1)));
  }

  CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result));

  size_t npoints = result.size();

  Rcpp::NumericMatrix chull(npoints, 2);
  for(size_t i = 0; i < npoints; i++) {
    chull(i, 0) = result[i].x();
    chull(i, 1) = result[i].y();
  }
  return chull;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cxhull3d(Rcpp::NumericMatrix pts) {
  std::vector<Point3> points;
  // define polyhedron to hold convex hull
  Polyhedron3 poly;

  for(int i = 0; i < pts.nrow(); i++) {
    points.push_back(Point3(pts(i, 0), pts(i, 1), pts(i, 2)));
  }

  CGAL::convex_hull_3(points.begin(), points.end(), poly);

  size_t npoints = poly.size_of_vertices();

  Rcpp::NumericMatrix chull(npoints, 3);
  int i = 0;
  for(Vertex_iterator v = poly.vertices_begin(); v != poly.vertices_end();
      ++v) {
    chull(i, 0) = v->point().x();
    chull(i, 1) = v->point().y();
    chull(i, 2) = v->point().z();
    i++;
  }
  return chull;
}

// [[Rcpp::export]]
Rcpp::List cxhull3d2(Rcpp::NumericMatrix pts) {
  size_t npoints = pts.nrow();
  std::vector<IPoint3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)), i + 1);
  }

  // define surface mesh to hold convex hull
  Mesh mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh, CHT(Pmap()));

  size_t nvertices = num_vertices(mesh);  // or number_of_vertices
  size_t nedges = num_edges(mesh);        // or number_of_edges
  size_t nfaces = num_faces(mesh);

  // Rcpp::NumericMatrix chull(nvertices, 3);

  //  Surface_mesh::Vertex_range::iterator  vb, ve;
  //  Surface_mesh::Vertex_range r = sm.vertices();
  // The iterators can be accessed through the C++ range API
  //  vb = r.begin();
  //  ve = r.end();
  // or the boost Range API
  // vb = boost::begin(r);
  // ve = boost::end(r);
  // or with boost::tie, as the CGAL range derives from std::pair
  //  for(boost::tie(vb, ve) = sm.vertices(); vb != ve; ++vb){
  //    Rcpp::Rcout << *vb << "\n";
  //  }

  Rcpp::List vertices(nvertices);
  std::vector<unsigned> ids(nvertices);
  {
    size_t i = 0;
    for(vertex_descriptor vd : mesh.vertices()) {
      IPoint3 ivertex = mesh.point(vd);
      Rcpp::NumericVector vertex(3);
      vertex[0] = ivertex.first.x();
      vertex[1] = ivertex.first.y();
      vertex[2] = ivertex.first.z();
      unsigned id = ivertex.second;
      ids[i] = id;
      Rcpp::List L = Rcpp::List::create(Rcpp::Named("id") = id,
                                        Rcpp::Named("point") = vertex);
      vertices[i] = L;
      i++;
    }
  }
  Rcpp::IntegerMatrix edges(nedges, 2);
  {
    size_t i = 0;
    for(edge_descriptor ed : mesh.edges()) {
      vertex_descriptor s = source(ed, mesh);
      vertex_descriptor t = target(ed, mesh);
      edges(i, 0) = ids[s];
      edges(i, 1) = ids[t];
      i++;
    }
  }
  //  int i = 0;
  //  Vertex_range vrange = sm.vertices();
  //  for(auto v = vrange.begin(); v != vrange.end();
  //      ++v) {
  //    chull(i, 0) = v.point().x();
  //    chull(i, 1) = v->point().y();
  //    chull(i, 2) = v->point().z();
  //    Rcpp::Rcout << "info : " << v(0)->info() << "\n";
  //    i++;
  //  }
  Rcpp::IntegerMatrix faces(nfaces, 3);
  {
    size_t i = 0;
    for(face_descriptor fa : mesh.faces()) {
      size_t j = 0
      for(vertex_descriptor vd : vertices_around_face(mesh.halfedge(fa), mesh)) {
        IPoint3 ivertex = mesh.point(vd);
        faces(i, j) = ivertex.second;
        j++;
      }
      i++;
    }
  }

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                                      Rcpp::Named("edges") = edges,
                                      Rcpp::Named("faces") = faces);

  return out;
}
