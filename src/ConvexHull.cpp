// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
#define CGAL_EIGEN3_ENABLED
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>

#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Kernel/global_functions.h>

#include <CGAL/Vector_3.h>

#include <CGAL/Object.h>

//#include <CGAL/Unique_hash_map.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <array>

// #include <CGAL/Polygon_mesh_processing/locate.h>
// #include <CGAL/AABB_traits.h>
// #include <CGAL/boost/graph/helpers.h>
// #include <CGAL/Dynamic_property_map.h>

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
//typedef CGAL::Surface_mesh<Point3> Mesh0;
// typedef Surface_mesh::Vertex_range Vertex_range;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Edge_index edge_descriptor;
typedef Mesh::Face_index face_descriptor;

typedef CGAL::Delaunay_triangulation_2<K> DT2;
typedef std::pair<Point2, unsigned> IPoint2;
typedef CGAL::First_of_pair_property_map<IPoint2> Pmap2;
// typedef CGAL::Extreme_points_traits_adapter_2<Pmap,
//                                               CGAL::Convex_hull_traits_3<K> >
//   CHT;
typedef CGAL::Surface_mesh<Point2> Mesh2;

// typedef K::FT FT;
// namespace PMP = CGAL::Polygon_mesh_processing;
// typedef PMP::Barycentric_coordinates<FT>                                Barycentric_coordinates;
// typedef PMP::Face_location<Mesh, FT>                                    Face_location;
// typedef typename boost::graph_traits<Mesh>::face_descriptor             face_descro;

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
Rcpp::NumericMatrix cxhull3d0(Rcpp::NumericMatrix pts) {
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
Rcpp::List cxhull3d(Rcpp::NumericMatrix pts) {

  size_t npoints = pts.nrow();
//  std::vector<Point3> points0(npoints);
  std::vector<IPoint3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
//    points0[i] = Point3(pts(i, 0), pts(i, 1), pts(i, 2));
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)), i + 1);
  }

  // define surface mesh to hold convex hull
  Mesh mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh, CHT(Pmap()));

  //Mesh0 mesh0;
  //CGAL::convex_hull_3(points0.begin(), points0.end(), mesh0);

  size_t nvertices = num_vertices(mesh);  // or number_of_vertices
  size_t nedges = num_edges(mesh);        // or number_of_edges
  size_t nfaces = num_faces(mesh);

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

  Rcpp::IntegerMatrix edges(nedges, 3);
  {
    size_t i = 0;
    for(edge_descriptor ed : mesh.edges()) {
      vertex_descriptor s = source(ed, mesh);
      vertex_descriptor t = target(ed, mesh);
      edges(i, 0) = ids[s];
      edges(i, 1) = ids[t];

      Mesh::Halfedge_index h0 = mesh.halfedge(ed, 0);
      Mesh::Face_index face0 = mesh.face(h0);
      // face_descro f0 = mesh.face(h0);
      //
      // Point3 source_pt = mesh.point(s).first;
      // Point3 target_pt = mesh.point(t).first;
      // const Point3 middle_edge = CGAL::midpoint(source_pt, target_pt);
      // Rcpp::Rcout << "SSSSSSSSSSSSS : " << source_pt.x() << "\n";
      // //face_descriptor fd = mesh.face(mesh.edge(ed));
      // Face_location query_location = PMP::locate(middle_edge, mesh0);
      // Rcpp::Rcout << "Is it on the face's border? " << (PMP::is_on_mesh_border(query_location, mesh0) ? "Yes" : "No") << "\n\n";


      std::array<unsigned, 3> vids0;
      std::array<Point3, 3> vpoints0;
      unsigned counter0 = 0;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(face0), mesh)) {
        IPoint3 ivertex = mesh.point(vd);
        vids0[counter0] = ivertex.second;
        vpoints0[counter0] = ivertex.first;
        counter0++;
        Rcpp::Rcout << "face0 : " << ivertex.second << "\n";
      }
      CGAL::Vector_3<K> normal0 = CGAL::normal(
        vpoints0[0], vpoints0[1], vpoints0[2]
      );
      Rcpp::Rcout << "normal0 : " << normal0 << "\n";

      Mesh::Halfedge_index h1 = mesh.halfedge(ed, 1);
      Mesh::Face_index face1 = mesh.face(h1);
      std::array<unsigned, 3> vids1;
      std::array<Point3, 3> vpoints1;
      unsigned counter1 = 0;
      for(vertex_descriptor vd :
            vertices_around_face(mesh.halfedge(face1), mesh)) {
        IPoint3 ivertex = mesh.point(vd);
        vids1[counter1] = ivertex.second;
        vpoints1[counter1] = ivertex.first;
        counter1++;
        Rcpp::Rcout << "face1 : " << ivertex.second << "\n";
      }
      CGAL::Vector_3<K> normal1 = CGAL::normal(
        vpoints1[0], vpoints1[1], vpoints1[2]
      );
      Rcpp::Rcout << "normal1 : " << normal1 << "\n";

      Rcpp::Rcout << "directionnormal1 : " << normal1.direction() << "\n";

      edges(i, 2) = normal0 == normal1 ? 0 : 1;

      i++;
    }
  }

  Rcpp::IntegerMatrix faces(nfaces, 3);
  {
    size_t i = 0;
    for(face_descriptor fa : mesh.faces()) {
      size_t j = 0;
      for(vertex_descriptor vd :
          vertices_around_face(mesh.halfedge(fa), mesh)) {
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


// [[Rcpp::export]]
Rcpp::List cxhull3dObj(Rcpp::NumericMatrix pts) {

  size_t npoints = pts.nrow();
  std::vector<Point3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = Point3(pts(i, 0), pts(i, 1), pts(i, 2));
  }

  // define object to hold convex hull
  CGAL::Object obj;
  CGAL::convex_hull_3(points.begin(), points.end(), obj);
  const Polyhedron3* mesh = CGAL::object_cast<Polyhedron3>(&obj);
  Rcpp::Rcout << "Polyhedron\n " << *mesh << "\n";

  Rcpp::IntegerMatrix faces(2, 3);

  Rcpp::List out = Rcpp::List::create(//Rcpp::Named("nedges") = edges.size(),
                                      Rcpp::Named("faces") = faces);

  return out;
}


// [[Rcpp::export]]
Rcpp::List del2d(Rcpp::NumericMatrix pts) {

  size_t npoints = pts.nrow();
  std::vector<Point2> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = Point2(pts(i, 0), pts(i, 1));
  }

  // define surface mesh to hold convex hull
  DT2 mesh(points.begin(), points.end());

  for (DT2::Finite_faces_iterator it = mesh.finite_faces_begin(); it != mesh.finite_faces_end(); it++)
  {
    Rcpp::Rcout << it->vertex(0)->point() << " -- ";
    Rcpp::Rcout << it->vertex(1)->point() << " -- ";
    Rcpp::Rcout << it->vertex(2)->point() << "\n";
  }


  Rcpp::IntegerMatrix faces(2, 3);

  Rcpp::List out = Rcpp::List::create(//Rcpp::Named("nedges") = edges.size(),
    Rcpp::Named("faces") = faces);

  return out;
}
