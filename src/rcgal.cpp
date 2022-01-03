// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
//#define CGAL_EIGEN3_ENABLED
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/property_map.h>

#include <CGAL/convex_hull_3.h>

#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Kernel/global_functions.h>

#include <CGAL/Vector_3.h>

//#include <CGAL/Object.h>

//#include <CGAL/Unique_hash_map.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/squared_distance_2.h>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <array>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2 Point2;
typedef K::Point_3 Point3;
typedef std::vector<Point2> Points2;
typedef std::pair<Point3, unsigned> IPoint3;

typedef CGAL::
    Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point2>::type>
        CHT2;

typedef CGAL::First_of_pair_property_map<IPoint3> Pmap;
typedef CGAL::Extreme_points_traits_adapter_3<Pmap,
                                              CGAL::Convex_hull_traits_3<K>>
    CHT;

typedef CGAL::Surface_mesh<IPoint3> Mesh;

typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Edge_index edge_descriptor;
typedef Mesh::Face_index face_descriptor;

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> DT2;
typedef std::pair<Point2, unsigned> IPoint2;

// [[Rcpp::export]]
Rcpp::List cxhull2d(Rcpp::NumericMatrix pts) {
  Points2 points;
  for(int i = 0; i < pts.nrow(); i++) {
    points.push_back(Point2(pts(i, 0), pts(i, 1)));
  }

  std::vector<std::size_t> indices(points.size()), out;
  std::iota(indices.begin(), indices.end(), 0);
  CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(out),
                      CHT2(CGAL::make_property_map(points)));

  const size_t npoints = out.size();

  Rcpp::IntegerMatrix ids(npoints, 2);
  Rcpp::NumericMatrix chull(npoints, 2);
  double perimeter = 0;
  Point2 pt0 = points[out[npoints - 1]];
  for(size_t i = 0; i < npoints; i++) {
    const size_t id = out[i];
    ids[i] = id;
    Point2 pt1 = points[id];
    perimeter += sqrt(CGAL::squared_distance(pt0, pt1));
    chull(i, 0) = pt1.x();
    chull(i, 1) = pt1.y();
    pt0 = pt1;
  }

  return Rcpp::List::create(Rcpp::Named("verticesIds") = ids,
                            Rcpp::Named("vertices") = chull,
                            Rcpp::Named("perimeter") = perimeter);
}

// [[Rcpp::export]]
Rcpp::List cxhull3d(Rcpp::NumericMatrix pts) {
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
      CGAL::Vector_3<K> normal0 =
          CGAL::normal(vpoints0[0], vpoints0[1], vpoints0[2]);
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
      CGAL::Vector_3<K> normal1 =
          CGAL::normal(vpoints1[0], vpoints1[1], vpoints1[2]);
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
Rcpp::List del2d(Rcpp::NumericMatrix pts) {
  size_t npoints = pts.nrow();
  std::vector<IPoint2> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point2(pts(i, 0), pts(i, 1)), i + 1);
  }

  // compute Delaunay mesh
  DT2 mesh(points.begin(), points.end());

  const size_t nfaces = mesh.number_of_faces();
  const size_t h =
      2 * npoints - 2 - nfaces;  // number of vertices of convex hull
  const size_t nedges = 3 * npoints - 3 - h;

  Rcpp::IntegerMatrix faces(nfaces, 3);
  {
    size_t i = 0;
    for(DT2::Finite_faces_iterator fit = mesh.finite_faces_begin();
        fit != mesh.finite_faces_end(); fit++) {
      faces(i, 0) = fit->vertex(0)->info();
      faces(i, 1) = fit->vertex(1)->info();
      faces(i, 2) = fit->vertex(2)->info();
      i++;
    }
  }

  Rcpp::Rcout << "**********************************"
              << "\n";

  DT2::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 2);
  {
    size_t i = 0;
    for(DT2::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
        eit++) {
      std::pair<DT2::Face_handle, int> edge = *eit;
      edges(i, 0) = edge.first->vertex((edge.second + 1) % 3)->info();
      edges(i, 1) = edge.first->vertex((edge.second + 2) % 3)->info();
      i++;
    }
  }

  Rcpp::List out = Rcpp::List::create(Rcpp::Named("faces") = faces,
                                      Rcpp::Named("edges") = edges);

  return out;
}
