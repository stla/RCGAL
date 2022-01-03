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

#include <CGAL/Vector_2.h>
#include <CGAL/Vector_3.h>

//#include <CGAL/Object.h>

//#include <CGAL/Unique_hash_map.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/squared_distance_2.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

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

typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb2;
typedef CGAL::Triangulation_data_structure_2<Vb2> Tds2;
typedef CGAL::Delaunay_triangulation_2<K, Tds2> DT2;
typedef std::pair<Point2, unsigned> IPoint2;

typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K> Vb3;
typedef CGAL::Triangulation_data_structure_3<Vb3> Tds3;
typedef CGAL::Delaunay_triangulation_3<K, Tds3> DT3;
typedef K::Tetrahedron_3 Tetrahedron3;

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

  Rcpp::IntegerVector ids(npoints);
  Rcpp::IntegerMatrix edges(npoints, 2);
  Rcpp::NumericMatrix chull(npoints, 2);
  // Rcpp::NumericMatrix midpoints(npoints, 2);
  Rcpp::NumericMatrix normals(npoints, 2);
  // Rcpp::NumericVector lengths(npoints);
  double perimeter = 0;
  double surface = 0;
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
Rcpp::List cxhull3d(Rcpp::NumericMatrix pts) {
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
        Rcpp::Rcout << "face0 : " << ivertex.second << "\n";
      }
      const CGAL::Vector_3<K> normal0 =
          CGAL::unit_normal(vpoints0[0], vpoints0[1], vpoints0[2]);
      Rcpp::Rcout << "normal0 : " << normal0 << "\n";

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
        Rcpp::Rcout << "face1 : " << ivertex.second << "\n";
      }
      const CGAL::Vector_3<K> normal1 =
          CGAL::unit_normal(vpoints1[0], vpoints1[1], vpoints1[2]);
      Rcpp::Rcout << "normal1 : " << normal1 << "\n";

      edges(i, 2) = normal0 == normal1 ? 0 : 1;

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

// [[Rcpp::export]]
Rcpp::List del2d(Rcpp::NumericMatrix pts) {
  const size_t npoints = pts.nrow();
  std::vector<IPoint2> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point2(pts(i, 0), pts(i, 1)), i + 1);
  }

  // compute Delaunay mesh
  const DT2 mesh(points.begin(), points.end());

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

  const DT2::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 2);
  {
    size_t i = 0;
    for(DT2::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
        eit++) {
      const std::pair<DT2::Face_handle, int> edge = *eit;
      edges(i, 0) = edge.first->vertex((edge.second + 1) % 3)->info();
      edges(i, 1) = edge.first->vertex((edge.second + 2) % 3)->info();
      i++;
    }
  }

  return Rcpp::List::create(Rcpp::Named("faces") = faces,
                            Rcpp::Named("edges") = edges);
}

// [[Rcpp::export]]
Rcpp::List del3d(Rcpp::NumericMatrix pts) {
  const size_t npoints = pts.nrow();
  std::vector<IPoint3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)), i + 1);
  }

  // compute Delaunay mesh
  const DT3 mesh(points.begin(), points.end());

  const size_t nfacets = mesh.number_of_finite_facets();
  const size_t ncells = mesh.number_of_finite_cells();

  Rcpp::IntegerMatrix facets(nfacets, 3);
  {
    size_t i = 0;
    for(DT3::Finite_facets_iterator fit = mesh.finite_facets_begin();
        fit != mesh.finite_facets_end(); fit++) {
      std::pair<DT3::Cell_handle, int> facet = *fit;
      DT3::Vertex_handle v0 = facet.first->vertex((facet.second + 1) % 4);
      DT3::Vertex_handle v1 = facet.first->vertex((facet.second + 2) % 4);
      DT3::Vertex_handle v2 = facet.first->vertex((facet.second + 3) % 4);
      facets(i, 0) = v0->info();
      facets(i, 1) = v1->info();
      facets(i, 2) = v2->info();
      i++;
    }
  }

  Rcpp::IntegerMatrix cells(ncells, 4);
  Rcpp::NumericVector volumes(ncells);
  {
    size_t i = 0;
    for(DT3::Finite_cells_iterator cit = mesh.finite_cells_begin();
        cit != mesh.finite_cells_end(); cit++) {
      cells(i, 0) = cit->vertex(0)->info();
      cells(i, 1) = cit->vertex(1)->info();
      cells(i, 2) = cit->vertex(2)->info();
      cells(i, 3) = cit->vertex(3)->info();
      Tetrahedron3 th = CGAL::Tetrahedron_3<K>(
          cit->vertex(0)->point(), cit->vertex(1)->point(),
          cit->vertex(2)->point(), cit->vertex(3)->point());
      volumes(i) = th.volume();
      i++;
    }
    cells.attr("volumes") = volumes;
  }

  return Rcpp::List::create(Rcpp::Named("cells") = cells,
                            Rcpp::Named("facets") = facets);
}
