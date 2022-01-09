// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
//#define CGAL_EIGEN3_ENABLED

// #include <CGAL/assertions.h>
// #undef CGAL_error
// #define CGAL_error
// #undef CGAL_error_msg
// #define CGAL_error_msg(msg)

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

#include <CGAL/utility.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/tuple.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/poisson_surface_reconstruction.h>

#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>

#include <CGAL/Projection_traits_xy_3.h>



#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <array>
#include <map>
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

typedef CGAL::Advancing_front_surface_reconstruction<> AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3 AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;
typedef K::Vector_3 Vector3;

typedef CGAL::Simple_cartesian<double> KSC;

typedef std::pair<Point3, Vector3> P3wn;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

typedef CGAL::Projection_traits_xy_3<K> Pxy;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Pxy> Vb_xy;
typedef CGAL::Triangulation_data_structure_2<Vb_xy> Tds_xy;
typedef CGAL::Delaunay_triangulation_2<Pxy, Tds_xy> Delaunay_xy;

typedef Rcpp::NumericVector Dvector;

Rcpp::String stringTriple(const unsigned i,
                          const unsigned j,
                          const unsigned k) {
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
    std::to_string(i), "-", std::to_string(j), "-", std::to_string(k));
  return Rcpp::collapse(stringids);
}

Rcpp::String stringPair(const unsigned i,
                          const unsigned j) {
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
    std::to_string(i), "-", std::to_string(j));
  return Rcpp::collapse(stringids);
}

// [[Rcpp::export]]
Rcpp::List del3d_cpp(Rcpp::NumericMatrix pts) {
  const size_t npoints = pts.nrow();
  std::vector<IPoint3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)), i + 1);
  }
  
  // compute Delaunay mesh
  const DT3 mesh(points.begin(), points.end());
  
  const size_t nfacets = mesh.number_of_finite_facets();
  const size_t ncells = mesh.number_of_finite_cells();
  const size_t nedges = mesh.number_of_finite_edges();
  
  //const DT3::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 3);
  std::map<Rcpp::String, size_t> edgesMap = {};
  {
    size_t i = 0;
    const DT3::Finite_edges itedges = mesh.finite_edges();
    for(DT3::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end(); eit++) {
      const DT3::Edge edge = *eit;
      const DT3::Vertex::Info v1_info = edge.first->vertex(edge.second)->info();
      const DT3::Vertex::Info v2_info = edge.first->vertex(edge.third)->info();
      const unsigned imin = v1_info < v2_info ? v1_info : v2_info;
      const unsigned imax = v1_info > v2_info ? v1_info : v2_info;
      //Rcpp::Rcout << i << "@@@@@" << imin << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ";
      edges(i, 0) = imin;
      edges(i, 1) = imax;
      const Rcpp::String edgeAsString = stringPair(imin, imax);
      i++;
      edgesMap[edgeAsString] = i;
    }
    // for(DT3::Finite_edges_iterator eit = mesh.finite_edges_begin(); eit != mesh.finite_edges_end();
    // eit++) {
    //   //const CGAL::Triple<DT3::Cell_handle, int, int> edge = *eit;
    //   unsigned i0 = eit->first->vertex(eit->second % 4)->info();
    //   unsigned i1 = eit->first->vertex(eit->third % 4)->info();
    //   unsigned imin = i0 < i1 ? i0 : i1;
    //   unsigned imax = i0 > i1 ? i0 : i1;
    //   Rcpp::Rcout << i << "@@@@@" << imin << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ";
    //   edges(i, 0) = imin;
    //   edges(i, 1) = imax;
    //   const Rcpp::String edgeAsString = stringPair(imin, imax);
    //   i++;
    //   edgesMap[edgeAsString] = i;
    // }
  }
  
  std::map<Rcpp::String, size_t> facetsMap = {};
  Rcpp::IntegerMatrix facets(nfacets, 4);
  //Imatrix facetsOnHull(0, 3);
  {
    size_t i = 0;
    for(DT3::Finite_facets_iterator fit = mesh.finite_facets_begin();
        fit != mesh.finite_facets_end(); fit++) {
      std::pair<DT3::Cell_handle, int> facet = *fit;
      DT3::Vertex_handle v0 = facet.first->vertex((facet.second + 1) % 4);
      DT3::Vertex_handle v1 = facet.first->vertex((facet.second + 2) % 4);
      DT3::Vertex_handle v2 = facet.first->vertex((facet.second + 3) % 4);
      const unsigned id0 = v0->info();
      const unsigned id1 = v1->info();
      const unsigned id2 = v2->info();
      bool onhull = mesh.is_infinite(facet.first) ||
        mesh.is_infinite(mesh.mirror_facet(facet).first);
      facets(i, 3) = onhull;
      std::array<unsigned, 3> ids = {id0, id1, id2};
      std::sort(ids.begin(), ids.end());
      facets(i, 0) = ids[0];
      facets(i, 1) = ids[1];
      facets(i, 2) = ids[2];
      // if(onhull){
      //   bool notcoplanar = true;
      //   unsigned j = 0;
      //   while(notcoplanar && j < nedges){
      //     std::vector<unsigned> v_intersection;
      //     std::array<unsigned, 2> edge = {edges(j, 0), edges(j, 1)};
      //     std::set_intersection(ids.begin(), ids.end(),
      //                           edge.begin(), edge.end(),
      //                           std::back_inserter(v_intersection));
      //     if(v_intersection.size() == 2){
      //
      //     }
      //   }
      // Imatrix.conservativeResize(Imatrix.rows() + 1, Eigen::NoChange);
      // Ivector ivec << id0, id1, id2;
      // Imatrix.bottomRows(1) = ivec;
      //}
      const Rcpp::String facetAsString = stringTriple(ids[0], ids[1], ids[2]);
      i++;
      facetsMap[facetAsString] = i;
    }
    Rcpp::CharacterVector columnNames =
      Rcpp::CharacterVector::create("i1", "i2", "i3", "onhull");
    Rcpp::colnames(facets) = columnNames;
    facets.attr("info") =
      "The `onhull` column indicates whether the face is on the convex hull.";
  }
//  {
    for(size_t i = 0; i < nfacets-1; i++){
      Rcpp::IntegerVector facet_i = facets(i, Rcpp::_);
      if(facet_i(3)){
        facet_i = Rcpp::IntegerVector::create(facet_i(0), facet_i(1), facet_i(2));
        for(size_t j = i+1; j < nfacets; j++){
          Rcpp::IntegerVector facet_j = facets(j, Rcpp::_);
          if(facet_j(3)){
            facet_j = Rcpp::IntegerVector::create(facet_j(0), facet_j(1), facet_j(2));
            std::vector<int> v_intersection;
            std::set_intersection(facet_i.begin(), facet_i.end(),
                                  facet_j.begin(), facet_j.end(),
                                  std::back_inserter(v_intersection));
            //Rcpp::Rcout << v_intersection.size() << "******\n";
            if(v_intersection.size() == 2){
              int v0 = v_intersection[0];
              int v1 = v_intersection[1];

              Rcpp::Rcout << v0 << " +++ ";
              Rcpp::Rcout << v1 << "\n";
              
              Rcpp::Rcout << facet_i << " --- ";
              Rcpp::Rcout << facet_j << "\n";
              std::vector<int> v_union;
              std::set_union(facet_i.begin(), facet_i.end(),
                                    facet_j.begin(), facet_j.end(),
                                    std::back_inserter(v_union));
              Point3 p0 = points[v_union[0]-1].first;
              Point3 p1 = points[v_union[1]-1].first;
              Point3 p2 = points[v_union[2]-1].first;
              Point3 p3 = points[v_union[3]-1].first;
              bool coplanar = CGAL::coplanar(p0, p1, p2, p3);
              Rcpp::Rcout << "coplanar" << coplanar << "\n";
              if(coplanar){
                Rcpp::Rcout << "COPLANAR" << "\n";
                edges(edgesMap[stringPair(v0, v1)], 2) = 0;
              }else{
                Rcpp::Rcout << "RRRAAAAAAAAAAAAAAAAAAAA" << "\n";
                Rcpp::Rcout << edgesMap[stringPair(v0, v1)] << "\n";
                edges(edgesMap[stringPair(v0, v1)], 2) = 1;
              }
            }
          }
        }
      }
    }
//  }
  
  Rcpp::List cells(ncells);
  double totalVolume = 0.0;
  {
    size_t i = 0;
    for(DT3::Finite_cells_iterator cit = mesh.finite_cells_begin();
        cit != mesh.finite_cells_end(); cit++) {
      const unsigned id0 = cit->vertex(0)->info();
      const unsigned id1 = cit->vertex(1)->info();
      const unsigned id2 = cit->vertex(2)->info();
      const unsigned id3 = cit->vertex(3)->info();
      Rcpp::IntegerVector cell(4);
      cell(0) = id0;
      cell(1) = id1;
      cell(2) = id2;
      cell(3) = id3;
      std::array<unsigned, 4> ids = {id0, id1, id2, id3};
      std::sort(ids.begin(), ids.end());
      const Rcpp::IntegerVector faces = Rcpp::IntegerVector::create(
        facetsMap[stringTriple(ids[0], ids[1], ids[2])],
                 facetsMap[stringTriple(ids[0], ids[1], ids[3])],
                          facetsMap[stringTriple(ids[0], ids[2], ids[3])],
                                   facetsMap[stringTriple(ids[1], ids[2], ids[3])]);
      Tetrahedron3 th = CGAL::Tetrahedron_3<K>(
        cit->vertex(0)->point(), cit->vertex(1)->point(),
        cit->vertex(2)->point(), cit->vertex(3)->point());
      const double volume = th.volume();
      totalVolume += volume;
      cells[i] = Rcpp::List::create(Rcpp::Named("cell") = cell,
                                    Rcpp::Named("faces") = faces,
                                    Rcpp::Named("volume") = volume);
      i++;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("cells") = cells, Rcpp::Named("facets") = facets,
    Rcpp::Named("edges") = edges, Rcpp::Named("volume") = totalVolume);
}

// double volume_under_triangle(Dvector v0, Dvector v1, Dvector v2){
//   double x0 = v0(0);
//   double y0 = v0(1);
//   double z0 = v0(2);
//   double x1 = v1(0);
//   double y1 = v1(1);
//   double z1 = v1(2);
//   double x2 = v2(0);
//   double y2 = v2(1);
//   double z2 = v2(2);
//   return (z0+z1+z2) * (x0*y1 - x1*y0 + x1*y2 - x2*y1 + x2*y0 - x0*y2) / 6.0;
// }
// 
// // vol <- function(trgl){ # calculates volume under a triangle
// //   with(
// //     trgl, 
// //     sum(z)*(x[1]*y[2]-x[2]*y[1]+x[2]*y[3]-x[3]*y[2]+x[3]*y[1]-x[1]*y[3])/6
// //   )
// // }
// 
// // [[Rcpp::export]]
// Rcpp::List del2d_xy_cpp(Rcpp::NumericMatrix pts) {
//   const size_t npoints = pts.nrow();
// 
//   // compute Delaunay mesh
//   Delaunay_xy mesh;
//   Delaunay_xy::Vertex_handle vh;
//   for(size_t i = 0; i < npoints; i++) {
//     vh = mesh.insert(Point3(pts(i, 0), pts(i, 1), pts(i, 2)));
//     vh->info() = i;
//   }
// 
//   const size_t nfaces = mesh.number_of_faces();
//   Rcpp::Rcout << nfaces << "\n";
//   const size_t h =
//     2 * npoints - 2 - nfaces;  // number of vertices of convex hull
//   const size_t nedges = 3 * npoints - 3 - h;
//   
//   Rcpp::IntegerMatrix faces(nfaces, 3);
//   Rcpp::NumericVector volumes(nfaces);
//   double totalVolume = 0.0;
//   {
//     size_t i = 0;
//     for(Delaunay_xy::Finite_faces_iterator fit = mesh.finite_faces_begin();
//         fit != mesh.finite_faces_end(); fit++) {
//       unsigned i0 = fit->vertex(0)->info();
//       unsigned i1 = fit->vertex(1)->info();
//       unsigned i2 = fit->vertex(2)->info();
//       faces(i, 0) = i0 + 1;
//       faces(i, 1) = i1 + 1;
//       faces(i, 2) = i2 + 1;
//       double volume = volume_under_triangle(
//         pts(i0, Rcpp::_), pts(i1, Rcpp::_), pts(i2, Rcpp::_)
//       );
//       volumes(i) = volume;
//       totalVolume += volume;
//       i++;
//     }
//     faces.attr("volumes") = volumes;
//   }
//   
//   const Delaunay_xy::Finite_edges itedges = mesh.finite_edges();
//   Rcpp::IntegerMatrix edges(nedges, 2);
//   {
//     size_t i = 0;
//     for(Delaunay_xy::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
//     eit++) {
//       const std::pair<Delaunay_xy::Face_handle, int> edge = *eit;
//       edges(i, 0) = edge.first->vertex((edge.second + 1) % 3)->info();
//       edges(i, 1) = edge.first->vertex((edge.second + 2) % 3)->info();
//       i++;
//     }
//   }
//   
//   return Rcpp::List::create(Rcpp::Named("faces") = faces,
//                             Rcpp::Named("edges") = edges,
//                             Rcpp::Named("volume") = totalVolume);
// }

// 
// // [[Rcpp::export]]
// Rcpp::NumericMatrix compute_normals_cpp(Rcpp::NumericMatrix pts,
//                                         unsigned nb_neighbors) {
//   const size_t npoints = pts.nrow();
//   std::vector<P3wn> points(npoints);
//   for(size_t i = 0; i < npoints; i++) {
//     points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)),
//                                Vector3(0.0, 0.0, 0.0));
//   }
// 
//   CGAL::jet_estimate_normals<Concurrency_tag>(
//       points, nb_neighbors,
//       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
//           .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));
// 
//   CGAL::mst_orient_normals(
//       points, nb_neighbors,
//       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
//           .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));
// 
//   Rcpp::NumericMatrix normals(npoints, 3);
//   for(size_t i = 0; i < npoints; i++) {
//     const Vector3 normal = points[i].second;
//     normals(i, 0) = normal.x();
//     normals(i, 1) = normal.y();
//     normals(i, 2) = normal.z();
//   }
// 
//   return normals;
// }
// 
// // [[Rcpp::export]]
// Rcpp::NumericMatrix compute_normals_cpp2(Rcpp::NumericMatrix pts,
//                                          unsigned nb_neighbors) {
//   const size_t npoints = pts.nrow();
//   std::vector<P3wn> points(npoints);
//   for(size_t i = 0; i < npoints; i++) {
//     points[i] = std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)),
//                                Vector3(0.0, 0.0, 0.0));
//   }
// 
//   CGAL::pca_estimate_normals<Concurrency_tag>(
//       points, nb_neighbors,
//       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
//           .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));
// 
//   CGAL::mst_orient_normals(
//       points, nb_neighbors,
//       CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
//           .normal_map(CGAL::Second_of_pair_property_map<P3wn>()));
// 
//   Rcpp::NumericMatrix normals(npoints, 3);
//   for(size_t i = 0; i < npoints; i++) {
//     const Vector3 normal = points[i].second;
//     normals(i, 0) = normal.x();
//     normals(i, 1) = normal.y();
//     normals(i, 2) = normal.z();
//   }
// 
//   return normals;
// }