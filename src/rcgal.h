#ifndef _RCGALHEADER_
#define _RCGALHEADER_
#endif

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

// #include <iostream>

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

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <array>
#include <limits>
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

struct FaceInfo2 {
  FaceInfo2() {}
  int nesting_level;
  bool in_domain() { return nesting_level % 2 == 1; }
};

typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Tfbi2;
typedef CGAL::Constrained_triangulation_face_base_2<K, Tfbi2> Ctfb2;
typedef CGAL::Triangulation_data_structure_2<Vb2, Ctfb2> Tds22;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds22, Itag> CDT;
typedef CDT::Point CDPoint;
typedef CDT::Face_handle CDFace_handle;

typedef Rcpp::NumericVector Dvector;

typedef Eigen::
  Matrix<unsigned, Eigen::Dynamic, 3, Eigen::RowMajor | Eigen::AutoAlign>
  Imatrix;
typedef Eigen::Matrix<unsigned, 1, 3, Eigen::RowMajor | Eigen::AutoAlign>
  Ivector;

// -------------------------------------------------------------------------- //
bool approxEqual(double, double, double);

bool approxEqualVectors(Vector3, Vector3, double);

Rcpp::String stringPair(const size_t, const size_t);

std::array<Rcpp::String, 3> triangleEdges(const size_t,
                                          const size_t,
                                          const size_t );

Rcpp::String stringTriple(const size_t, const size_t, const size_t);

double volume_under_triangle(Dvector, Dvector, Dvector);

void mark_domains0(CDT&,
                   CDFace_handle,
                   int,
                   std::list<CDT::Edge>&);

void mark_domains(CDT&);
