#ifndef _RCGALHEADER_
#define _RCGALHEADER_
#endif

// [[Rcpp::depends(RcppCGAL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp14)]]
#define CGAL_EIGEN3_ENABLED 1

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

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

// #include <CGAL/Nef_polyhedron_3.h>
// #include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/number_utils.h>

// #include <CGAL/Triangle_3.h>

#include <boost/multiprecision/gmp.hpp>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>


#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <array>
//#include <limits>
#include <map>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;

typedef K::Point_2 Point2;
typedef K::Point_3 Point3;
typedef EK::Point_3 EPoint3;
typedef std::vector<Point2> Points2;
typedef std::vector<Point3> Points3;
typedef std::pair<Point3, unsigned> IPoint3;

typedef CGAL::
    Convex_hull_traits_adapter_2<K, CGAL::Pointer_property_map<Point2>::type>
        CHT2;

typedef CGAL::First_of_pair_property_map<IPoint3> Pmap;
typedef CGAL::Extreme_points_traits_adapter_3<Pmap,
                                              CGAL::Convex_hull_traits_3<K>>
    CHT;

typedef CGAL::Surface_mesh<IPoint3> Mesh;
typedef CGAL::Surface_mesh<Point3> Mesh3;
typedef CGAL::Surface_mesh<EPoint3> EMesh3;

typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Edge_index edge_descriptor;
typedef Mesh::Face_index face_descriptor;

typedef Mesh3::Vertex_index m3_vertex_descriptor;
typedef Mesh3::Face_index m3_face_descriptor;
typedef Mesh3::Edge_index m3_edge_descriptor;

typedef boost::graph_traits<Mesh3>::vertex_descriptor boost_vertex_descriptor;
typedef boost::graph_traits<Mesh3>::face_descriptor boost_face_descriptor;

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
typedef EK::Vector_3 EVector3;

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

//typedef CGAL::Nef_polyhedron_3<EK> Nef;

typedef CGAL::Cartesian<boost::multiprecision::mpq_rational> QK;
typedef CGAL::Surface_mesh<QK::Point_3> QMesh3;
typedef QK::Point_3 QPoint3;
typedef QK::Vector_3 QVector3;

typedef Rcpp::NumericVector Dvector;

typedef Eigen::
    Matrix<unsigned, Eigen::Dynamic, 3, Eigen::RowMajor | Eigen::AutoAlign>
        Imatrix;
typedef Eigen::Matrix<unsigned, 1, 3, Eigen::RowMajor | Eigen::AutoAlign>
    Ivector;

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;
namespace mp = boost::multiprecision;

// -------------------------------------------------------------------------- //
bool approxEqual(double, double, double);

bool approxEqualVectors(Vector3, Vector3, double);

Rcpp::String stringPair(const size_t, const size_t);

std::array<Rcpp::String, 3> triangleEdges(const size_t,
                                          const size_t,
                                          const size_t);

Rcpp::String stringTriple(const size_t, const size_t, const size_t);

double volume_under_triangle(Dvector, Dvector, Dvector);

void mark_domains0(CDT&, CDFace_handle, int, std::list<CDT::Edge>&);

void mark_domains(CDT&);

//template <typename PointT>
//std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix);

// std::vector<std::vector<size_t>> matrix_to_faces(const Rcpp::IntegerMatrix);

//std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List);

// Polyhedron makePolyMesh(const Rcpp::NumericMatrix, const Rcpp::IntegerMatrix);

// Rcpp::List RPolyMesh(Polyhedron);

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(const Rcpp::List, const bool);

QMesh3 makeSurfQMesh(const Rcpp::List, const bool);

// template <typename MeshT>
// Rcpp::IntegerMatrix getEdges1(MeshT);
template <typename KernelT, typename MeshT, typename PointT>
Rcpp::IntegerMatrix getEdges2(MeshT, const double);

Rcpp::NumericMatrix getKNormals(Mesh3);
Rcpp::NumericMatrix getEKNormals(EMesh3);
Rcpp::NumericMatrix getQNormals(QMesh3);

Rcpp::List RSurfKMesh(Mesh3, const bool, const double);
Rcpp::List RSurfEKMesh(EMesh3, const bool, const double);
Rcpp::List RSurfQMesh(QMesh3, const bool, const double);


// template <typename KernelT, typename MeshT, typename PointT>
// Rcpp::List RSurfMesh(MeshT, const bool, const double, const bool);

// Mesh3 Poly2Mesh3(Polyhedron);
