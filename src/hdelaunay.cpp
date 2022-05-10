#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>
//#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
//#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<K> HDtt;
typedef HDtt::Point_2 HPoint;
typedef CGAL::Triangulation_data_structure_2<
    CGAL::Triangulation_vertex_base_with_id_2<HDtt>,
    CGAL::Hyperbolic_triangulation_face_base_2<HDtt>>
    HTds;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<HDtt, HTds> HDt;
//typedef boost::graph_traits<CGAL::Delaunay_triangulation_2<K, HTds>>::vertex_descriptor hvertex_descriptor;

// [[Rcpp::export]]
Rcpp::List htest(const Rcpp::NumericMatrix points) {
  //std::vector<std::pair<HPoint, unsigned>> hpts;

  std::vector<HPoint> hpts;
  const unsigned npoints = points.ncol();
  hpts.reserve(npoints);
  for(unsigned i = 0; i != npoints; i++) {
    const Rcpp::NumericVector pt = points(Rcpp::_, i);
    //hpts.emplace_back(std::make_pair(HPoint(pt[0], pt[1]), i));
    hpts.emplace_back(HPoint(pt(0), pt(1)));
  }
  HDt hdt;
  //CGAL::set_triangulation_ids(hdt);
  hdt.insert(hpts.begin(), hpts.end());
  Rcpp::NumericMatrix Vertices(2, hdt.number_of_vertices());
  int index = 0;
  for(HDt::All_vertices_iterator vd = hdt.all_vertices_begin();
      vd != hdt.all_vertices_end(); ++vd){
    vd->id() = index;
    HPoint pt = vd->point();
    Vertices(0, index) = pt.x();
    Vertices(1, index) = pt.y();
    Rcpp::Rcout << index << "\n";
    index++;
  }
  const size_t nedges = hdt.number_of_hyperbolic_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  size_t i = 0;
  for(HDt::All_edges_iterator ed = hdt.all_edges_begin();
      ed != hdt.all_edges_end(); ++ed) {
    Rcpp::IntegerVector edge_i(2);
    HDt::Vertex_handle sVertex = ed->first->vertex(HDt::cw(ed->second));
    //Rcpp::Rcout << sVertex->id();
    edge_i(0) = sVertex->id();
    HDt::Vertex_handle tVertex = ed->first->vertex(HDt::ccw(ed->second));
    edge_i(1) = tVertex->id();
    Edges(Rcpp::_, i) = edge_i;
    Rcpp::Rcout << i << "\n";
    i++;
  }
  return Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                            Rcpp::Named("edges") = Edges);
}
