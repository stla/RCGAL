#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  HDtt;
typedef HDtt::Point_2                                         HPoint;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<HDtt, Tds2>       HDt;

size_t htest(std::vector<std::vector<double>> points){

	std::vector<std::pair<HPoint, unsigned>> hpts;

	//std::vector<HPoint> hpts;
	const unsigned npoints = points.size();
	hpts.reserve(npoints);
	for(unsigned i=0; i != npoints; i++){
		const std::vector<double> pt = points[i];
		hpts.emplace_back(std::make_pair(HPoint(pt[0], pt[1]), i));
	}
	HDt hdt;
	hdt.insert(hpts.begin(), hpts.end());
	for(HDt::All_edges_iterator ed = hdt.all_edges_begin(); ed != hdt.all_edges_end(); ++ed){
		HDt::Vertex_handle sVertex = ed->first->vertex(HDt::cw(ed->second));
		Rcpp::Rcout << sVertex->info();
	}
	return hdt.number_of_vertices();
}