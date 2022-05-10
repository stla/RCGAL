#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

#include <CGAL/Hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_traits_2.h>

typedef CGAL::Hyperbolic_Delaunay_triangulation_traits_2<>  HDtt;
typedef HDtt::Point_2                                         HPoint;
typedef CGAL::Hyperbolic_Delaunay_triangulation_2<HDtt>       HDt;

size_t htest(std::vector<std::vector<double>> points){
	std::vector<HPoint> hpts;
	const size_t npoints = points.size();
	hpts.reserve(npoints);
	for(size_t i=0; i != npoints; i++){
		const std::vector<double> pt = points[i];
		hpts[i] = HPoint(pt[1], pt[2]);
	}
	HDt hdt;
	hdt.insert(hpts.begin(), hpts.end());
	for(HDt::Edge_iterator ed = hdt.all_edges_begin(); ed != hdt.all_edges_end(); ++ed){
		Rcpp::Rcout << source(ed, hdt);
	}
	return hdt.number_of_vertices();
}