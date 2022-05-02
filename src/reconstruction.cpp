#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// [[Rcpp::export]]
Rcpp::List AFSreconstruction_cpp(Rcpp::NumericMatrix pts) {
  const size_t npoints = pts.nrow();
  std::vector<Point3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = Point3(pts(i, 0), pts(i, 1), pts(i, 2));
  }

  AFS_triangulation3 dt(points.begin(), points.end());
  AFS_reconstruction reconstruction(dt);
  reconstruction.run();
  const AFS_Tds2& tds = reconstruction.triangulation_data_structure_2();

  // Eigen::MatrixXd normals(3, 0);
  Eigen::MatrixXd vertices(4, 0);
  unsigned counter = 0;
  for(AFS_Tds2::Face_iterator fit = tds.faces_begin(); fit != tds.faces_end();
  ++fit) {
    if(reconstruction.has_on_surface(fit)) {
      counter++;
      AFS_triangulation3::Facet f = fit->facet();
      AFS_triangulation3::Cell_handle ch = f.first;
      int ci = f.second;
      Point3 points[3];
      for(int i = 0, j = 0; i < 4; i++) {
        if(ci != i) {
          points[j] = ch->vertex(i)->point();
          j++;
        }
      }
      // Vector3 normal = CGAL::unit_normal(points[0], points[1], points[2]);
      // Eigen::VectorXd v(3);
      // v << normal.x(), normal.y(), normal.z();
      // Eigen::MatrixXd M(3, 3);
      // M << v, v, v;
      // normals.conservativeResize(Eigen::NoChange, normals.cols() + 3);
      // normals.rightCols(3) = M;
      for(size_t k = 0; k < 3; k++) {
        const Point3 p = points[k];
        Eigen::VectorXd w(4);
        w << p.x(), p.y(), p.z(), 1.0;
        vertices.conservativeResize(Eigen::NoChange, vertices.cols() + 1);
        vertices.rightCols(1) = w;
      }
    }
  }
  Rcpp::IntegerVector vtriangles(3 * counter);
  for(size_t i = 0; i < 3 * counter; i++) {
    vtriangles(i) = i + 1;
  }
  vtriangles.attr("dim") = Rcpp::Dimension(3, counter);
  Rcpp::IntegerMatrix triangles = Rcpp::as<Rcpp::IntegerMatrix>(vtriangles);

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                            // Rcpp::Named("normals") = normals,
                            Rcpp::Named("triangles") = triangles);
}

struct Perimeter {
  double bound;
  Perimeter(double bound) : bound(bound) {}
  template <typename AdvancingFront, typename Cell_handle>
  double operator()(const AdvancingFront& adv,
                  Cell_handle& c,
                  const int& index) const {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0) {
      return adv.smallest_radius_delaunay_sphere(c, index);
    }
    // If perimeter > bound, return infinity so that facet is not used
    double d = 0;
    d = sqrt(squared_distance(c->vertex((index + 1) % 4)->point(),
                              c->vertex((index + 2) % 4)->point()));
    if(d > bound)
      return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index + 2) % 4)->point(),
                               c->vertex((index + 3) % 4)->point()));
    if(d > bound)
      return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index + 1) % 4)->point(),
                               c->vertex((index + 3) % 4)->point()));
    if(d > bound)
      return adv.infinity();
    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere(c, index);
  }
};

// [[Rcpp::export]]
Rcpp::List AFSreconstruction_perimeter_cpp(Rcpp::NumericMatrix pts,
                                           double per) {
  const size_t npoints = pts.nrow();
  std::vector<KSC::Point_3> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] = KSC::Point_3(pts(i, 0), pts(i, 1), pts(i, 2));
  }

  // double per = 0;
  double radius_ratio_bound = 5.0;

  std::vector<std::array<size_t, 3>> facets;

  Perimeter perimeter(per);
  CGAL::advancing_front_surface_reconstruction(points.begin(), points.end(),
                                               std::back_inserter(facets),
                                               perimeter, radius_ratio_bound);

  size_t nfacets = facets.size();
  Rcpp::IntegerMatrix triangles(3, nfacets);
  for(size_t j = 0; j < nfacets; j++) {
    std::array<size_t, 3> facet = facets[j];
    triangles(0, j) = facet[0] + 1;
    triangles(1, j) = facet[1] + 1;
    triangles(2, j) = facet[2] + 1;
  }

  return Rcpp::List::create(Rcpp::Named("vertices") = pts,
                            Rcpp::Named("triangles") = triangles);
}

// [[Rcpp::export]]
Rcpp::List Poisson_reconstruction_cpp(Rcpp::NumericMatrix pts,
                                      Rcpp::NumericMatrix normals,
                                      double spacing,
                                      double sm_angle,
                                      double sm_radius,
                                      double sm_distance) {
  const size_t npoints = pts.nrow();
  std::vector<P3wn> points(npoints);
  for(size_t i = 0; i < npoints; i++) {
    points[i] =
      std::make_pair(Point3(pts(i, 0), pts(i, 1), pts(i, 2)),
                     Vector3(normals(i, 0), normals(i, 1), normals(i, 2)));
  }

  Polyhedron mesh;
  if(spacing == -1.0) {
    spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, 6,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>()));
  }

  bool psr = CGAL::poisson_surface_reconstruction_delaunay(
    points.begin(), points.end(), CGAL::First_of_pair_property_map<P3wn>(),
    CGAL::Second_of_pair_property_map<P3wn>(), mesh, spacing, sm_angle,
    sm_radius, sm_distance);

  if(!psr) {
    throw Rcpp::exception("Poisson surface reconstruction has failed.");
  }

  int id = 1;
  for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
      vit != mesh.vertices_end(); ++vit) {
    vit->id() = id;
    id++;
  }

  const size_t nfacets = mesh.size_of_facets();
  const size_t nvertices = mesh.size_of_vertices();

  Rcpp::IntegerMatrix facets(nfacets, 3);
  {
    size_t i = 0;
    for(Polyhedron::Facet_iterator fit = mesh.facets_begin();
        fit != mesh.facets_end(); fit++) {
      facets(i, 0) = fit->halfedge()->vertex()->id();
      facets(i, 1) = fit->halfedge()->next()->vertex()->id();
      facets(i, 2) = fit->halfedge()->opposite()->vertex()->id();
      i++;
    }
  }

  Rcpp::NumericMatrix vertices(nvertices, 3);
  {
    size_t i = 0;
    for(Polyhedron::Vertex_iterator vit = mesh.vertices_begin();
        vit != mesh.vertices_end(); vit++) {
      vertices(i, 0) = vit->point().x();
      vertices(i, 1) = vit->point().y();
      vertices(i, 2) = vit->point().z();
      i++;
    }
  }

  // std::ofstream("out.off") << std::setprecision(17) << mesh;

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,
                            Rcpp::Named("facets") = facets,
                            Rcpp::Named("spacing") = spacing);
}
