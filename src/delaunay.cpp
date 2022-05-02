#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// [[Rcpp::export]]
Rcpp::List del2d_cpp(Rcpp::NumericMatrix pts) {
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

  const DT2::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 3);
  {
    size_t i = 0;
    for(DT2::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
    eit++) {
      const std::pair<DT2::Face_handle, int> edge = *eit;
      const size_t i0 = edge.first->vertex((edge.second + 1) % 3)->info();
      const size_t i1 = edge.first->vertex((edge.second + 2) % 3)->info();
      edges(i, 0) = i0;
      edges(i, 1) = i1;
      const Rcpp::String i0i1 = stringPair(i0, i1);
      unsigned flag = 0;
      size_t j = 0;
      while(flag < 2 && j <= nfaces) {
        Rcpp::String pair01 =
          stringPair((size_t)faces(j, 0), (size_t)faces(j, 1));
        Rcpp::String pair02 =
          stringPair((size_t)faces(j, 0), (size_t)faces(j, 2));
        Rcpp::String pair12 =
          stringPair((size_t)faces(j, 1), (size_t)faces(j, 2));
        if(i0i1 == pair01 || i0i1 == pair02 || i0i1 == pair12) {
          flag++;
        }
        j++;
      }
      if(flag == 1) {
        edges(i, 2) = 1;
      }
      i++;
    }
    Rcpp::CharacterVector columnNames =
      Rcpp::CharacterVector::create("i1", "i2", "border");
    Rcpp::colnames(edges) = columnNames;
    edges.attr("info") = "The `border` column indicates border edges.";
  }

  return Rcpp::List::create(Rcpp::Named("faces") = faces,
                            Rcpp::Named("edges") = edges);
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

  std::map<Rcpp::String, size_t> facetsMap = {};
  Rcpp::IntegerMatrix facets(nfacets, 4);
  // Imatrix facetsOnHull(0, 3);
  std::vector<Rcpp::String> edgesOnHull(0);
  {
    size_t i = 0;
    for(DT3::Finite_facets_iterator fit = mesh.finite_facets_begin();
        fit != mesh.finite_facets_end(); fit++) {
      std::pair<DT3::Cell_handle, int> facet = *fit;
      DT3::Vertex_handle v0 = facet.first->vertex((facet.second + 1) % 4);
      DT3::Vertex_handle v1 = facet.first->vertex((facet.second + 2) % 4);
      DT3::Vertex_handle v2 = facet.first->vertex((facet.second + 3) % 4);
      const size_t id0 = v0->info();
      const size_t id1 = v1->info();
      const size_t id2 = v2->info();
      facets(i, 0) = id0;
      facets(i, 1) = id1;
      facets(i, 2) = id2;
      bool onhull = mesh.is_infinite(facet.first) ||
        mesh.is_infinite(mesh.mirror_facet(facet).first);
      if(onhull) {
        facets(i, 3) = onhull;
        std::array<Rcpp::String, 3> triangle = triangleEdges(id0, id1, id2);
        edgesOnHull.insert(edgesOnHull.end(), triangle.begin(), triangle.end());
      }
      std::array<size_t, 3> ids = {id0, id1, id2};
      std::sort(ids.begin(), ids.end());
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

  const DT3::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 3);
  {
    size_t i = 0;
    for(DT3::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
    eit++) {
      const CGAL::Triple<DT3::Cell_handle, int, int> edge = *eit;
      size_t i0 = edge.first->vertex(edge.second % 4)->info();
      size_t i1 = edge.first->vertex(edge.third % 4)->info();
      edges(i, 0) = i0;  // std::min(i0, i1);
      edges(i, 1) = i1;  // std::max(i0, i1);
      Rcpp::String i0i1 = stringPair(i0, i1);
      if(std::find(edgesOnHull.begin(), edgesOnHull.end(), i0i1) !=
         edgesOnHull.end()) {
        edges(i, 2) = 1;
      }
      i++;
    }
    Rcpp::CharacterVector columnNames =
      Rcpp::CharacterVector::create("i1", "i2", "onhull");
    Rcpp::colnames(edges) = columnNames;
    edges.attr("info") =
      "The `onhull` column indicates whether the edge is on the convex hull.";
  }

  // {
  //   for(size_t i = 0; i < nfacets-1; i++){
  //     Rcpp::IntegerVector facet_i = facets(i, Rcpp::_);
  //     if(facet_i(3)){
  //       for(size_t j = i+1; j < nfacets; j++){
  //         Rcpp::IntegerVector facet_j = facets(j, Rcpp::_);
  //         if(facet_j(3)){
  //           std::vector<int> v_intersection;
  //           std::set_intersection(facet_i.begin(), facet_i.end(),
  //                                 facet_j.begin(), facet_j.end(),
  //                                 std::back_inserter(v_intersection));
  //           if(v_intersection.size() == 2){
  //             Rcpp::Rcout << v_intersection(0) << " --- "
  //             Rcpp::Rcout << v_intersection(1) << "\n"
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

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


// [[Rcpp::export]]
Rcpp::List del2d_xy_cpp(Rcpp::NumericMatrix pts) {
  const unsigned npoints = pts.nrow();

  // compute Delaunay mesh
  Delaunay_xy mesh;
  Delaunay_xy::Vertex_handle vh;
  for(unsigned i = 0; i < npoints; i++) {
    vh = mesh.insert(Point3(pts(i, 0), pts(i, 1), pts(i, 2)));
    vh->info() = i;
  }

  const size_t nfaces = mesh.number_of_faces();
  const size_t h =
    2 * npoints - 2 - nfaces;  // number of vertices of convex hull
  const size_t nedges = 3 * npoints - 3 - h;

  Rcpp::IntegerMatrix faces(nfaces, 3);
  Rcpp::NumericVector volumes(nfaces);
  double totalVolume = 0.0;
  {
    size_t i = 0;
    for(Delaunay_xy::Finite_faces_iterator fit = mesh.finite_faces_begin();
        fit != mesh.finite_faces_end(); fit++) {
      unsigned i0 = fit->vertex(0)->info();
      unsigned i1 = fit->vertex(1)->info();
      unsigned i2 = fit->vertex(2)->info();
      faces(i, 0) = i0 + 1;
      faces(i, 1) = i1 + 1;
      faces(i, 2) = i2 + 1;
      double volume = volume_under_triangle(pts(i0, Rcpp::_), pts(i1, Rcpp::_),
                                            pts(i2, Rcpp::_));
      volumes(i) = volume;
      totalVolume += volume;
      i++;
    }
    faces.attr("volumes") = volumes;
  }

  const Delaunay_xy::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix edges(nedges, 2);
  {
    size_t i = 0;
    for(Delaunay_xy::Finite_edges_iterator eit = itedges.begin();
        eit != itedges.end(); eit++) {
      const std::pair<Delaunay_xy::Face_handle, int> edge = *eit;
      edges(i, 0) = edge.first->vertex((edge.second + 1) % 3)->info() + 1;
      edges(i, 1) = edge.first->vertex((edge.second + 2) % 3)->info() + 1;
      i++;
    }
  }

  return Rcpp::List::create(Rcpp::Named("faces") = faces,
                            Rcpp::Named("edges") = edges,
                            Rcpp::Named("volume") = totalVolume);
}

//// ---- Constrained Delaunay -------------------------------------------- ////

// [[Rcpp::export]]
Rcpp::IntegerMatrix del2d_constrained_cpp(Rcpp::NumericMatrix pts,
                                          Rcpp::IntegerMatrix edges) {
  const unsigned npoints = pts.nrow();
  std::vector<std::pair<CDT::Point, unsigned>> points(npoints);
  for(unsigned i = 0; i < npoints; ++i) {
    points[i] = std::make_pair(CDT::Point(pts(i, 0), pts(i, 1)), i);
  }
  CDT cdt;
  {
    const size_t nedges = edges.nrow();
    for(size_t k = 0; k < nedges; ++k) {
      cdt.insert_constraint(points[edges(k, 0) - 1].first,
                            points[edges(k, 1) - 1].first);
    }
  }
  cdt.insert(points.begin(), points.end());
  const size_t nfaces = cdt.number_of_faces();
  Rcpp::IntegerMatrix faces(nfaces, 3);
  mark_domains(cdt);
  size_t nfaces_out;
  {
    size_t i = 0;
    for(CDFace_handle f : cdt.finite_face_handles()){
      if(f->info().in_domain()){
        unsigned id0 = f->vertex(0)->info();
        unsigned id1 = f->vertex(1)->info();
        unsigned id2 = f->vertex(2)->info();
        std::array<unsigned, 3> ids = {id0 + 1, id1 + 1, id2 + 1};
        std::sort(ids.begin(), ids.end());
        faces(i, 0) = ids[0];
        faces(i, 1) = ids[1];
        faces(i, 2) = ids[2];
        ++i;
      }
    }
    nfaces_out = i;
  }
  return faces(Rcpp::Range(0, nfaces_out-1), Rcpp::_);
}
