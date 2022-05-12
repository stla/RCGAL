#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

#include <boost/multiprecision/gmp.hpp>
namespace mp = boost::multiprecision;
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>


typedef CGAL::Cartesian<mp::mpq_rational> QK;
typedef CGAL::Surface_mesh<QK::Point_3> QMesh;

// [[Rcpp::export]]
int testgmp() {
  mp::mpq_rational phi(33137336090491, 20480000000000);
  mp::mpq_rational a(1099511627776, 1904410004001);
  mp::mpq_rational b = a / phi;
  mp::mpq_rational c = a * phi;

  QMesh mesh1, mesh2, mesh3, mesh4, mesh5;

  QMesh::Vertex_index v0 = mesh1.add_vertex(QK::Point_3(0, b, c));
  QMesh::Vertex_index v1 = mesh1.add_vertex(QK::Point_3(b, -c, 0));
  QMesh::Vertex_index v2 = mesh1.add_vertex(QK::Point_3(a, a, -a));
  QMesh::Vertex_index v3 = mesh1.add_vertex(QK::Point_3(-c, 0, b));
  QMesh::Vertex_index v4 = mesh2.add_vertex(QK::Point_3(a, -a, -a));
  QMesh::Vertex_index v5 = mesh2.add_vertex(QK::Point_3(a, a, a));
  QMesh::Vertex_index v6 = mesh2.add_vertex(QK::Point_3(-a, -a, a));
  QMesh::Vertex_index v7 = mesh2.add_vertex(QK::Point_3(-a, a, -a));
  QMesh::Vertex_index v8 = mesh3.add_vertex(QK::Point_3(c, 0, b));
  QMesh::Vertex_index v9 = mesh3.add_vertex(QK::Point_3(-a, a, a));
  QMesh::Vertex_index v10 = mesh3.add_vertex(QK::Point_3(-b, -c, 0));
  QMesh::Vertex_index v11 = mesh3.add_vertex(QK::Point_3(0, b, -c));
  QMesh::Vertex_index v12 = mesh4.add_vertex(QK::Point_3(a, -a, a));
  QMesh::Vertex_index v13 = mesh4.add_vertex(QK::Point_3(b, c, 0));
  QMesh::Vertex_index v14 = mesh4.add_vertex(QK::Point_3(-c, 0, b));
  QMesh::Vertex_index v15 = mesh4.add_vertex(QK::Point_3(0, -b, -c));
  QMesh::Vertex_index v16 = mesh5.add_vertex(QK::Point_3(-a, -a, -a));
  QMesh::Vertex_index v17 = mesh5.add_vertex(QK::Point_3(-b, c, 0));
  QMesh::Vertex_index v18 = mesh5.add_vertex(QK::Point_3(c, 0, -b));
  QMesh::Vertex_index v19 = mesh5.add_vertex(QK::Point_3(0, -b, c));

  mesh1.add_face(v0, v1, v2);
  mesh1.add_face(v2, v1, v3);
  mesh1.add_face(v3, v1, v0);
  mesh1.add_face(v0, v2, v3);
  mesh2.add_face(v4, v5, v6);
  mesh2.add_face(v6, v5, v7);
  mesh2.add_face(v7, v5, v4);
  mesh2.add_face(v4, v6, v7);
  mesh3.add_face(v8, v9, v10);
  mesh3.add_face(v10, v9, v11);
  mesh3.add_face(v11, v9, v8);
  mesh3.add_face(v8, v10, v11);
  mesh4.add_face(v12, v13, v14);
  mesh4.add_face(v14, v13, v15);
  mesh4.add_face(v15, v13, v12);
  mesh4.add_face(v12, v14, v15);
  mesh5.add_face(v16, v17, v18);
  mesh5.add_face(v18, v17, v19);
  mesh5.add_face(v19, v17, v16);
  mesh5.add_face(v16, v18, v19);

  QMesh inter12;
  if(PMP::corefine_and_compute_intersection(mesh1, mesh2, inter12)) {
    std::cout << "Intersection th1-th2 successfully computed.\n";
    if(PMP::does_self_intersect(inter12)) {
      std::cout << "Intersection th1-th2 self-intersects.\n";
      return 0;
    }
    if(!PMP::does_bound_a_volume(inter12)) {
      std::cout << "Intersection th1-th2 does not bound a volume.\n";
      return 0;
    }
    QMesh inter34;
    if(PMP::corefine_and_compute_intersection(mesh3, mesh4, inter34)) {
      std::cout << "Intersection th3-th4 successfully computed.\n";
      if(PMP::does_self_intersect(inter34)) {
        std::cout << "Intersection th3-th4 self-intersects.\n";
        return 0;
      }
      if(!PMP::does_bound_a_volume(inter34)) {
        std::cout << "Intersection th3-th4 does not bound a volume.\n";
        return 0;
      }
      QMesh inter1234;
      if(PMP::corefine_and_compute_intersection(inter12, inter34, inter1234)) {
        std::cout << "Intersection th1-th2-th3-th4 successfully computed.\n";
        if(PMP::does_self_intersect(inter1234)) {
          std::cout << "Intersection th1-th2-th3-th4 self-intersects.\n";
          return 0;
        }
        if(!PMP::does_bound_a_volume(inter1234)) {
          std::cout
              << "Intersection th1-th2-th3-th4 does not bound a volume.\n";
          return 0;
        }
        QMesh inter12345;
        if(PMP::corefine_and_compute_intersection(inter1234, mesh5,
                                                  inter12345)) {
          std::cout << "Final intersection successfully computed.\n";
          CGAL::IO::write_polygon_mesh("inter_tetrahedra.off", inter12345,
                                       CGAL::parameters::stream_precision(17));
          return 0;
        } else {
          std::cout << "Final intersection failed.\n";
          return 1;
        }
      } else {
        std::cout << "Intersection th1-th2-th3-th4 failed.\n";
        return 1;
      }
    } else {
      std::cout << "Intersection th3-th4 failed.\n";
      return 1;
    }
  } else {
    std::cout << "Intersection th1-th2 failed.\n";
    return 1;
  }
}