// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cxhull2d_cpp
Rcpp::List cxhull2d_cpp(Rcpp::NumericMatrix pts);
RcppExport SEXP _RCGAL_cxhull2d_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(cxhull2d_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// cxhull3d_cpp
Rcpp::List cxhull3d_cpp(Rcpp::NumericMatrix pts, double epsilon);
RcppExport SEXP _RCGAL_cxhull3d_cpp(SEXP ptsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(cxhull3d_cpp(pts, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// del2d_cpp
Rcpp::List del2d_cpp(Rcpp::NumericMatrix pts);
RcppExport SEXP _RCGAL_del2d_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(del2d_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// del3d_cpp
Rcpp::List del3d_cpp(Rcpp::NumericMatrix pts);
RcppExport SEXP _RCGAL_del3d_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(del3d_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// del2d_xy_cpp
Rcpp::List del2d_xy_cpp(Rcpp::NumericMatrix pts);
RcppExport SEXP _RCGAL_del2d_xy_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(del2d_xy_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// del2d_constrained_cpp
Rcpp::IntegerMatrix del2d_constrained_cpp(Rcpp::NumericMatrix pts, Rcpp::IntegerMatrix edges);
RcppExport SEXP _RCGAL_del2d_constrained_cpp(SEXP ptsSEXP, SEXP edgesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type edges(edgesSEXP);
    rcpp_result_gen = Rcpp::wrap(del2d_constrained_cpp(pts, edges));
    return rcpp_result_gen;
END_RCPP
}
// hdelaunay_K
Rcpp::List hdelaunay_K(const Rcpp::NumericMatrix points);
RcppExport SEXP _RCGAL_hdelaunay_K(SEXP pointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type points(pointsSEXP);
    rcpp_result_gen = Rcpp::wrap(hdelaunay_K(points));
    return rcpp_result_gen;
END_RCPP
}
// hdelaunay_EK
Rcpp::List hdelaunay_EK(const Rcpp::NumericMatrix points);
RcppExport SEXP _RCGAL_hdelaunay_EK(SEXP pointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type points(pointsSEXP);
    rcpp_result_gen = Rcpp::wrap(hdelaunay_EK(points));
    return rcpp_result_gen;
END_RCPP
}
// jet_normals_cpp
Rcpp::NumericMatrix jet_normals_cpp(Rcpp::NumericMatrix pts, unsigned nb_neighbors);
RcppExport SEXP _RCGAL_jet_normals_cpp(SEXP ptsSEXP, SEXP nb_neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nb_neighbors(nb_neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(jet_normals_cpp(pts, nb_neighbors));
    return rcpp_result_gen;
END_RCPP
}
// pca_normals_cpp
Rcpp::NumericMatrix pca_normals_cpp(Rcpp::NumericMatrix pts, unsigned nb_neighbors);
RcppExport SEXP _RCGAL_pca_normals_cpp(SEXP ptsSEXP, SEXP nb_neighborsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< unsigned >::type nb_neighbors(nb_neighborsSEXP);
    rcpp_result_gen = Rcpp::wrap(pca_normals_cpp(pts, nb_neighbors));
    return rcpp_result_gen;
END_RCPP
}
// AFSreconstruction_cpp
Rcpp::List AFSreconstruction_cpp(Rcpp::NumericMatrix pts);
RcppExport SEXP _RCGAL_AFSreconstruction_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(AFSreconstruction_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// AFSreconstruction_perimeter_cpp
Rcpp::List AFSreconstruction_perimeter_cpp(Rcpp::NumericMatrix pts, double per);
RcppExport SEXP _RCGAL_AFSreconstruction_perimeter_cpp(SEXP ptsSEXP, SEXP perSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< double >::type per(perSEXP);
    rcpp_result_gen = Rcpp::wrap(AFSreconstruction_perimeter_cpp(pts, per));
    return rcpp_result_gen;
END_RCPP
}
// Poisson_reconstruction_cpp
Rcpp::List Poisson_reconstruction_cpp(Rcpp::NumericMatrix pts, Rcpp::NumericMatrix normals, double spacing, double sm_angle, double sm_radius, double sm_distance);
RcppExport SEXP _RCGAL_Poisson_reconstruction_cpp(SEXP ptsSEXP, SEXP normalsSEXP, SEXP spacingSEXP, SEXP sm_angleSEXP, SEXP sm_radiusSEXP, SEXP sm_distanceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< double >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< double >::type sm_angle(sm_angleSEXP);
    Rcpp::traits::input_parameter< double >::type sm_radius(sm_radiusSEXP);
    Rcpp::traits::input_parameter< double >::type sm_distance(sm_distanceSEXP);
    rcpp_result_gen = Rcpp::wrap(Poisson_reconstruction_cpp(pts, normals, spacing, sm_angle, sm_radius, sm_distance));
    return rcpp_result_gen;
END_RCPP
}
// SurfMesh
Rcpp::List SurfMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool merge, const bool normals, const double epsilon);
RcppExport SEXP _RCGAL_SurfMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP mergeSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfMesh(rmesh, isTriangle, triangulate, merge, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// SurfEMesh
Rcpp::List SurfEMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool merge, const bool normals, const double epsilon);
RcppExport SEXP _RCGAL_SurfEMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP mergeSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfEMesh(rmesh, isTriangle, triangulate, merge, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// SurfQMesh
Rcpp::List SurfQMesh(const Rcpp::List rmesh, const bool isTriangle, const bool triangulate, const bool merge, const bool normals, const double epsilon);
RcppExport SEXP _RCGAL_SurfQMesh(SEXP rmeshSEXP, SEXP isTriangleSEXP, SEXP triangulateSEXP, SEXP mergeSEXP, SEXP normalsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh(rmeshSEXP);
    Rcpp::traits::input_parameter< const bool >::type isTriangle(isTriangleSEXP);
    Rcpp::traits::input_parameter< const bool >::type triangulate(triangulateSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfQMesh(rmesh, isTriangle, triangulate, merge, normals, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// Intersection2_K
Rcpp::List Intersection2_K(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Intersection2_K(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection2_K(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Intersection2_EK
Rcpp::List Intersection2_EK(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Intersection2_EK(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection2_EK(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Intersection_Q
Rcpp::List Intersection_Q(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Intersection_Q(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Intersection_Q(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Difference_K
Rcpp::List Difference_K(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Difference_K(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_K(rmesh1, rmesh2, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Difference_EK
Rcpp::List Difference_EK(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Difference_EK(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_EK(rmesh1, rmesh2, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Difference_Q
Rcpp::List Difference_Q(const Rcpp::List rmesh1, const Rcpp::List rmesh2, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Difference_Q(SEXP rmesh1SEXP, SEXP rmesh2SEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh1(rmesh1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmesh2(rmesh2SEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Difference_Q(rmesh1, rmesh2, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Union_K
Rcpp::List Union_K(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Union_K(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_K(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Union_EK
Rcpp::List Union_EK(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Union_EK(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_EK(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}
// Union_Q
Rcpp::List Union_Q(const Rcpp::List rmeshes, const bool merge, const bool normals);
RcppExport SEXP _RCGAL_Union_Q(SEXP rmeshesSEXP, SEXP mergeSEXP, SEXP normalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type rmeshes(rmeshesSEXP);
    Rcpp::traits::input_parameter< const bool >::type merge(mergeSEXP);
    Rcpp::traits::input_parameter< const bool >::type normals(normalsSEXP);
    rcpp_result_gen = Rcpp::wrap(Union_Q(rmeshes, merge, normals));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RCGAL_cxhull2d_cpp", (DL_FUNC) &_RCGAL_cxhull2d_cpp, 1},
    {"_RCGAL_cxhull3d_cpp", (DL_FUNC) &_RCGAL_cxhull3d_cpp, 2},
    {"_RCGAL_del2d_cpp", (DL_FUNC) &_RCGAL_del2d_cpp, 1},
    {"_RCGAL_del3d_cpp", (DL_FUNC) &_RCGAL_del3d_cpp, 1},
    {"_RCGAL_del2d_xy_cpp", (DL_FUNC) &_RCGAL_del2d_xy_cpp, 1},
    {"_RCGAL_del2d_constrained_cpp", (DL_FUNC) &_RCGAL_del2d_constrained_cpp, 2},
    {"_RCGAL_hdelaunay_K", (DL_FUNC) &_RCGAL_hdelaunay_K, 1},
    {"_RCGAL_hdelaunay_EK", (DL_FUNC) &_RCGAL_hdelaunay_EK, 1},
    {"_RCGAL_jet_normals_cpp", (DL_FUNC) &_RCGAL_jet_normals_cpp, 2},
    {"_RCGAL_pca_normals_cpp", (DL_FUNC) &_RCGAL_pca_normals_cpp, 2},
    {"_RCGAL_AFSreconstruction_cpp", (DL_FUNC) &_RCGAL_AFSreconstruction_cpp, 1},
    {"_RCGAL_AFSreconstruction_perimeter_cpp", (DL_FUNC) &_RCGAL_AFSreconstruction_perimeter_cpp, 2},
    {"_RCGAL_Poisson_reconstruction_cpp", (DL_FUNC) &_RCGAL_Poisson_reconstruction_cpp, 6},
    {"_RCGAL_SurfMesh", (DL_FUNC) &_RCGAL_SurfMesh, 6},
    {"_RCGAL_SurfEMesh", (DL_FUNC) &_RCGAL_SurfEMesh, 6},
    {"_RCGAL_SurfQMesh", (DL_FUNC) &_RCGAL_SurfQMesh, 6},
    {"_RCGAL_Intersection2_K", (DL_FUNC) &_RCGAL_Intersection2_K, 3},
    {"_RCGAL_Intersection2_EK", (DL_FUNC) &_RCGAL_Intersection2_EK, 3},
    {"_RCGAL_Intersection_Q", (DL_FUNC) &_RCGAL_Intersection_Q, 3},
    {"_RCGAL_Difference_K", (DL_FUNC) &_RCGAL_Difference_K, 4},
    {"_RCGAL_Difference_EK", (DL_FUNC) &_RCGAL_Difference_EK, 4},
    {"_RCGAL_Difference_Q", (DL_FUNC) &_RCGAL_Difference_Q, 4},
    {"_RCGAL_Union_K", (DL_FUNC) &_RCGAL_Union_K, 3},
    {"_RCGAL_Union_EK", (DL_FUNC) &_RCGAL_Union_EK, 3},
    {"_RCGAL_Union_Q", (DL_FUNC) &_RCGAL_Union_Q, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RCGAL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
