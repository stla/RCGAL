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
// PolyMesh
Rcpp::List PolyMesh(const Rcpp::NumericMatrix points, const Rcpp::IntegerMatrix faces);
RcppExport SEXP _RCGAL_PolyMesh(SEXP pointsSEXP, SEXP facesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type faces(facesSEXP);
    rcpp_result_gen = Rcpp::wrap(PolyMesh(points, faces));
    return rcpp_result_gen;
END_RCPP
}
// SurfMesh
Rcpp::List SurfMesh(const Rcpp::NumericMatrix points, const Rcpp::IntegerMatrix faces);
RcppExport SEXP _RCGAL_SurfMesh(SEXP pointsSEXP, SEXP facesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type faces(facesSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfMesh(points, faces));
    return rcpp_result_gen;
END_RCPP
}
// SurfMeshWithNormals
Rcpp::List SurfMeshWithNormals(const Rcpp::NumericMatrix points, const Rcpp::IntegerMatrix faces);
RcppExport SEXP _RCGAL_SurfMeshWithNormals(SEXP pointsSEXP, SEXP facesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type faces(facesSEXP);
    rcpp_result_gen = Rcpp::wrap(SurfMeshWithNormals(points, faces));
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
    {"_RCGAL_jet_normals_cpp", (DL_FUNC) &_RCGAL_jet_normals_cpp, 2},
    {"_RCGAL_pca_normals_cpp", (DL_FUNC) &_RCGAL_pca_normals_cpp, 2},
    {"_RCGAL_AFSreconstruction_cpp", (DL_FUNC) &_RCGAL_AFSreconstruction_cpp, 1},
    {"_RCGAL_AFSreconstruction_perimeter_cpp", (DL_FUNC) &_RCGAL_AFSreconstruction_perimeter_cpp, 2},
    {"_RCGAL_Poisson_reconstruction_cpp", (DL_FUNC) &_RCGAL_Poisson_reconstruction_cpp, 6},
    {"_RCGAL_PolyMesh", (DL_FUNC) &_RCGAL_PolyMesh, 2},
    {"_RCGAL_SurfMesh", (DL_FUNC) &_RCGAL_SurfMesh, 2},
    {"_RCGAL_SurfMeshWithNormals", (DL_FUNC) &_RCGAL_SurfMeshWithNormals, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_RCGAL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
