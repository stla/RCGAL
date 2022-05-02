#ifndef _RCGALHEADER_
#include "rcgal.h"
#endif

// double sepsilon = sqrt(std::numeric_limits<double>::epsilon());

bool approxEqual(double x, double y, double epsilon) {
  return fabs(x - y) <= epsilon;
}

bool approxEqualVectors(Vector3 v, Vector3 w, double epsilon) {
  return approxEqual(v.x(), w.x(), epsilon) &&
    approxEqual(v.y(), w.y(), epsilon) &&
    approxEqual(v.z(), w.z(), epsilon);
}

Rcpp::String stringPair(const size_t i, const size_t j) {
  const size_t i0 = std::min(i, j);
  const size_t i1 = std::max(i, j);
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
    std::to_string(i0), "-", std::to_string(i1));
  return Rcpp::collapse(stringids);
}

std::array<Rcpp::String, 3> triangleEdges(const size_t i0,
                                          const size_t i1,
                                          const size_t i2) {
  return {stringPair(i0, i1), stringPair(i0, i2), stringPair(i1, i2)};
}

Rcpp::String stringTriple(const size_t i, const size_t j, const size_t k) {
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
    std::to_string(i), "-", std::to_string(j), "-", std::to_string(k));
  return Rcpp::collapse(stringids);
}

double volume_under_triangle(Dvector v0, Dvector v1, Dvector v2) {
  double x0 = v0(0);
  double y0 = v0(1);
  double z0 = v0(2);
  double x1 = v1(0);
  double y1 = v1(1);
  double z1 = v1(2);
  double x2 = v2(0);
  double y2 = v2(1);
  double z2 = v2(2);
  return (z0 + z1 + z2) *
    (x0 * y1 - x1 * y0 + x1 * y2 - x2 * y1 + x2 * y0 - x0 * y2) / 6.0;
}

void mark_domains0(CDT& ct,
                   CDFace_handle start,
                   int index,
                   std::list<CDT::Edge>& border) {
  if(start->info().nesting_level != -1) {
    return;
  }
  std::list<CDFace_handle> queue;
  queue.push_back(start);
  while(!queue.empty()) {
    CDFace_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++) {
        CDT::Edge e(fh, i);
        CDFace_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1) {
          if(ct.is_constrained(e))
            border.push_back(e);
          else
            queue.push_back(n);
        }
      }
    }
  }
}

void mark_domains(CDT& cdt)
{
  for(CDFace_handle f : cdt.all_face_handles()){
    f->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains0(cdt, cdt.infinite_face(), 0, border);
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDFace_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains0(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}


