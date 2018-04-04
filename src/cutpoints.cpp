#include <Rcpp.h>
#include "CollectorList.h"
using namespace Rcpp;

void FindCutPoints(double low, double high,
                   double x1, double y1, double z1,
                   double x2, double y2, double z2,
                   double *x, double *y, double *z,
                   int *npt)
{
  
  double c;
  
  if (z1 > z2 ) {
    if (z2 > high || z1 < low) return;
    if (z1 < high) {
      x[*npt] = x1;
      y[*npt] = y1;
      z[*npt] = z1;
      ++*npt;
    } else if (z1 == R_PosInf) {
      x[*npt] = x2;
      y[*npt] = y1;
      z[*npt] = z2;
      ++*npt;
    } else { /* z1 >= high, z2 in range */
      c = (z1 - high) / (z1 - z2);
      x[*npt] = x1 + c * (x2 - x1);
      y[*npt] = y1;
      z[*npt] = z1 + c * (z2 - z1);
      ++*npt;
    }
    if (z2 == R_NegInf) {
      x[*npt] = x1;
      y[*npt] = y1;
      z[*npt] = z1;
      ++*npt;
    } else if (z2 <= low) { /* and z1 in range */
      c = (z2 -low) / (z2 - z1);
      x[*npt] = x2 - c * (x2 - x1);
      y[*npt] = y1;
      z[*npt] = z2 - c * (z2 - z1);
      ++*npt;
    }
  } else if (z1 < z2) {
    if (z2 < low || z1 > high) return;
    if (z1 > low) {
      x[*npt] = x1;
      y[*npt] = y1;
      z[*npt] = z1;
      ++*npt;
    } else if (z1 == R_NegInf) {
      x[*npt] = x2;
      y[*npt] = y1;
      z[*npt] = z2;;
      ++*npt;
    } else { /* and z2 in range */
      c = (z1 - low) / (z1 - z2);
      x[*npt] = x1 + c * (x2 - x1);
      y[*npt] = y1;
      z[*npt] = z1 + c * (z2 - z1);
      ++*npt;
    }
    if (z2 < high) {
#ifdef OMIT
      /* Don't repeat corner vertices */
      x[*npt] = x2;
      y[*npt] = y2;
      z[*npt] = z2;
      ++*npt;
#endif
    } else if (z2 == R_PosInf) {
      x[*npt] = x1;
      y[*npt] = y1;
      z[*npt] = z1;
      ++*npt;
    } else { /* z2 high, z1 in range */
      c = (z2 - high) / (z2 - z1);
      x[*npt] = x2 - c * (x2 - x1);
      y[*npt] = y1;
      z[*npt] = z2 - c * (z2 - z1);
      ++*npt;
    }
  } else {
    if(low <= z1 && z1 <= high) {
      x[*npt] = x1;
      y[*npt] = y1;
      z[*npt] = z1;
      ++*npt;
#ifdef OMIT
      /* Don't repeat corner vertices */
      x[*npt] = x2;
      y[*npt] = y2;
      z[*npt] = z2;
      ++*npt;
#endif
    }
  }
}


void FindPolygonVertices(double low, double high,
                         double x1, double x2, double y1, double y2,
                         double z11, double z21, double z12, double z22,
                         double *x, double *y, double *z, int *npt)
{
  *npt = 0;
  FindCutPoints(low, high, x1,  y1,  z11, x2,  y1,  z21, x, y, z, npt);
  FindCutPoints(low, high, y1,  x2,  z21, y2,  x2,  z22, y, x, z, npt);
  FindCutPoints(low, high, x2,  y2,  z22, x1,  y2,  z12, x, y, z, npt);
  FindCutPoints(low, high, y2,  x1,  z12, y1,  x1,  z11, y, x, z, npt);
}

NumericVector CreateNvector(int n, double *vec) {
  NumericVector aa(n); 
  for (int ii = 0; ii < n; ii++) {
    aa[ii] = vec[ii]; 
  }
  return aa; 
}
//' Filled contour
//'
//' @param x vector x coordinate
//' @param y vector x coordinate
//' @param z vector z matrix
//' @export
// [[Rcpp::export]]
List fcontour(NumericVector x, NumericVector y, NumericMatrix z, NumericVector c) {
  int i, j, k, npt;
  int ai; 
  int nx = x.length();
  int ny = y.length();
  int nc = c.length();
  double px[8], py[8], pz[8];
  IntegerVector nn(nx * ny * nc);
  int ncount = 0; 
  NumericVector rx(8); 
  CollectorList outX; 
  CollectorList outY; 
  
  for (i = 1; i < nx; i++) {
    for (j = 1; j < ny; j++) {
      for (k = 1; k < nc ; k++) {
        FindPolygonVertices(c[k - 1], c[k],
                            x[i - 1], x[i],
                                       y[j - 1], y[j],
                                                  z[i - 1 + (j - 1) * nx],
                                                   z[i + (j - 1) * nx],
                                                    z[i - 1 + j * nx],
                                                     z[i + j * nx],
                                                      px, py, pz, &npt);
        
        
        if (npt > 2) {
          outX.push_back(CreateNvector(npt, px)); 
          outY.push_back(CreateNvector(npt, py)); 
        }

      }
    }}
  
  Rcpp::List out(3);
  // Rcpp::CharacterVector names(3);
// pass out the raw npts-length polygons 
// these need to be capture per nc above, but for now just bundled together
  out[0] =outX.vector();
  out[1] =outY.vector();

  return out; 
}
