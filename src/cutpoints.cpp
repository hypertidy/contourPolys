#include <Rcpp.h>
#include "CollectorList.h"
#include <stdlib.h>     /* abs */
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
  CollectorList outL;
  CollectorList outU; 
  
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
          outL.push_back(IntegerVector::create(c[k-1]));
          outU.push_back(IntegerVector::create(c[k])); 
        }

      }
    }}
  
  Rcpp::List out(4);
  // Rcpp::CharacterVector names(3);
// pass out the raw npts-length polygons 
// these need to be capture per nc above, but for now just bundled together
  out[0] =outX.vector();
  out[1] =outY.vector();
  out[2] = outL.vector();
  out[3] = outU.vector();
  return out; 
}

NumericVector CreateN1vector(int n, double *vec) {
  NumericVector aa(n + 1); 
  for (int ii = 0; ii < n; ii++) {
    aa[ii] = vec[ii]; 
  }
  aa[n] = aa[0]; 
  return aa; 
}
// calculate area from raw cut points, n is the length of xvec and yvec and neither are closed
// so we link the last to the first
// area is positive/negative or zero
double calculate_N1_area(int n, double *xvec, double *yvec) {
  double xsum = 0.0;
  double ysum = 0.0;
  for (int ii = 0; ii < (n-1); ii++) {
    xsum = xsum + xvec[ii] * yvec[ii+1];
  }
  for (int ii = 0; ii < (n-1); ii++) {
    ysum = ysum + yvec[ii] * xvec[ii+1];
  }
  // final sum to close the ring
  xsum = xsum + xvec[n-1] * yvec[0];
  ysum = ysum + yvec[n-1] * xvec[0];
  
  
  return (xsum - ysum)/2.0;
}
//' Filled contour
//'
//' @param x vector x coordinate
//' @param y vector x coordinate
//' @param z vector z matrix
//' @export
//' @examples
//' library(raster)
//'   data("topo", package = "contourPolys")
//'   levels <- c(-6000, -4000, -2000, 0, 2000, 4000)
//'   fc <- fcontour_sf(xFromCol(topo), rev(yFromRow(topo)), t(as.matrix(flip(topo, "y"))), c = levels)
//'   g <- purrr::map(fc[[1]], ~sf::st_polygon(list(.x)))
//'   ik <- unlist(fc[[2]])
//'   library(dplyr)
//'   x <- st_sf(geometry = st_sfc(g), kk = ik) %>% group_by(kk) %>% summarize() %>% st_cast("MULTIPOLYGON")
//'   ramp2 <- grDevices::colorRampPalette(c("#54A3D1", "#60B3EB", 
//'                                          "#78C8F0", "#98D1F5", "#B5DCFF", "#BDE1F0", "#CDEBFA", 
//'                                          "#D6EFFF", "#EBFAFF", "grey92", "grey94", "grey96", "white"))
//'     plot(x, col = ramp2(nrow(x)))
// [[Rcpp::export]]
List fcontour_sf(NumericVector x, NumericVector y, NumericMatrix z, NumericVector c) {
  int i, j, k, npt;
  int ai; 
  int nx = x.length();
  int ny = y.length();
  int nc = c.length();
  double px[8], py[8], pz[8];
  IntegerVector nn(nx * ny * nc);
  int ncount = 0; 
  NumericVector rx(8); 
  CollectorList outI; 
//  CollectorList outY;
//  CollectorList outL;
//  CollectorList outU; 
CollectorList outA;
CollectorList out_area; 
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
          //DONE 1.  calculate area and ignore any that are 0  
          //2. build sfc down here, one for each k (so put k in outer)
          //3. call union directly from here?  (or send out a decent GC with no invalid geoms?)
         double area = calculate_N1_area(npt, px, py); 
         if (!(abs(area) > 0))  {
           
          // Rprintf("area: %f\n", area);
         } else {
        
          NumericMatrix mat =  Rcpp::cbind(CreateN1vector(npt, px), CreateN1vector(npt, py)); 
          outA.push_back(mat);
          NumericVector ik(1); 
          ik[0] = k;
          outI.push_back(ik); 
          NumericVector ak(1); 
          ak[0] = area; 
          out_area.push_back(ak); 
         }
        }
        
      }
    }}
  
  Rcpp::List out(3);
  // Rcpp::CharacterVector names(3);
  // pass out the raw npts-length polygons 
  // these need to be capture per nc above, but for now just bundled together
   out[0] =outA.vector();
  out[1] =outI.vector();
  out[2] = out_area.vector();
//  out[2] = outL.vector();
//  out[3] = outU.vector();
  return out; 
}
