#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _contourPolys_fcontour(SEXP, SEXP, SEXP, SEXP);
extern SEXP _contourPolys_fcontour_sf(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_contourPolys_fcontour",    (DL_FUNC) &_contourPolys_fcontour,    4},
    {"_contourPolys_fcontour_sf", (DL_FUNC) &_contourPolys_fcontour_sf, 4},
    {NULL, NULL, 0}
};

void R_init_contourPolys(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
