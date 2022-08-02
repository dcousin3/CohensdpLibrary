#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(subhyg0f1)(void *, void *, void *);
extern void F77_NAME(subhyg1f1)(void *, void *, void *, void *);
extern void F77_NAME(subhyg2f1)(void *, void *, void *, void *, void *);
extern void F77_NAME(subkprimecdf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(subkprimeidf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(subkprimepdf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublprimecdf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublprimeidf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublprimepdf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublsecondcdf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublsecondidf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sublsecondpdf)(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"subhyg0f1",     (DL_FUNC) &F77_NAME(subhyg0f1),     3},
    {"subhyg1f1",     (DL_FUNC) &F77_NAME(subhyg1f1),     4},
    {"subhyg2f1",     (DL_FUNC) &F77_NAME(subhyg2f1),     5},
    {"subkprimecdf",  (DL_FUNC) &F77_NAME(subkprimecdf),  8},
    {"subkprimeidf",  (DL_FUNC) &F77_NAME(subkprimeidf),  8},
    {"subkprimepdf",  (DL_FUNC) &F77_NAME(subkprimepdf),  8},
    {"sublprimecdf",  (DL_FUNC) &F77_NAME(sublprimecdf),  7},
    {"sublprimeidf",  (DL_FUNC) &F77_NAME(sublprimeidf),  7},
    {"sublprimepdf",  (DL_FUNC) &F77_NAME(sublprimepdf),  7},
    {"sublsecondcdf", (DL_FUNC) &F77_NAME(sublsecondcdf), 8},
    {"sublsecondidf", (DL_FUNC) &F77_NAME(sublsecondidf), 8},
    {"sublsecondpdf", (DL_FUNC) &F77_NAME(sublsecondpdf), 8},
    {NULL, NULL, 0}
};

void R_init_CohensdpLibrary(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
