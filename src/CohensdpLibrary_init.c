#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(subhyg0f1)    (double *, double *, double *);
extern void F77_NAME(subhyg1f1)    (double *, double *, double *, double *);
extern void F77_NAME(subhyg2f1)    (double *, double *, double *, double *, double *);
extern void F77_NAME(sublprimepdf) (double *, double *, double *, double *, short *,  short *, double *);
extern void F77_NAME(sublprimecdf) (double *, double *, double *, double *, short *,  short *, double *);
extern void F77_NAME(sublprimeidf) (double *, double *, double *, double *, short *,  short *, double *);
extern void F77_NAME(subkprimepdf) (double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(subkprimecdf) (double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(subkprimeidf) (double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(sublsecondpdf)(double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(sublsecondcdf)(double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(sublsecondidf)(double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(subfbdeltafromobsdpobsrpdf)(double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(subfbdeltafromobsdpobsrcdf)(double *, double *, double *, double *, double *, short *, short *, double *);
extern void F77_NAME(subfbdeltafromobsdpobsridf)(double *, double *, double *, double *, double *, short *, short *, double *);


static const R_FortranMethodDef FortranEntries[] = {
    {"subhyg0f1",                  (DL_FUNC) &F77_NAME(subhyg0f1),                  3},
    {"subhyg1f1",                  (DL_FUNC) &F77_NAME(subhyg1f1),                  4},
    {"subhyg2f1",                  (DL_FUNC) &F77_NAME(subhyg2f1),                  5},
    {"sublprimepdf",               (DL_FUNC) &F77_NAME(sublprimepdf),               7},
    {"sublprimecdf",               (DL_FUNC) &F77_NAME(sublprimecdf),               7},
    {"sublprimeidf",               (DL_FUNC) &F77_NAME(sublprimeidf),               7},
    {"subkprimepdf",               (DL_FUNC) &F77_NAME(subkprimepdf),               8},
    {"subkprimecdf",               (DL_FUNC) &F77_NAME(subkprimecdf),               8},
    {"subkprimeidf",               (DL_FUNC) &F77_NAME(subkprimeidf),               8},
    {"sublsecondpdf",              (DL_FUNC) &F77_NAME(sublsecondpdf),              8},
    {"sublsecondcdf",              (DL_FUNC) &F77_NAME(sublsecondcdf),              8},
    {"sublsecondidf",              (DL_FUNC) &F77_NAME(sublsecondidf),              8},
    {"subfbdeltafromobsdpobsrpdf", (DL_FUNC) &F77_NAME(subfbdeltafromobsdpobsrpdf), 8},
    {"subfbdeltafromobsdpobsrcdf", (DL_FUNC) &F77_NAME(subfbdeltafromobsdpobsrcdf), 8},
    {"subfbdeltafromobsdpobsridf", (DL_FUNC) &F77_NAME(subfbdeltafromobsdpobsridf), 8},
    {NULL, NULL, 0} 
};

void R_init_CohensdpLibrary(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE); 
}
