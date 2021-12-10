#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "JLPM.h"

static R_FortranMethodDef FortRout[] = {
  {"irtsre", (DL_FUNC) &F77_SUB(irtsre), 65},
  {"proba_irtsre", (DL_FUNC) &F77_SUB(proba_irtsre), 41},
  {"loglik1", (DL_FUNC) &F77_SUB(loglik1), 48},
  {"loglik2", (DL_FUNC) &F77_SUB(loglik2), 60},
  {NULL, NULL, 0}
};


void R_init_JLPM(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
