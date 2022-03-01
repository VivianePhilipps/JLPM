#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "JLPM.h"

static R_FortranMethodDef FortRout[] = {
  {"loglik", (DL_FUNC) &F77_SUB(loglik), 47},
  {NULL, NULL, 0}
};


void R_init_JLPM(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, NULL, FortRout, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
