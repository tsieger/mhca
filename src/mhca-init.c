#include "mhclust.h"
#include <R_ext/Rdynload.h>

void R_init_mhca(DllInfo *info) {
    // This is the easier way described at http://r-pkgs.had.co.nz/src.html (Hadley Wickham):
    //R_RegisterCCallable("mhca", "mhclust",  (DL_FUNC) &mhclust_);
    // However, "R CMD check --as-cran" is not happy with it, complaining about
    // missing calls to: `R_registerRoutines`, `R_useDynamicSymbols`.

    // Therefore, we must go the more complex way:
    static const R_CallMethodDef callMethods[]  = {
        {"mhclust_", (DL_FUNC) &mhclust_,20},
        {NULL, NULL, 0}
    };
    R_registerRoutines(info,NULL,callMethods,NULL,NULL);
    R_useDynamicSymbols(info,FALSE);
    R_forceSymbols(info,TRUE);
}
