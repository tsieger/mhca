#include "mhclust.h"
#include <R_ext/Rdynload.h>

void R_init_mhca(DllInfo *info) {
    R_RegisterCCallable("mhca", "mhclust",  (DL_FUNC) &mhclust_);
}