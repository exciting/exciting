#include "util.h"

extern XC(func_info_type) XC(func_info_lca_omc);
extern XC(func_info_type) XC(func_info_lca_lch);


const XC(func_info_type) *XC(lca_known_funct)[] = {
  &XC(func_info_lca_omc),
  &XC(func_info_lca_lch),
  NULL
};
