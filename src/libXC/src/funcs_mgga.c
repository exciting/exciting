#include "util.h"

extern XC(func_info_type) XC(func_info_mgga_x_lta);
extern XC(func_info_type) XC(func_info_mgga_x_tpss);
extern XC(func_info_type) XC(func_info_mgga_x_m06l);
extern XC(func_info_type) XC(func_info_mgga_x_gvt4);
extern XC(func_info_type) XC(func_info_mgga_x_tau_hcth);
extern XC(func_info_type) XC(func_info_mgga_x_br89);
extern XC(func_info_type) XC(func_info_mgga_x_bj06);
extern XC(func_info_type) XC(func_info_mgga_x_tb09);
extern XC(func_info_type) XC(func_info_mgga_x_rpp09);
extern XC(func_info_type) XC(func_info_mgga_c_tpss);
extern XC(func_info_type) XC(func_info_mgga_c_vsxc);


const XC(func_info_type) *XC(mgga_known_funct)[] = {
  &XC(func_info_mgga_x_lta),
  &XC(func_info_mgga_x_tpss),
  &XC(func_info_mgga_x_m06l),
  &XC(func_info_mgga_x_gvt4),
  &XC(func_info_mgga_x_tau_hcth),
  &XC(func_info_mgga_x_br89),
  &XC(func_info_mgga_x_bj06),
  &XC(func_info_mgga_x_tb09),
  &XC(func_info_mgga_x_rpp09),
  &XC(func_info_mgga_c_tpss),
  &XC(func_info_mgga_c_vsxc),
  NULL
};
