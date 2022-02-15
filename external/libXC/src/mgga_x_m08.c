/*
 Copyright (C) 2006-2007 M.A.L. Marques

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_MGGA_X_M08_HX       219 /* M08-HX functional of Minnesota */
#define XC_MGGA_X_M08_SO       220 /* M08-SO functional of Minnesota */
/*#define XC_MGGA_X_M11          221 *//* M11 functional of Minnesota */
/*#define XC_MGGA_X_M11_L        222 *//* M11-L functional of Minnesota */

static const FLOAT a_m08_hx[12] = {
   1.3340172e+00, -9.4751087e+00, -1.2541893e+01,  9.1369974e+00,  3.4717204e+01,  5.8831807e+01,
   7.1369574e+01,  2.3312961e+01,  4.8314679e+00, -6.5044167e+00, -1.4058265e+01,  1.2880570e+01
};
static const FLOAT b_m08_hx[12] = {
  -8.5631823e-01,  9.2810354e+00,  1.2260749e+01, -5.5189665e+00, -3.5534989e+01, -8.2049996e+01,
  -6.8586558e+01,  3.6085694e+01, -9.3740983e+00, -5.9731688e+01,  1.6587868e+01,  1.3993203e+01
};

static const FLOAT a_m08_so[12] = {
  -3.4888428e-01, -5.8157416e+00,  3.7550810e+01,  6.3727406e+01, -5.3742313e+01, -9.8595529e+01,
   1.6282216e+01,  1.7513468e+01, -6.7627553e+00,  1.1106658e+01,  1.5663545e+00,  8.7603470e+00
};
static const FLOAT b_m08_so[12] = {
   7.8098428e-01,  5.4538178e+00, -3.7853348e+01, -6.2295080e+01,  4.6713254e+01,  8.7321376e+01,
   1.6053446e+01,  2.0126920e+01, -4.0343695e+01, -5.8577565e+01,  2.0890272e+01,  1.0946903e+01
};

static const FLOAT a_m11[12] = {
  -0.18399900e+00, -1.39046703e+01,  1.18206837e+01,  3.10098465e+01, -5.19625696e+01,  1.55750312e+01,
  -6.94775730e+00, -1.58465014e+02, -1.48447565e+00,  5.51042124e+01, -1.34714184e+01,  0.00000000e+00
};
static const FLOAT b_m11[12] = {
   0.75599900e+00,  1.37137944e+01, -1.27998304e+01, -2.93428814e+01,  5.91075674e+01, -2.27604866e+01,
  -1.02769340e+01,  1.64752731e+02,  1.85349258e+01, -5.56825639e+01,  7.47980859e+00,  0.00000000e+00
};

static const FLOAT a_m11_l[12] = {
   8.121131e-01,  1.738124e+01,  1.154007e+00,  6.869556e+01,  1.016864e+02, -5.887467e+00, 
   4.517409e+01, -2.773149e+00, -2.617211e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00
};
static const FLOAT b_m11_l[12] = {
   1.878869e-01, -1.653877e+01,  6.755753e-01, -7.567572e+01, -1.040272e+02,  1.831853e+01,
  -5.573352e+01, -3.520210e+00,  3.724276e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00
};
static const FLOAT c_m11_l[12] = {
  -4.386615e-01, -1.214016e+02, -1.393573e+02, -2.046649e+00,  2.804098e+01, -1.312258e+01,
  -6.361819e+00, -8.055758e-01,  3.736551e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00
};
static const FLOAT d_m11_l[12] = {
   1.438662e+00,  1.209465e+02,  1.328252e+02,  1.296355e+01,  5.854866e+00, -3.378162e+00,
  -4.423393e+01,  6.844475e+00,  1.949541e+01,  0.000000e+00,  0.000000e+00,  0.000000e+00
};

typedef struct{
  const FLOAT *a, *b;
} mgga_x_m08_params;


static void
mgga_x_m08_init(XC(func_type) *p)
{
  mgga_x_m08_params *params;

  assert(p != NULL);

  p->n_func_aux  = 2;
  p->func_aux    = (XC(func_type) **) malloc(2*sizeof(XC(func_type) *));
  p->func_aux[0] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));
  p->func_aux[1] = (XC(func_type) *)  malloc(  sizeof(XC(func_type)));

  XC(func_init)(p->func_aux[0], XC_GGA_X_PBE,  p->nspin);
  XC(func_init)(p->func_aux[1], XC_GGA_X_RPBE, p->nspin);

  XC(gga_x_pbe_set_params) (p->func_aux[1], 0.552, 10.0/81.0);

  assert(p->params == NULL);
  p->params = malloc(sizeof(mgga_x_m08_params));
  params = (mgga_x_m08_params *) (p->params);

  switch(p->info->number){
  case XC_MGGA_X_M08_HX: 
    params->a = a_m08_hx;
    params->b = b_m08_hx;
    break;
  case XC_MGGA_X_M08_SO:
    params->a = a_m08_so;
    params->b = b_m08_so;
    break;
    /*
  case XC_MGGA_X_M11:
    params->a = a_m11;
    params->b = b_m11;
    break;
  case XC_MGGA_X_M11_L:
    params->a = a_m11_l;
    params->b = b_m11_l;
    break;
    */
  default:
    fprintf(stderr, "Internal error in mgga_x_m08\n");
    exit(1);
  }
}


static void 
func(const XC(func_type) *pt, XC(mgga_work_x_t) *r)
{
  mgga_x_m08_params *params;

  FLOAT ep_f, ep_dfdx, ep_d2fdx2, er_f, er_dfdx, er_d2fdx2;
  FLOAT fw1, fw2, dfw1dt, dfw2dt;

  assert(pt != NULL && pt->params != NULL);
  params = (mgga_x_m08_params *) (pt->params);
  
  XC(gga_x_pbe_enhance)(pt->func_aux[0], r->order, r->x, &ep_f, &ep_dfdx, &ep_d2fdx2);
  XC(gga_x_pbe_enhance)(pt->func_aux[1], r->order, r->x, &er_f, &er_dfdx, &er_d2fdx2);
  
  XC(mgga_series_w)(r->order, 12, params->a, r->t, &fw1, &dfw1dt);
  XC(mgga_series_w)(r->order, 12, params->b, r->t, &fw2, &dfw2dt);

  r->f = ep_f*fw1 + er_f*fw2 ;

  if(r->order < 1) return;

  r->dfdx = ep_dfdx*fw1 + er_dfdx*fw2;
  r->dfdt = ep_f*dfw1dt + er_f*dfw2dt;
  r->dfdu = 0.0;

  if(r->order < 2) return;

}


#include "work_mgga_x.c"


XC(func_info_type) XC(func_info_mgga_x_m08_hx) = {
  XC_MGGA_X_M08_HX,
  XC_EXCHANGE,
  "M08-HX functional of Minnesota",
  XC_FAMILY_MGGA,
  "Y Zhao and DG Truhlar, J. Chem. Theory Comput. 4, 1849-1868 (2008)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

XC(func_info_type) XC(func_info_mgga_x_m08_so) = {
  XC_MGGA_X_M08_SO,
  XC_EXCHANGE,
  "M08-SO functional of Minnesota",
  XC_FAMILY_MGGA,
  "Y Zhao and DG Truhlar, J. Chem. Theory Comput. 4, 1849-1868 (2008)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

/* Warning: to be added in the future
XC(func_info_type) XC(func_info_mgga_x_m11) = {
  XC_MGGA_X_M11,
  XC_EXCHANGE,
  "M11 functional of Minnesota",
  XC_FAMILY_MGGA,
  "R Peverati, and DG Truhlar, J. Phys. Chem. Lett. 2, 2810 (2011)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

XC(func_info_type) XC(func_info_mgga_x_m11_l) = {
  XC_MGGA_X_M11_L,
  XC_EXCHANGE,
  "M11-L functional of Minnesota",
  XC_FAMILY_MGGA,
  "R Peverati, and DG Truhlar, J. Phys. Chem. Lett. 3, 117 (2012)",
  XC_FLAGS_3D | XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC,
  1e-32, 1e-32, 1e-32, 1e-32,
  mgga_x_m08_init,
  NULL, NULL, NULL,
  work_mgga_x,
};

*/
