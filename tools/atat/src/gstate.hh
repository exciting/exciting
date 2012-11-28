#ifndef _GSTATE_H_
#define _GSTATE_H_

#include "arraylist.h"
#include "integer.h"

void calc_convex_hull(Array<int> *p_gs_str, const Array<Real> &concentration, const Array<Real> &energy);
void calc_convex_hull(Array<Real> *p_gs_concentration, Array<Real> *p_gs_energy,
			const Array<Real> &concentration, const Array<Real> &energy);
void convex_hull_from_gs(Array<Real> *p_gs_concentration, Array<Real> *p_gs_energy, const Array<int> &gs_str,
			 const Array<Real> &concentration, const Array<Real> &energy);
int ground_state_problem(Array<int> *p_problem, const Array<int> &ground_state, const Array<Real> &concentration, const Array<Real> &fitted_energy, Real minc=-MAXFLOAT, Real maxc=MAXFLOAT);
int ground_state_problem(Array<int> *p_problem, const Array<Real> &gs_concentration, const Array<Real> &gs_energy, const Array<Real> &concentration, const Array<Real> &energy, Real minc=-MAXFLOAT, Real maxc=MAXFLOAT);
void calc_formation(Array<Real> *p_form_value, const Array<Real> &concentration, const Array<Real> &value, Real *e_c0=NULL, Real *e_c1=NULL);
void calc_hull_energy(Array<Real> *hull_energy, const Array<int> &ground_state, const Array<Real> &concentration, const Array<Real> &energy);

#endif
