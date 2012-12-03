#ifndef __LATTYPE_H__
#define __LATTYPE_H__

#include "vectmac.h"
#include "stringo.h"

enum BravaisType {BT_P=1, BT_F=2, BT_I=4, BT_C=8, BT_R=16, BT_c=32, BT_t=64, BT_o=128, BT_m=256, BT_a=512, BT_h=1024,
                  BT_cP=BT_c|BT_P,
		  BT_cF=BT_c|BT_F,
		  BT_cI=BT_c|BT_I,
		  BT_tP=BT_t|BT_P,
		  BT_tI=BT_t|BT_I,
		  BT_oP=BT_o|BT_P,
		  BT_oF=BT_o|BT_F,
		  BT_oI=BT_o|BT_I,
		  BT_oC=BT_o|BT_C,
		  BT_mP=BT_m|BT_P,
		  BT_mC=BT_m|BT_C,
		  BT_mI=BT_m|BT_I, // oddity in the int'l tables of xtal;
		  BT_aP=BT_a|BT_P,
		  BT_hP=BT_h|BT_P,
		  BT_hR=BT_h|BT_R
};

rMatrix3d find_reduced_cell(const rMatrix3d &cell);
int is_reduced_cell(const rMatrix3d &cell);
BravaisType find_bravais_type(rMatrix3d *p_conventional, const rMatrix3d &cell);
void print_bravais_type(AutoString *s, BravaisType bt);
rMatrix3d find_symmetric_cell(const rMatrix3d &cell);
void find_smallest_supercell_enclosing_sphere(rMatrix3d *psupercell, const rMatrix3d cell, Real r);

#endif
