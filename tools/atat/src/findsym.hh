#ifndef _FINDSYM_H_
#define _FINDSYM_H_

#include "xtalutil.h"

void apply_symmetry(Array<rVector3d> *b, const rMatrix3d &point_op, const rVector3d &trans, const Array<rVector3d> &a);
void apply_symmetry(MultiCluster *b, const rMatrix3d &point_op, const rVector3d &trans, const MultiCluster &a);

int equivalent_by_symmetry(const rVector3d &a, const rVector3d &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
int equivalent_by_symmetry(const Array<rVector3d> &a, const Array<rVector3d> &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
int equivalent_by_symmetry(const MultiCluster &a, const MultiCluster &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

void calc_concentration(int *pmin_type, int *pmax_type, Array<Real> *pconc, const Array<int> &atom_type);

// assumes that the structure with smaller cell is given with the smallest;
// possible cell;
int equivalent_by_symmetry(const Structure &a, const Structure &b,
                           const rMatrix3d &cell,
                           const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

class SpaceGroup {
 public:
  rMatrix3d cell;
  Array<rMatrix3d> point_op;
  Array<rVector3d> trans;
  SpaceGroup(void): cell(), point_op(), trans() {}
  SpaceGroup(const SpaceGroup &a): cell(a.cell), point_op(a.point_op), trans(a.trans) {}
  SpaceGroup(const rMatrix3d &_cell,
                       const Array<rMatrix3d> &_point_op,
                       const Array<rVector3d> &_trans):
    cell(_cell), point_op(_point_op), trans(_trans) {}
  int operator() (const rVector3d &a, const rVector3d &b) const {
    return equivalent_by_symmetry(a,b,cell,point_op,trans);
  }
  int operator() (const Array<rVector3d> &a, const Array<rVector3d> &b) const {
    return equivalent_by_symmetry(a,b,cell,point_op,trans);
  }
  int operator() (const MultiCluster &a, const MultiCluster &b) const {
    return equivalent_by_symmetry(a,b,cell,point_op,trans);
  }
  int operator() (const Structure &a, const Structure &b) const {
    return equivalent_by_symmetry(a,b,cell,point_op,trans);
  }
};

void find_pointgroup(Array<rMatrix3d> *point_op, const rMatrix3d &cell);

void pointgroup_from_spacegroup(Array<rMatrix3d> *pointgroup, const Array<rMatrix3d> &point_op);

void pointgroup_from_spacegroup(const Structure &str, Array<rMatrix3d> *pointgroup,  const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

void find_primitive_unit_cell(rMatrix3d *pcell, const Structure &lat);

void find_spacegroup(Array<rMatrix3d> *point_op, Array<rVector3d> *trans,
                      const rMatrix3d &cell,
                      const Array<rVector3d> &atom_pos, const Array<int> &atom_type);

rMatrix3d find_almost_reduced_cell(const rMatrix3d &cell);

void find_all_equivalent_cell(Array<rMatrix3d> *psupercell, const rMatrix3d &cell, Real radius);

void find_supercells(Array<rMatrix3d> *supercell, int min_volume, int max_volume, const rMatrix3d &unitcell, const Array<rMatrix3d> &pointgroup);

void find_supercells_2D(Array<rMatrix3d> *supercell, int min_volume, int max_volume, const rMatrix3d &unitcell, const Array<rMatrix3d> &pointgroup);

void sort_supercells(Array<rMatrix3d> *p_supercell, const Array<rMatrix3d> &pointgroup);

int find_point_op_type(rVector3d *pspecial_dir, const rMatrix3d &point_op);

void find_atom_point_group(Array<rMatrix3d> *patom_point_group,
			   const rVector3d &atom_pos,
			   const rMatrix3d &cell, const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);

enum PointGroupType {
  pt_grp_none,
  pt_grp_1,
  pt_grp_1_,
  pt_grp_2,
  pt_grp_m,
  pt_grp_2om,
  pt_grp_222,
  pt_grp_2mm,
  pt_grp_mmm,
  pt_grp_4,
  pt_grp_4_,
  pt_grp_4om,
  pt_grp_422,
  pt_grp_4mm,
  pt_grp_4_2m,
  pt_grp_4ommm,
  pt_grp_3,
  pt_grp_3_,
  pt_grp_32,
  pt_grp_3m,
  pt_grp_3_m,
  pt_grp_6,
  pt_grp_6_,
  pt_grp_6om,
  pt_grp_622,
  pt_grp_6mm,
  pt_grp_6_m2,
  pt_grp_6ommm,
  pt_grp_23,
  pt_grp_m3_,
  pt_grp_432,
  pt_grp_4_3m,
  pt_grp_m3_m
};

PointGroupType find_point_group_type(const Array<rMatrix3d> &point_op);

int contains_pure_translations(const Array<rMatrix3d> &point_op, const Array<rVector3d> &trans);
//int contains_pure_translations(const Structure &str);
int contains_pure_translations(const Structure &str, rMatrix3d unitcell);
int contains_pure_translations_or_lexico_successor(const Structure &str, rMatrix3d unitcell);

void generate_space_group(SpaceGroup *p_fullgroup, const SpaceGroup &generator);

#endif
