#include "tensor.h"
#include "findsym.h"

void apply_symmetry(rTensor *ptt, const rMatrix3d &op, const rTensor &t);
int symmetrize(rTensor *psymt, const rTensor &t, const SpaceGroup &space_group);
//void symmetrize(rTensor *psymt, const Array<Array<int> > &flip);
void transpose(rTensor *psymt, const rTensor &t, int row, int col);
void transpose(rTensor *psymt, int row, int col);
void transpose(rTensor *psymt, const rTensor &t, const Array<int> &map);
int symmetrize(rTensor *psymt, const rTensor &t, const SpaceGroup &space_group, const Array<Array<int> > &flip);
void calc_flip_group(Array<Array<int> > *pgroup, int rank, const Array<Array<int> > &flip);
void symmetrize_flip(rTensor *psymt, const Array<Array<int> > &flipgroup);
void find_symmetric_basis(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip);
void find_symmetry_breaking_basis(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip);

void find_symmetric_basis_projection(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group, const Array<Array<int> > &flip, const Array<rTensor> &proj, int paral );
void all_flip(Array<Array<int> > *pflip, int rank);
void calc_sym_harmonics(Array<rTensor> *pbasis, int rank, const SpaceGroup &space_group);
void print_harm(const rTensor &t);
