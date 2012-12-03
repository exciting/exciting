#ifndef _PARSE_H_
#define _PARSE_H_

#include <iostream.h>
#include <strstream.h>
#include <ctype.h>
#include "xtalutil.h"
#include "stringo.h"
#include "clus_str.h"
#include "linalg.h"
#include "gceutil.h"
#include "getvalue.h"

typedef Array<AutoString> ArrayAutoString;
typedef Array<int> Arrayint;

int is_eol(istream &file);

void read_cell(rMatrix3d *pcell, istream &file);

// reads in a lattice from a file;
// IN:
//  file: any input stream (including cin). the format is:
//
//   First, the coordinate system is specified, either as
//    a b c alpha beta gamma
//   or as:
//    ax ay az
//    bx by bz
//    cx cy cz
//   Then the lattice vectors are listed, expressed in the coordinate system just defined:
//    ua ub uc
//    va vb vc
//    wa wb wc
//   Finally, atom positions and atomic symbols are given, expressed in the same coordinate system
//   as the lattice vectors:
//    atom1a atom1b atom1c atom1type1,atom1type2,etc.
//    atom2a atom2b atom2c atom2type2,atom2type2,etc.
//    etc.
// OUT:
//  pcell: Each column of this matrix contains a lattice vector in cartesian;
//  patom_pos: cartesian positions of each atom;
//  psite_type: for each atom, index into *psite_site_type;
//  psite_type_list: an array of indices into *patom_label;
//  patom_label: an array of atomic symbols describing each type of atom;
//  paxes: (optional) the coordinate system specified at the beginning of the file
//   (one lattice vector per column);
//
// An example:
//  The fcc lattice for the MgO-CaO system:
//   3.2 3.2 3.2 90 90 90
//   0   0.5 0.5
//   0.5 0   0.5
//   0.5 0.5 0
//   0 0 0 Ca,Mg
//   0.5 0.5 0.5 O
//  then the variables are initialized as follows:
//  *pcell:
//    0   1.6 1.6
//    1.6 0   1.6
//    1.6 1.6 0
//  *patom_pos:
//    0   0   0
//    1.6 1.6 1.6
//  *psite_type:
//    0 1
//  *psite_type_list
//    0 1
//    2
//  *patom_label:
//   Ca Mg O
// End of example;

void parse_lattice_file(rMatrix3d *pcell, Array<rVector3d> *patom_pos, Array<int> *psite_type, Array<Arrayint> *psite_type_list, Array<AutoString> *patom_label, istream &file, rMatrix3d *paxes=NULL);

// reads in a structure from a file;
// IN:
//  atom_label: An array of atomic symbols describing each type of atom that will appear in the file;
//  file: any input stream (including cin). the format is:
//
//   First, the coordinate system is specified, either as
//    a b c alpha beta gamma
//   or as:
//    ax ay az
//    bx by bz
//    cx cy cz
//   Then the lattice vectors are listed, expressed in the coordinate system just defined:
//    ua ub uc
//    va vb vc
//    wa wb wc
//   Finally, atom positions and atomic symbol are given, expressed in the same coordinate system
//   as the lattice vectors:
//    atom1a atom1b atom1c atom1type
//    atom2a atom2b atom2c atom2type
//    etc.
// OUT:
//  pcell: Each column of this matrix contains a lattice vector in cartesian;
//  patom_pos: cartesian positions of each atom;
//  patom_type: an array of indices into atom_label;
//  paxes: (optional) the coordinate system specified at the beginning of the file
//   (one lattice vector per column);
//  it return 1 upon successeful completion and 0 otherwise;
//
// An example (continued from above):
//  Pure CaO in the MgO-CaO system:
//   3.2 3.2 3.2 90 90 90
//   0   0.5 0.5
//   0.5 0   0.5
//   0.5 0.5 0
//   0 0 0 Ca
//   0.5 0.5 0.5 O
//  then the variables are initialized as follows:
//  *pcell:
//    0   1.6 1.6
//    1.6 0   1.6
//    1.6 1.6 0
//  *patom_pos:
//    0   0   0
//    1.6 1.6 1.6
//  *patom_type:
//    0
//  (assuming the atom_label contains):
//   Ca Mg O
// End of example;
int parse_structure_file(rMatrix3d *pcell, Array<rVector3d> *patom_pos, Array<int> *patom_type, const Array<AutoString> &atom_label, istream &file, rMatrix3d *paxes=NULL);

// takes a structure as read by parse_structure_file and makes it conform;
// to the same standard as parse_lattice_file.;
// It converts a pstr->atom_type that contains indices into atom_label;
// into a  pstr->atom_type that contains indices into site_type_list.;
// ARG:
//  pstr: (IN/OUT) structure to fix (as read by parse_structure_file);
//  lat: (IN) a lattice (as read by parse_lattice_file);
//  site_type_list: (IN) (as read by parse_lattice_file);
//  it return 1 upon successeful completion and 0 otherwise;
int fix_atom_type(Structure *pstr, const Structure &lat, const Array<Arrayint> &site_type_list, int drop_fixed_site);

int fix_atom_type(Structure *pstr, const Array<Arrayint> &site_type_list);

void write_axes(const rMatrix3d &axes, ostream &file, int doabc);

void write_structure(const Structure &str, const Structure &lat,
		     const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc=0);
void write_structure(const Structure &str, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc=0);

void skip_to_next_structure(istream &strfile);

void write_lattice(const Structure &lat, const Array<Arrayint> &site_type_list, const Array<AutoString> &atom_label, const rMatrix3d &axes, ostream &file, int doabc=0);

Real read_clusters_and_eci(LinkedList<Cluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes);
Real read_clusters_and_eci(LinkedList<MultiCluster> *clusterlist, LinkedList<Real> *ecilist,
                           istream &clusterfile, istream &ecifile, const rMatrix3d &axes);

int file_exists(const char *filename);

#endif
