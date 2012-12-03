#include <fstream.h>
#include <strstream.h>
#include "mrefine.h"
#include "getvalue.h"
#include "mpiinterf.h"

// first do a "find and replace" to change myown to a tag of your choice everywhere. This will be the tag specified on the command line with -fa=tag
// rename this file, for instance, to mrefine_myown.c++ and add mrefine_myown.o in the PLUGINMMAPS section of the makefile

// define a derived class of ClusterExpansion where member functions can be overridded;
class myown_ClusterExpansion: public ClusterExpansion {
 public:
  myown_ClusterExpansion(const Structure &_parent_lattice,
		      const Array<Arrayint> &_site_type_list,
		      const Array<AutoString> &_atom_label,
		      const SpaceGroup &_spacegroup, CorrFuncTable *_pcorrfunc);
  // uncomment if want to override and add function definitio below;
  // void find_best_cluster_choice(CEFitInfo *pfitinfo);  
  // uncomment if want to override and add function definitio below;
  // StructureInfo *find_best_structure(void);
};

// register plug-in name;
SpecificPlugIn<ClusterExpansionCreator, CustomClusterExpansionCreator<myown_ClusterExpansion> > myown_ClusterExpansionPlugIn("myown");

myown_ClusterExpansion::myown_ClusterExpansion(const Structure &_parent_lattice,
  const Array<Arrayint> &_site_type_list,
  const Array<AutoString> &_atom_label,
  const SpaceGroup &_spacegroup, CorrFuncTable *_pcorrfunc):
  ClusterExpansion(_parent_lattice,_site_type_list,_atom_label,_spacegroup,_pcorrfunc) {
  // add initialization code here if needed. Do not delete this function even if not needed;
}

/* uncomment if needed. You may copy in the code from the corresponding function in mrefine.c++ as a startiung point;
void myown_ClusterExpansion::find_best_cluster_choice(CEFitInfo *pfitinfo) {
}
*/

/* uncomment if needed. You may copy in the code from the corresponding function in mrefine.c++ as a startiung point;
StructureInfo * myown_ClusterExpansion::find_best_structure(void) {
}
*/
