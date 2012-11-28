#include "mmclib.h"
#include "getvalue.h"
#include "parse.h"
#include "lstsqr.h"

typedef int *pint;
typedef int **ppint;
typedef Real *pReal;
typedef Real **ppReal;
typedef Real ***pppReal;
typedef Complex *pComplex;

MultiMonteCarlo::MultiMonteCarlo(const Structure &_lattice, const Array<Array<int> > &_site_type_list, const iVector3d &_supercell, 
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list, 
		const Array<Array<Array<Real> > > &_corrfunc) :
               lattice(_lattice), site_type_list(_site_type_list), supercell(_supercell), corrfunc(_corrfunc), cur_rho(), cur_conc(), mu(), allowed_flip_site(),allowed_flip_before(),allowed_flip_after(), flip_span(-1,-1,-1) {
    Real max_clus_len=0.;
    {
      LinkedListIterator<MultiCluster> c(cluster_list);
      for ( ; c; c++) {
        Real l=get_cluster_length(c->clus);
        if (l>max_clus_len) max_clus_len=l;
      }
    }
    margin=find_sphere_bounding_box(lattice.cell,max_clus_len);
    iVector3d supercell_orig=supercell;
    for (int mult=1; 1 ; mult++) {
      supercell=mult*supercell_orig;
      int i;
      for (i=0; i<3; i++) {
        if (supercell(i) < margin(i)) break;
      }
      if (i==3) break;
    }
    total_box=supercell+2*margin;
    int spin_total_size=total_box(0)*total_box(1)*total_box(2)*lattice.atom_pos.get_size();
    spin=new SPIN_TYPE[spin_total_size];
    for (int i=0; i<spin_total_size; i++) {
      spin[i]=0;
    }
    cur_energy=0;
    int spin_size=supercell(0)*supercell(1)*supercell(2)*lattice.atom_pos.get_size();
    spin_orig=new SPIN_TYPE[spin_size];
    for (int i=0; i<spin_size; i++) {
      spin_orig[i]=0;
    }
    cur_disorder_param=0.;

    rspin_size=(Real)spin_size;

    site_in_cell=lattice.atom_pos.get_size();
    {
      Array<int> frozen_site(site_in_cell);
      zero_array(&frozen_site);
      active_site_in_cell=site_in_cell;
      if (file_exists("frozensite.in")) {
	ifstream file("frozensite.in");
	while (skip_delim(file)) {
	  int s=-1;
	  file >> s;
	  if (!frozen_site(s)) {active_site_in_cell--;}
	  frozen_site(s)=1;
	}
      }
      which_site=new int[active_site_in_cell];
      int s=0;
      for (int i=0; i<active_site_in_cell; i++,s++) {
	while (frozen_site(s)==1) {s++;}
	which_site[i]=s;
      }
    }
    nb_spin_val=new int[site_in_cell];
    for (int i=0; i<site_in_cell; i++) {
	nb_spin_val[i]=site_type_list(lattice.atom_type(i)).get_size();
    }
    nb_clusters=new int[site_in_cell];
    cluster_size=new pint[site_in_cell];
    site_offset=new ppint[site_in_cell];
    spin_val_clus=new pppReal[site_in_cell];
    eci=new pReal[site_in_cell];
    which_cluster=new pint[site_in_cell];
    total_clusters=cluster_list.get_size();
    cur_rho.resize(total_clusters);
    zero_array(&cur_rho);
    pcur_rho=cur_rho.get_buf();
    new_rho=new Real[total_clusters];
    rcluster_mult_per_atom=new Real[total_clusters];
    {
      LinkedListIterator<MultiCluster> c(cluster_list);
      for (int i=0; c; c++,i++) {
	rcluster_mult_per_atom[i]=calc_multiplicity(*c,lattice.cell,space_group.point_op,space_group.trans)/(Real)lattice.atom_pos.get_size();
      }
    }

    E_ref=0.;
    which_is_empty=-1;
    nb_point=0;
    {
      LinkedListIterator<MultiCluster> c(cluster_list);
      int cluster_index=0;
      for (; c; c++, cluster_index++) {
	  if (c->clus.get_size()==0) {which_is_empty=cluster_index;}
	  if (c->clus.get_size()==1) {nb_point++;}
	}
      which_is_point=new int[nb_point];
      point_mult=new Real[nb_point];
	c.init(cluster_list);
      cluster_index=0;
      int point_index=0;
      for (; c; c++, cluster_index++) {
            if (c->clus.get_size()==1) {
              which_is_point[point_index]=cluster_index;
              point_mult[point_index]=rcluster_mult_per_atom[cluster_index];
              point_index++;
            }
	}
    }
    cur_conc.resize(nb_point);
    mu.resize(nb_point);
    zero_array(&mu);
    pmu=mu.get_buf();

    rMatrix3d inv_cell=!lattice.cell;
    for (int s=0; s<lattice.atom_pos.get_size(); s++) {
      LinkedList<MultiCluster> long_cluster_list;
      LinkedList<int> long_cluster_i_list;
      LinkedListIterator<MultiCluster> c(cluster_list);
      int cluster_index=0;
      
      for ( ; c; c++,cluster_index++) {
	if (c->clus.get_size()!=0) {
            Array<MultiCluster> clusters;
            find_equivalent_clusters(&clusters, *c,space_group.cell,space_group.point_op,space_group.trans);
            for (int ec=0; ec<clusters.get_size(); ec++) {
	      for (int center_atom=0; center_atom<clusters(ec).clus.get_size(); center_atom++) {
		rVector3d lat_shift=lattice.atom_pos(s)-clusters(ec).clus(center_atom);
		if (is_int(inv_cell*lat_shift)) {
                  MultiCluster *to_add=new MultiCluster(clusters(ec).clus.get_size());
                  to_add->clus(0)=clusters(ec).clus(center_atom)+lat_shift;
		  to_add->site_type(0)=clusters(ec).site_type(center_atom);
		  to_add->func(0)=clusters(ec).func(center_atom);
                  int j=1;
                  for (int i=0; i<clusters(ec).clus.get_size(); i++) {
                    if (i!=center_atom) {
                      to_add->clus(j)=clusters(ec).clus(i)+lat_shift;
                      to_add->site_type(j)=clusters(ec).site_type(i);
                      to_add->func(j)=clusters(ec).func(i);
                      j++;
                    }
                  }
                  long_cluster_list << to_add;
                  long_cluster_i_list << new int(cluster_index);
                }
              }
            }
        }
      }
      nb_clusters[s]=long_cluster_list.get_size();
      cluster_size[s]=new int[nb_clusters[s]];
      eci[s]=new Real[nb_clusters[s]];
      which_cluster[s]=new int[nb_clusters[s]];
      site_offset[s]=new pint[nb_clusters[s]];
      spin_val_clus[s]=new ppReal[nb_clusters[s]];
      c.init(long_cluster_list);
      LinkedListIterator<int> c_i(long_cluster_i_list);
      for (int ic=0; ic<nb_clusters[s]; ic++, c++, c_i++) {
        cluster_size[s][ic]=c->clus.get_size();
	eci[s][ic]=0.;
        which_cluster[s][ic]=*c_i;
        site_offset[s][ic]=new int[c->clus.get_size()];
        spin_val_clus[s][ic]=new pReal[c->clus.get_size()];
        for (int j=0; j<c->clus.get_size(); j++) {
          int offset_in_cell=which_atom(lattice.atom_pos,c->clus(j),inv_cell);
          iVector3d offset_cell=to_int(inv_cell*(c->clus(j)-lattice.atom_pos(offset_in_cell)));
          site_offset[s][ic][j]=((offset_cell(0)*total_box(1) + offset_cell(1))*total_box(2) + offset_cell(2))*site_in_cell + offset_in_cell - s;
	  spin_val_clus[s][ic][j]=new Real[nb_spin_val[offset_in_cell]];
	  for (int l=0; l<nb_spin_val[offset_in_cell]; l++) {
            spin_val_clus[s][ic][j][l]=_corrfunc(c->site_type(j))(c->func(j))(l);
          }
        }
      }
    }
  }

MultiMonteCarlo::~MultiMonteCarlo(void) {
    for (int s=0; s<site_in_cell; s++) {
      for (int ic=0; ic<nb_clusters[s]; ic++) {
        for (int j=0; j<cluster_size[s][ic]; j++) {
	    delete[] spin_val_clus[s][ic][j];
        }
        delete[] spin_val_clus[s][ic];
        delete[] site_offset[s][ic];
      }
      delete[] spin_val_clus[s];
      delete[] site_offset[s];
      delete[] cluster_size[s];
      delete[] eci[s];
      delete[] which_cluster[s];
    }
    delete[] which_site;
    delete[] spin_val_clus;
    delete[] site_offset;
    delete[] cluster_size;
    delete[] eci;
    delete[] which_cluster;
    delete[] nb_clusters;

    delete[] nb_spin_val;
    delete[] new_rho;
    delete[] rcluster_mult_per_atom;
    delete[] which_is_point;
    delete[] point_mult;

    delete[] spin;
    delete[] spin_orig;
}

void MultiMonteCarlo::calc_delta_point_corr(Array<Real> *pcorr, int site, int type) {
  Array<Real> allcorr(total_clusters);
  zero_array(&allcorr);
  for (int c=0; c<nb_clusters[site]; c++) {
    if (cluster_size[site][c]==1) {
      allcorr(which_cluster[site][c])=spin_val_clus[site][c][0][type];
    }
  }
  for (int p=0; p<nb_point; p++) {
    (*pcorr)(p)+=allcorr(which_is_point[p]); // point_mult[p];
  }
}

void MultiMonteCarlo::calc_delta_point_corr(Array<Real> *pcorr, const Array<int> &sites,const Array<int> &types) {
  pcorr->resize(nb_point);
  zero_array(pcorr);
  for (int i=0; i<sites.get_size(); i++) {
    calc_delta_point_corr(pcorr,sites(i),types(i));
  }
}

void MultiMonteCarlo::find_all_allowed_flips(const Array2d<Real> &corr_constraints, const iVector3d &_flip_span) {
  flip_span=_flip_span;
  LinkedList<Array<int> > site_list;
  LinkedList<Array<int> > before_flip_list;
  LinkedList<Array<Array<int> > > after_flip_list;
  Array<int> site_done(site_in_cell);
  zero_array(&site_done);
  int nb_site_done=0;
  int nb_in_basis=0;
  Array<Array<Real> > corr_basis;
  int nbflip=1;
  while (1) {
    Array<int> maxsite(nbflip);
    for (int i=0; i<nbflip; i++) {maxsite(i)=site_in_cell;}
    MultiDimIterator<Array<int> > sites(maxsite);
    for ( ; sites; sites++) {
      Array<int> maxtype(maxsite.get_size());
      for (int i=0; i<maxsite.get_size(); i++) {maxtype(i)=nb_spin_val[((Array<int> &)sites)(i)];}
      MultiDimIterator<Array<int> > before_types(maxtype);
      int nb_tries=1;
      for (int i=0; i<maxtype.get_size(); i++) {nb_tries*=maxtype(i);}
      Array<int> before_flip(nb_tries);
      for (int i=0; i<before_flip.get_size(); i++) {before_flip(i)=-1;}
      int nb_flip=0;
      LinkedList<Array<int> > after_flip;
      for (; before_types; before_types++) {
        Array<Real> before_corr;
        calc_delta_point_corr(&before_corr, (const Array<int> &)sites,(const Array<int> &)before_types);
        MultiDimIterator<Array<int> > after_types(maxtype);
        for (; after_types; after_types++) {
          int same;
          for (same=0; same<nbflip; same++) {
            if (((Array<int> &)before_types)(same)==((Array<int> &)after_types)(same)) break;
	  }
	  if (same==nbflip) {
	    Array<Real> after_corr;
	    calc_delta_point_corr(&after_corr, (const Array<int> &)sites,(const Array<int> &)after_types);
	    Array<Real> delta_corr;
	    diff(&delta_corr,after_corr,before_corr);
	    Array<Real> c;
	    if (corr_constraints.get_size()(1)!=delta_corr.get_size()) {
	      ERRORQUIT("Incorrect number of point correlations to use the hybrid grandcanonical/canonical mode.\nPlease make sure clusters.out includes all point correlations");
	    }
	    product(&c,corr_constraints,delta_corr);
	    if (inner_product(c,c)<sqr(zero_tolerance)) {
	      int index=0;
	      for (int i=0; i<maxtype.get_size(); i++) {index=index*maxtype(i)+((Array<int> &)before_types)(i);}
	      before_flip(index)=nb_flip;
	      after_flip << new Array<int>(after_types);
	      nb_flip++;
	      for (int l=0; l<((Array<int> &)sites).get_size();  l++) {
	       cerr << ((Array<int> &)sites)(l) << ": " << ((Array<int> &)before_types)(l) << "->" << ((Array<int> &)after_types)(l) << endl;
	      }
	      cerr << endl;
	      build_basis(&corr_basis,&nb_in_basis,delta_corr);
	    }
	  }
	}
      }
      if (nb_flip>0) {
	for (int i=0; i<((Array<int> &)sites).get_size(); i++) {
	  if (site_done(((Array<int> &)sites)(i))==0) {
	    site_done(((Array<int> &)sites)(i))=1;
	    nb_site_done++;
	  }
	}
	site_list << new Array<int>(sites);
	before_flip_list << new Array<int>(before_flip);
	Array<Array<int> > *pafter_flip_array=new Array<Array<int> >();
	LinkedList_to_Array(pafter_flip_array,after_flip);
	after_flip_list << pafter_flip_array;
      }
    }
    if (nb_in_basis >= nb_point-corr_constraints.get_size()(0) && (nb_site_done==site_in_cell)) break;
    nbflip++;
  }
  LinkedList_to_Array(&allowed_flip_site,site_list);
  LinkedList_to_Array(&allowed_flip_before,before_flip_list);
  LinkedList_to_Array(&allowed_flip_after,after_flip_list);
}

int MultiMonteCarlo::force_spin_flip(int *cell, int incell, int newspin) {
  int i;
  int mcell[3],mcellscan[3],offset;
  int oldspin;
  Real denergy,dconc;
  int cluster_count,site_count;
  int **pcluster;
  Real ***pppspin_val_clus,**ppspin_val_clus;
  int *pwhich_cluster;
  int *poffset,*psize;
  Real *peci, *prho, *prho2, *prcluster_mult_per_atom;
  Real rho;
  for (i=0; i<3; i++) {
    mcell[i]=cell[i]+margin(i);
  }

  offset=((mcell[0]*total_box(1) + mcell[1])*total_box(2) + mcell[2])*site_in_cell + incell;
  oldspin=spin[offset];
  denergy=0.;
  for (i=total_clusters, prho=new_rho; i>0; i--, prho++) {
    *prho=0.;
  }
  for (cluster_count=nb_clusters[incell], pcluster=site_offset[incell], pppspin_val_clus=spin_val_clus[incell], psize=cluster_size[incell], pwhich_cluster=which_cluster[incell], peci=eci[incell]; cluster_count>0; cluster_count--, pcluster++, pppspin_val_clus++, psize++, pwhich_cluster++, peci++) {
    poffset=*pcluster;
    ppspin_val_clus=*pppspin_val_clus;
    rho=(*ppspin_val_clus)[newspin]-(*ppspin_val_clus)[oldspin];
    poffset++;
    ppspin_val_clus++;
    for (site_count=(*psize)-1; site_count>0; site_count--, poffset++, ppspin_val_clus++) {
      rho*=(*ppspin_val_clus)[spin[offset+(*poffset)]];
    }
    new_rho[*pwhich_cluster]+=rho;
    denergy+=(*peci)*rho;
  }
  for (i=0; i<nb_point; i++) {
    denergy-=pmu[i]*new_rho[which_is_point[i]];
  }

  cur_energy+=denergy/rspin_size;
  for (i=0, prho=new_rho, prho2=pcur_rho, prcluster_mult_per_atom=rcluster_mult_per_atom; i<total_clusters; i++, prho++, prho2++, prcluster_mult_per_atom++) {
    *prho2+=(*prho)/(*prcluster_mult_per_atom)/rspin_size;
  }

  for (i=0; i<3; i++) {
    if (mcell[i]>=supercell(i)) {mcell[i]-=supercell(i);}
  }
  for (mcellscan[0]=mcell[0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
    for (mcellscan[1]=mcell[1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
      for (mcellscan[2]=mcell[2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
	spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell]=newspin;
      }
    }
  }
  offset=((cell[0]*supercell(1) + cell[1])*supercell(2) + cell[2])*site_in_cell + incell;
  cur_disorder_param+=((newspin!=spin_orig[offset])-(oldspin!=spin_orig[offset]))/rspin_size;
  extension_update_spin_flip(cell,incell,oldspin,newspin);
  return oldspin;
}

void MultiMonteCarlo::save_state(MultiMonteCarloState *pstate) {
  pstate->cur_energy=cur_energy;
  pstate->cur_disorder_param=cur_disorder_param;
  pstate->cur_rho.resize(total_clusters);
  for (int i=0; i<total_clusters; i++) {pstate->cur_rho(i)=pcur_rho[i];}
  pstate->saved_spins.delete_all();
}

void MultiMonteCarlo::save_spin(MultiMonteCarloState *pstate, int *cell, int incell) {
  FlipInfo *pinfo=new FlipInfo();
  for (int i=0; i<3; i++) {
    pinfo->cell[i]=cell[i]+margin(i);
  }
  pinfo->incell=incell;
  int offset=((pinfo->cell[0]*total_box(1) + pinfo->cell[1])*total_box(2) + pinfo->cell[2])*site_in_cell + incell;
  pinfo->spin=spin[offset];
  pstate->saved_spins << pinfo;
}

void MultiMonteCarlo::restore_state(const MultiMonteCarloState &state) {
  cur_energy=state.cur_energy;
  cur_disorder_param=state.cur_disorder_param;
  for (int i=0; i<total_clusters; i++) {pcur_rho[i]=state.cur_rho(i);}
  LinkedListIterator<FlipInfo> it(state.saved_spins);
  for ( ; it; it++) {
    int mcell[3];
    for (int i=0; i<3; i++) {
      mcell[i]=it->cell[i];
      if (mcell[i]>=supercell(i)) {mcell[i]-=supercell(i);}
    }
    int mcellscan[3];
    for (mcellscan[0]=mcell[0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
      for (mcellscan[1]=mcell[1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
	for (mcellscan[2]=mcell[2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
	  spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + it->incell]=it->spin;
	}
      }
    }
  }
}

int MultiMonteCarlo::access(int *cell, int incell) {
  int mcell[3];
  for (int i=0; i<3; i++) {
    mcell[i]=cell[i]+margin(i);
  }
  int offset=((mcell[0]*total_box(1) + mcell[1])*total_box(2) + mcell[2])*site_in_cell + incell;
  return spin[offset];
}

void MultiMonteCarlo::multi_spin_flip(void) {
  int which_flip,before_index,after_index;
  Array<iVector3d> cells;
  Array<int> oldspin;
  do {
    which_flip=random(allowed_flip_site.get_size());
    cells.resize(allowed_flip_site(which_flip).get_size());
    iVector3d center;
    if (flip_span(0)>=0) {
      for (int j=0; j<3; j++) {
	center(j)=random(supercell(j));
      }
    }
    for (int i=0; i<cells.get_size(); i++) {
      int ii;
      do {
	for (int j=0; j<3; j++) {
	  if (flip_span(0)>=0) {
	    cells(i)(j)=(supercell(j)+center(j)+random(2*flip_span(j)+1)-flip_span(j)) % supercell(j);
	  }
	  else {
	    cells(i)(j)=random(supercell(j));
	  }
	}
	for (ii=0; ii<i; ii++) {
	  if (allowed_flip_site(which_flip)(ii)==allowed_flip_site(which_flip)(i) && cells(ii)==cells(i)) break;
	}
      } while (ii!=i);
    }
    before_index=0;
    oldspin.resize(cells.get_size());
    for (int i=0; i<cells.get_size(); i++) {
      before_index=before_index*nb_spin_val[allowed_flip_site(which_flip)(i)]+access(cells(i).get_buf(),allowed_flip_site(which_flip)(i));
      oldspin(i)=access(cells(i).get_buf(),allowed_flip_site(which_flip)(i));
    }
    after_index=allowed_flip_before(which_flip)(before_index);
  } while (after_index==-1);
  MultiMonteCarloState saved_state;
  save_state(&saved_state);
  Real before_energy=cur_energy;
  Real before_extension_energy=extension_get_energy();
  extension_save_state();
  //  cerr << "---<" << endl;
  //  cerr << get_cur_energy() << " " << (-0.5*get_cur_concentration()(0)+1.0*get_cur_concentration()(1)) << endl;
  //  cerr << which_flip << endl << before_index << endl << after_index << endl;
  for (int i=0; i<cells.get_size(); i++) {
    save_spin(&saved_state,cells(i).get_buf(),allowed_flip_site(which_flip)(i));
    force_spin_flip(cells(i).get_buf(),allowed_flip_site(which_flip)(i),allowed_flip_after(which_flip)(after_index)(i));
    //extension_update_spin_flip(cells(i).get_buf(),allowed_flip_site(which_flip)(i),oldspin(i),allowed_flip_after(which_flip)(after_index)(i));

    // cerr << i << ") " << allowed_flip_site(which_flip)(i) << allowed_flip_after(which_flip)(after_index)(i) << endl;
  }
  Real denergy=(cur_energy-before_energy+extension_get_energy()-before_extension_energy)*rspin_size;
  // cerr << get_cur_energy() << " " << (-0.5*get_cur_concentration()(0)+1.0*get_cur_concentration()(1)) << endl;
  int accept=0;
  if (denergy<0) {
    accept=1;
  }
  else {
    if (uniform01()<exp(-denergy/T)) accept=1;
  }
  if (!accept) {
    restore_state(saved_state);
    extension_undo_spin_flip();
  }
  else {
    extension_forget_state();
  }
  // cerr << (accept ? 'A' : 'R') << endl;
  // cerr << get_cur_energy() << " " << (-0.5*get_cur_concentration()(0)+1.0*get_cur_concentration()(1)) << endl;
  // cerr << ">---" << endl;
}

void MultiMonteCarlo::set_eci(const Array<Real> &new_eci) {
  if (which_is_empty>=0) {
    E_ref=new_eci(which_is_empty)/lattice.atom_pos.get_size();
  } else {
    E_ref=0.;
  }
  for (int s=0; s<site_in_cell; s++) {
    for (int c=0; c<nb_clusters[s]; c++) {
      eci[s][c]=new_eci(which_cluster[s][c]);
    }
  }
  calc_from_scratch();
}

void MultiMonteCarlo::init_random(const Array<Array<Real> > &conc) {
    Array<Array<Real> > sconc;
    sconc=conc;
    for (int i=0; i<sconc.get_size(); i++) {
      Real s=0.;
      for (int j=0; j<sconc(i).get_size(); j++) {
         s+=conc(i)(j);
         sconc(i)(j)=s;
      }
      for (int j=0; j<sconc(i).get_size(); j++) {
         sconc(i)(j)/=s;
      }     
    }
    iMatrix3d msupercell;
    msupercell.diag(supercell);
    BoundingBox<int,3> bb(iVector3d(0,0,0),total_box-iVector3d(1,1,1));
    MultiDimIterator<iVector3d> cell(supercell);
    for ( ; cell; cell++) {
      for (int offset_in_cell=0; offset_in_cell<site_in_cell; offset_in_cell++) {
        Real r=uniform01();
        int the_spin=0;
        for (; the_spin<sconc(offset_in_cell).get_size()-1; the_spin++) {
          if (r<sconc(offset_in_cell)(the_spin)) break;
        }
        MultiDimIterator<iVector3d> image(iVector3d(-1,-1,-1),iVector3d(1,1,1));
        for ( ; image; image++) {
          iVector3d image_cell=(iVector3d &)cell+msupercell*(iVector3d &)image+margin;
          if (bb.is_in(image_cell)) {
            int offset=((image_cell(0)*total_box(1) + image_cell(1))*total_box(2) + image_cell(2))*site_in_cell + offset_in_cell;
            spin[offset]=the_spin;
          }
        }
      }
    }
    int spin_size=supercell(0)*supercell(1)*supercell(2)*site_in_cell;
    for (int i=0; i<spin_size; i++) {
      spin_orig[i]=0;
    }
    calc_from_scratch();
    extension_calc_from_scratch();
  }

void MultiMonteCarlo::init_structure(const Structure &str) {
  cerr << "begin init" << endl;
  rMatrix3d inv_cell=!str.cell;
  iMatrix3d msupercell;
  msupercell.diag(supercell);
  BoundingBox<int,3> bb(iVector3d(0,0,0),total_box-iVector3d(1,1,1));
  MultiDimIterator<iVector3d> cell(supercell);
  int quickload=0;
  if (file_exists("quickload")) {quickload=1;}
  cerr << "quickload=" << quickload << endl;;
  int at=-1;
  for (; cell; cell++) {
    rVector3d rcell=lattice.cell*to_real((iVector3d &)cell);
    for (int s=0; s<site_in_cell; s++) {
      if (quickload) {
	at++;
      }
      else {
	at=which_atom(str.atom_pos,rcell+lattice.atom_pos(s),inv_cell);
      }
      SPIN_TYPE newspin=(SPIN_TYPE)( str.atom_type(at)==-1 ? random(nb_spin_val[s]) : str.atom_type(at));
      MultiDimIterator<iVector3d> image(iVector3d(-1,-1,-1),iVector3d(1,1,1));
      for ( ; image; image++) {
        iVector3d image_cell=cell+msupercell*(iVector3d &)image+margin;
        if (bb.is_in(image_cell)) {
          int offset=((image_cell(0)*total_box(1) + image_cell(1))*total_box(2) + image_cell(2))*site_in_cell + s;
          spin[offset]=newspin;
        }
      }
      iVector3d thecell=(iVector3d &)cell;
      int offset=((thecell(0)*supercell(1) + thecell(1))*supercell(2) + thecell(2))*site_in_cell + s;
      spin_orig[offset]=newspin;
    }
  }
  cerr << "end init" << endl;
  calc_from_scratch();
  extension_calc_from_scratch();
  cerr << "end calc" << endl;
}

void MultiMonteCarlo::calc_from_scratch(void) {
    cur_disorder_param=0.;
    cur_energy=0;
    for (int i=0; i<total_clusters; i++) {cur_rho(i)=0.;}
    MultiDimIterator<iVector3d> cur_cell(supercell);
    for ( ; cur_cell; cur_cell++) {
      iVector3d m_cur_cell=(iVector3d &)cur_cell+margin;
      for (int s=0; s<site_in_cell; s++) {
        int moffset=((m_cur_cell(0)*total_box(1) + m_cur_cell(1))*total_box(2) + m_cur_cell(2))*site_in_cell + s;
        iVector3d cell=cur_cell;
        int offset=((cell(0)*supercell(1) + cell(1))*supercell(2) + cell(2))*site_in_cell + s;
        cur_disorder_param+=(spin[moffset]!=spin_orig[offset]);
        for (int c=0; c<nb_clusters[s]; c++) {
          Real rho=1;
          for (int i=0; i<cluster_size[s][c]; i++) {
            rho*=spin_val_clus[s][c][i][spin[moffset+site_offset[s][c][i]]];
          }
	  rho/=(Real)(cluster_size[s][c]);
          cur_energy+=rho*eci[s][c];
	  pcur_rho[which_cluster[s][c]]+=rho/rcluster_mult_per_atom[which_cluster[s][c]];
        }
      }
    }
    cur_disorder_param/=rspin_size;
    cur_energy/=rspin_size;
    product(&cur_rho,cur_rho,1./rspin_size);
    if (which_is_empty>=0) {cur_rho(which_is_empty)=1.;}
    for (int i=0; i<nb_point; i++) {
      cur_energy-=mu(i)*cur_rho(which_is_point[i])*point_mult[i];
    }
  }

const Array<Real> & MultiMonteCarlo::get_cur_concentration(void) {
  for (int i=0; i<nb_point; i++) {
    cur_conc(i)=cur_rho(which_is_point[i])*point_mult[i];
  }
  return cur_conc;
}

void MultiMonteCarlo::set_T_mu(Real _T, const Array<Real> &_mu) {
    T=_T;
    for (int i=0; i<nb_point; i++) {
      cur_energy-=(_mu(i)-mu(i))*cur_rho(which_is_point[i])*point_mult[i];
    }
    mu=_mu;
    pmu=mu.get_buf();
}

void MultiMonteCarlo::spin_flip(void) {
    int i;
    int cell[3],mcell[3],mcellscan[3],incell,inactivecell,offset;
    int oldspin,newspin;
    Real denergy,dconc;
    int cluster_count,site_count;
    int **pcluster;
    Real ***pppspin_val_clus,**ppspin_val_clus;
    int *pwhich_cluster;
    int *poffset,*psize;
    Real *peci, *prho, *prho2, *prcluster_mult_per_atom;
    Real rho;
    int accept;
    for (i=0; i<3; i++) {
      cell[i]=random(supercell(i));
      mcell[i]=cell[i]+margin(i);
    }
    // incell=random(site_in_cell);
    inactivecell=random(active_site_in_cell);
    incell=which_site[inactivecell];
    offset=((mcell[0]*total_box(1) + mcell[1])*total_box(2) + mcell[2])*site_in_cell + incell;
    oldspin=spin[offset];
    newspin=(oldspin+1+random(nb_spin_val[incell]-1)) % nb_spin_val[incell];
    denergy=0.;
    for (i=total_clusters, prho=new_rho; i>0; i--, prho++) {
      *prho=0.;
    }
    for (cluster_count=nb_clusters[incell], pcluster=site_offset[incell], pppspin_val_clus=spin_val_clus[incell], psize=cluster_size[incell], pwhich_cluster=which_cluster[incell], peci=eci[incell]; cluster_count>0; cluster_count--, pcluster++, pppspin_val_clus++, psize++, pwhich_cluster++, peci++) {
      poffset=*pcluster;
      ppspin_val_clus=*pppspin_val_clus;
      rho=(*ppspin_val_clus)[newspin]-(*ppspin_val_clus)[oldspin];
      poffset++;
      ppspin_val_clus++;
      for (site_count=(*psize)-1; site_count>0; site_count--, poffset++, ppspin_val_clus++) {
        rho*=(*ppspin_val_clus)[spin[offset+(*poffset)]];
      }
      new_rho[*pwhich_cluster]+=rho;
      denergy+=(*peci)*rho;
    }
    for (i=0; i<nb_point; i++) {
      denergy-=pmu[i]*new_rho[which_is_point[i]];
    }
    extension_save_state();
    Real save_extension_energy=extension_get_energy();
    extension_update_spin_flip(cell,incell,oldspin,newspin);
    Real dextension_energy=(extension_get_energy()-save_extension_energy)*rspin_size;
    accept=0;
    if (denergy+dextension_energy<0) {
      accept=1;
    }
    else {
      if (uniform01()<exp(-(denergy+dextension_energy)/T)) accept=1;
    }
    if (accept) {
      cur_energy+=denergy/rspin_size;
      for (i=0, prho=new_rho, prho2=pcur_rho, prcluster_mult_per_atom=rcluster_mult_per_atom; i<total_clusters; i++, prho++, prho2++, prcluster_mult_per_atom++) {
         *prho2+=(*prho)/(*prcluster_mult_per_atom)/rspin_size;
      }

      for (i=0; i<3; i++) {
        if (mcell[i]>=supercell(i)) {mcell[i]-=supercell(i);}
      }
      for (mcellscan[0]=mcell[0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
        for (mcellscan[1]=mcell[1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
          for (mcellscan[2]=mcell[2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
            spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell]=newspin;
          }
        }
      }
      offset=((cell[0]*supercell(1) + cell[1])*supercell(2) + cell[2])*site_in_cell + incell;
      cur_disorder_param+=((newspin!=spin_orig[offset])-(oldspin!=spin_orig[offset]))/rspin_size;
    }
    if (accept) {
      extension_forget_state();
    }
    else {
      extension_undo_spin_flip();
    }
}

void MultiMonteCarlo::run(int mc_passes, int mode) {
  int maxn=mc_passes*supercell(0)*supercell(1)*supercell(2)*site_in_cell;
  if (mode==1) {
    for (int n=0; n<maxn; n++) {
      spin_flip();
    }
  }
  else {
    for (int n=0; n<maxn; n++) {
      multi_spin_flip();
    }
  }
}

void MultiMonteCarlo::view(const Array<Arrayint> &labellookup, const Array<AutoString> &atom_label, ofstream &file, const rMatrix3d &axes) {
  for (int i=0; i<3; i++) {
    file << axes.get_column(i) << endl;
  }
  for (int i=0; i<3; i++) {
    file << (!axes)*(supercell(i)*lattice.cell.get_column(i)) << endl;
  }
  rMatrix3d iaxes=!axes;
  MultiDimIterator<iVector3d> cur_cell(supercell);
  for ( ; cur_cell; cur_cell++) {
    iVector3d m_cur_cell=(iVector3d &)cur_cell+margin;
    for (int s=0; s<site_in_cell; s++) {
      int moffset=((m_cur_cell(0)*total_box(1) + m_cur_cell(1))*total_box(2) + m_cur_cell(2))*site_in_cell + s;
      file << iaxes*(lattice.cell*to_real(cur_cell)+lattice.atom_pos(s)) << " " << atom_label(labellookup(lattice.atom_type(s))(spin[moffset])) << endl;
    }
  }
}

void MultiMonteCarlo::get_thermo_data(Array<Real> *pdata) {
  pdata->resize(2+nb_point+total_clusters);
  (*pdata)(0)=get_cur_energy();
  (*pdata)(1)=get_cur_disorder_param();
  int i=2;
  Array<Real> conc;
  conc=get_cur_concentration();
  for (int j=0; j<conc.get_size(); j++) {
    (*pdata)(i)=conc(j);
    i++;
  }
  Array<Real> rho;
  rho=get_cur_corr();
  for (int j=0; j<rho.get_size(); j++) {
    (*pdata)(i)=rho(j);
    i++;
  }
}

void run_mc(GenericAccumulator *accum, MultiMonteCarlo *pmc, int mode) {
  Array<Real> data;
  pmc->get_thermo_data(&data);
  while (accum->new_data(data)) {
    pmc->run(1,mode);
    pmc->get_thermo_data(&data);
  }
}

#include "anyfft.h"

KSpaceMultiMonteCarlo::KSpaceMultiMonteCarlo(const Structure &_lattice, const Array<Array<int> > &_site_type_list,
		const iVector3d &_supercell,
		const SpaceGroup &space_group, const LinkedList<MultiCluster> &cluster_list,
		const Array<Array<Array<Real> > > &_corrfunc, KSpaceECI *_p_kspace_eci):
  MultiMonteCarlo(_lattice,_site_type_list,_supercell,space_group,cluster_list,_corrfunc), flipped_spins(), ft_spin(), ft_eci(), dir_eci(), convol(), ref_x() {
  p_kspace_eci=_p_kspace_eci;
  nsite=_lattice.atom_pos.get_size();
  size=supercell(0)*supercell(1)*supercell(2);
  rsize=(Real)size;

  ft_spin.resize(nsite);
  for (int s=0; s<nsite; s++) {
    ft_spin(s).resize(nb_spin_val[s]);
    for (int ss=0; ss<ft_spin(s).get_size(); ss++) {
      ft_spin(s)(ss).resize(size);
    }
  }
  convol=ft_spin;

  ft_eci.resize(nsite);
  for (int s=0; s<nsite; s++) {
    ft_eci(s).resize(s+1);
    for (int t=0; t<=s; t++) {
      ft_eci(s)(t).resize(nb_spin_val[s]);
      for (int ss=0; ss<ft_eci(s)(t).get_size(); ss++) {
        ft_eci(s)(t)(ss).resize(nb_spin_val[t]);
        for (int tt=0; tt<ft_eci(s)(t)(ss).get_size(); tt++) {
          ft_eci(s)(t)(ss)(tt).resize(size);
        }
      }
    }
  }
  dir_eci=ft_eci;

  fft_time=0;
  flip_time=0;
  ref_x=get_cur_concentration();
  threshold_dx=1e-2;
  cur_recip_E=0.;
  psave_spins=NULL;
}

void plot_3d_surf(const rMatrix3d &k_mesh, const rMatrix3d &rec_lat, const iVector3d &supercell, Array<Complex> a) {
  {
    ofstream file("tmp.gnu");
    file << "set hidden" << endl 
	 << "splot 'tmp.out' u 1:2:3 w l" << endl
	 << "pause -1" << endl;
  }
    {
      ofstream file("tmp.out");
  for (int z=0; z<supercell(2); z++) {
      for (int y=0; y<supercell(1); y++) {
	for (int x=0; x<supercell(0); x++) {
	  int offset=((z*supercell(1) + y)*supercell(0) + x);
	      rVector3d k;
	      if (k_mesh==rec_lat) {
		 k=k_mesh*to_real(iVector3d(x,y,z));
	      } else {
		 k=flip_into_brillouin_1(k_mesh*to_real(iVector3d(x,y,z)),rec_lat);
	      }
	  file << k(0) << " " << k(1) << " " << k(2) << " " << real(a(offset)) << " " << imag(a(offset)) << endl;
	}
	file << endl;
      }
    }
    //    system("gnuplot tmp.gnu");
  }
}

void KSpaceMultiMonteCarlo::set_k_space_eci(void) {
  ref_x=get_cur_concentration();
  for (int s=0; s<nsite; s++) {
    for (int t=0; t<=s; t++) {
      for (int ss=0; ss<ft_eci(s)(t).get_size(); ss++) {
        for (int tt=0; tt<ft_eci(s)(t)(ss).get_size(); tt++) {
	  zero_array(&(ft_eci(s)(t)(ss)(tt)));
	}
      }
    }
  }
  p_kspace_eci->get_k_space_eci(&ft_eci, ref_x);
  dir_eci=ft_eci;
  for (int s=0; s<nsite; s++) {
    for (int t=0; t<=s; t++) {
      for (int ss=0; ss<dir_eci(s)(t).get_size(); ss++) {
        for (int tt=0; tt<dir_eci(s)(t)(ss).get_size(); tt++) {
          fftnd(dir_eci(s)(t)(ss)(tt).get_buf(),3,supercell.get_buf(),-1);
/*
rMatrix3d id;
id.identity();
plot_3d_surf(id,id,supercell,dir_eci(s)(t)(ss)(tt));
cerr << s << " " << t << " " << ss << " " << tt << " " << "eci" << endl;
*/
        }
      }
    }
  }
}

void KSpaceMultiMonteCarlo::extension_calc_from_scratch(void) {
//cerr << "b_total=" << E_ref+cur_energy+cur_recip_E+mu*cur_conc << endl;
  Array<Real> diffx;
  diff(&diffx,get_cur_concentration(),ref_x);
  if (max_norm(diffx)>threshold_dx) {
    set_k_space_eci();
  }
  fft_time=clock();
  flipped_spins.delete_all();
  psave_spins=NULL;
  flip_time=0;
  SPIN_TYPE *spin_corner=spin+((margin(0)*total_box(1) + margin(1))*total_box(2) + margin(2))*nsite;
  MultiDimIterator<iVector3d> cell(supercell);
  for ( ; cell; cell++) {
    iVector3d &rcell=(iVector3d &)cell;
    int offset=((rcell(0)*total_box(1) + rcell(1))*total_box(2) + rcell(2))*nsite;
    //int i=(rcell(0)*supercell(1) + rcell(1))*supercell(2) + rcell(2); //permbug
    int i=(rcell(2)*supercell(1) + rcell(1))*supercell(0) + rcell(0);
    for (int s=0; s<nsite; s++) {
      ft_spin(s)(0)(i)=Complex(1.,0.);
      for (int ss=1; ss<ft_spin(s).get_size(); ss++) {
        ft_spin(s)(ss)(i)=Complex(corrfunc(nb_spin_val[s]-2)(ss-1)(spin_corner[offset+s]),0.);
      }
    }
  }
//Array<Array<Array<Complex> > > savesp(ft_spin);
  for (int s=0; s<nsite; s++) {
    for (int ss=0; ss<ft_spin(s).get_size(); ss++) {
      fftnd(ft_spin(s)(ss).get_buf(),3,supercell.get_buf(),1);
//rMatrix3d id;
//id.identity();
//plot_3d_surf(id,id,supercell,ft_spin(s)(ss));
// cerr << "p" << endl;
    }
  }
  cur_recip_E=0.;
  for (int i=0; i<size; i++) {
    for (int s=0; s<nsite; s++) {
      for (int ss=0; ss<nb_spin_val[s]; ss++) {
        convol(s)(ss)(i)=0.;
        for (int t=0; t<nsite; t++) {
	  for (int tt=0; tt<nb_spin_val[t]; tt++) {
	      if (t<=s) {
		  convol(s)(ss)(i)+=ft_eci(s)(t)(ss)(tt)(i)*ft_spin(t)(tt)(i);
	      }
	      else {
		  convol(s)(ss)(i)+=conj(ft_eci(t)(s)(tt)(ss)(i))*ft_spin(t)(tt)(i);
	      }
          }
        }
      }
    }
  }
  for (int i=0; i<size; i++) {
    for (int s=0; s<nsite; s++) {
      for (int ss=0; ss<ft_spin(s).get_size(); ss++) {
	cur_recip_E+=real(conj(ft_spin(s)(ss)(i))*convol(s)(ss)(i));
//if (fabs(cur_recip_E)>0.) {
//cerr << cur_recip_E << endl;
//}
      }
    }
  }
//cerr << "end" << endl;
  for (int s=0; s<nsite; s++) {
    for (int ss=0; ss<convol(s).get_size(); ss++) {
      fftnd(convol(s)(ss).get_buf(),3,supercell.get_buf(),-1);
    }
  }
/*
  Real debE=0.;
  for (int s=0; s<nsite; s++) {
    for (int ss=0; ss<convol(s).get_size(); ss++) {
	for (int i=0; i<size; i++) {
	    debE+=real(convol(s)(ss)(i)*savesp(s)(ss)(i));
	}
    }
  }
  cerr << "debcE=" << debE/rspin_size << endl;

  {
      Real debE=0.;
      MultiDimIterator<iVector3d> cella(supercell);
      for ( ; cella; cella++) {
	  MultiDimIterator<iVector3d> cellb(supercell);
	  for ( ; cellb; cellb++) {
	      iVector3d &rcella=(iVector3d &)cella;
	      iVector3d &rcellb=(iVector3d &)cellb;
	      iVector3d dr=-(cellb-cella);
	      for (int j=0; j<3; j++) {
		  dr(j)=(dr(j)+supercell(j)) % supercell(j);
	      }
	      int offseta=((rcella(2)*supercell(1) + rcella(1))*supercell(0) + rcella(0));
	      int offsetb=((rcellb(2)*supercell(1) + rcellb(1))*supercell(0) + rcellb(0));
	      int offsetd=((dr(2)*supercell(1) + dr(1))*supercell(0) + dr(0));
	      int offsetp=((((supercell(2)-dr(2))%supercell(2))*supercell(1) + ((supercell(1)-dr(1))%supercell(1)))*supercell(0) + ((supercell(0)-dr(0))%supercell(0)));
	      for (int s=0; s<nsite; s++) {
		  for (int t=0; t<nsite; t++) {
		      for (int ss=0; ss<nb_spin_val[s]; ss++) {
			  for (int tt=0; tt<nb_spin_val[t]; tt++) {
			      if (t<=s) {
				  debE+=real(dir_eci(s)(t)(ss)(tt)(offsetd)*savesp(s)(ss)(offseta)*savesp(t)(tt)(offsetb));
			      } else {
				  debE+=real(dir_eci(t)(s)(tt)(ss)(offsetp)*savesp(s)(ss)(offseta)*savesp(t)(tt)(offsetb));
			      }
			  }
		      }
		  }
	      }
	  }
      }
      cerr << "debE=" << debE/rspin_size << endl;
  }
*/

  cur_recip_E/=(rsize*rspin_size);
cerr << "rE=" << cur_recip_E << endl;
  fft_time=clock()-fft_time;
}

void KSpaceMultiMonteCarlo::extension_save_state(void) {
  LinkedListIterator<FlippedSpin> tmp(flipped_spins);
  psave_spins=(FlippedSpin *)tmp;
  save_recip_E=cur_recip_E;
//cerr << "save" << flipped_spins.get_size() << endl;
}

void KSpaceMultiMonteCarlo::extension_forget_state(void) {
  psave_spins=NULL;
  int nbflip=flipped_spins.get_size();
  Array<Real> diffx;
  diff(&diffx,get_cur_concentration(),ref_x);
//cerr << "forget" << endl;
  if (max_norm(diffx)>threshold_dx) {
    extension_calc_from_scratch();
  } else {
    if (flip_time>0) {
//cerr << "time: " << nbflip << " " << flip_time << " " << 2*fft_time/flip_time << endl;
      if (nbflip > 2*fft_time/flip_time) {
//cerr << "before= " << extension_get_energy() << endl;
	  extension_calc_from_scratch();
//cerr << "after=  " << extension_get_energy() << endl;
      }
    }
  }
}

void KSpaceMultiMonteCarlo::extension_update_spin_flip(int *cell, int incell, int oldspin, int newspin) {

  flip_time=clock();
  int dr[3];
  int mdr[3];
  int *psupercell=supercell.get_buf();

  Array<Real> rspinchange(nb_spin_val[incell]);
  rspinchange(0)=0.;
  for (int ss=1; ss<rspinchange.get_size(); ss++) {
    rspinchange(ss)=corrfunc(rspinchange.get_size()-2)(ss-1)(newspin)-corrfunc(rspinchange.get_size()-2)(ss-1)(oldspin);
  }


  Real dE=0.;
  for (int ss=1; ss<nb_spin_val[incell]; ss++) {
    for (int tt=1; tt<nb_spin_val[incell]; tt++) {
      dE+=real(dir_eci(incell)(incell)(ss)(tt)(0)*rspinchange(ss)*rspinchange(tt));
    }
  }

//dE=0.;

  //int offset=((cell[0]*psupercell[1] + cell[1])*psupercell[2] + cell[2]); //permbug
  int offset=((cell[2]*psupercell[1] + cell[1])*psupercell[0] + cell[0]);

  for (int ss=1; ss<rspinchange.get_size(); ss++) {
    dE+=2.*real(convol(incell)(ss)(offset)*rspinchange(ss));
  }

  LinkedListIterator<FlippedSpin> i(flipped_spins);
  for (; i; i++) {
    for (int j=0; j<3; j++) {
      dr[j] =(cell[j]-(i->cell[j])+psupercell[j]) % psupercell[j];
      mdr[j]=(i->cell[j]-cell[j]+psupercell[j]) % psupercell[j];
    }
    //int offsetdr=((dr[0]*psupercell[1] + dr[1])*psupercell[2] + dr[2]); //permbug
    int offsetdr =((dr[2]*psupercell[1] + dr[1])*psupercell[0] + dr[0]);
    int offsetmdr=((mdr[2]*psupercell[1] + mdr[1])*psupercell[0] + mdr[0]);
    for (int ss=1; ss<nb_spin_val[incell]; ss++) {
      for (int tt=1; tt<nb_spin_val[i->incell]; tt++) {
        if (incell>=i->incell) {
          dE+=2.*real(dir_eci(incell)(i->incell)(ss)(tt)(offsetdr)*rspinchange(ss)*i->dspin(tt));
        } else {
          dE+=2.*real(dir_eci(i->incell)(incell)(tt)(ss)(offsetmdr)*rspinchange(ss)*i->dspin(tt));
	}
      }
    }
  }
  cur_recip_E+=dE/(rspin_size);

/*
 cerr << "Eafterflip=" << cur_recip_E << endl;
 int ck=0;
 if (ck==1) {
     extension_calc_from_scratch();
     cerr << cur_recip_E << endl;
 }
*/
  flipped_spins.push_front(new FlippedSpin(cell,incell,rspinchange));
  flip_time=clock()-flip_time;
}

void KSpaceMultiMonteCarlo::extension_undo_spin_flip(void) {
  LinkedListIterator<FlippedSpin> flip(flipped_spins);
  while ((FlippedSpin *)flip != psave_spins) {
    delete flipped_spins.detach(flip);
  }
  cur_recip_E=save_recip_E;
  psave_spins=NULL;
//cerr << "undo" << flipped_spins.get_size() << endl;
}

KSpaceMultiMonteCarlo::~KSpaceMultiMonteCarlo() {
  // free all dynamically allocated memory;
}
