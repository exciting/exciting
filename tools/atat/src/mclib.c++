#include "mclib.h"

typedef int *pint;
typedef int **ppint;
typedef Real *pReal;
typedef Complex *pComplex;

MonteCarlo::MonteCarlo(const Structure &_lattice, const iVector3d &_supercell, const SpaceGroup &space_group,
             const LinkedList<Cluster> &cluster_list):
               lattice(_lattice), supercell(_supercell), cur_rho() {
    Real max_clus_len=0.;
    {
      LinkedListIterator<Cluster> c(cluster_list);
      for ( ; c; c++) {
        Real l=get_cluster_length(*c);
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
    spin_changed=new SPIN_TYPE[spin_size];
    for (int i=0; i<spin_size; i++) {
      spin_changed[i]=0;
    }
    cur_disorder_param=0.;
    cur_conc=0.;

    rspin_size=(Real)spin_size;

    site_in_cell=lattice.atom_pos.get_size();
    nb_clusters=new int[site_in_cell];
    nb_x_clusters=new int[site_in_cell];
    cluster_size=new pint[site_in_cell];
    site_offset=new ppint[site_in_cell];
    eci=new pReal[site_in_cell];
    which_cluster=new pint[site_in_cell];
    total_clusters=cluster_list.get_size();
    cur_rho.resize(total_clusters);
    zero_array(&cur_rho);
    pcur_rho=cur_rho.get_buf();
    rcluster_mult_per_atom=new Real[total_clusters];
    {
      LinkedListIterator<Cluster> c(cluster_list);
      for (int i=0; c; c++,i++) {
	rcluster_mult_per_atom[i]=calc_multiplicity(*c,lattice.cell,space_group.point_op,space_group.trans)/(Real)lattice.atom_pos.get_size();
      }
    }

    E_ref=0.;
    which_is_empty=-1;
    rMatrix3d inv_cell=!lattice.cell;
    for (int s=0; s<lattice.atom_pos.get_size(); s++) {
      LinkedList<Cluster> long_cluster_list;
      LinkedList<int> long_cluster_i_list;
      LinkedListIterator<Cluster> c(cluster_list);
      int cluster_index=0;
      
      for ( ; c; c++,cluster_index++) {
	if (c->get_size()==0) {which_is_empty=cluster_index;}
        else {
            Array<Cluster> clusters;
            find_equivalent_clusters(&clusters, *c,space_group.cell,space_group.point_op,space_group.trans);
            for (int ec=0; ec<clusters.get_size(); ec++) {
	      for (int center_atom=0; center_atom<clusters(ec).get_size(); center_atom++) {
		rVector3d lat_shift=lattice.atom_pos(s)-clusters(ec)(center_atom);
		if (is_int(inv_cell*lat_shift)) {
                  Cluster *to_add=new Cluster(clusters(ec).get_size()-1);
                  int j=0;
                  for (int i=0; i<clusters(ec).get_size(); i++) {
                    if (i!=center_atom) {
                      (*to_add)(j)=clusters(ec)(i)+lat_shift;
                      j=j+1;
                    }
                  }
                  long_cluster_list << to_add;
                  long_cluster_i_list << new int(cluster_index);
                }
              }
            }
        }
      }
      nb_x_clusters[s]=long_cluster_list.get_size();
      cluster_size[s]=new int[nb_x_clusters[s]];
      eci[s]=new Real[nb_x_clusters[s]];
      which_cluster[s]=new int[nb_x_clusters[s]];
      site_offset[s]=new pint[nb_x_clusters[s]];
      c.init(long_cluster_list);
      LinkedListIterator<int> c_i(long_cluster_i_list);
      for (int ic=0; ic<nb_x_clusters[s]; ic++, c++, c_i++) {
        cluster_size[s][ic]=c->get_size();
	nb_clusters[s]=ic+1;
	eci[s][ic]=0.;
        which_cluster[s][ic]=*c_i;
        site_offset[s][ic]=new int[c->get_size()];
        for (int j=0; j<c->get_size(); j++) {
          int offset_in_cell=which_atom(lattice.atom_pos,(*c)(j),inv_cell);
          iVector3d offset_cell=to_int(inv_cell*((*c)(j)-lattice.atom_pos(offset_in_cell)));
          site_offset[s][ic][j]=((offset_cell(0)*total_box(1) + offset_cell(1))*total_box(2) + offset_cell(2))*site_in_cell + offset_in_cell - s;
        }
      }
    }
    mu=0.;
  }

MonteCarlo::~MonteCarlo(void) {
  //to write: delete all dynamically allocated arrays;
}

void MonteCarlo::set_eci(const Array<Real> &new_eci) {
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

void MonteCarlo::init_random(Real concentration) {
    Real c=(1+concentration)/2;
    iMatrix3d msupercell;
    msupercell.diag(supercell);
    BoundingBox<int,3> bb(iVector3d(0,0,0),total_box-iVector3d(1,1,1));
    MultiDimIterator<iVector3d> cell(supercell);
    for ( ; cell; cell++) {
      for (int offset_in_cell=0; offset_in_cell<site_in_cell; offset_in_cell++) {
        int the_spin=( uniform01()<c ? +1 : -1 );
        MultiDimIterator<iVector3d> image(iVector3d(-1,-1,-1),iVector3d(1,1,1));
        for ( ; image; image++) {
          iVector3d image_offset_cell=(iVector3d &)cell+msupercell*(iVector3d &)image+margin;
          if (bb.is_in(image_offset_cell)) {
            int offset=((image_offset_cell(0)*total_box(1) + image_offset_cell(1))*total_box(2) + image_offset_cell(2))*site_in_cell + offset_in_cell;
            spin[offset]=the_spin;
          }
        }
      }
    }
    int spin_size=supercell(0)*supercell(1)*supercell(2)*site_in_cell;
    for (int i=0; i<spin_size; i++) {
      spin_changed[i]=0;
    }
    calc_from_scratch();
    extension_calc_from_scratch();
  }

void MonteCarlo::init_structure(const Structure &str) {
  if (str.atom_type(0)==0) {init_random(0);}
  else {
    Structure super;
    iMatrix3d msupercell;
    msupercell.diag(supercell);
    super.cell=lattice.cell*to_real(msupercell);
    find_all_atom_in_supercell(&super.atom_pos,&super.atom_type, str.atom_pos, str.atom_type, str.cell, super.cell);
    rMatrix3d inv_cell=!lattice.cell;
    BoundingBox<int,3> bb(iVector3d(0,0,0),total_box-iVector3d(1,1,1));
    for (int i=0; i<super.atom_pos.get_size(); i++) {
      int offset_in_cell=which_atom(lattice.atom_pos,super.atom_pos(i),inv_cell);
      iVector3d offset_cell=to_int(inv_cell*(super.atom_pos(i)-lattice.atom_pos(offset_in_cell)));
      MultiDimIterator<iVector3d> image(iVector3d(-1,-1,-1),iVector3d(1,1,1));
      for ( ; image; image++) {
        iVector3d image_offset_cell=offset_cell+msupercell*(iVector3d &)image+margin;
        if (bb.is_in(image_offset_cell)) {
          int offset=((image_offset_cell(0)*total_box(1) + image_offset_cell(1))*total_box(2) + image_offset_cell(2))*site_in_cell + offset_in_cell;
          spin[offset]=(SPIN_TYPE)(super.atom_type(i));
        }
      }
    }
    int spin_size=supercell(0)*supercell(1)*supercell(2)*site_in_cell;
    for (int i=0; i<spin_size; i++) {
      spin_changed[i]=0;
    }
    calc_from_scratch();
    extension_calc_from_scratch();
  }
}

void MonteCarlo::set_concentration(Real concentration) {
    int nbspin=(int)(rspin_size*(concentration-get_cur_concentration())/2);
    int toflip=(nbspin>0 ? -1 : 1);
    nbspin=abs(nbspin);
    int cell[3],mcell[3],mcellscan[3],incell,offset;
    
    for (int s=0; s<nbspin; s++) {
      do {
	for (int i=0; i<3; i++) {
	  cell[i]=random(supercell(i));
	  mcell[i]=cell[i]+margin(i);
	}
	incell=random(site_in_cell);
	offset=((mcell[0]*total_box(1) + mcell[1])*total_box(2) + mcell[2])*site_in_cell + incell;
      } while (spin[offset]!=toflip);
      for (int i=0; i<3; i++) {
        if (mcell[i]>=supercell(i)) {mcell[i]-=supercell(i);}
      }
      for (mcellscan[0]=mcell[0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
        for (mcellscan[1]=mcell[1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
          for (mcellscan[2]=mcell[2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
            spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell]*=-1;
          }
        }
      }
      offset=((cell[0]*supercell(1) + cell[1])*supercell(2) + cell[2])*site_in_cell + incell;
      spin_changed[offset]^=1;
    }
    calc_from_scratch();
    extension_calc_from_scratch();
}

void MonteCarlo::calc_from_scratch(void) {
    cur_conc=0.;
    cur_disorder_param=0.;
    cur_energy=0;
    for (int i=0; i<total_clusters; i++) {cur_rho(i)=0.;}
    MultiDimIterator<iVector3d> cur_cell(supercell);
    for ( ; cur_cell; cur_cell++) {
      iVector3d m_cur_cell=(iVector3d &)cur_cell+margin;
      for (int s=0; s<site_in_cell; s++) {
        int moffset=((m_cur_cell(0)*total_box(1) + m_cur_cell(1))*total_box(2) + m_cur_cell(2))*site_in_cell + s;
        cur_conc+=(Real)spin[moffset];
        iVector3d cell=cur_cell;
        int offset=((cell(0)*supercell(1) + cell(1))*supercell(2) + cell(2))*site_in_cell + s;
        cur_disorder_param+=(Real)spin_changed[offset];
        for (int c=0; c<nb_clusters[s]; c++) {
          int rho=spin[moffset];
          for (int i=0; i<cluster_size[s][c]; i++) {
            rho*=spin[moffset+site_offset[s][c][i]];
          }
	  Real rrho=(Real)rho/(Real)(cluster_size[s][c]+1);
          cur_energy+=rrho*eci[s][c];
	  pcur_rho[which_cluster[s][c]]+=rrho/rcluster_mult_per_atom[which_cluster[s][c]];
        }
      }
    }
    cur_conc/=rspin_size;
    cur_disorder_param/=rspin_size;
    cur_energy/=rspin_size;
    product(&cur_rho,cur_rho,1./rspin_size);
    if (which_is_empty>=0) {cur_rho(which_is_empty)=1.;}
    cur_energy-=mu*cur_conc;
  }

void MonteCarlo::init_run(Real _T, Real _mu) {
    T=_T;
    cur_energy-=(_mu-mu)*cur_conc;
    mu=_mu;
}

void MonteCarlo::spin_flip(void) {
    int i;
    int cell[3],mcell[3],mcellscan[3],incell,offset;
    Real denergy,d_recip_energy,dconc;
    int cluster_count,site_count;
    int **pcluster;
    int *pwhich_cluster;
    int *poffset,*psize;
    Real *peci;
    int rho;
    int accept;
    for (i=0; i<3; i++) {
      cell[i]=random(supercell(i));
      mcell[i]=cell[i]+margin(i);
    }
    incell=random(site_in_cell);
    offset=((mcell[0]*total_box(1) + mcell[1])*total_box(2) + mcell[2])*site_in_cell + incell;
    denergy=0.;
    for (cluster_count=nb_clusters[incell], pcluster=site_offset[incell], psize=cluster_size[incell], peci=eci[incell]; cluster_count>0; cluster_count--, pcluster++, psize++, peci++) {
      rho=spin[offset];
      for (site_count=*psize, poffset=*pcluster; site_count>0; site_count--, poffset++) {
        rho*=spin[offset+(*poffset)];
      }
      denergy+=-2.*(*peci)*(Real)rho;
    }
    dconc=-(Real)(2*spin[offset]);
    denergy+=-mu*dconc;
    d_recip_energy=extension_spin_flip_energy(cell,incell,-spin[offset]);
    accept=0;
    if (denergy+d_recip_energy<0) {
      accept=1;
    }
    else {
      if (uniform01()<exp(-(denergy+d_recip_energy)/T)) accept=1;
    }
    if (accept) {
      cur_energy+=denergy/rspin_size;
      cur_conc+=dconc/rspin_size;

      for (cluster_count=nb_clusters[incell], pcluster=site_offset[incell], psize=cluster_size[incell], pwhich_cluster=which_cluster[incell]; cluster_count>0; cluster_count--, pcluster++, psize++, pwhich_cluster++) {
	rho=spin[offset];
	for (site_count=*psize, poffset=*pcluster; site_count>0; site_count--, poffset++) {
	  rho*=spin[offset+(*poffset)];
	}
	pcur_rho[*pwhich_cluster]+=-2.*(Real)rho/rcluster_mult_per_atom[*pwhich_cluster]/rspin_size;
      }

      for (i=0; i<3; i++) {
        if (mcell[i]>=supercell(i)) {mcell[i]-=supercell(i);}
      }
      for (mcellscan[0]=mcell[0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
        for (mcellscan[1]=mcell[1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
          for (mcellscan[2]=mcell[2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
            spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell]*=-1;
          }
        }
      }
      offset=((cell[0]*supercell(1) + cell[1])*supercell(2) + cell[2])*site_in_cell + incell;
      cur_disorder_param+=(Real)(1-2*spin_changed[offset])/rspin_size;
      spin_changed[offset]^=1;
      extension_update_spin_flip(cell,incell,spin[offset],d_recip_energy);
    }
}

void MonteCarlo::spin_double_flip(void) {
    int i,f;
    int array_cell[6],array_mcell[6],array_bmcell[6];
    int *cell[2]={array_cell,array_cell+3};
    int *mcell[2]={array_mcell,array_mcell+3};
    int *bmcell[2]={array_bmcell,array_bmcell+3};
    int mcellscan[3],incell[2],offset[2],nomoffset;
    Real denergy,d_total_energy;
    Real d_recip_energy[2];
    int cluster_count,site_count;
    int **pcluster;
    int *pwhich_cluster;
    int *poffset,*psize;
    Real *peci;
    int rho;
    int accept;

    offset[1]=0;
    for (f=0; f<2; f++) {
      do {
        for (i=0; i<3; i++) {
          cell[f][i]=random(supercell(i));
          mcell[f][i]=cell[f][i]+margin(i);
        }
        incell[f]=random(site_in_cell);
        offset[f]=((mcell[f][0]*total_box(1) + mcell[f][1])*total_box(2) + mcell[f][2])*site_in_cell + incell[f];
      } while (f==1 && spin[offset[0]]==spin[offset[1]]);
      for (i=0; i<3; i++) {
	bmcell[f][i]=mcell[f][i];
	if (bmcell[f][i]>=supercell(i)) {bmcell[f][i]-=supercell(i);}
      }
    }

    denergy=0.;
    for (f=0; f<2; f++) {
      for (cluster_count=nb_clusters[incell[f]], pcluster=site_offset[incell[f]], psize=cluster_size[incell[f]], peci=eci[incell[f]]; cluster_count>0; cluster_count--, pcluster++, psize++, peci++) {
        rho=spin[offset[f]];
        for (site_count=*psize, poffset=*pcluster; site_count>0; site_count--, poffset++) {
          rho*=spin[offset[f]+(*poffset)];
        }
        denergy+=-2.*(*peci)*(Real)rho;
      }
      d_recip_energy[f]=extension_spin_flip_energy(cell[f],incell[f],-spin[offset[f]]);
      if (f==0) {
        for (mcellscan[0]=bmcell[0][0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
          for (mcellscan[1]=bmcell[0][1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
            for (mcellscan[2]=bmcell[0][2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
              spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell[0]]*=-1;
            }
          }
        }
        extension_update_spin_flip(cell[0],incell[0],spin[offset[0]],d_recip_energy[f]);
      }
    }

    d_total_energy=denergy+d_recip_energy[0]+d_recip_energy[1];
    accept=0;
    if (d_total_energy<0) {
      accept=1;
    }
    else {
      if (uniform01()<exp(-d_total_energy/T)) accept=1;
    }
 
    if (accept) {
      for (f=0; f<2; f++) {
        if (f==1) {
          for (mcellscan[0]=bmcell[f][0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
            for (mcellscan[1]=bmcell[f][1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
              for (mcellscan[2]=bmcell[f][2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
                spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell[f]]*=-1;
              }
            }
          }
        }
        for (cluster_count=nb_clusters[incell[f]], pcluster=site_offset[incell[f]], psize=cluster_size[incell[f]], pwhich_cluster=which_cluster[incell[f]]; cluster_count>0; cluster_count--, pcluster++, psize++, pwhich_cluster++) {
          rho=-spin[offset[f]];
          for (site_count=*psize, poffset=*pcluster; site_count>0; site_count--, poffset++) {
	      rho*=spin[offset[f]+(*poffset)];
          }
	    pcur_rho[*pwhich_cluster]+=-2.*(Real)rho/rcluster_mult_per_atom[*pwhich_cluster]/rspin_size;
        }
        nomoffset=((cell[f][0]*supercell(1) + cell[f][1])*supercell(2) + cell[f][2])*site_in_cell + incell[f];
        cur_disorder_param+=(Real)(1-2*spin_changed[nomoffset])/rspin_size;
        spin_changed[nomoffset]^=1;
      }
      cur_energy+=denergy/rspin_size;
      extension_update_spin_flip(cell[1],incell[1],spin[offset[1]],d_recip_energy[1]);
    }
    else {
      for (mcellscan[0]=bmcell[0][0]; mcellscan[0]<total_box(0); mcellscan[0]+=supercell(0)) {
        for (mcellscan[1]=bmcell[0][1]; mcellscan[1]<total_box(1); mcellscan[1]+=supercell(1)) {
          for (mcellscan[2]=bmcell[0][2]; mcellscan[2]<total_box(2); mcellscan[2]+=supercell(2)) {
            spin[((mcellscan[0]*total_box(1) + mcellscan[1])*total_box(2) + mcellscan[2])*site_in_cell + incell[0]]*=-1;
          }
        }
      }
      extension_undo_spin_flip();
    }
}

void MonteCarlo::run(int mc_passes, int mode) {
  int maxn=mc_passes*supercell(0)*supercell(1)*supercell(2)*site_in_cell;
  if (mode==1) {
    for (int n=0; n<maxn; n++) {
      spin_flip();
    }
  }
  else {
    for (int n=0; n<maxn; n++) {
      spin_double_flip();
    }
  }
}

void MonteCarlo::view(const Array<Arrayint> &labellookup, const Array<AutoString> &atom_label, ofstream &file, const rMatrix3d &axes) {
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
      file << iaxes*(lattice.cell*to_real(cur_cell)+lattice.atom_pos(s)) << " " << atom_label(labellookup(lattice.atom_type(s))((1+spin[moffset])/2)) << endl;
    }
  }
}

void run_mc_until_converged(MonteCarlo *pmc, int mode, Equilibrator *eq, Real prec, int which_elem) {
  eq->init(prec,nb_value_accum+pmc->get_total_clusters(),which_elem);
  Array<Real> data(nb_value_accum+pmc->get_total_clusters());
  do {
    pmc->run(1,mode);
    data(0)=pmc->get_cur_energy();
    data(1)=pmc->get_cur_concentration();
    data(2)=1.-pmc->get_cur_disorder_param();
    const Array<Real> &cur_corr=pmc->get_cur_corr();
    for (int i=0; i<cur_corr.get_size(); i++) {data(nb_value_accum+i)=cur_corr(i);}
  } while (eq->new_data(data));
}

void run_mc(MonteCarlo *pmc, int mode, int n_step, Accumulator *eq) {
  Array<Real> data(nb_value_accum+pmc->get_total_clusters());
  if (n_step==0) {
    data(0)=pmc->get_cur_energy();
    data(1)=pmc->get_cur_concentration();
    data(2)=1.-pmc->get_cur_disorder_param();
    const Array<Real> &cur_corr=pmc->get_cur_corr();
    for (int i=0; i<cur_corr.get_size(); i++) {data(nb_value_accum+i)=cur_corr(i);}
    eq->new_data(data);
  }
  else {
    for (int n=0; n<n_step; n++) {
      pmc->run(1,mode);
      data(0)=pmc->get_cur_energy();
      data(1)=pmc->get_cur_concentration();
      data(2)=1.-pmc->get_cur_disorder_param();
      const Array<Real> &cur_corr=pmc->get_cur_corr();
      for (int i=0; i<cur_corr.get_size(); i++) {data(nb_value_accum+i)=cur_corr(i);}
      eq->new_data(data);
    }
  }
}

void run_mc(MCOutputData *pmcdata, MonteCarlo *pmc, int mode, int n_step, Real prec, int which_elem ) {
  pmcdata->T=pmc->get_T();
  pmcdata->mu=pmc->get_mu();
  GenericAccumulator *pgacc;
  Accumulator acc;
  Equilibrator eq;
  if (prec==0.) {
    run_mc(pmc,mode,n_step,&acc);
    pgacc=&acc;
    pmcdata->n_equil=0;
    pmcdata->n_step=n_step;
  }
  else {
    run_mc_until_converged(pmc,mode,&eq,prec,which_elem);
    pgacc=&eq;
    eq.get_step(&(pmcdata->n_equil),&(pmcdata->n_step));
  }
  const Array<Real> &mean=pgacc->get_mean();
  const Array<Real> &var=pgacc->get_var();
  pmcdata->E=mean(0);
  pmcdata->x=mean(1);
  pmcdata->lro=mean(2);
  pmcdata->heatcap=var(0);
  pmcdata->suscept=var(1);
  pmcdata->corr.resize(pmc->get_total_clusters());
  for (int i=0; i<pmcdata->corr.get_size(); i++) {
    (pmcdata->corr)(i)=mean(nb_value_accum+i);
  }
  pmc->get_cluster_mult(&(pmcdata->mult));
}

Equilibrator::Equilibrator(void): cur_sum(), cur_sum2(), bin_sum(), bin_sum2(), good_val(), buf_val() {}

Equilibrator::Equilibrator(Real _prec, int data_size, int _which_elem, int init_granularity, int nb_bin):
  cur_sum(), cur_sum2(), bin_sum(), bin_sum2(), good_val(), buf_val() {
  init(_prec,data_size,_which_elem,init_granularity,nb_bin);
}

void Equilibrator::init(Real _prec, int data_size, int _which_elem, int init_granularity, int nb_bin) {
  prec=_prec;
  which_elem=_which_elem;
  if (nb_bin%2==1) ERRORQUIT("nb_bin must be even");
  granularity=init_granularity;
  cur_corr_len=granularity/2;
  cur_bin=1;
  cur_cnt=0;
  cur_sum.resize(data_size);
  cur_sum2.resize(data_size);
  zero_array(&cur_sum);
  zero_array(&cur_sum2);
  bin_sum.resize(nb_bin+1);
  bin_sum2.resize(nb_bin+1);
  for (int i=0; i<bin_sum.get_size(); i++) {
    bin_sum(i).resize(data_size);
    bin_sum2(i).resize(data_size);
    zero_array(&bin_sum(i));
    zero_array(&bin_sum2(i));
  }
  buf_corr.resize(nb_bin+1);
  buf_val.resize(granularity);
}

int Equilibrator::new_data(const Array<Real> &data) {
  buf_val(cur_cnt)=data(which_elem);
  sum(&cur_sum,cur_sum,data);
  Array<Real> data2;
  product_diag(&data2,data,data);
  sum(&cur_sum2,cur_sum2,data2);
  cur_cnt++;
  if (cur_cnt==granularity) {
    bin_sum(cur_bin)=cur_sum;
    bin_sum2(cur_bin)=cur_sum2;
    
    Real corr;
    Real last_corr=0.5;
    int done=0;
    while (!done) {
      corr=0.;
      for (int i=0; i<granularity-cur_corr_len; i++) {
	corr+=buf_val(i)*buf_val(i+cur_corr_len);
      }
      corr/=(Real)(granularity-cur_corr_len);
      corr-=sqr((bin_sum(cur_bin)(which_elem)-bin_sum(cur_bin-1)(which_elem))/(Real)granularity);
      Real var=(bin_sum2(cur_bin)(which_elem)-bin_sum2(cur_bin-1)(which_elem))/(Real)granularity
	-sqr((bin_sum(cur_bin)(which_elem)-bin_sum(cur_bin-1)(which_elem))/(Real)granularity);
      if (var>0) {
	corr/=var;
      }
      else {
	corr=1.;
      }
      done=1;
      if (corr<0.) corr=0.;
      if (corr<0.25 && last_corr<=0.75 && cur_corr_len>=2) {cur_corr_len/=2; done=0;}
      if (corr>0.75 && last_corr>=0.25 && cur_corr_len<granularity/4) {cur_corr_len*=2; done=0;}
      last_corr=corr;
    }
    corr=pow(corr,1./(Real)cur_corr_len);
    buf_corr(cur_bin)=corr;
    
    int b;
    for (b=cur_bin/2; b>=1; b--) {
      Real block1=(bin_sum(cur_bin  )(which_elem)-bin_sum(cur_bin-  b)(which_elem))/(Real)(granularity*b);
      Real block2=(bin_sum(cur_bin-b)(which_elem)-bin_sum(cur_bin-2*b)(which_elem))/(Real)(granularity*b);
      if (fabs(block1-block2)<prec) {
	int bin1=cur_bin-2*b;
	Real n=(Real)(granularity*2*b);
	Real avg_corr=0.;
	for (int i=bin1+1; i<=cur_bin; i++) {
	  avg_corr+=corr;
	}
	avg_corr/=(Real)(cur_bin-bin1);
	Real var=((bin_sum2(cur_bin)(which_elem)-bin_sum2(bin1)(which_elem))/n-sqr((bin_sum(cur_bin)(which_elem)-bin_sum(bin1)(which_elem))/n))/n;
	if ((1.-avg_corr)!=0.) {
	  var*=(1.+avg_corr)/(1.-avg_corr);
	  if (var<sqr(prec) && var>0) {
	    diff(&good_val,bin_sum(cur_bin),bin_sum(bin1));
	    product(&good_val,good_val,1./n);
	    diff(&good_val2,bin_sum2(cur_bin),bin_sum2(bin1));
	    product(&good_val2,good_val2,1./n);
	    Array<Real> tmp;
	    product_diag(&tmp,good_val,good_val);
	    diff(&good_val2,good_val2,tmp);
	    equil=bin1*granularity;
	    step=(cur_bin-bin1)*granularity;
	    break;
	  }
	}
      }
    }
    
    cur_cnt=0;
    cur_bin++;
    if (cur_bin==bin_sum.get_size()) {
      int i,j;
      for (i=1,j=2; j<bin_sum.get_size(); i++,j+=2) {
	bin_sum(i)=bin_sum(j);
	bin_sum2(i)=bin_sum2(j);
	buf_corr(i)=(buf_corr(j-1)+buf_corr(j))/2.;
      }
      for (; i<bin_sum.get_size(); i++) {
	zero_array(&bin_sum(i));
	zero_array(&bin_sum2(i));
      }
      cur_bin=(bin_sum.get_size()-1)/2+1;
      granularity*=2;
      buf_val.resize(granularity);
    }
  }
  return (good_val.get_size()==0);
}

int Accumulator::new_data(const Array<Real> &data) {
  if (cur_sum.get_size()==0) {
    cur_sum.resize(data.get_size());
    zero_array(&cur_sum);
    cur_sum2.resize(data.get_size());
    zero_array(&cur_sum2);
    cur_n=0;
  }
  sum(&cur_sum,cur_sum,data);
  Array<Real> data2;
  product_diag(&data2,data,data);
  sum(&cur_sum2,cur_sum2,data2);
  cur_n++;
  return 0;
}

void Accumulator::accum(void) {
  product(&good_val,cur_sum,1./(Real)cur_n);
  product(&good_val2,cur_sum2,1./(Real)cur_n);
  Array<Real> tmp;
  product_diag(&tmp,good_val,good_val);
  diff(&good_val2,good_val2,tmp);
}

#include "anyfft.h"

KSpaceMonteCarlo::KSpaceMonteCarlo(const Structure &_lattice, const iVector3d &_supercell, const SpaceGroup &space_group,
				   const LinkedList<Cluster> &cluster_list, KSpaceECI *_p_kspace_eci):
  MonteCarlo(_lattice,_supercell,space_group,cluster_list),flipped_spins() {
  p_kspace_eci=_p_kspace_eci;
  nsite=_lattice.atom_pos.get_size();
  size=supercell(0)*supercell(1)*supercell(2);
  rsize=(Real)size;
  ft_spin=new pComplex[nsite];
  convol=new pComplex[nsite];
  ft_eci=new pComplex[nsite*(nsite+1)/2];
  dir_eci=new pComplex[nsite*(nsite+1)/2];
  for (int s=0; s<nsite; s++) {
    ft_spin[s]=new Complex[size];
    convol[s]=new Complex[size];
  }
  for (int s=0; s<nsite; s++) {
    for (int t=0; t<=s; t++) {
      ft_eci[unrollsym(s,t)]=new Complex[size];
      dir_eci[unrollsym(s,t)]=new Complex[size];
    }
  }
  fft_time=0;
  flip_time=0;
  ref_x=MAXFLOAT;
  threshold_dx=1e-2;
  cur_recip_E=0.;
}

void plot_3d_surf(const rMatrix3d &k_mesh, const rMatrix3d &rec_lat, const iVector3d &supercell, Complex *pa) {
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
	  int offset=((x*supercell(1) + y)*supercell(2) + z);
	  rVector3d k=flip_into_brillouin_1(k_mesh*to_real(iVector3d(x,y,z)),rec_lat);
	  file << k(0) << " " << k(1) << " " << k(2) << " " << real(pa[offset]) << " " << imag(pa[offset]) << endl;
	}
	file << endl;
      }
    }
    //    system("gnuplot tmp.gnu");
  }
}

void KSpaceMonteCarlo::set_k_space_eci(void) {
  ref_x=get_cur_concentration();
  rMatrix3d mat_supercell;
  mat_supercell.diag(to_real(supercell));
  rMatrix3d recip_lat=!(~(lattice.cell));
  rMatrix3d k_mesh=!(~(lattice.cell*mat_supercell));
  Real x=(1.+get_cur_concentration())/2.;
  //cerr << "x=" << x << endl;
  MultiDimIterator<iVector3d> cell(supercell);
  iVector3d &vcell=(iVector3d &)cell;
  Array2d<Complex> cur_ft_eci;
  for (; cell; cell++) {
    //int offset=((vcell(0)*supercell(1) + vcell(1))*supercell(2) + vcell(2)); //permbug
    int offset=((vcell(2)*supercell(1) + vcell(1))*supercell(0) + vcell(0));
    rVector3d k=flip_into_brillouin_1(k_mesh*to_real(cell),recip_lat);
    p_kspace_eci->get_k_space_eci(&cur_ft_eci,k,x);
    for (int s=0; s<nsite; s++) {
      for (int t=0; t<=s; t++) {
	int st=unrollsym(s,t);
	ft_eci[st][offset]=cur_ft_eci(s,t);
	dir_eci[st][offset]=ft_eci[st][offset];
      }
    }
  }
#ifdef DEBUG
  plot_3d_surf(k_mesh,recip_lat,supercell,dir_eci[0]);
#endif
  for (int s=0; s<nsite; s++) {
    for (int t=0; t<=s; t++) {
      fftnd(dir_eci[unrollsym(s,t)],3,supercell.get_buf(),-1);
    }
  }
#ifdef DEBUG
  plot_3d_surf(k_mesh,recip_lat,supercell,dir_eci[0]);
#endif
}

void KSpaceMonteCarlo::extension_calc_from_scratch(void) {
//cerr << "b_total=" << E_ref+cur_energy+cur_recip_E+mu*cur_conc << endl;
  if (fabs(get_cur_concentration()-ref_x)>threshold_dx) {
    set_k_space_eci();
  }
  fft_time=clock();
  flipped_spins.delete_all();
  flip_time=0;
  SPIN_TYPE *spin_corner=spin+((margin(0)*total_box(1) + margin(1))*total_box(2) + margin(2))*nsite;
  MultiDimIterator<iVector3d> cell(supercell);
  for (; cell; cell++) {
    iVector3d &rcell=(iVector3d &)cell;
    int offset=((rcell(0)*total_box(1) + rcell(1))*total_box(2) + rcell(2))*nsite;
    int i=(rcell(2)*supercell(1) + rcell(1))*supercell(0) + rcell(0); //permbugfix
    for (int s=0; s<nsite; s++) {
      ft_spin[s][i]=(Complex)(spin_corner[offset+s]);
    }
  }
  for (int s=0; s<nsite; s++) {
    fftnd(ft_spin[s],3,supercell.get_buf(),1);
  }
  cur_recip_E=0.;
//rMatrix3d id;
//id.identity();
//plot_3d_surf(id,id,supercell,ft_spin[0]);
  for (int i=0; i<size; i++) {
    for (int s=0; s<nsite; s++) {
      convol[s][i]=0.;
      for (int t=0; t<nsite; t++) {
	convol[s][i]+=ft_eci[unrollsym(s,t)][i]*ft_spin[t][i];
      }
    }
  }
  for (int i=0; i<size; i++) {
    for (int s=0; s<nsite; s++) {
      cur_recip_E+=real(conj(ft_spin[s][i])*convol[s][i]);
    }
  }
  for (int s=0; s<nsite; s++) {
    fftnd(convol[s],3,supercell.get_buf(),-1);
  }
  cur_recip_E/=(rsize*rspin_size);
//  cerr << "a_total=" << E_ref+cur_energy+cur_recip_E+mu*cur_conc << endl;
  fft_time=clock()-fft_time;
}

void KSpaceMonteCarlo::extension_update_spin_flip(int *cell, int incell, int newspin, Real d_recip_energy) {
  flipped_spins.push_front(new FlippedSpin(cell,incell,newspin,d_recip_energy));
  cur_recip_E+=d_recip_energy/rspin_size;
  int nbflip=flipped_spins.get_size();
  if (fabs(get_cur_concentration()-ref_x)>threshold_dx && (nbflip % 2)==0) {
    extension_calc_from_scratch();
  } else {
    if (flip_time>0) {
//cerr << "time: " << nbflip << " " << flip_time << " " << 2*fft_time/flip_time << endl;
      if (nbflip > 2*fft_time/flip_time && (nbflip % 2)==0) {
	extension_calc_from_scratch();
      }
    }
  }
}

void KSpaceMonteCarlo::extension_undo_spin_flip(void) {
  LinkedListIterator<FlippedSpin> last_flip(flipped_spins);
  cur_recip_E-=(last_flip->d_recip_energy)/rspin_size;
  delete flipped_spins.detach(last_flip);
}

Real KSpaceMonteCarlo::extension_spin_flip_energy(int *cell, int incell, int newspin) {
//  cerr << "f_total=" << E_ref+cur_energy+cur_recip_E+mu*cur_conc << endl;
  flip_time=clock();
  int dr[3];
  int *psupercell=supercell.get_buf();
  Real rspinchange=2.*(Real)newspin;
  //int offset=((cell[0]*psupercell[1] + cell[1])*psupercell[2] + cell[2]); //permbug
  int offset=((cell[2]*psupercell[1] + cell[1])*psupercell[0] + cell[0]);
  Real dE=4.*real(dir_eci[unrollsym(incell,incell)][0]);
  dE+=2.*real(convol[incell][offset]*rspinchange);
  LinkedListIterator<FlippedSpin> i(flipped_spins);
  for (; i; i++) {
    for (int j=0; j<3; j++) {
      dr[j]=(i->cell[j]-cell[j]+psupercell[j]) % psupercell[j];
    }
    //int offsetdr=((dr[0]*psupercell[1] + dr[1])*psupercell[2] + dr[2]); //permbug
    int offsetdr=((dr[2]*psupercell[1] + dr[1])*psupercell[0] + dr[0]);
    dE+=2.*real(dir_eci[unrollsym(i->incell,incell)][offsetdr])*2.*(i->newspin)*rspinchange;
  }
  flip_time=clock()-flip_time;
  return dE;
}

KSpaceMonteCarlo::~KSpaceMonteCarlo() {
  // free all dynamically allocated memory;
}

Array<Real> mclibdummyarray;
