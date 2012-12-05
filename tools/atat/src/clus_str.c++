#include "clus_str.h"
#include "arraylist.h"

ClusterBank::ClusterBank(const rMatrix3d &cell,
              const Array<rVector3d> &atom_pos,
              int ntuple,
              const SpaceGroup &_equivalent_by_symmetry):
  cluster_list(), current_cluster(),
  current_index(0), previous_length(0.), multiplet(cell,atom_pos,ntuple),
  equivalent_by_symmetry(_equivalent_by_symmetry) {
    if (ntuple>1) {
      cluster_list << new Cluster(multiplet);
    }
    else {
       for (int at=0; at<atom_pos.get_size(); at++) {
         Cluster point(1);
         point(0)=atom_pos(at);
         add_unique(&cluster_list,point,equivalent_by_symmetry);
       }
    }
    current_cluster.init(cluster_list);
}

ClusterBank::ClusterBank(const ClusterBank &clusterbank,int ntuple):
  cluster_list(), current_cluster(),
  current_index(0),
  previous_length(0.),
  multiplet(clusterbank.multiplet,ntuple),
  equivalent_by_symmetry(clusterbank.equivalent_by_symmetry) {
    cluster_list << new Cluster(multiplet);
    current_cluster.init(cluster_list);
}

void ClusterBank::reset(void) {
  current_index=0;
  current_cluster.init(cluster_list);
}

void ClusterBank::operator++(int) {
  previous_length=get_length_quick(*current_cluster);
  LinkedListIterator<Cluster> save=current_cluster;
  current_cluster++;
  current_index++;
  if (!current_cluster && save->get_size()>1) {
    do {multiplet++;} while (!add_unique(&cluster_list,Cluster(multiplet),equivalent_by_symmetry));
    current_cluster=save;
    current_cluster++;
  }
}

MultiClusterBank::MultiClusterBank(const Structure &_lat,
				   int ntuple,
				   const SpaceGroup &_equivalent_by_symmetry):
  lat(_lat), cluster_bank(_lat.cell,_lat.atom_pos,ntuple,_equivalent_by_symmetry),
    cluster_list(), current_cluster(),
    current_index(0), previous_length(0.), 
    equivalent_by_symmetry(_equivalent_by_symmetry) {
  if (ntuple>1) {
    make_new_clusters();
  }
  else {
    for (int at=0; at<lat.atom_pos.get_size(); at++) {
      MultiCluster point(1);
      point.clus(0)=lat.atom_pos(at);
      point.site_type(0)=(lat.atom_type(at)-2);
      for (int t=0; t<lat.atom_type(at)-1; t++) {
	point.func(0)=t;
	add_unique(&cluster_list,point,equivalent_by_symmetry);
      }
    }
  }
  current_cluster.init(cluster_list);
}

MultiClusterBank::MultiClusterBank(const MultiClusterBank &bank, int ntuple):
  lat(bank.lat), cluster_bank(bank.lat.cell,bank.lat.atom_pos,ntuple,bank.equivalent_by_symmetry),
    cluster_list(), current_cluster(),
    current_index(0), previous_length(0.), 
    equivalent_by_symmetry(bank.equivalent_by_symmetry) {
  make_new_clusters();
  current_cluster.init(cluster_list);
}

void MultiClusterBank::reset(void) {
  current_index=0;
  current_cluster.init(cluster_list);
}

void MultiClusterBank::make_new_clusters(void) {
  MultiCluster new_clus;
  new_clus.clus=cluster_bank;
  new_clus.site_type.resize(new_clus.clus.get_size());
  new_clus.func.resize(new_clus.clus.get_size());
  Array<int> minconfig(new_clus.clus.get_size());
  zero_array(&minconfig);
  rMatrix3d inv_cell=!(lat.cell);
  for (int at=0; at<new_clus.clus.get_size(); at++) {
    new_clus.site_type(at)=lat.atom_type(which_atom(lat.atom_pos,new_clus.clus(at),inv_cell))-2;
  }
  LinkedList<MultiCluster> new_clus_list;
  MultiDimIterator<Array<int> > config(minconfig,new_clus.site_type);
  for (; config; config++) {
    new_clus.func=config;
    add_unique(&new_clus_list,new_clus,equivalent_by_symmetry);
  }
  transfer_list(&cluster_list,&new_clus_list);
}

void MultiClusterBank::operator++(int) {
  previous_length=get_length_quick(current_cluster->clus);
  current_cluster++;
  current_index++;
  if (!current_cluster && previous_length!=0) {
    cluster_bank++;
    make_new_clusters();
  }
}


