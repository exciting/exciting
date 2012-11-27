#include <fstream.h>
#include "integer.h"
#include "arraylist.h"
#include "linalg.h"
#include "getvalue.h"

void output_medit(ostream &file, const LinkedList<Array<int> > &triangles, const LinkedList<int> &tritype, const Array<Array<Real> > &pts) {
  Real s=sqrt(3.)/2.;
  Real c=0.5;
  Real Ts=1e-3;
  file << "MeshVersionFormatted 0" << endl << "Dimension"  << endl << 3 << endl;
  file << "Vertices" << endl;
  file << pts.get_size() << endl;
  for (int i=0; i<pts.get_size(); i++) {
    file << pts(i)(1)+c*pts(i)(2) << " " << s*pts(i)(2) << " " << Ts*pts(i)(0) << " " << 1 << endl;
  }
  file << "Triangles" << endl;
  file << triangles.get_size() << endl;
  LinkedListIterator<Array<int> > it(triangles);
  LinkedListIterator<int> itp(tritype);
  for ( ; it; it++, itp++) {
    file << 1+(*it)(0) << " " << 1+(*it)(1) << " " << 1+(*it)(2) << " " << (*itp) << endl;
  }
}

void triangulate(LinkedList<Array<int> > *ptriangles, LinkedList<int> *ptritype, const Array<Array<Real> > &pts, int offset=0, int type=1) {
  int ze=0;
  LinkedList<LinkedList<int> > rings;
  for (int i=0; i<pts.get_size(); i++) {
    LinkedListIterator<LinkedList<int> > j(rings);
    for ( ; j; j++) {
      LinkedListIterator<int> k(*j);
      if (pts(*k)(ze) > pts(i)(ze) - zero_tolerance) break;
    }
    int done=0;
    if (j) {
      LinkedListIterator<int> k(*j);
      if (near_zero(pts(*k)(ze)-pts(i)(ze))) {
        (*j) << new int(i);
        done=1;
      }
    }
    if (!done) {
      LinkedList<int> *pll=new LinkedList<int>;
      (*pll) << new int(i);
      rings.add(pll,j);
    }
  }
  /*  {
    LinkedListIterator<LinkedList<int> > j(rings);
    for (; j; j++) {
      LinkedListIterator<int> k(*j);
      for ( ; k; k++) {
	cout << pts(*k)(0) << " " << pts(*k)(1) << " " << pts(*k)(2) << endl;
      }
      cout << endl;
    }
  }
  */
  LinkedListIterator<LinkedList<int> > j(rings);
  Array<Real> start;
  Array<Real> lastpt;
  for ( ; j; j++) {
    if (start.get_size()==0) {
      LinkedListIterator<int> k(*j);
      start=pts(*k);
    }
    int ifirstpt;
    LinkedList<int> sorted;
    lastpt=start;
    while (j->get_size()>0) {
      LinkedListIterator<int> k(*j);
      LinkedListIterator<int> bestk;
      Real bestd=MAXFLOAT;
      for (; k; k++) {
        Array<Real> tmp;
        diff(&tmp,lastpt,pts(*k));
        Real d=inner_product(tmp,tmp);
        if (d<bestd) {
	    bestd=d;
	    bestk=k;
	  }
      }
      lastpt=pts(*bestk);
      if (sorted.get_size()==0) {ifirstpt=*bestk;}
      sorted << j->detach(bestk);
    }
    //sorted << new int(ifirstpt);
    start=pts(ifirstpt);
    transfer_list(&(*j),&sorted);
  }
  /*  {
    LinkedListIterator<LinkedList<int> > j(rings);
    for (; j; j++) {
      LinkedListIterator<int> k(*j);
      for ( ; k; k++) {
	cout << pts(*k)(0) << " " << pts(*k)(1) << " " << pts(*k)(2) << endl;
      }
      cout << endl;
    }
    }*/

  LinkedListIterator<LinkedList<int> > j1(rings);
  LinkedListIterator<LinkedList<int> > j2(rings);
  if (!j2) {ERRORQUIT("Not enough z values.");}
  j2++;
  for ( ; j2; j1++, j2++) {
    LinkedListIterator<int> k1(*j1);
    LinkedListIterator<int> k2(*j2);
    while (k1 && k2) {
      LinkedListIterator<int> k1n(k1);
      LinkedListIterator<int> k2n(k2);
      k1n++;
      k2n++;
      if (!k1n && !k2n) {break;}
      int pick;
      if (!k1n) {pick=2;}
      else if (!k2n) {pick=1;}
      else {
	Array<Real> m,d;
	sum(&m,pts(*k1),pts(*k2));
	product(&m,m,0.5);
	diff(&d,pts(*k1),pts(*k2));
	Real ld2=inner_product(d,d);
	Array<Real> r1,r2;
	diff(&r1,pts(*k1n),m);
	diff(&r2,pts(*k2n),m);
	// Real p1=fabs(inner_product(r1,d)/sqrt(ld2*inner_product(r1,r1)));
	// Real p2=fabs(inner_product(r2,d)/sqrt(ld2*inner_product(r2,r2)));
	Real p1=inner_product(r1,r1);
	Real p2=inner_product(r2,r2);
	if (p1>p2) {pick=2;} else {pick=1;}
      }
      Array<int> *ptri=new Array<int>(3);
      (*ptri)(0)=*k1;
      (*ptri)(1)=*k2;
      if (pick==2) {
	(*ptri)(2)=offset+*k2n;
	k2=k2n;
      }
      else {
	(*ptri)(2)=offset+*k1n;
	k1=k1n;
      }
      (*ptriangles) << ptri;
      if (ptritype) {*ptritype << new int(type);}
    }
  }
}

int main(int argc, char *argv[]) {
  LinkedList<Array<Real> > ptslist[2];
  AutoString oldline[2];
  int nbl=0;
  while (!cin.eof()) {
    AutoString line;
    get_string(&line,cin,"\n");
    cin.get();
    if (strlen(line)==0) {
      if (nbl>=2) {
	for (int p=0; p<2; p++) {
	  istrstream buf(oldline[p]);
	  Array<Real> pt(3);
	  Real tmp;
	  buf >> pt(0) >> pt(1) >> pt(2);
	  ptslist[p] << new Array<Real>(pt);
	}
      }
      nbl=0;
    }
    else {
      oldline[0]=oldline[1];
      oldline[1]=line;
      nbl++;
    }
  }
  LinkedList<Array<int> > tri_list;
  LinkedList<int> tritype_list;
  Array<Array<Real> > ptsarray;
  LinkedList_to_Array(&ptsarray,ptslist[0]);
  triangulate(&tri_list,&tritype_list,ptsarray,0,1);
  if (0) {
    LinkedList_to_Array(&ptsarray,ptslist[1]);
    triangulate(&tri_list,&tritype_list,ptsarray,ptslist[0].get_size(),2);
    transfer_list(&(ptslist[0]),&(ptslist[1]));
    LinkedList_to_Array(&ptsarray,ptslist[0]);
  }
  output_medit(cout,tri_list,tritype_list,ptsarray); 
}
