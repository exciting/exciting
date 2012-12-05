#include "chull.h"
#include "lstsqr.h"
#include "getvalue.h"

class EdgeFlag {
public:
  int done;
  Array<int> edge;
  EdgeFlag(void): done(0),edge(0) {}
  //EdgeFlag(const EdgeFlag& e): done(e.done),edge(e.edge) {}
  EdgeFlag(int d, const Array<int>& e): done(d),edge(e) {}
};

void recurse_convex_hull(LinkedList<PolytopeFace> *hull, LinkedList<EdgeFlag> *edges, const Array<Array<Real> > &x, Array<int> edge, Array<Real> fwrd, Array<Real> up) {
  int dim=x(0).get_size();

  //  cerr << "RECURSE" << endl;
  //  cerr << up << endl;
  //  cerr << fwrd << endl;
  //  cerr << edge << endl;
  //  LinkedListIterator<EdgeFlag> ii(*edges);
  //  for (; ii; ii++) {
    //    cerr << ii->done << endl;
    //    cerr << ii->edge << endl;
  //  }

  // look for the point that requires the least tilting to be touched by the plane perpendicular to "up";
  Array<Real> midedge(x(edge(0)).get_size());
  zero_array(&midedge);
  for (int j=0; j<edge.get_size(); j++) {
    sum(&midedge,midedge,x(edge(j)));
  }
  product(&midedge,midedge,1./(Real)(edge.get_size()));
  Real best_theta=-M_PI;
  Real best_dist=0.;
  int best_i=-1;
  for (int i=0; i<x.get_size(); i++) {
    Array<Real> dx;
    diff(&dx,x(i),midedge);
    Real f=inner_product(fwrd,dx);
    Real u=inner_product(up,dx);
    Real theta=atan2(f,-u);
    Real distp=sqrt(u*u+f*f);
    Real dist=sqrt(inner_product(dx,dx));
    //    cerr << i << " " << theta << " " << dist << endl;
    if (distp > zero_tolerance && (theta>best_theta+zero_tolerance || (fabs(theta-best_theta)<zero_tolerance && dist>best_dist))) {
      best_theta=theta;
      best_dist=dist;
      best_i=i;
    }
  }
  //  cerr << best_i << endl;
  if (best_i==-1) return;
  // the initial edge and this new point make up a new face;
  PolytopeFace *pface=new PolytopeFace();
  pface->pts.resize(dim);
  for (int j=0; j<dim-1; j++) {
    pface->pts(j)=edge(j);
  }
  pface->pts(dim-1)=best_i;
  //  cerr << pface->pts << endl;

  // which way is up on this new face;
  Array2d<Real> eqs(dim,dim);
  Array<Real> rhs(dim);
  one_array(&rhs);
  for (int j=0; j<dim; j++) {
    //    cerr << x(pface->pts(j)) << endl;
    for (int k=0; k<dim; k++) {
      eqs(j,k)=x(pface->pts(j))(k);
    }
  }
  solve_linsys(&(pface->up),eqs,rhs);
  normalize(&(pface->up));

  (*hull) << pface;

  // create edges;
  Array<Array<int> > newedge(dim);
  for (int d=0; d<dim; d++) {
    newedge(d).resize(dim-1);
    int l=0;
    for (int k=0; k<dim; k++) {
      if (k!=d) {
	newedge(d)(l)=pface->pts(k);
	l++;
      }
    }
    sort_array(&(newedge(d)));
    // check for overlapping plane;
    LinkedListIterator<EdgeFlag> i(*edges);
    for ( ; i; i++) {
      if (i->edge == newedge(d)) break;
    }
    if (i) {
      if (i->done==2) {return;}
    }
  }

  // now we keep track of the new edges created so not to do them again later;
  for (int d=0; d<dim; d++) {
    LinkedListIterator<EdgeFlag> i(*edges);
    for ( ; i; i++) {
      if (i->edge == newedge(d)) break;
    }
    if (i) {
      i->done++;
    }
    else {
      (*edges) << new EdgeFlag(1,newedge(d));
    }
  }

  // now recurse for each new edges;
  for (int d=0; d<dim; d++) {
    int newpt=pface->pts(d);
    LinkedListIterator<EdgeFlag> i(*edges);
    for ( ; i; i++) {
      if (i->edge == newedge(d)) break;
    }
    if (i->done<=1) {
      // creates the new forward vector for the current new edge;
      Array2d<Real> proj_onto(dim,dim-2);
      for (int j=0; j<dim-2; j++) {
	for (int k=0; k<dim; k++) {
	  proj_onto(k,j)=x(newedge(d)(j+1))(k)-x(newedge(d)(j))(k);
	}
      }
      Array<Real> dx;
      diff(&dx,x(newedge(d)(0)),x(newpt));
      Array<Real> newfwrd;
      predict_ols(&newfwrd,proj_onto,dx);
      // cerr << proj_onto << endl;
      diff(&newfwrd,dx,newfwrd);
      product(&newfwrd,newfwrd,1./sqrt(inner_product(newfwrd,newfwrd)));
      //      cerr << "newfwrd " << newfwrd << endl;
      // now lookup for the face that touches this new edge;
      recurse_convex_hull(hull,edges,x, newedge(d), newfwrd, pface->up);
    }
  }
}

void de_mean(Array<Array<Real> > *px, const Array<Array<Real> > &org_x) {
  int dim=org_x(0).get_size();
  Array<Real> meanx(dim);
  zero_array(&meanx);
  for (int i=0; i<org_x.get_size(); i++) {
    sum(&meanx,meanx,org_x(i));
  }
  product(&meanx,meanx,1./(Real)org_x.get_size());
  px->resize(org_x.get_size());
  for (int i=0; i<org_x.get_size(); i++) {
    diff(&((*px)(i)),org_x(i),meanx);
  }
}

void calc_convex_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x, const Array<Real> &ground) {
  int dim=org_x(0).get_size();
  // calculate mean of x and substract it out, so that origin is in convex hull;
  if (org_x.get_size()==0) return;
  Array<Array<Real> > x;
  de_mean(&x,org_x);

  // find arbitrary starting point, the i one with largest x(i)(0);
  int ixmax;
  Real xmax=-MAXFLOAT;
  for (int i=0; i<x.get_size(); i++) {
    if (x(i)(0)>xmax) {
      ixmax=i;
      xmax=x(i)(0);
    }
  }

  // find starting edge by tilting starting plane until it touches dim-1 points;
  Array<Real> fwrd;
  Array<int> edge(dim-1);
  edge(0)=ixmax;
  Array2d<Real> proj_onto;
  for (int d=1; d<dim; d++) {
    int ibest=-1;
    Real best_proj=-1.;
    Real best_ldx=0.;
    for (int i=0; i<x.get_size(); i++) {
      Array<Real> dx,pdx;
      Real ldx;
      diff(&dx,x(i),x(ixmax));
      if (d>1) {
        predict_ols(&pdx,proj_onto,dx);
	diff(&dx,dx,pdx);
      }
      ldx=normalize(&dx);
      if (ldx>zero_tolerance && (dx(0)>best_proj+zero_tolerance || (fabs(dx(0)-best_proj)<zero_tolerance && (ldx>best_ldx)) ) ) {
        ibest=i;
        best_proj=dx(0);
        best_ldx=ldx;
        if (d==dim-1) {
          fwrd=dx;
        }
      }
    }
    if (ibest==-1) return;
    if (d<dim-1) {
      edge(d)=ibest;
      proj_onto.resize(iVector2d(dim,d));
      for (int j=1; j<=d; j++) {
        for (int k=0; k<dim; k++) {
          proj_onto(k,j-1)=x(edge(j))(k)-x(ixmax)(k);
        }
      }
    }
  }
  sort_array(&edge);
  Array<Real> up(dim);
  zero_array(&up);
  up(0)=1.;
  if (dim<=2) {
    fwrd.resize(dim);
    zero_array(&fwrd);
    fwrd(1)=1.;
  }
  else {
    Array<Real> pup;
    predict_ols(&pup,proj_onto,up);
    diff(&up,up,pup);
    normalize(&up);
    normalize(&fwrd);
  }
  LinkedList<EdgeFlag> edges;
  recurse_convex_hull(hull, &edges, x, edge, fwrd, up);
  LinkedListIterator<PolytopeFace> i(*hull);
  for (; i; i++) {
    i->c=inner_product(i->up,org_x(i->pts(0)));
  }
  if (inner_product(ground,ground)>sqr(zero_tolerance)) {
    LinkedListIterator<PolytopeFace> i(*hull);
    while (i) {
      if (inner_product(ground,i->up)>zero_tolerance) {
	i++;
      }
      else {
	delete hull->detach(i);
      }
    }
  }
}

void calc_convex_hull_degenerate(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x, const Array<Real> &ground) {
  if (org_x.get_size()==0) return;
  int dmax=org_x(0).get_size();
  Array<Array<Real> > basis(1);
  basis(0)=ground;
  int d=0;
  for (int i=0; i<org_x.get_size(); i++) {
    if (d==dmax) break;
    build_basis(&basis,&d,org_x(i));
  }
  Array2d<Real> bb(d,d),ibb;
  for (int i=0; i<d; i++) {
    for (int j=0; j<d; j++) {
      bb(i,j)=inner_product(basis(i),basis(j));
    }
  }
  invert_matrix(&ibb,bb);
  Array2d<Real> b(d,dmax);
  for (int i=0; i<d; i++) {
    for (int j=0; j<dmax; j++) {
      b(i,j)=basis(i)(j);
    }
  }
  Array2d<Real> proj(d,dmax);
  product(&proj,ibb,b);
  Array<Array<Real> > x(org_x.get_size());
  for (int i=0; i<x.get_size(); i++) {
    product(&(x(i)),proj,org_x(i));
  }
  Array<Real> tground;
  product(&tground,proj,ground);
  calc_convex_hull(hull, x, tground);
  update_normals(hull, org_x);
}

void update_normals(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x) {
  int dim=org_x(0).get_size();
  Array<Array<Real> > x;
  de_mean(&x,org_x);
  LinkedListIterator<PolytopeFace> i(*hull);
  for (; i; i++) {
    Array2d<Real> eqs(i->pts.get_size(),dim);
    Array<Real> rhs(i->pts.get_size());
    one_array(&rhs);
    for (int j=0; j<i->pts.get_size(); j++) {
      for (int k=0; k<dim; k++) {
	eqs(j,k)=x(i->pts(j))(k);
      }
    }
#ifdef DEGE_GS
    Array2d<Real> ieqs;
    invert_matrix_wide(&ieqs,eqs);
    product(&(i->up),ieqs,rhs);
#else
    solve_linsys(&(i->up),eqs,rhs);
#endif
    normalize(&(i->up));
    i->c=inner_product(i->up,org_x(i->pts(0)));
  }
}

int flag_outside_hull(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, int flag_plane_too) {
  int oneoutside=0;
  poutside->resize(x.get_size());
  zero_array(poutside);
  LinkedListIterator<PolytopeFace> i(*hull);
  for (; i; i++) {
    int badplane=0;
    for (int j=0; j<x.get_size(); j++) {
      if (inner_product(i->up,x(j))>i->c+zero_tolerance) {
	(*poutside)(j)=1;
	oneoutside=1;
	badplane=1;
      }
    }
    if (flag_plane_too && badplane) {
      for (int j=0; j<i->pts.get_size(); j++) {
	(*poutside)(i->pts(j))=1;
      }
    }
  }
  return oneoutside;
}

void clip_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const LinkedList<LinearInequality> &ineq) {
  LinkedListIterator<PolytopeFace> i(*hull);
  while (i) {
    int donotcut=0;
    for (int k=0; k<i->pts.get_size(); k++) {
      int cutit=0;
      LinkedListIterator<LinearInequality> j(ineq);
      for (; j; j++) {
	if (inner_product(x(i->pts(k)),j->v)>j->c+zero_tolerance) {
	  cutit=1;
	  break;
	}
      }
      if (!cutit) {
	donotcut=1;
	break;
      }
    }
    if (donotcut) {
      i++;
    }
    else {
      delete hull->detach(i);
    }
  }
}

void calc_formation(Array<Real> *pfe, Array<Real> *ppure, const Array<Array<Real> > &x, const Array<Real> &e) {
  if (ppure->get_size()==0) {
    LinkedList<Array<Real> > pts;
    LinkedList<Real> es;
    for (int i=0; i<x.get_size(); i++) {
      LinkedListIterator<Array<Real> > insertp(pts);
      LinkedListIterator<Real> inserte(es);
      Real l2=inner_product(x(i),x(i));
      while (insertp) {
	Real lcurp=inner_product(*insertp,*insertp);
	if (l2>lcurp) break;
	if (fabs(l2-lcurp)<zero_tolerance && e(i)<(*inserte)) break;
	insertp++;
	inserte++;
      }
      pts.add(new Array<Real>(x(i)),insertp);
      es.add(new Real(e(i)),inserte);
    }
    
    Array<Array<Real> > basis;
    int nb_in_basis=0;
    {
      LinkedListIterator<Array<Real> > ip(pts);
      LinkedListIterator<Real> ie(es);
      while (ip) {
	if (!build_basis(&basis,&nb_in_basis,*ip)) {
	  delete pts.detach(ip);
	  delete es.detach(ie);
	}
	else {
	  ip++;
	  ie++;
	}
      }
    }
    if (nb_in_basis==0) {
      *pfe=e;
      return;
    }
    int dim=basis(0).get_size();
    Array2d<Real> X(nb_in_basis,dim);
    LinkedListIterator<Array<Real> > ip(pts);
    for (int j=0; j<nb_in_basis; j++, ip++) {
      for (int k=0; k<dim; k++) {
	X(j,k)=(*ip)(k);
      }
    }
    // cerr << X << endl;
    Array<Real> r;
    LinkedList_to_Array(&r,es);
    // cerr << r << endl;
    Array2d<Real> inv,M;
    outer_product(&inv,X,X);
    invert_matrix(&inv,inv);
    inner_product(&M,X,inv);
    product(ppure,M,r);
  }
  pfe->resize(e.get_size());
  for (int i=0; i<e.get_size(); i++) {
    (*pfe)(i)=e(i)-inner_product(*ppure,x(i));
  }
}

int is_point_in_hull(const Array<Real> &x, const LinkedList<LinearInequality> &ineq_list) {
  LinkedListIterator<LinearInequality> i(ineq_list);
  for (; i; i++) {
    Real a=0.;
    for (int j=0; j<x.get_size(); j++) {
      a+=x(j)*(i->v(j));
    }
    if (a > i->c+zero_tolerance) break;
  }
  return (!i);
}

void read_inequalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s) {
  skip_delim(s," \t\n");
  while (!s.eof()) {
    LinearInequality ineq(label.get_size());
    char dir;
    zero_array(&(ineq.v));
    while (1) {
      skip_delim(s," \t+");
      Real a=MAXFLOAT;
      s >> a;
      if (a==MAXFLOAT || !s) {
	ERRORQUIT("Error parsing inequality");
      }
      skip_delim(s,"* \t\n");
      AutoString w;
      get_string(&w,s,"\t\n +-<>");
      if (!is_in_array(label,w)) {
        cerr << "unknown atom label " << w << endl;
        ERRORQUIT("Aborting.");
      }
      ineq.v(index_in_array(label,w))+=a;
      skip_delim(s," \t\n");
      dir=s.get();
      if (strchr("<>",dir)) {break;}
      s.putback(dir);
    }
    skip_delim(s,"=");
    s >> ineq.c;
    if (dir=='>') {
      product(&(ineq.v),ineq.v,-1.);
      ineq.c=-ineq.c;
    }
    (*ineq_list) << new LinearInequality(ineq);
    skip_delim(s,", \t\n");
  }
}

void read_equalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s) {
  skip_delim(s," \t\n");
  while (!s.eof()) {
    LinearInequality ineq(label.get_size());
    char dir;
    zero_array(&(ineq.v));
    while (1) {
      skip_delim(s," \t+");
      Real a=MAXFLOAT;
      s >> a;
      if (a==MAXFLOAT || !s) {
	ERRORQUIT("Error parsing equality");
      }
      skip_delim(s,"* \t\n");
      AutoString w;
      get_string(&w,s,"\t\n +-=");
      if (!is_in_array(label,w)) {
        cerr << "unknown atom label " << w << endl;
        ERRORQUIT("Aborting.");
      }
      ineq.v(index_in_array(label,w))+=a;
      skip_delim(s," \t\n");
      dir=s.get();
      if (dir=='=') {break;}
      s.putback(dir);
    }
    skip_delim(s,"=");
    s >> ineq.c;
    (*ineq_list) << new LinearInequality(ineq);
    skip_delim(s,", \t\n");
  }
}

void ineq_list_to_Array(Array2d<Real> *a, Array<Real> *c, const LinkedList<LinearInequality> &l) {
  LinkedListIterator<LinearInequality> it(l);
  iVector2d size(l.get_size(),it->v.get_size());
  if (a) {
    a->resize(size);
    for (int i=0; i<size(0); i++) {
      for (int j=0; j<size(1); j++) {
	(*a)(i,j)=(it->v(j));
      }
      it++;
    }
  }
  it.init(l);
  if (c) {
    c->resize(size(0));
    for (int i=0; i<size(0); i++) {
      (*c)(i)=it->c;
    }
  }
}

void calc_convex_hull_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e) {
  Array<Array<Real> > allx;
  paste_col_vector(&allx,x,e);
  Array<Real> ground(allx(0).get_size());
  zero_array(&ground);
  ground(ground.get_size()-1)=-1.;
  calc_convex_hull(hull,allx,ground);
}

void update_normals_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e) {
  Array<Array<Real> > allx;
  paste_col_vector(&allx,x,e);
  update_normals(hull,allx);
}

int flag_outside_hull_p1(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e, int flag_plane_too) {
  Array<Array<Real> > allx;
  paste_col_vector(&allx,x,e);
  return flag_outside_hull(poutside,hull,allx,flag_plane_too);
}
