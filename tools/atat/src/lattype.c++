#include "lattype.h"
#include "findsym.h"

rMatrix3d find_reduced_cell(const rMatrix3d &cell) {
  if (is_reduced_cell(cell)) {return cell;}
  Real d=fabs(det(cell));
  LatticePointIterator lp[3];
  rMatrix3d curcell;
  int done=0;
  lp[2].init(cell,1);

  rVector3d a=(rVector3d)lp[2];
  while (near_zero(norm(a^(rVector3d)lp[2]))) {lp[2]++;}
  rVector3d b=(rVector3d)lp[2];
  while (near_zero((b^a)*(rVector3d)lp[2])) {lp[2]++;}
  Real lc=norm((rVector3d)lp[2]);
  lp[2].init(cell,1);
  while (norm((rVector3d)lp[2])<lc-zero_tolerance) {lp[2]++;}

  while (!done) {
    lp[1].init(cell,1);
    while (norm2(lp[1])<=norm2(lp[2])+zero_tolerance && !done) {
      lp[0].init(cell,1);
      while (norm2(lp[0])<=norm2(lp[1])+zero_tolerance && !done) {
	for  (int i=0; i<3; i++) {
	  curcell.set_column(i,lp[i]);
	}
	if (near_zero(fabs(det(curcell))-d)) {
	  if (is_reduced_cell(curcell)) done=1;
	}
	lp[0]++;
      }
      lp[1]++;
    }
    lp[2]++;
  }
  return curcell;
}

int is_reduced_cell(const rMatrix3d &cell) {
  Real tiny=zero_tolerance;
  if (det(cell)<tiny) return 0;
  rVector3d a=cell.get_column(0);
  rVector3d b=cell.get_column(1);
  rVector3d c=cell.get_column(2);
  Real aa=a*a;
  Real bb=b*b;
  Real cc=c*c;
  Real ab=a*b;
  Real bc=b*c;
  Real ac=a*c;

  if (! (aa<=bb+tiny && bb<=cc+tiny) ) return 0;
  if (! (fabs(bc)<=bb/2.+tiny && fabs(ac)<=aa/2.+tiny && fabs(ab)<=aa/2.+tiny)) return 0;
  if (bc>tiny && ac>tiny && ab>tiny) { // Type I cell;
    if (near_zero(aa-bb)) {
      if (! (bc <= ac+tiny)) return 0;
    }
    if (near_zero(bb-cc)) {
      if (! (ac<=ab+tiny)) return 0;
    }
    if (near_zero(bc-bb/2.)) {
      if (! (ab<=2.*ac+tiny)) return 0;
    }
    if (near_zero(ac-aa/2.)) {
      if (! (ab<=2.*bc+tiny)) return 0;
    }
    if (near_zero(ab-aa/2.)) {
      if (! (ac<=2.*bc+tiny)) return 0;
    }
    return 1;
  }
  else if ((bc<=tiny && ac<=tiny && ab<=tiny)) { // Type II cell;
    if (! (fabs(bc)+fabs(ac)+fabs(ab)<=(aa+bb)/2.+tiny)) return 0;
    if (near_zero(aa-bb)) {
      if (! (fabs(bc) <= fabs(ac)+tiny)) return 0;
    }
    if (near_zero(bb-cc)) {
      if (! (fabs(ac)<=fabs(ab)+tiny)) return 0;
    }
    if (near_zero(fabs(bc)-bb/2.)) {
      if (! (near_zero(ab))) return 0;
    }
    if (near_zero(fabs(ac)-aa/2.)) {
      if (! (near_zero(ab))) return 0;
    }
    if (near_zero(fabs(ab)-aa/2.)) {
      if (! (near_zero(ac))) return 0;
    }
    if (near_zero(fabs(bc)+fabs(ac)+fabs(ab)-(aa+bb)/2.)) {
      if (! (aa <= 2.*fabs(ac)+fabs(ab)+tiny) ) return 0;
    }
    return 2;
  }
  else return 0;
}

int tconv[][9]={
  { 1,-1, 1,   1, 1,-1,   -1, 1, 1}, // 1;
  { 1,-1, 0,  -1, 0, 1,   -1,-1,-1}, // 2;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 3;
  { 1, 0, 1,   1, 1, 0,    0, 1, 1}, // 5;
  { 1,-1, 0,  -1, 0, 1,   -1,-1,-1}, // 4;
  { 0, 1, 1,   1, 0, 1,    1, 1, 0}, // 6;
  { 1, 0, 1,   1, 1, 0,    0, 1, 1}, // 7;
  {-1,-1, 0,  -1, 0,-1,    0,-1,-1}, // 8;
  { 1, 0, 0,  -1, 1, 0,   -1,-1, 3}, // 9;
  { 1, 1, 0,   1,-1, 0,    0, 0,-1}, // 10;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 11;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 12;
  { 1, 1, 0,  -1, 1, 0,    0, 0, 1}, // 13;
  { 1, 0, 0,   0, 1, 0,    1, 1, 2}, // 15;
  {-1,-1, 0,   1,-1, 0,    1, 1, 2}, // 16;
  { 1, 1, 0,  -1, 1, 0,    0, 0, 1}, // 14;
  { 1,-1, 0,   1, 1, 0,   -1, 0,-1}, // 17;
  { 0,-1, 1,   1,-1,-1,    1, 0, 0}, // 18;
  {-1, 0, 0,   0,-1, 1,   -1, 1, 1}, // 19;
  { 0, 1, 1,   0, 1,-1,   -1, 0, 0}, // 20;
  { 0, 1, 0,   0, 0, 1,    1, 0, 0}, // 21;
  { 0, 1, 0,   0, 0, 1,    1, 0, 0}, // 22;
  { 0, 1, 1,   0,-1, 1,    1, 0, 0}, // 23;
  { 1, 2, 1,   0,-1, 1,    1, 0, 0}, // 24;
  { 0, 1, 1,   0,-1, 1,    1, 0, 0}, // 25;
  { 1, 0, 0,  -1, 2, 0,   -1, 0, 2}, // 26;
  {-1, 2, 0,  -1, 0, 0,    0,-1, 1}, // 27;
  {-1, 0, 0,  -1, 0, 2,    0, 1, 0}, // 28;
  { 1, 0, 0,   1,-2, 0,    0, 0,-1}, // 29;
  { 0, 1, 0,   0, 1,-2,   -1, 0, 0}, // 30;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 31;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 32;
  { 0,-1, 0,   0, 1, 2,   -1, 0, 0}, // 40;
  { 0,-1, 0,  -1, 0, 0,    0, 0,-1}, // 35;
  { 1, 0, 0,  -1, 0,-2,    0, 1, 0}, // 36;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}, // 33;
  {-1, 0, 0,   1, 2, 0,    0, 0,-1}, // 38;
  {-1, 0, 0,   0, 0,-1,    0,-1, 0}, // 34;
  {-1, 0, 0,   0,-1, 0,    1, 1, 2}, // 42;
  { 0,-1,-2,   0,-1, 0,   -1, 0, 0}, // 41;
  { 1, 0, 2,   1, 0, 0,    0, 1, 0}, // 37;
  {-1,-2, 0,  -1, 0, 0,    0, 0,-1}, // 39;
  {-1, 0, 0,  -1,-1,-2,    0,-1, 0}, // 43;
  { 1, 0, 0,   0, 1, 0,    0, 0, 1}  // 44;
};

BravaisType bt_table[]={
  BT_cF, // 1;
  BT_hR, // 2;
  BT_cP, // 3;
  BT_cI, // 5;
  BT_hR, // 4;
  BT_tI, // 6;
  BT_tI, // 7;
  BT_oI, // 8;
  BT_hR, // 9;
  BT_mC, // 10;
  BT_tP, // 11;
  BT_hP, // 12;
  BT_oC, // 13;
  BT_tI, // 15;
  BT_oF, // 16;
  BT_mC, // 14;
  BT_mC, // 17;
  BT_tI, // 18;
  BT_oI, // 19;
  BT_mC, // 20;
  BT_tP, // 21;
  BT_hP, // 22;
  BT_oC, // 23;
  BT_hR, // 24;
  BT_mC, // 25;
  BT_oF, // 26;
  BT_mC, // 27;
  BT_mC, // 28;
  BT_mC, // 29;
  BT_mC, // 30;
  BT_aP, // 31;
  BT_oP, // 32;
  BT_oC, // 40;
  BT_mP, // 35;
  BT_oC, // 36;
  BT_mP, // 33;
  BT_oC, // 38;
  BT_mP, // 34;
  BT_oI, // 42;
  BT_mC, // 41;
  BT_mC, // 37;
  BT_mC, // 39;
  BT_mI, // 43;  // oddity in the International tables of xtallography;
  BT_aP  // 44;
};

int ct_table[]={
  1, // 1;
  1, // 2;
  2, // 3;
  2, // 5;
  2, // 4;
  2, // 6;
  2, // 7;
  2, // 8;
  1, // 9;
  1, // 10;
  2, // 11;
  2, // 12;
  2, // 13;
  2, // 15;
  2, // 16;
  2, // 14;
  2, // 17;
  1, // 18;
  1, // 19;
  1, // 20;
  2, // 21;
  2, // 22;
  2, // 23;
  2, // 24;
  2, // 25;
  1, // 26;
  1, // 27;
  1, // 28;
  1, // 29;
  1, // 30;
  1, // 31;
  2, // 32;
  2, // 40;
  2, // 35;
  2, // 36;
  2, // 33;
  2, // 38;
  2, // 34;
  2, // 42;
  2, // 41;
  2, // 37;
  2, // 39;
  2, // 43;
  2  // 44;
};

BravaisType find_bravais_type(rMatrix3d *p_conventional, const rMatrix3d &cell) {
  rVector3d a=cell.get_column(0);
  rVector3d b=cell.get_column(1);
  rVector3d c=cell.get_column(2);
  Real A=a*a;
  Real B=b*b;
  Real C=c*c;
  Real D=b*c;
  Real E=a*c;
  Real F=a*b;
  int celltype=(D>zero_tolerance ? 1 : 2);
  rVector3d v[44];
  BravaisType bt;
  int i_min;
  int i_max;
  rVector3d DEF(D,E,F);

  v[ 0]=rVector3d( A/2., A/2., A/2.);
  v[ 1]=rVector3d( D   , D   , D   );
  v[ 2]=rVector3d( 0   , 0   , 0   );
  v[ 3]=rVector3d(-A/3.,-A/3.,-A/3.);
  v[ 4]=rVector3d( D   , D   , D   );
  v[ 5]=rVector3d( D   , D   , F   );
  v[ 6]=rVector3d( D   , E   , E   );
  v[ 7]=rVector3d( D   , E   , F   );
  
  v[ 8]=rVector3d( A/2., A/2., A/2.);
  v[ 9]=rVector3d( D   , D   , F   );
  v[10]=rVector3d( 0   , 0   , 0   );
  v[11]=rVector3d( 0   , 0   ,-A/2.);
  v[12]=rVector3d( 0   , 0   , F   );
  v[13]=rVector3d(-A/2.,-A/2., 0.  );
  v[14]=rVector3d( D   , D   , F   );
  v[15]=rVector3d( D   , D   , F   );
  v[16]=rVector3d( D   , E   , F   );
  
  v[17]=rVector3d( A/4., A/2., A/2.);
  v[18]=rVector3d( D   , A/2., A/2.);
  v[19]=rVector3d( D   , E   , E   );
  v[20]=rVector3d( 0   , 0   , 0   );
  v[21]=rVector3d(-B/2., 0   , 0   );
  v[22]=rVector3d( D   , 0   , 0   );
  v[23]=rVector3d( D   ,-A/3.,-A/3.);
  v[24]=rVector3d( D   , E   , E   );
  
  v[25]=rVector3d( A/4., A/2., A/2.);
  v[26]=rVector3d( D   , A/2., A/2.);
  v[27]=rVector3d( D   , A/2., D*2.);
  v[28]=rVector3d( D   , D*2., A/2.);
  v[29]=rVector3d( B/2., E   , E*2.);
  v[30]=rVector3d( D   , E   , F   );
  v[31]=rVector3d( 0   , 0   , 0   );
  v[32]=rVector3d(-B/2., 0   , 0   );
  v[33]=rVector3d( D   , 0   , 0   );
  v[34]=rVector3d( 0   ,-A/2., 0.  );
  v[35]=rVector3d( 0   , E   , 0.  );
  v[36]=rVector3d( 0   , 0.  ,-A/2.);
  v[37]=rVector3d( 0   , 0.  , F   );
  v[38]=rVector3d(-B/2.,-A/2., 0.  );
  v[39]=rVector3d(-B/2., E   , 0.  );
  v[40]=rVector3d( D   ,-A/2., 0.  );
  v[41]=rVector3d( D   , 0.  ,-A/2.);
  v[42]=rVector3d( D   , E   , F   );
  v[43]=rVector3d( D   , E   , F   );
  
  int i;
  for (i=0; i<44; i++) {
    if (i>=0 && i<=7 && !(near_zero(A-B) && near_zero(B-C))) continue;
    if (i>=8 && i<=16 && !near_zero(A-B)) continue;
    if (i>=17 && i<=24 && !near_zero(B-C)) continue;
    if (near_zero(norm(v[i]-DEF)) && celltype==ct_table[i]) {
      if (i==5 || i==6 || i==7 || i==14 || i==16 || i==23 || i==42) {
	if (near_zero(2*fabs(D+E+F)-(A+B))) {
	  if (i==42) {
	    if (near_zero(fabs(2*D+F)-B)) break;
	  }
	  else {
	    break;
	  }
	}
      }
      else {
	break;
      }
    }
  }
  iMatrix3d T;
  T.set(tconv[i]);
  (*p_conventional)=rMatrix3d(cell)*(~to_real(T));
  return bt_table[i];
}

rMatrix3d find_symmetric_cell(const rMatrix3d &cell) {
  Real transfo[][9]={
    { 1. , 0. , 0. ,   0. , 1. , 0. ,   0. , 0. , 1.  },
    { 0. , 0.5, 0.5,   0.5, 0. , 0.5,   0.5, 0.5, 0.  },
    {-0.5, 0.5, 0.5,   0.5,-0.5, 0.5,   0.5, 0.5,-0.5 },
    { 0.5,-0.5, 0. ,   0.5, 0.5, 0. ,   0. , 0. , 1.  },
    { 2. , 1. , 1. ,  -1. , 1. , 1. ,  -1. ,-2. , 1.  }
  };
  rMatrix3d reduced_cell=find_reduced_cell(cell);
  rMatrix3d conv_cell;
  BravaisType bt=find_bravais_type(&conv_cell, reduced_cell);
  int t;
  if (bt & BT_P) {
    t=0;
  }
  else if (bt & BT_F) {
    t=1;
  }
  else if (bt & BT_I) {
    t=2;
  }
  else if (bt & BT_C) {
    t=3;
  }
  else if (bt & BT_R) {
    t=4;
  }
  rMatrix3d T;
  T.set(transfo[t]);
  if (t==4) T=T*(1./3.);
  return conv_cell*(~T);
}

void print_bravais_type(AutoString *s, BravaisType bt) {
  char codes[]="PFICRctomah";
  for (int i=strlen(codes)-1; i>=0; i--) {
    if (bt & (1<<i)) (*s)+=codes[i];
  }
}

void find_smallest_supercell_enclosing_sphere(rMatrix3d *psupercell, const rMatrix3d cell, Real r) {
  rVector3d v[3];
  LatticePointIterator l(cell);
  while (norm(l)<r-zero_tolerance) l++;
  v[0]=l;
  while ( norm(v[0]^l)<zero_tolerance
	  || fabs((l*v[0])/norm2(v[0]))>(0.5+zero_tolerance) ) l++;
  v[1]=l;
  while (1) {
    while ( (v[0]^v[1])*l<zero_tolerance || norm(l)<r-zero_tolerance ) l++;
    MultiDimIterator<iVector2d> p(iVector2d(-1,-1),iVector2d(1,1));
    for (; p; p++) {
      rVector2d r=to_real(p);
      rVector3d u=r(0)*v[0]+r(1)*v[1];
      if (fabs(l*u/norm2(u))>0.5+zero_tolerance) break;
    }
    if (!p) break;
    l++;
  }
  v[2]=l;
  for (int i=0; i<3; i++) {psupercell->set_column(i,v[i]);}
  *psupercell=find_symmetric_cell(*psupercell);
}
