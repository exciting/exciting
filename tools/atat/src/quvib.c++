#include "array.h"
#include <fstream.h>

Real evolve_imag_time(int n, Real stiff, Real freq, int umesh, Real tscale) {
  Real A=stiff/2;
  Real B=freq*freq/stiff/2.;
  Array<Real> psi0,psi1;
  Array<Real> org_psi(umesh);
  Real du=2./(Real)(umesh-1);

  Real u=-1.;
  for (int i=0; i<umesh; i++, u+=du) {
    org_psi(i)=sin((Real)n*M_PI/2 * (u+1))/sqrt(2);
  }
  psi1=org_psi;
  Real dt=tscale * du*du/B;
  for (Real t=0.; t<=1.; t+=dt) {
    Real u=-1.+du;
    psi0=psi1;
    for (int i=1; i<umesh-1; i++, u+=du) {
      psi1(i) -= ( A*u*u*psi0(i) - B*(psi0(i-1)+psi0(i+1)-2*psi0(i))/(2*du*du) ) * dt;
    }
/*
    {
      ofstream debug("qbh.tmp");
      Real u=-1.;
      for (int i=0; i<umesh; i++, u+=du) {
	debug << u << " " << psi1(i) << endl;
      }
      //      system("gnuplot qbh.gnu");
    }
*/
  }
  Real innerprod=0.;
  for (int i=0; i<umesh; i++, u+=du) {
    innerprod+=psi1(i)*org_psi(i)*du;
  }
  return innerprod;
}


Real calc_quantum_part_func_boxharm(Real stiff, Real freq, Real tol, int umesh, Real tscale) {
  Real Z=0.;
  for (int n=1; ; n++) {
    Real dZ=evolve_imag_time(n,stiff,freq,umesh,tscale);
    cerr << dZ << endl;
    Z+=dZ;
    if (fabs(dZ)<tol) break;
  }
  return Z;
}

int main(void) {
  for (Real s=1; s<=10.01 ; s+=1.0) {
    for (Real w=1; w<=10.01; w+=1.0) {
      cout << s << " " << w << " " << calc_quantum_part_func_boxharm(s,w,1e-3,50,0.95) << endl;
    }
    cout << endl;
  }
}
