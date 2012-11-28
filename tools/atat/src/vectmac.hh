#ifndef __VECTMAC_H__
#define __VECTMAC_H__

#define STREAM_VECTOR
#include "fxvector.h"

#define rVector2d Vector2d<Real>
#define iVector2d Vector2d<int>
#define rVector2d Vector2d<Real>
#define iVector2d Vector2d<int>
#define rMatrix2d Matrix2d<Real>
#define iMatrix2d Matrix2d<int>
#define iBoundingBox2d BoundingBox<int,2>

#define rVector3d Vector3d<Real>
#define iVector3d Vector3d<int>
#define rVector3d Vector3d<Real>
#define iVector3d Vector3d<int>
#define rMatrix3d Matrix3d<Real>
#define iMatrix3d Matrix3d<int>
#define iBoundingBox3d BoundingBox<int,3>
#define rBoundingBox3d BoundingBox<Real,3>

inline rVector2d to_real(const iVector2d &v) {
  return rVector2d((Real)v(0),(Real)v(1));
}

inline rVector3d to_real(const iVector3d &v) {
  return rVector3d((Real)v(0),(Real)v(1),(Real)v(2));
}

inline rMatrix3d to_real(const iMatrix3d &a) {
  rMatrix3d b;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b(i,j)=(Real)a(i,j);
    }
  }
  return rMatrix3d(b);
}

inline int iround(Real x) {
  return (int)rint(x);
}

inline iVector2d to_int(const rVector2d &v) {
  return iVector2d(iround(v(0)),iround(v(1)));
}

inline iVector3d to_int(const rVector3d &v) {
  return iVector3d(iround(v(0)),iround(v(1)),iround(v(2)));
}

inline iMatrix3d to_int(const rMatrix3d &a) {
  iMatrix3d b;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b(i,j)=iround(a(i,j));
    }
  }
  return b;
}

#endif
