#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <assert.h>

template <class T>
class Vector {
public:
  inline Vector(int n=0);
  inline Vector(const T *, int);
  inline Vector(const Vector<T> &);
  inline ~Vector();

  inline Vector<T> &operator = (const Vector<T> &);
  inline Vector<T> &operator = (const T &);

  inline void load(const T *, int);

  inline void resize(int);
  inline int size() const;
  inline void clear();
  
  inline const T &operator()(int) const;
  inline T &operator()(int);
  
private:
  int n_;
  T *v_;
};


template <class T>
class Matrix2D {
public:
  inline Matrix2D(int nr=0, int nc=0);
  inline Matrix2D(const T *, int nr, int nc);
  inline Matrix2D(const Matrix2D<T> &);
  inline ~Matrix2D();

  inline Matrix2D<T> &operator = (const Matrix2D<T> &);
  inline Matrix2D<T> &operator = (const T &);

  inline void load(const T *, int, int);

  inline void resize(int, int);
  inline int rows() const;
  inline int cols() const;
  inline void clear();
  
  inline const T &operator()(int, int) const;
  inline T &operator()(int, int);
  
private:
  int nr_, nc_;
  T *v_;
};


template <class T>
class Matrix3D {
public:
  inline Matrix3D(int nr=0, int nc=0, int nd=0);
  inline Matrix3D(const T *, int nr, int nc, int nd);
  inline Matrix3D(const Matrix3D<T> &);
  inline ~Matrix3D();

  inline Matrix3D<T> &operator = (const Matrix3D<T> &);
  inline Matrix3D<T> &operator = (const T &);

  inline void load(const T *, int, int, int);

  inline void resize(int, int, int);
  inline int rows() const;
  inline int cols() const;
  inline int depth() const;
  inline void clear();
  
  inline const T &operator()(int, int, int) const;
  inline T &operator()(int, int, int);
  
private:
  int nr_, nc_, nd_;
  T *v_;
};


typedef Vector<double> DVector;
typedef Matrix2D<double> DMatrix2D;
typedef Matrix3D<double> DMatrix3D;

typedef Vector<int> IVector;
typedef Matrix2D<int> IMatrix2D;
typedef Matrix3D<int> IMatrix3D;

typedef Matrix2D<unsigned char> BMatrix2D;

typedef Matrix2D<char> CMatrix2D;


template <class T>
Vector<T>::Vector(int n)
  : n_(n),
    v_(NULL)
{
  if(n_ > 0) {
    v_ = new T[n_];
  }
}

template <class T>
Vector<T>::Vector(const T *d, int n)
  : n_(0),
    v_(NULL)
{
  load(d, n);
}

template <class T>
Vector<T>::Vector(const Vector<T> &a)
  : n_(a.n_),
    v_(NULL)
{
  if(n_ > 0) {
    v_ = new T[n_];
  }
  for(int i=0; i<n_; i++) {
    v_[i] = a.v_[i];
  }
}

template <class T>
Vector<T>::~Vector() {
  if(v_) delete[] v_;
}

template <class T>
Vector<T> &
Vector<T>::operator = (const Vector<T> &a) {
  if(v_) delete[] v_;
  v_ = NULL;
  n_ = a.n_;
  if(n_ > 0) {
    v_ = new T[n_];
  }
  for(int i=0; i<n_; i++) {
    v_[i] = a.v_[i];
  }
  return *this;
}

template <class T>
Vector<T> &
Vector<T>::operator = (const T &a) {
  for(int i=0; i<n_; i++) {
    v_[i] = a;
  }
  return *this;
}

template <class T>
void
Vector<T>::load(const T *d, int n) {
  if(v_) delete[] v_;
  v_ = NULL;
  n_ = n;
  if(n_ > 0) {
    v_ = new T[n_];
  }
  for(int i=0; i<n_; i++) {
    v_[i] = d[i];
  }  
}

template <class T>
void
Vector<T>::resize(int n) {
  if(n != n_) {
    if(v_) delete[] v_;
    v_ = NULL;
    n_ = n;
    if(n_ > 0) {
      v_ = new T[n_];
    }    
  }
}

template <class T>
int
Vector<T>::size() const {
  return n_;
}

template <class T>
void
Vector<T>::clear() {
  if(v_) delete[] v_;
  v_ = NULL;
  n_ = 0;
}

template <class T>
const T &
Vector<T>::operator()(int i) const {
  assert(0 <= i && i < n_);
  return v_[i];
}

template <class T>
T &
Vector<T>::operator()(int i) {
  assert(0 <= i && i < n_);
  return v_[i];
}






template <class T>
Matrix2D<T>::Matrix2D(int nr, int nc)
  : nr_(nr),
    nc_(nc),
    v_(NULL)
{
  if(nr_ > 0 && nc_ > 0) {
    v_ = new T[nr_*nc_];
  }
}

template <class T>
Matrix2D<T>::Matrix2D(const T *d, int nr, int nc)
  : nr_(0),
    nc_(0),
    v_(NULL)
{
  load(d, nr, nc);
}

template <class T>
Matrix2D<T>::Matrix2D(const Matrix2D<T> &a)
  : nr_(a.nr_),
    nc_(a.nc_),
    v_(NULL)
{
  if(nr_ > 0 && nc_ > 0) {
    v_ = new T[nr_*nc_];
  }
  for(int i=0; i<nr_*nc_; i++) {
    v_[i] = a.v_[i];
  }
}

template <class T>
Matrix2D<T>::~Matrix2D() {
  if(v_) delete[] v_;
}

template <class T>
Matrix2D<T> &
Matrix2D<T>::operator = (const Matrix2D<T> &a) {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = a.nr_;
  nc_ = a.nc_;
  if(nr_ > 0 && nc_ > 0) {
    v_ = new T[nr_*nc_];
  }
  for(int i=0; i<nr_*nc_; i++) {
    v_[i] = a.v_[i];
  }
  return *this;
}

template <class T>
Matrix2D<T> &
Matrix2D<T>::operator = (const T &a) {
  for(int i=0; i<nr_*nc_; i++) {
    v_[i] = a;
  }
  return *this;
}

template <class T>
void
Matrix2D<T>::load(const T *d, int nr, int nc) {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = nr;
  nc_ = nc;
  if(nr_ > 0 && nc_ > 0) {
    v_ = new T[nr_*nc_];
  }
  for(int i=0; i<nr_*nc_; i++) {
    v_[i] = d[i];
  }  
}

template <class T>
void
Matrix2D<T>::resize(int nr, int nc) {
  if(nr*nc != nr_*nc_) {
    if(v_) delete[] v_;
    v_ = NULL;
    nr_ = nr;
    nc_ = nc;
    if(nr_ > 0 && nc_ > 0) {
      v_ = new T[nr_*nc_];
    }
  } else {
    nr_ = nr;
    nc_ = nc;    
  }
}

template <class T>
int
Matrix2D<T>::rows() const {
  return nr_;
}

template <class T>
int
Matrix2D<T>::cols() const {
  return nc_;
}

template <class T>
void
Matrix2D<T>::clear() {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = 0;
  nc_ = 0;
}

template <class T>
const T &
Matrix2D<T>::operator()(int r, int c) const {
  assert(0 <= r && r < nr_);
  assert(0 <= c && c < nc_);
  return v_[r*nc_+c];
}

template <class T>
T &
Matrix2D<T>::operator()(int r, int c) {
  assert(0 <= r && r < nr_);
  assert(0 <= c && c < nc_);
  return v_[r*nc_+c];
}








template <class T>
Matrix3D<T>::Matrix3D(int nr, int nc, int nd)
  : nr_(nr),
    nc_(nc),
    nd_(nd),
    v_(NULL)
{
  if(nr_ > 0 && nc_ > 0 && nd_ > 0) {
    v_ = new T[nr_*nc_*nd_];
  }
}

template <class T>
Matrix3D<T>::Matrix3D(const T *d, int nr, int nc, int nd)
  : nr_(0),
    nc_(0),
    nd_(0),
    v_(NULL)
{
  load(d, nr, nc, nd);
}

template <class T>
Matrix3D<T>::Matrix3D(const Matrix3D<T> &a)
  : nr_(a.nr_),
    nc_(a.nc_),
    nd_(a.nd_),
    v_(NULL)
{
  if(nr_ > 0 && nc_ > 0 && nd_ > 0) {
    v_ = new T[nr_*nc_*nd_];
  }
  for(int i=0; i<nr_*nc_*nd_; i++) {
    v_[i] = a.v_[i];
  }
}

template <class T>
Matrix3D<T>::~Matrix3D() {
  if(v_) delete[] v_;
}

template <class T>
Matrix3D<T> &
Matrix3D<T>::operator = (const Matrix3D<T> &a) {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = a.nr_;
  nc_ = a.nc_;
  nd_ = a.nd_;
  if(nr_ > 0 && nc_ > 0 && nd_ > 0) {
    v_ = new T[nr_*nc_*nd_];
  }
  for(int i=0; i<nr_*nc_*nd_; i++) {
    v_[i] = a.v_[i];
  }
  return *this;
}

template <class T>
Matrix3D<T> &
Matrix3D<T>::operator = (const T &a) {
  for(int i=0; i<nr_*nc_*nd_; i++) {
    v_[i] = a;
  }
  return *this;
}

template <class T>
void
Matrix3D<T>::load(const T *d, int nr, int nc, int nd) {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = nr;
  nc_ = nc;
  nd_ = nd;
  if(nr_ > 0 && nc_ > 0 && nd_ > 0) {
    v_ = new T[nr_*nc_*nd_];
  }
  for(int i=0; i<nr_*nc_*nd_; i++) {
    v_[i] = d[i];
  }  
}

template <class T>
void
Matrix3D<T>::resize(int nr, int nc, int nd) {
  if(nr*nc != nr_*nc_*nd_) {
    if(v_) delete[] v_;
    v_ = NULL;
    nr_ = nr;
    nc_ = nc;
    nd_ = nd;
    if(nr_ > 0 && nc_ > 0 && nd_ > 0) {
      v_ = new T[nr_*nc_*nd_];
    }
  } else {
    nr_ = nr;
    nc_ = nc;
    nd_ = nd;
  }
}

template <class T>
int
Matrix3D<T>::rows() const {
  return nr_;
}

template <class T>
int
Matrix3D<T>::cols() const {
  return nc_;
}

template <class T>
int
Matrix3D<T>::depth() const {
  return nd_;
}

template <class T>
void
Matrix3D<T>::clear() {
  if(v_) delete[] v_;
  v_ = NULL;
  nr_ = 0;
  nc_ = 0;
  nd_ = 0;
}

template <class T>
const T &
Matrix3D<T>::operator()(int r, int c, int d) const {
  assert(0 <= r && r < nr_);
  assert(0 <= c && c < nc_);
  assert(0 <= d && d < nd_);
  return v_[r*nc_*nd_+c*nd_+d];
}

template <class T>
T &
Matrix3D<T>::operator()(int r, int c, int d) {
  assert(0 <= r && r < nr_);
  assert(0 <= c && c < nc_);
  assert(0 <= d && d < nd_);
  return v_[r*nc_*nd_+c*nd_+d];
}

#endif
