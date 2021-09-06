#ifndef H_VECTOR
#define H_VECTOR

#include <iomanip>
#include <ostream>
#include <cmath>

template< int N, typename T>
class Vector{
 private:
  T data[N];
 public:
  Vector(){ for(int i=0;i<N;i++) data[i] = (T)0; }
  Vector(const T &v){ for(int i=0;i<N;i++) data[i] = v; }

  Vector(const Vector &v){ for(int i=0;i<N;i++) data[i] = v[i]; }
  ~Vector(){}

  //for scalar
  const Vector &operator=(const T v) {
    for(int i=0;i<N;i++) data[i] = v;
    return *this;
  }
  const Vector &operator+=(const T v) {
    for(int i=0;i<N;i++) data[i] += v;
    return *this;
  }
  const Vector &operator-=(const T v) {
    for(int i=0;i<N;i++) data[i] -= v;
    return *this;
  }
  const Vector &operator*=(const T v) {
    for(int i=0;i<N;i++) data[i] *= v;
    return *this;
  }
  const Vector &operator/=(const T v) {
    for(int i=0;i<N;i++) data[i] /= v;
    return *this;
  }
  Vector operator+(const T v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] + v;
    return tmp;
  }
  Vector operator-(const T v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] - v;
    return tmp;
  }
  Vector operator*(const T v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] * v;
    return tmp;
  }
  Vector operator/(const T v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] / v;
    return tmp;
  }

  //for vector
  const Vector &operator=(const Vector v) {
    for(int i=0;i<N;i++) data[i] = v[i];
    return *this;
  }
  const Vector &operator+=(const Vector v) {
    for(int i=0;i<N;i++) data[i] += v[i];
    return *this;
  }
  const Vector &operator-=(const Vector v) {
    for(int i=0;i<N;i++) data[i] -= v[i];
    return *this;
  }
  const Vector &operator*=(const Vector v) {
    for(int i=0;i<N;i++) data[i] *= v[i];
    return *this;
  }
  const Vector &operator/=(const Vector v) {
    for(int i=0;i<N;i++) data[i] /= v[i];
    return *this;
  }
  Vector operator+(const Vector v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] + v[i];
    return tmp;
  }
  Vector operator-(const Vector v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] - v[i];
    return tmp;
  }
  Vector operator*(const Vector v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] * v[i];
    return tmp;
  }
  Vector operator/(const Vector v) const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = data[i] / v[i];
    return tmp;
  }

  //others
  Vector operator-() const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = -data[i];
    return tmp;
  }
  T& operator[](int i){ return data[i];}
  const T &operator[](int i) const { return data[i];}
  operator       T* ()      { return data;}
  operator const T* ()const { return data;}

  //cast
  template<class T2>
  operator Vector<N,T2> () {
    Vector<N,T2> tmp;
    for(int i=0;i<N;i++){
      tmp[i] = (T2)data[i];
    }
    return tmp;
  }
  template<class T2>
  operator Vector<N,T2> () const {
    Vector<N,T2> tmp;
    for(int i=0;i<N;i++){
      tmp[i] = (T2)data[i];
    }
    return tmp;
  }


  //inverse
  const Vector inv() const {
    Vector tmp;
    for(int i=0;i<N;i++) tmp[i] = (T)1 / data[i];
    return tmp;
  }

  // ostream
  friend std::ostream &operator<<(std::ostream &s, const Vector &v){
    for(int i=0;i<N;i++){
      //s << std::setw(10);
      //s << ' ' << std::fixed << v[i];
      s << ' ' << std::scientific << v[i];
    }
    return s;
  }

  // math lib
  friend T sum(const Vector &v){
    T tmp = (T)0;
    for(int i=0;i<N;i++) tmp += v[i];
    return tmp;
  }
  friend T norm(const Vector &v){
    T tmp = (T)0;
    for(int i=0;i<N;i++) tmp += v[i]*v[i];
    return sqrt(tmp);
  }
  friend T norm2(const Vector &v){
    T tmp = (T)0;
    for(int i=0;i<N;i++) tmp += v[i]*v[i];
    return tmp;
  }
};

typedef Vector<3,double>  dvec2;
typedef Vector<3,double>  dvec3;
typedef Vector<4,double>  dvec4;
typedef Vector<3,float>	  fvec2;
typedef Vector<3,float>	  fvec3;
typedef Vector<4,float>	  fvec4;
typedef Vector<3,int>	  ivec2;
typedef Vector<3,int>	  ivec3;
typedef Vector<4,int>	  ivec4;

typedef Vector< 3, Vector<3,double> > dvec33;
typedef Vector< 4, Vector<4,double> > dvec44;

#endif
