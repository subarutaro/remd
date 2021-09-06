#include <iostream>
#include <iomanip>

#include <cmath>
#include <cstring>

class Vec3{
public:
  union{
    double data[3];
    struct{
      double x;
      double y;
      double z;
    };
  };

  Vec3(){}
  Vec3(const double _x,const double _y,const double _z) : x(_x),y(_y),z(_z) {}
  ~Vec3(){}

  //vector arithmetics
  const Vec3 operator+(const Vec3& v) const {
    Vec3 tmp;
    for(int i=0;i<3;i++) tmp[i] = data[i] + v[i];
    return tmp;
  }
  const Vec3 operator-(const Vec3& v) const {
    Vec3 tmp;
    for(int i=0;i<3;i++) tmp[i] = data[i] - v[i];
    return tmp;
  }
  const Vec3 operator*(const Vec3& v) const {
    Vec3 tmp;
    for(int i=0;i<3;i++) tmp[i] = data[i] * v[i];
    return tmp;
  }
  const Vec3 operator/(const Vec3& v) const {
    Vec3 tmp;
    for(int i=0;i<3;i++) tmp[i] = data[i] / v[i];
    return tmp;
  }
  // scalar arithmetics
  const Vec3 operator+(const double& s) const {
    Vec3 v;
    for(int i=0;i<3;i++) v[i] = data[i] + s;
    return v;
  }
  const Vec3 operator-(const double& s) const {
    Vec3 v;
    for(int i=0;i<3;i++) v[i] = data[i] - s;
    return v;
  }
  const Vec3 operator*(const double& s) const {
    Vec3 v;
    for(int i=0;i<3;i++) v[i] = data[i] * s;
    return v;
  }
  const Vec3 operator/(const double& s) const {
    Vec3 v;
    for(int i=0;i<3;i++) v[i] = data[i] / s;
    return v;
  }

  const bool operator==(const Vec3 v) const {
    bool isEqual = true;
    for(int i=0;i<3;i++) isEqual &= (data[i] == v[i]);
    return isEqual;
  }
  const bool operator!=(const Vec3 v) const {
    return !(*this == v);
  }
  const bool operator<(const Vec3 v) const {
    return z < v.z;
  }
  const bool operator>(const Vec3 v) const {
    return z > v.z;
  }
  const bool operator<=(const Vec3 v) const {
    return z <= v.z;
  }
  const bool operator>=(const Vec3 v) const {
    return z >= v.z;
  }

  double& operator[](const int i){return data[i];}
  const double& operator[](const int i) const {return data[i];}

  friend const double sum(const Vec3& v){
    double s = 0.0;
    for(int i=0;i<3;i++) s += v[i];
    return s;
  }
  friend const double dot(const Vec3& v0,const Vec3& v1){
    return sum(v0 * v1);
  }
  friend const Vec3 cross(const Vec3& v0,const Vec3& v1){
     return Vec3(v0[1]*v1[2] - v0[2]*v1[1],
		 v0[2]*v1[0] - v0[0]*v1[2],
		 v0[0]*v1[1] - v0[1]*v1[0]);
  }
  friend const double norm(const Vec3 v){
    return sqrt(sum(v*v));
  }

  friend std::istream& operator>>(std::istream &is,Vec3 &v){
    for(int i=0;i<3;i++) is >> v[i];
    return is;
  }
  friend std::ostream& operator<<(std::ostream &os,const Vec3 &v){
    for(int i=0;i<3;i++) os << ' ' << v[i];
    return os;
  }
};
