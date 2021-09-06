#ifndef CUH_VECTOR
#define CUH_VECTOR

/*--- for float2 ---*/
inline __host__ __device__
float2 make_float2(float s){
  return make_float2(s,s);
}
inline __host__ __device__
float2 make_float2(float3 v){
  return make_float2(v.x,v.y);
}
inline __host__ __device__
float2 operator-(float2 &v){
  return make_float2(-v.x,-v.y);
}

inline __host__ __device__
float2 operator+(float2 v1,float2 v2){
  return make_float2(v1.x + v2.x,v1.y + v2.y);
}

inline __host__ __device__
float2 operator+(float2 v1,float3 v2){
  return make_float2(v1.x + v2.x,v1.y + v2.y);
}

inline __host__ __device__
float2 operator+(float2 v,float s){
  return make_float2(v.x + s,v.y + s);
}

inline __host__ __device__
float2 operator-(float2 v1,float2 v2){
  return make_float2(v1.x - v2.x,v1.y - v2.y);
}
inline __host__ __device__
float2 operator-(float2 v,float s){
  return make_float2(v.x - s,v.y - s);
}

inline __host__ __device__
float2 operator*(float2 v1,float2 v2){
  return make_float2(v1.x * v2.x,v1.y * v2.y);
}
inline __host__ __device__
float2 operator*(float2 v,float s){
  return make_float2(v.x * s,v.y * s);
}
inline __host__ __device__
float2 operator*(float s,float2 v){
  return make_float2(v.x * s,v.y * s);
}

inline __host__ __device__
void operator+=(float2 &v1,float2 v2){
  v1.x += v2.x;
  v1.y += v2.y;
}
inline __host__ __device__
void operator+=(float2 &v1,float3 v2){
  v1.x += v2.x;
  v1.y += v2.y;
}
inline __host__ __device__
void operator+=(float2 &v,float s){
  v.x += s;
  v.y += s;
}

inline __host__ __device__
void operator-=(float2 &v1,float2 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
}
inline __host__ __device__
void operator-=(float2 &v1,float3 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
}
inline __host__ __device__
void operator-=(float2 &v,float s){
  v.x -= s;
  v.y -= s;
}

inline __host__ __device__
void operator*=(float2 &v1,float2 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
}
inline __host__ __device__
void operator*=(float2 &v,float s){
  v.x *= s;
  v.y *= s;
}

inline __host__ __device__
void operator/=(float2 &v1,float2 v2){
  v1.x /= v2.x;
  v1.y /= v2.y;
}
inline __host__ __device__
void operator/=(float2 &v,float s){
  v.x /= s;
  v.y /= s;
}

inline __host__ __device__
float norm(float2 v){
  return sqrt(v.x*v.x + v.y*v.y);
}
inline __host__ __device__
float norm2(float2 v){
  return (v.x*v.x + v.y*v.y);
}
inline __host__ __device__
float dist(float2 v1,float2 v2){
  return norm(v1-v2);
}

inline __host__ __device__
float sum(float2 v){
  return v.x + v.y;
}
/*--- for float3 ---*/
inline __host__ __device__
float3 make_float3(float s){
  return make_float3(s,s,s);
}
inline __host__ __device__
float3 make_float3(float4 v){
  return make_float3(v.x,v.y,v.z);
}
inline __host__ __device__
float3 operator-(float3 &v){
  return make_float3(-v.x,-v.y,-v.z);
}

inline __host__ __device__
float3 operator+(float3 v1,float3 v2){
  return make_float3(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z);
}

inline __host__ __device__
float3 operator+(float3 v1,float4 v2){
  return make_float3(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z);
}

inline __host__ __device__
float3 operator+(float3 v,float s){
  return make_float3(v.x + s,v.y + s,v.z + s);
}

inline __host__ __device__
float3 operator-(float3 v1,float3 v2){
  return make_float3(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z);
}
inline __host__ __device__
float3 operator-(float3 v,float s){
  return make_float3(v.x - s,v.y - s,v.z - s);
}

inline __host__ __device__
float3 operator*(float3 v1,float3 v2){
  return make_float3(v1.x * v2.x,v1.y * v2.y,v1.z * v2.z);
}
inline __host__ __device__
float3 operator*(float3 v,float s){
  return make_float3(v.x * s,v.y * s,v.z * s);
}
inline __host__ __device__
float3 operator*(float s,float3 v){
  return make_float3(v.x * s,v.y * s,v.z * s);
}

inline __host__ __device__
void operator+=(float3 &v1,float3 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
}
inline __host__ __device__
void operator+=(float3 &v1,float4 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
}
inline __host__ __device__
void operator+=(float3 &v,float s){
  v.x += s;
  v.y += s;
  v.z += s;
}

inline __host__ __device__
void operator-=(float3 &v1,float3 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
}
inline __host__ __device__
void operator-=(float3 &v1,float4 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
}
inline __host__ __device__
void operator-=(float3 &v,float s){
  v.x -= s;
  v.y -= s;
  v.z -= s;
}

inline __host__ __device__
void operator*=(float3 &v1,float3 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
}
inline __host__ __device__
void operator*=(float3 &v,float s){
  v.x *= s;
  v.y *= s;
  v.z *= s;
}

inline __host__ __device__
void operator/=(float3 &v1,float3 v2){
  v1.x /= v2.x;
  v1.y /= v2.y;
  v1.z /= v2.z;
}
inline __host__ __device__
void operator/=(float3 &v,float s){
  v.x /= s;
  v.y /= s;
  v.z /= s;
}

inline __host__ __device__
float norm(float3 v){
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
inline __host__ __device__
float norm2(float3 v){
  return (v.x*v.x + v.y*v.y + v.z*v.z);
}
inline __host__ __device__
float dist(float3 v1,float3 v2){
  return norm(v1-v2);
}

inline __host__ __device__
float sum(float3 v){
  return v.x + v.y + v.z;
}

/*--- for float4 ---*/
inline __host__ __device__
float4 make_float4(float s){
  return make_float4(s,s,s,s);
}
inline __host__ __device__
float4 make_float4(float3 v){
  return make_float4(v.x,v.y,v.z,0.f);
}
inline __host__ __device__
float4 make_float4(double4 v){
  return make_float4((float)v.x,(float)v.y,(float)v.z,(float)v.w);
}

inline __host__ __device__
float4 operator-(float4 &v){
  return make_float4(-v.x,-v.y,-v.z,-v.w);
}

inline __host__ __device__
float4 operator+(float4 v1,float4 v2){
  return make_float4(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z, v1.w + v2.w);
}
inline __host__ __device__
float4 operator+(float4 v1,float3 v2){
  return make_float4(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z, v1.w);
}
inline __host__ __device__
float4 operator+(float4 v,float s){
  return make_float4(v.x + s,v.y + s,v.z + s,v.w + s);
}

inline __host__ __device__
float4 operator-(float4 v1,float4 v2){
  return make_float4(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z,v1.w - v2.w);
}
inline __host__ __device__
float4 operator-(float4 v1,float3 v2){
  return make_float4(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z,v1.w);
}
inline __host__ __device__
float4 operator-(float4 v,float s){
  return make_float4(v.x - s,v.y - s,v.z - s,v.w - s);
}

inline __host__ __device__
float4 operator*(float4 v1,float4 v2){
  return make_float4(v1.x * v2.x,v1.y * v2.y,v1.z * v2.z,v1.w * v2.w);
}
inline __host__ __device__
float4 operator*(float4 v,float s){
  return make_float4(v.x * s,v.y * s,v.z * s,v.w * s);
}
inline __host__ __device__
float4 operator*(float s,float4 v){
  return make_float4(v.x * s,v.y * s,v.z * s,v.w * s);
}

inline __host__ __device__
void operator+=(float4 &v1,float4 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
  v1.w += v2.w;
}
inline __host__ __device__
void operator+=(float4 &v1,float3 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
}
inline __host__ __device__
void operator+=(float4 &v,float s){
  v.x += s;
  v.y += s;
  v.z += s;
  v.w += s;
}

inline __host__ __device__
void operator-=(float4 &v1,float4 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
  v1.w -= v2.w;
}
inline __host__ __device__
void operator-=(float4 &v1,float3 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
}
inline __host__ __device__
void operator-=(float4 &v,float s){
  v.x -= s;
  v.y -= s;
  v.z -= s;
  v.w -= s;
}

inline __host__ __device__
void operator*=(float4 &v1,float4 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
  v1.w *= v2.w;
}

inline __host__ __device__
void operator*=(float4 &v,float s){
  v.x *= s;
  v.y *= s;
  v.z *= s;
  v.w *= s;
}

inline __host__ __device__
void operator/=(float4 &v1,float4 v2){
  v1.x /= v2.x;
  v1.y /= v2.y;
  v1.z /= v2.z;
  v1.w /= v2.w;
}

inline __host__ __device__
void operator/=(float4 &v,float s){
  v.x /= s;
  v.y /= s;
  v.z /= s;
  v.w /= s;
}

inline __host__ __device__
float norm(float4 v){
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}
inline __host__ __device__
float norm2(float4 v){
  return (v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}
inline __host__ __device__
float dist(float4 v1,float4 v2){
  return norm(v1-v2);
}

inline __host__ __device__
float sum(float4 v){
  return v.x + v.y + v.z + v.w;
}


/*--- for double4 ---*/
inline __host__ __device__
double4 make_double4(double s){
  return make_double4(s,s,s,s);
}
inline __host__ __device__
double4 make_double4(double3 v){
  return make_double4(v.x,v.y,v.z,0.f);
}
inline __host__ __device__
double4 make_double4(float4 v){
  return make_double4((double)v.x,(double)v.y,(double)v.z,(double)v.w);
}

inline __host__ __device__
double4 operator-(double4 &v){
  return make_double4(-v.x,-v.y,-v.z,-v.w);
}

inline __host__ __device__
double4 operator+(double4 v1,double4 v2){
  return make_double4(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z, v1.w + v2.w);
}
inline __host__ __device__
double4 operator+(double4 v1,double3 v2){
  return make_double4(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z, v1.w);
}
inline __host__ __device__
double4 operator+(double4 v,double s){
  return make_double4(v.x + s,v.y + s,v.z + s,v.w + s);
}

inline __host__ __device__
double4 operator-(double4 v1,double4 v2){
  return make_double4(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z,v1.w - v2.w);
}
inline __host__ __device__
double4 operator-(double4 v1,double3 v2){
  return make_double4(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z,v1.w);
}
inline __host__ __device__
double4 operator-(double4 v,double s){
  return make_double4(v.x - s,v.y - s,v.z - s,v.w - s);
}

inline __host__ __device__
double4 operator*(double4 v1,double4 v2){
  return make_double4(v1.x * v2.x,v1.y * v2.y,v1.z * v2.z,v1.w * v2.w);
}
inline __host__ __device__
double4 operator*(double4 v,double s){
  return make_double4(v.x * s,v.y * s,v.z * s,v.w * s);
}
inline __host__ __device__
double4 operator*(double s,double4 v){
  return make_double4(v.x * s,v.y * s,v.z * s,v.w * s);
}

inline __host__ __device__
void operator+=(double4 &v1,double4 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
  v1.w += v2.w;
}
inline __host__ __device__
void operator+=(double4 &v1,double3 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
}
inline __host__ __device__
void operator+=(double4 &v,double s){
  v.x += s;
  v.y += s;
  v.z += s;
  v.w += s;
}

inline __host__ __device__
void operator-=(double4 &v1,double4 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
  v1.w -= v2.w;
}
inline __host__ __device__
void operator-=(double4 &v1,double3 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
}
inline __host__ __device__
void operator-=(double4 &v,double s){
  v.x -= s;
  v.y -= s;
  v.z -= s;
  v.w -= s;
}

inline __host__ __device__
void operator*=(double4 &v1,double4 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
  v1.w *= v2.w;
}

inline __host__ __device__
void operator*=(double4 &v,double s){
  v.x *= s;
  v.y *= s;
  v.z *= s;
  v.w *= s;
}

inline __host__ __device__
void operator/=(double4 &v1,double4 v2){
  v1.x /= v2.x;
  v1.y /= v2.y;
  v1.z /= v2.z;
  v1.w /= v2.w;
}

inline __host__ __device__
void operator/=(double4 &v,double s){
  v.x /= s;
  v.y /= s;
  v.z /= s;
  v.w /= s;
}

inline __host__ __device__
double norm(double4 v){
  return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}
inline __host__ __device__
double norm2(double4 v){
  return (v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w);
}
inline __host__ __device__
double dist(double4 v1,double4 v2){
  return norm(v1-v2);
}

inline __host__ __device__
double sum(double4 v){
  return v.x + v.y + v.z + v.w;
}

/*--- for int3 ---*/
inline __host__ __device__
int3 make_int3(int s){
  return make_int3(s,s,s);
}
inline __host__ __device__
int3 make_int3(int3 v){
  return make_int3(v.x,v.y,v.z);
}
inline __host__ __device__
int3 operator-(int3 &v){
  return make_int3(-v.x,-v.y,-v.z);
}

inline __host__ __device__
int3 operator+(int3 v1,int3 v2){
  return make_int3(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z);
}
inline __host__ __device__
int3 operator+(int3 v1,int s){
  return make_int3(v1.x + s,v1.y + s,v1.z + s);
}

inline __host__ __device__
int3 operator-(int3 v1,int3 v2){
  return make_int3(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z);
}
inline __host__ __device__
int3 operator-(int3 v1,int s){
  return make_int3(v1.x - s,v1.y - s,v1.z - s);
}

inline __host__ __device__
int3 operator*(int3 v1,int3 v2){
  return make_int3(v1.x * v2.x,v1.y * v2.y,v1.z * v2.z);
}
inline __host__ __device__
int3 operator*(int3 v1,int s){
  return make_int3(v1.x * s,v1.y * s,v1.z * s);
}
inline __host__ __device__
int3 operator*(int s,int3 v1){
  return make_int3(v1.x * s,v1.y * s,v1.z * s);
}

inline __host__ __device__
void operator+=(int3 &v1,int3 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
}
inline __host__ __device__
void operator+=(int3 &v1,int s){
  v1.x += s;
  v1.y += s;
  v1.z += s;
}

inline __host__ __device__
void operator-=(int3 &v1,int3 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
}
inline __host__ __device__
void operator-=(int3 &v1,int s){
  v1.x -= s;
  v1.y -= s;
  v1.z -= s;
}

inline __host__ __device__
void operator*=(int3 &v1,int3 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
}
inline __host__ __device__
void operator*=(int3 &v1,int s){
  v1.x *= s;
  v1.y *= s;
  v1.z *= s;
}

/*--- for int4 ---*/
inline __host__ __device__
int4 make_int4(int s){
  return make_int4(s,s,s,s);
}
inline __host__ __device__
int4 make_int4(int3 v){
  return make_int4(v.x,v.y,v.z,0);
}
inline __host__ __device__
int4 operator-(int4 &v){
  return make_int4(-v.x,-v.y,-v.z,-v.w);
}

inline __host__ __device__
int4 operator+(int4 v1,int4 v2){
  return make_int4(v1.x + v2.x,v1.y + v2.y,v1.z + v2.z, v1.w + v2.w);
}
inline __host__ __device__
int4 operator+(int4 v1,int s){
  return make_int4(v1.x + s,v1.y + s,v1.z + s,v1.w + s);
}

inline __host__ __device__
int4 operator-(int4 v1,int4 v2){
  return make_int4(v1.x - v2.x,v1.y - v2.y,v1.z - v2.z,v1.w - v2.w);
}
inline __host__ __device__
int4 operator-(int4 v1,int s){
  return make_int4(v1.x - s,v1.y - s,v1.z - s,v1.w - s);
}

inline __host__ __device__
int4 operator*(int4 v1,int4 v2){
  return make_int4(v1.x * v2.x,v1.y * v2.y,v1.z * v2.z,v1.w * v2.w);
}
inline __host__ __device__
int4 operator*(int4 v1,int s){
  return make_int4(v1.x * s,v1.y * s,v1.z * s,v1.w * s);
}
inline __host__ __device__
int4 operator*(int s,int4 v1){
  return make_int4(v1.x * s,v1.y * s,v1.z * s,v1.w * s);
}

inline __host__ __device__
void operator+=(int4 &v1,int4 v2){
  v1.x += v2.x;
  v1.y += v2.y;
  v1.z += v2.z;
  v1.w += v2.w;
}
inline __host__ __device__
void operator+=(int4 &v1,int s){
  v1.x += s;
  v1.y += s;
  v1.z += s;
  v1.w += s;
}

inline __host__ __device__
void operator-=(int4 &v1,int4 v2){
  v1.x -= v2.x;
  v1.y -= v2.y;
  v1.z -= v2.z;
  v1.w -= v2.w;
}
inline __host__ __device__
void operator-=(int4 &v1,int s){
  v1.x -= s;
  v1.y -= s;
  v1.z -= s;
  v1.w -= s;
}

inline __host__ __device__
void operator*=(int4 &v1,int4 v2){
  v1.x *= v2.x;
  v1.y *= v2.y;
  v1.z *= v2.z;
  v1.w *= v2.w;
}
inline __host__ __device__
void operator*=(int4 &v1,int s){
  v1.x *= s;
  v1.y *= s;
  v1.z *= s;
  v1.w *= s;
}



#endif
