#ifndef TYPEDEF_H
#define TYPEDEF_H

#define IntMax    (1LL<<32)
#define IntMaxInv (1./(D)IntMax)
#define IntMaxInvF (1.f/(float)IntMax)

#define I int
#define D double
/*
typedef struct{
  I x;
  I y;
  I z;
}I3;

typedef struct{
  D x;
  D y;
  D z;
}D3;
//*/
/*
typedef struct{
  int x;
  int y;
  int z;
  int w;
}int4;

typedef struct{
  float x;
  float y;
  float z;
  float w;
}float4;
//*/

typedef struct{
  D pot;
  D vir;
  D cnst_p;
  D cnst_v;
}LRC;

#endif







