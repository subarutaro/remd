extern "C" {
  void getbitpattern_float(float);
  void get_cputime(double*, double*);
  int getIntegerMulOf(int, int);
  double getbspline(double, int);
  float getbsplinef(float, int);
  void getcoefpme(double*, int, int);
  void getcoefpmef(float*, int, int);
  void getInterpolationCurve(double, int, double*, int, int);
  void getInterpolationCurvef(float, int, float*, int, int);
  void myfft1d(int, double, double*);
  void myfft1df(int, float, float*);
  void myfft3d(int, int, int, double*, double);
  void myfft3df(int, int, int, float*, float);
  double minvald(double*, int);
  void initCellindex(double*, double, int*, int*, int, int*, int*);
  void showComptime(double*, int, int, const char*);
  void ffttest(void);
  void getNeighborCell(double*, double, int*, int*, int*, int);
  void initrandomuni(unsigned long, int);
  double getrandomuni(void);
}
