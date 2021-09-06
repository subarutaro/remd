#include <cstdlib>
#include <iostream>
#include <fstream>

#include <string>

using namespace std;

int main(int argc, char **argv){
  const int nsite = atoi(argv[1]);

  int nmol,natom;
  cin  >> nmol >> natom;
  cout << nmol << " " << nmol * nsite << endl;

  double r[3],v[3],a[4],b[4],i[4],type,index;

  for(int n=0;n<nmol;n++){
    cin >> r[0] >> r[1] >> r[2];
    cin >> v[0] >> v[1] >> v[2];
    cin >> a[0] >> a[1] >> a[2] >> a[3];
    cin >> b[0] >> b[1] >> b[2] >> b[3];
    cin >> i[0] >> i[1] >> i[2] >> i[3];
    cin >> type >> index;

    cout << r[0] << " " << r[1] << " " << r[2] << " ";
    cout << v[0] << " " << v[1] << " " << v[2] << " ";
    cout << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " ";
    cout << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " ";
    cout << i[0] << " " << i[1] << " " << i[2] << " " << i[3] << " ";
    cout << type << " " << n * nsite << endl;
  }

  double tst[3],bst[4];
  cin  >> tst[0] >> tst[1] >> tst[2];
  cout << tst[0] << " " << tst[1] << " " << tst[2]<< endl;

  cin  >> bst[0] >> bst[1] >> bst[2] >> bst[3];
  cout << bst[0] << " " << bst[1] << " " << bst[2] << " " << bst[3] << endl;

  double l[3];
  cin  >> l[0] >> l[1] >> l[2];
  cout << l[0] << " " << l[1] << " " << l[2] << endl;
}
