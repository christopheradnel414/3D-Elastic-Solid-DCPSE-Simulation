#ifndef READGEOM_HPP
#define READGEOM_HPP

#include <vector>
#include <string>

using namespace std;

void readgeom  (double& E,
                double& v,
                double& h,
                vector<double>& x,
                vector<double>& y,
                vector<double>& z,
                vector<double>& nx,
                vector<double>& ny,
                vector<double>& nz,
                vector<int>& nxi,
                vector<int>& nyi,
                vector<int>& nzi,
                vector<double>& dispX,
                vector<double>& dispY,
                vector<double>& dispZ,
                vector<int>& dispXi,
                vector<int>& dispYi,
                vector<int>& dispZi,
                vector<double>& forceX,
                vector<double>& forceY,
                vector<double>& forceZ,
                vector<int>& forceXi,
                vector<int>& forceYi,
                vector<int>& forceZi,
                vector<double>& bx,
                vector<double>& by,
                vector<double>& bz,
                string geomfilename,
                string geomparamfilename);

#endif